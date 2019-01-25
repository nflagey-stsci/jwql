"""Various utility functions for the ``jwql`` project.

Authors
-------

    - Bryan Hilbert

Use
---

    This module can be imported as such:

    >>> from jwql.instrument_monitors import pipeline_tools
    pipeline_steps = pipeline_tools.completed_pipeline_steps(filename)
 """

from collections import OrderedDict
import copy
import numpy as np

from astropy.io import fits
from jwst import datamodels
from jwst.dq_init import DQInitStep
from jwst.dark_current import DarkCurrentStep
from jwst.firstframe import FirstFrameStep
from jwst.group_scale import GroupScaleStep
from jwst.ipc import IPCStep
from jwst.jump import JumpStep
from jwst.lastframe import LastFrameStep
from jwst.linearity import LinearityStep
from jwst.persistence import PersistenceStep
from jwst.ramp_fitting import RampFitStep
from jwst.refpix import RefPixStep
from jwst.rscd import RSCD_Step
from jwst.saturation import SaturationStep
from jwst.superbias import SuperBiasStep

from jwql.utils.utils import JWST_INSTRUMENTS


# Define the fits header keyword that accompanies each step
PIPE_KEYWORDS = {'S_GRPSCL': 'group_scale', 'S_DQINIT': 'dq_init', 'S_SATURA': 'saturation',
                 'S_IPC': 'ipc', 'S_REFPIX': 'refpix', 'S_SUPERB': 'superbias',
                 'S_PERSIS': 'persistence', 'S_DARK': 'dark_current', 'S_LINEAR': 'linearity',
                 'S_FRSTFR': 'firstframe', 'S_LASTFR': 'lastframe', 'S_RSCD': 'rscd',
                 'S_JUMP': 'jump', 'S_RAMP': 'rate'}

PIPELINE_STEP_MAPPING = {'dq_init': DQInitStep, 'dark_current': DarkCurrentStep,
                         'firstframe': FirstFrameStep, 'group_scale': GroupScaleStep,
                         'ipc': IPCStep, 'jump': JumpStep, 'lastframe': LastFrameStep,
                         'linearity': LinearityStep, 'persistence': PersistenceStep,
                         'rate': RampFitStep, 'refpix': RefPixStep, 'rscd': RSCD_Step,
                         'saturation': SaturationStep, 'superbias': SuperBiasStep}

print('Remove line below before merging. Uppercase list of inst. is in a PR already')
JWST_INSTRUMENTS = [entry.upper() for entry in JWST_INSTRUMENTS]


def completed_pipeline_steps(filename):
    """
    Return a list of the completed pipeline steps for a given file.

    Parameters
    ----------
    filename : str
        File to examine

    Returns
    -------
    completed : collections.OrderedDict
        Dictionary with boolean entry for each pipeline step,
        indicating which pipeline steps have been run on filename
    """
    # Initialize using PIPE_KEYWORDS so that entries are guaranteed to
    # be in the correct order
    completed = OrderedDict({})
    for key in PIPE_KEYWORDS.values():
        completed[key] = False

    header = fits.getheader(filename)
    for key in PIPE_KEYWORDS.keys():
        value = header.get(key)
        if value == 'COMPLETE':
            completed[PIPE_KEYWORDS[key]] = True
    return completed


def get_pipeline_steps(instrument):
    """Get the names and order of the calwebb_detector1
    pipeline steps for a given instrument. Use values that match up with the values in the
    PIPE_STEP defintion in definitions.py

    Parameters
    ----------

    instrument : str
        Name of JWST instrument

    Returns
    -------

    steps : collections.OrderedDict
        Dictionary of step names (and modules? do we care?)
    """
    instrument = instrument.upper()
    if instrument not in JWST_INSTRUMENTS:
        raise ValueError("WARNING: {} is not a valid instrument name.".format(instrument))
    # all_steps = Detector1Pipeline.step_defs

    # Order is important in 'steps' lists below!!
    if instrument == 'MIRI':
        steps = ['group_scale', 'dq_init', 'saturation', 'ipc', 'firstframe', 'lastframe',
                 'linearity', 'rscd', 'dark_current', 'refpix', 'persistence', 'jump', 'rate']
        # No persistence correction for MIRI
        steps.remove('persistence')
    else:
        steps = ['group_scale', 'dq_init', 'saturation', 'ipc', 'superbias', 'refpix', 'linearity',
                 'persistence', 'dark_current', 'jump', 'rate']

        # No persistence correction for NIRSpec
        if instrument == 'NIRSPEC':
            steps.remove('persistence')

    # IPC correction currently not done for any instrument
    steps.remove('ipc')

    # Initialize using PIPE_KEYWORDS so the steps will be in the right
    # order
    req = OrderedDict({})
    for key in steps:
        req[key] = True
    for key in PIPE_KEYWORDS.values():
        if key not in req.keys():
            req[key] = False

    return req


def image_stack(file_list):
    """Given a list of fits files containing 2D images, read in all data
    and place into a 3D stack

    Parameters
    ----------
    file_list : list
        List of fits file names

    Returns
    -------
    cube : numpy.ndarray
        3D stack of the 2D images
    """
    for i, input_file in enumerate(file_list):
        model = datamodels.open(input_file)
        image = model.data

        # Stack all inputs together into a single 3D image cube
        if i == 0:
            ndim_base = image.shape
            if len(ndim_base) == 3:
                cube = copy.deepcopy(image)
            elif len(ndim_base) == 2:
                cube = np.expand_dims(image, 0)

        ndim = image.shape
        if ndim_base[-2:] == ndim[-2:]:
            if len(ndim) == 2:
                image = np.expand_dims(image, 0)
            elif len(ndim) > 3:
                #raise ValueError("4-dimensional input images not supported.")
                print('')
                print('using the initial frame for early testing!!!!!')
                print('remove line below before merging!!')
                image = np.expand_dims(image[0, :, :], 0)
            cube = np.vstack((cube, image))
        else:
            raise ValueError("Input images are of inconsistent size in x/y dimension.")
    return cube


def run_calwebb_detector1_steps(input_file, steps):
    """Run the steps of calwebb_detector1 specified in the steps
    dictionary on the input file

    Parameters
    ----------
    input_file : str
        File on which to run the pipeline steps

    steps : collections.OrderedDict
        Keys are the individual pipeline steps (as seen in the
        PIPE_KEYWORDS values above). Boolean values indicate
        whether a step should be run or not. Steps are run in the
        official calwebb_detector1 order.
    """
    first_step_to_be_run = True
    for step_name in steps.keys():
        if steps[step_name]:
            if first_step_to_be_run:
                model = PIPELINE_STEP_MAPPING[step_name].call(input_file)
                first_step_to_be_run = False
            else:
                model = PIPELINE_STEP_MAPPING[step_name].call(model)
            suffix = step_name
    output_filename = input_file.replace('.fits', '_{}.fits'.format(suffix))
    print('in run_calwebb_detector1_steps, input and output filenames are:')
    print(input_file)
    print(output_filename)
    print('STEPS RUN:', steps)
    if suffix != 'rate':
        model.save(output_filename)
    else:
        model[0].save(output_filename)


def steps_to_run(input_file, all_steps, finished_steps):
    """Given a list of pipeline steps that need to be completed as
    well as a list of steps that have already been completed, return
    a list of steps remaining to be done.

    Parameters
    ----------
    all_steps : collections.OrderedDict
        A dictionary of all steps that need to be completed

    finished_steps : collections.OrderedDict
        A dictionary with keys equal to the pipeline steps and boolean
        values indicating whether a particular step has been completed
        or not (i.e. output from completed_pipeline_steps)

    Returns
    -------
    steps_to_run : collections.OrderedDict
        A dictionaru with keys equal to the pipeline steps and boolean
        values indicating whether a particular step has yet to be run.
    """
    torun = copy.deepcopy(finished_steps)
    for key in all_steps:
        if all_steps[key] == finished_steps[key]:
            torun[key] = False
        elif ((all_steps[key] is True) & (finished_steps[key] is False)):
            torun[key] = True
        elif ((all_steps[key] is False) & (finished_steps[key] is True)):
            print(("WARNING: Input file {} has had {} step run, "
                   "but the requirements say that it should not "
                   "be. Need a new input file.".format(input_file, key)))
    return torun
