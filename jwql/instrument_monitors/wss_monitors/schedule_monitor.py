#! /usr/bin/env python

"""This module contains code for the WSS schedule monitor.

Author
------
    - Ben Sunnquist

Use
---
    This module can be used from the command line as such:

    ::

        python schedule_monitor.py

    To only update the background trending plots:

    ::

        m = ClawMonitor()
        m.make_background_plots()
"""

import datetime
import logging
import os
import warnings

from astropy.convolution import Gaussian2DKernel, convolve
from astropy.io import fits
from astropy.stats import gaussian_fwhm_to_sigma, sigma_clipped_stats
from astropy.table import Table
from astropy.time import Time
from astroquery.mast import Mast
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from photutils.segmentation import detect_sources, detect_threshold
from scipy.ndimage import binary_dilation

from jwql.utils import monitor_utils
from jwql.utils.constants import ON_GITHUB_ACTIONS, ON_READTHEDOCS
from jwql.utils.logging_functions import log_info, log_fail
from jwql.utils.utils import ensure_dir_exists, filesystem_path, get_config
from jwst_backgrounds import jbt

if not ON_GITHUB_ACTIONS and not ON_READTHEDOCS:
    # Need to set up django apps before we can access the models
    import django  # noqa: E402 (module level import not at top of file)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "jwql.website.jwql_proj.settings")
    django.setup()

    # Import * is okay here because this module specifically only contains database models
    # for this monitor
    from jwql.website.apps.jwql.monitor_models.schedule import *  # noqa: E402 (module level import not at top of file)

matplotlib.use('Agg')
warnings.filterwarnings('ignore', message="nan_treatment='interpolate', however, NaN values detected post convolution*")
warnings.filterwarnings('ignore', message='Input data contains invalid values (NaNs or infs)*')


class ScheduleMonitor():
    """Class for executing the wss-schedule monitor.

    TBD (based on claw_monitor.py)
    """

    def __init__(self):
        """Initialize an instance of the ``ClawMonitor`` class.
        """

        # Define and setup the output directories for the claw and background plots.
        self.output_dir_schedule = os.path.join(get_config()['outputs'], 'wss', 'schedule')
        ensure_dir_exists(self.output_dir_schedule)

        # Get the claw monitor database tables
        self.query_table = WSSScheduleQueryHistory

    def make_some_plots(self, plot_type='bkg'):
        """Makes plots of ...

        TBD (see claw_monitor.py for instance)
        """

    def process(self):
        """The main method for processing.  See module docstrings for further details.

        TBD (see claw_monitor.py for instance)
        """

    def query_mast(self):
        """Query MAST for new imagind data, any instrument

        Returns
        -------
        t : astropy.table.table.Table
            A table summarizing the new nircam imaging data.
        """

        server = "https://mast.stsci.edu"
        JwstObs = Mast()
        JwstObs._portal_api_connection.MAST_REQUEST_URL = server + "/portal_jwst/Mashup/Mashup.asmx/invoke"
        JwstObs._portal_api_connection.MAST_DOWNLOAD_URL = server + "/jwst/api/v0.1/download/file"
        JwstObs._portal_api_connection.COLUMNS_CONFIG_URL = server + "/portal_jwst/Mashup/Mashup.asmx/columnsconfig"
        JwstObs._portal_api_connection.MAST_BUNDLE_URL = server + "/jwst/api/v0.1/download/bundle"
        service = 'Mast.Jwst.Filtered.Nircam'
        FIELDS = ['filename', 'program', 'observtn', 'category', 'instrume', 'productLevel', 'filter',
                  'pupil', 'subarray', 'detector', 'datamodl', 'date_beg_mjd', 'effexptm']
        params = {"columns": ",".join(FIELDS),
                  "filters": [{"paramName": "pupil", "values": ['CLEAR', 'F162M', 'F164N', 'F323N', 'F405N', 'F466N', 'F470N']},
                              {"paramName": "exp_type", "values": ['NRC_IMAGE']},
                              {"paramName": "datamodl", "values": ['ImageModel']},  # exclude calints, which are cubemodel
                              {"paramName": "productLevel", "values": ['2b']},  # i.e. cal.fits
                              {"paramName": "subarray", "values": ['FULL']}, ]
                  }
        t = JwstObs.service_request(service, params)
        t = t[(t['date_beg_mjd'] > self.query_start_mjd) & (t['date_beg_mjd'] < self.query_end_mjd)]
        t.sort('date_beg_mjd')
        filetypes = np.array([row['filename'].split('_')[-1].replace('.fits', '') for row in t])
        t = t[filetypes == 'cal']  # only want cal.fits files, no e.g. i2d.fits

        return t

    @log_fail
    @log_info
    def run(self):
        """The main method.  See module docstrings for further details."""

        logging.info('Begin logging for schedule_monitor')

        # Query MAST for new NIRCam full-frame imaging data from the last 2 days
        self.query_end_mjd = Time.now().mjd
        self.query_start_mjd = self.query_end_mjd - 2
        mast_table = self.query_mast()
        logging.info('{} files found between {} and {}.'.format(len(mast_table), self.query_start_mjd, self.query_end_mjd))

        # Define pivot wavelengths
        self.filter_wave = {'F070W': 0.704, 'F090W': 0.902, 'F115W': 1.154, 'F150W': 1.501, 'F150W2': 1.659,
                            'F200W': 1.989, 'F212N': 2.121, 'F250M': 2.503, 'F277W': 2.762, 'F300M': 2.989,
                            'F322W2': 3.232, 'F356W': 3.568, 'F410M': 4.082, 'F430M': 4.281, 'F444W': 4.408,
                            'F480M': 4.874}

        # Create observation-level median stacks for each filter/pupil combo, in pixel-space
        combos = np.array(['{}_{}_{}_{}'.format(str(row['program']), row['observtn'], row['filter'], row['pupil']).lower() for row in mast_table])
        mast_table['combos'] = combos
        monitor_run = False
        for combo in np.unique(combos):
            mast_table_combo = mast_table[mast_table['combos'] == combo]
            if 'long' in mast_table_combo['filename'][0]:
                self.channel = 'LW'
            else:
                self.channel = 'SW'
            self.proposal, self.obs, self.fltr, self.pupil = combo.split('_')
            self.outfile = os.path.join(self.output_dir_claws, 'prop{}_obs{}_{}_{}_cal_norm_skyflat.png'.format(str(self.proposal).zfill(5),
                                        self.obs, self.fltr, self.pupil).lower())
            existing_files = []
            for row in mast_table_combo:
                try:
                    existing_files.append(filesystem_path(row['filename']))
                except Exception as e:
                    pass
            self.files = np.array(existing_files)
            self.detectors = np.array(mast_table_combo['detector'])
            if (not os.path.exists(self.outfile)) & (len(existing_files) == len(mast_table_combo)):
                logging.info('Working on {}'.format(self.outfile))
                self.process()
                monitor_run = True
            else:
                logging.info('{} already exists or is missing cal files ({}/{} files found).'.format(self.outfile, len(existing_files), len(mast_table_combo)))

        # Update the background trending plots, if any new data exists
        if len(mast_table) > 0:
            logging.info('Making background trending plots.')
            self.make_background_plots(plot_type='bkg')
            self.make_background_plots(plot_type='bkg_rms')
            self.make_background_plots(plot_type='model')

        # Update the query history
        new_entry = {'instrument': 'nircam',
                     'start_time_mjd': self.query_start_mjd,
                     'end_time_mjd': self.query_end_mjd,
                     'run_monitor': monitor_run,
                     'entry_date': datetime.datetime.now()}
        entry = self.query_table(**new_entry)
        entry.save()

        logging.info('Claw Monitor completed successfully.')


if __name__ == '__main__':

    module = os.path.basename(__file__).strip('.py')
    start_time, log_file = monitor_utils.initialize_instrument_monitor(module)

    monitor = ScheduleMonitor()
    monitor.run()

    monitor_utils.update_monitor_table(module, start_time, log_file)
