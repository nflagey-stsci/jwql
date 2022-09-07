"""Script that can be used to query MAST and return basic info
about all proposals. This information is used to help populate
the instrument archive pages

Authors
-------

    - Bryan Hilbert

Use
---

    This module is called as follows:
    ::

        from jwql.websites.apps.jwql.archvie_database_update import get_updates
        instrument = 'nircam'
        get_updates(insturument)

Dependencies
------------
    The user must have a configuration file named ``config.json``
    placed in the ``jwql`` directory.
"""

from jwql.utils.constants import JWST_INSTRUMENT_NAMES_MIXEDCASE
from jwql.utils.utils import filename_parser, filesystem_path, get_config

from jwql.website.apps.jwql.data_containers import get_instrument_proposals, get_filenames_by_instrument
from jwql.website.apps.jwql.data_containers import get_proposal_info, mast_query_filenames_by_instrument


def get_updates(inst):
    """Generate the page listing all archived proposals in the database

    Parameters
    ----------
    inst : str
        Name of JWST instrument

    Returns
    -------

    """
    # Ensure the instrument is correctly capitalized
    inst = JWST_INSTRUMENT_NAMES_MIXEDCASE[inst.lower()]
    filesystem = get_config()['filesystem']

    # Dictionary to hold summary information for all proposals
    all_proposals = get_instrument_proposals(inst)
    all_proposal_info = {'num_proposals': 0,
                         'proposals': [],
                         'min_obsnum': [],
                         'thumbnail_paths': [],
                         'num_files': []}

    # Dictionary to hold summary observation of exp_types
    exp_types = {}

    # Get list of all files for the given instrument
    for proposal in all_proposals:
        # Get lists of all public and proprietary files for the program
        filenames_public, metadata_public, filenames_proprietary, metadata_proprietary = get_all_possible_filenames_for_proposal(inst, proposal)

        # Find the location in the filesystem for all files
        filepaths_public = files_in_filesystem(filenames_public, 'public')
        filepaths_proprietary = files_in_filesystem(filenames_proprietary, 'proprietary')
        filenames = filepaths_public + filepaths_proprietary

        # Get set of unique rootnames
        num_files = 0
        rootnames = set(['_'.join(f.split('/')[-1].split('_')[:-1]) for f in filenames])
        for rootname in rootnames:
            filename_dict = filename_parser(rootname)

            # Weed out file types that are not supported by generate_preview_images
            if 'stage_3' not in filename_dict['filename_type']:
                num_files += 1

        if len(filenames) > 0:

            # Gather information about the proposals for the given instrument
            proposal_info = get_proposal_info(filenames)

            # Each observation number in each proposal can have a list of exp_types (e.g. NRC_TACQ, NRC_IMAGE)
            for obsnum in proposal_info['observation_nums']:
                match_pub = metadata_public['observtn'] == int(obsnum)
                exp_types = list(set(metadata_public['observtn'][match_pub]))
                match_prop = metadata_proprietary['observtn'] == int(obsnum)
                exp_types = list(set(exp_types + metadata_proprietary['observtn'][match_prop]))
                all_start_dates = np.array(metadata_public['expstart'][match_pub] + metadata_proprietary['expstart'][match_prop])
                starting_date = np.min(all_start_dates)
                all_end_dates = np.array(metadata_public['expend'][match_pub] + metadata_proprietary['expend'][match_prop])
                latest_date = np.max(all_end_dates)

                # Update the appropriate database table
                update_database_table(inst, proposal, obsnum, proposal_info['thumbnail_paths'][0], num_files,
                                      exp_types, starting_date, latest_date)


def get_all_possible_filenames_for_proposal(instrument, proposal_num):
    """
    """
    filename_query = mast_query_filenames_by_instrument(instrument, proposal_num,
                                                        other_columns=['exp_type', 'observtn', 'expstart', 'expend'])
    public, public_meta = get_filenames_by_instrument(instrument, proposal_num, restriction='public',
                                                      query_response=filename_query,
                                                      other_columns=['exp_type', 'observtn', 'expstart', 'expend'])
    proprietary, proprietary_meta = get_filenames_by_instrument(instrument, proposal_num, restriction='proprietary',
                                                                query_response=filename_query,
                                                                other_columns=['exp_type', 'observtn', 'expstart', 'expend'])
    return public, public_meta, proprietary, proprietary_meta


def files_in_filesystem(files, permission_type):
    """Determine locations in the filesystem for the input files

    Parameters
    ----------
    files : list
        List of filenames from MAST query

    permission_type : str
        Permission level of the input files: 'public' or 'proprietary'

    Return
    ------
    filenames : list
        List of full paths within the filesystem for the input files
    """
    if permission_type not in ['public', 'proprietary']:
        raise ValueError('permission type needs to be either "public" or "proprietary"')

    filenames = []
    for filename in files:
        try:
            relative_filepath = filesystem_path(filename, check_existence=False)
            full_filepath = os.path.join(filesystem, permission_type, relative_filepath)
            filenames.append(full_filepath)
        except ValueError:
            print('Unable to determine filepath for {}'.format(filename))
    return filenames


def update_database_table(instrument, prop, obs, thumbnail, files, types, startdate, enddate):
    """

    Parameters
    ----------
    instrument : str
        Instrument name

    prop : str
        Proposal ID. 5-digit string

    obs : str
        Observation number. 3-digit string

    thumbnail : str
        Full path to the thumbnail image for the proposal

    files : int
        Number of files in the observation

    types : list
        List of exposure types of the data in the observation

    startdate : float
        Date of the beginning of the observation in MJD

    enddate : float
        Date of the ending of the observation in MJD
    """
    logging.info('')

    # First, check to see if ExposureType instances exist for all of the elements
    # in types. For any that are missing, create and save an entry to the db.
    for etype in types:
        type_instance = ExposureType.object.get(exp_type=etype)
        if len(type_instance) == 0:
            type_instance = ExposureType(exp_type=etype)
            type_instance.save()

    # Check to see if there is an extry for the instrument/proposal/observation
    existing = Archive.objects.get(instrument=instrument,
                                   proposal__prop_id=prop,
                                   proposal__observation__obsnum=obs)

    # If the entry exists, check the exp_type list of the entry and add the current
    # exptype if necessary.
    if len(existing) > 0:
        # exp_type is not allowed to be null, so no need to check for that condition here
        #existing_exps = [entry.instrument.proposal.observation.exposure_type.exp_type async for entry in existing]
        existing_exps = existing.instrument.proposal.observation.exposure_type.exp_type

        # Loop over the list of new types
        for etype in types:

            # If the new type is not in the list of exp_types already in the DB entry, then we need to add it
            if etype not in existing_exps:

                # Get the ExposureType instance for this new exp_type to be added
                type_instance = ExposureType.object.get(exp_type=etype)

                # Add the new exptype to the instrument/proposal/observation entry
                #existing.instrument.proposal.observation.exposure_type.exp_type.append(type_instance)
                existing_exps.append(type_instance)

        # Update the entry with the new list of exp_types
        existing.update(exp_type=existing_exps)

    # If the instrument/proposal/observation entry does not yet exist, we need to create it
    else:
        # Generate a list of the query sets for the exposure types to be added
        exp_instances = [ExposureType.object.get(exp_type=etype) for etype in types]

        # Create the instance for the observation
        obs_instance = Observation(obsnum=obs, number_of_files=files, exposure_type=exp_instances,
                                   obsstart=startdate, obsend=enddate)

        # Check to see if the instrument/proposal combination exists in the db. It may be that
        # there are entries for the proposal already, but for other observation numbers. If that's the case,
        # then we just need to update the observation list. There should only be one entry if there are any
        # at all, so we can use get()
        prop_instance = Archive.objects.get(instrument=instrument, proposal__prop_id=prop)

        # If the proposal entry exists, update the observation list
        if len(prop_instance) != 0:
            obs_list = prop_instance.observation
            obs_list.append(obs_instance)
            prop_instance.update(observation=obs_list)

        # If the proposal entry does not exist, we need to create it
        else:
            prop_instance = Archive(instrument=instrument, proposal__prop_id=prop,
                                    proposal__observation=obs_instance,
                                    proposal__observation__exposure_type=exp_instances)
            prop_instance.save()



if __name__ == '__main__':
    instruments = ['nircam', 'miri', 'nirspec', 'niriss', 'fgs']
    for instrument in instruments:
        get_updates(instrument)
