import os, shutil

from myutils.logger import get_logger
from myutils.decorators import register_function, REGISTERED_FUNCTIONS

from ..utils.mollie import get_physical_props, write_setup

AVAILABLE_DUST_RT = set(['hyperion'])
AVAILABLE_LINE_RT = set(['mollie'])
THIS = os.path.dirname(os.path.abspath(__file__))

def rt_pipe(model, logger=get_logger(__name__)):
    """Configure and run the required RTs.

    Paramters:
        model (model): input model.
    """
    # Get RTs to run
    rts = model.get_rt_list()

    # Run RTs in the correct order
    # Dust first
    continuum = AVAILABLE_DUST_RT.intersection(rts)
    file_track = {}
    for cont in continuum:
        logger.info('Running continuum RT: %s', cont)
        pipe = rt.lower()+'_pipe'
        file_track.update(REGISTERED_FUNCTIONS[pipe](model))
    else:
        logger.info('No continuum RT requested')

    # Line
    line = AVAILABLE_LINE_RT.intersection(rts)
    for rt in line:
        logger.info('Running line RT: %s', rt)
        pipe = rt.lower()+'_pipe'
        file_track.update(REGISTERED_FUNCTIONS[pipe](model, logger=logger))
    else:
        logger.info('No line RT requested')

    return file_track

@register_function
def hyperion_pipe(model):
    raise NotImplementedError

@register_function
def mollie_pipe(model, from_file=False, logger=get_logger(__name__)):
    """Run the line RT.

    """
    MOLLIE = os.path.abspath(os.path.join(THIS, '../mollie'))
    rt = 'mollie'

    # Load grid 
    grids, cell_sizes, oversample = model.load_grids(rt)

    # Configure and write model
    fitslist = get_physical_props(model.params, grids, cell_sizes,
            MOLLIE, oversample=oversample, logger=logger)
    exit()

    # Other configurations
    nproc = model.setup.getint(rt, 'nproc')
    server = model.setup.get(rt, 'server', fallback=None)
    model_dir = os.path.expanduser(model.config.get('DEFAULT', 'model_dir',
        fallback='./'))
    model_dir = os.path.join(model_dir, 'Model_%s' % model.name)

    # Check the model directory tree exist or create it
    exit()

    # Run model
    mollie_file = os.path.join(MOLLIE, 'ModelCube0')
    setup_template = os.path.join(MOLLIE, 'src/setup_template.c')
    original_file = mollie_file
    file_track = {}
    for sect in model.images.sections():
        print('='*80)
        logger.info('Running model for: %s', sect)
        model_name = rt.upper()+'_'+sect.upper()

        # Setup RT for line
        write_setup(sect, model, setup_template)

        # line RT
        if from_file:
            mollie_file = os.path.join(mollie_dir, model_name)

            # Run if it does not exist
            if not os.path.isfile(mollie_file):
                run_mollie(mollie_dir, logger=logger, np=nproc, server=server,
                        from_email='folguin@phys.nthu.edu.tw',
                        to_email='f.olguin.ch@gmail.com')

                # Move to model directory
                logger.info('Moving model %s -> %s', os.path.basename(original_file),
                        mollie_file)
                shutil.move(original_file, mollie_file)
            file_track[section] = mollie_file
        else:
            run_mollie(MOLLIE, logger=logger, np=nproc, server=server, 
                    from_email='folguin@phys.nthu.edu.tw', 
                    to_email='f.olguin.ch@gmail.com')

            # Move to model directory
            logger.info('Moving model %s -> %s', os.path.basename(mollie_file),
                    model_name)
            shutil.move(mollie_file, os.path.join(model_dir, model_name))
            file_track[section] = os.path.join(model_dir, model_name)


    # Move FITS files
    for ffile in fitslist:
        shutil.move(ffile, model_dir)

    return file_track

