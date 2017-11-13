import os, argparse
from configparser import ConfigParser

from myutils.logger import get_logger

from MySoRadfit.objects.source import LoadSource
from MySoRadfit.rt.hyperion_model import Hyperion
#from utils.parsers.args import LoadConfig

def main():
    """Main function for the RT fitter"""
    # Important standard files and directories
    DIRECTORIES = 'directories.cfg'
    LOOKUP = ['~/.config/MySoRadfit', './MySoRadfit', './MySoRadfit/config']

    # Logger
    logger = get_logger(__name__)

    # Look for configuration file
    logger.debug('Loading main configuration file')
    for path in LOOKUP:
        file_name = os.path.join(path, DIRECTORIES)
        realpath = os.path.expanduser(path)
        if os.path.isdir(realpath) and os.path.isfile(file_name):
            config = ConfigParser()
            config.read(file_name)
            break

    # Read command line arguments
    logger.debug('Loading line arguments')
    parser = argparse.ArgumentParser()
    #config = parser.add_argument('--config', action=LoadConfig, default=config)
    parser.add_argument('--stats', action='store_true', 
            help='Compute stats only')
    parser.add_argument('--use', nargs='*', type=list, default=['all'], 
            help='Use these data')
    parser.add_argument('source', action=LoadSource, 
            source_dir=config.get('DIRECTORIES','sources'),
            help='Source name')
    args = parser.parse_args()

    # Workflow
    model = Hyperion('Model.0009')
    model.load_data('model',
            os.path.expanduser('~/Sources/AFGL2591/Models/hyperion/Model.0009.rtout'))
    model.load_data('sed', source=args.source, incl=0)
    print model['sed'].interpolate(args.source['sed']['wlg'])
    
if __name__=='__main__':
    main()
