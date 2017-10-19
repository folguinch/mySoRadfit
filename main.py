import os, argparse
from configparser import ConfigParser

from myutils.logger import get_logger

from objects.source import LoadSource
from utils.parsers.args import LoadConfig

def main():
    """Main function for the RT fitter"""
    # Important standard files and directories
    DIRECTORIES = 'directories.cfg'
    LOOKUP = ['~/.config/rt_fitter', './', './config']

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

if __name__=='__main__':
    main()
