import os, argparse
from configparser import ConfigParser

from objects.source import LoadSource

def main():
    """Main function for the RT fitter"""
    # Look for configuration file
    LOOKUP = ['~/.config/rt_fitter', './', './config']
    for path in LOOKUP:
        file_name = os.path.join(path,'directories.cfg')
        realpath = os.path.expanduser(path)
        if os.path.isdir(realpath) and os.path.isfile(file_name):
            config = ConfigParser()
            config.read(file_name)
            break

    # Read command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('source', action=LoadSource, config=config)
    parser.add_argument('--stats', action='store_true', 
            help='Compute stats only')
    parser.add_argument('--use', nargs='*', type=list, default=['all'], 
            help='Use these data')
    par = parser.parse_args()

    # Workflow

if __name__=='__main__':
    main()
