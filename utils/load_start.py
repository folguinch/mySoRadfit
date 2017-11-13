import os
from configparser import ConfigParser

from myutils.logger import get_logger

def load_config(config='settings.cfg', config_dir=[]):
    """Load a configuration file.

    The configuration file is searched in the directories in ``config_dir``.
    The default directories for looking the configuration file are:
    *~/.config/mySoRadfit*, *mySoRadfit*, *mySoRadfit/config*.
    The directories in ``config_dir`` take precedence.

    Parameters:
        config (str, default=settings.cfg): file name of the configuration
            file.
        config_dir (list): list of directories where ``config`` is searched.

    Returns:
        parser (ConfigParser): the configuration parser.
    """
    # Directories
    radfit_dir = os.path.join(os.path.dirname(__file__),'..')
    radfit_dir = os.path.realpath(radfit_dir)
    LOOKUP = config_dir + ['~/.config/mySoRadfit', radfit_dir, 
            os.path.join(radfit_dir, 'config')]

    # Look for configuration file
    parser = ConfigParser()
    for path in LOOKUP:
        file_name = os.path.join(path, config)
        realpath = os.path.expanduser(path)
        if os.path.isdir(realpath) and os.path.isfile(file_name):
            parser.read(file_name)
            break
    else:
        raise IOError('Could not load file %s' % config)

    return parser

def load_logger(name):
    return get_logger(name)

def load_all(name, config='settings.cfg', config_dir=[]):
    """Load all the regular start files and logger.

    Parameters:
        name (str): name of the logger.
        config (str, default=settings.cfg): file name of the configuration.
        config_dir (list): diretories where ``config`` may be located.

    Returns:
        logger (logging.Logger): logging handler.
        parser (ConfigParser): the configuration parser.
    """
    return load_logger(name), load_config(config, config_dir)

if __name__=='__main__':
    # Test
    load_config()
