import os
#from abc import ABCMeta, abstractmethod
from configparser import ExtendedInterpolation

from myutils.logger import get_logger
from myutils.myconfigparser import myConfigParser

class Model(object):
    """Model class.

    Attributes:
        config (myConfigParser): model configuration file name.
        params (myConfigParser): model parameters file name.
        setup (myConfigParser): model rt and output setup file name.
        logger: logging system.
    """

    logger = get_logger(__name__)

    def __init__(self, name='model', config=None, params=None, setup=None):
        """Initialize the model.

        A model can be created with just a name. A generic name is provided,
        thus if none of the parameters are given no exception is rised.

        If a model *name* and a configuration file are given, the name in the
        configuration file will be replaced by *name*.

        The configuration file is standard for each model and contains all the
        information to load/recover a model. Therefore, it includes the params
        and setup file names. If all of the latter and model name are given the
        parameter *config* is ignored.

        Parameters:
            name: model name.
            config: model configuration file.
            params: model parameter file.
            setup: model setup file.
        """
        # Initialize 
        self.config = myConfigParser(interpolation=ExtendedInterpolation())
        self.params = None
        self.setup = None

        # Load configuration
        if config and not params and not setup:
            self.config = load_config(config, name)
        else:
            self.config['DEFAULT'] = {'name': name, 'params': params,
                    'setup': setup}
            if params:
                self.load_params(params)
            if setup:
                self.load_setup(setup)

    @property
    def name(self):
        return self.config['DEFAULT']['name']

    @name.setter
    def name(self, value):
        self.config['DEFAULT']['name'] = value

    def _load_configparser(self, parser, filename):
        # Verify file
        name = os.path.realpath(os.path.expanduser(filename))
        assert os.path.isfile(name)

        # Load file
        parser.read(name)

        return parser

    def _load_parser(self, filename):
        parser = myConfigParser(interpolation=ExtendedInterpolation())
        parser = _load_configparser(parser, filename)

        return parser

    def load_config(self, config, name=None):
        """Load configuration file.

        Parameters:
            config: configuration file name.
            name: model name.
        """
        if name is None and config.get('DEFAULT', 'name'):
            name = config['DEFAULT']['name']

        self.config = self._load_configparser(self.config, config)
        self.name = name

        # Load setup and params
        if self.config.get('DEFAULT','setup'):
            self.load_setup(self.config.get('DEFAULT','setup'))
        if self.config.get('DEFAULT','params'):
            self.load_params(self.config.get('DEFAULT','params'))

    def load_setup(self, filename):
        """Load setup file.

        Parameters:
            filename: name of the setup file.
        """
        self.setup = self._load_parser(filename)

    def load_params(self, filename):
        """Load parameters file.

        Parameters:
            filename: name of the parameter file.
        """
        self.params = self._load_parser(filename)
