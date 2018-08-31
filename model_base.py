import os
from collections import OrderedDict
from abc import ABCMeta, abstractmethod
from configparser import ExtendedInterpolation

from myutils.logger import get_logger
from myutils.myconfigparser import myConfigParser
from myutils.array_utils import load_struct_array


class BaseModel(object):
    """Model class.

    Attributes:
        config (myConfigParser): model configuration file name.
        params (OrderedDict): model physical parameters.
        setup (myConfigParser): model rt setup.
        images (myConfigParser): model images setup.
        logger: logging system.
    """

    __metaclass__ = ABCMeta
    logger = get_logger(__name__)

    def __init__(self, name=None, config=None, params=None, setup=None,
            source_names=['source'], locs=[(0,0,0)]):
        """Initialize the model.

        A model can be created with just a name. A generic name is provided,
        thus if none of the parameters are given no exception is raised.

        If a model *name* and a configuration file are given, the name in the
        configuration file will be replaced by *name*.

        The configuration file is standard for each model and contains all the
        information to load/recover a model. Therefore, it includes the params
        and setup file names. If all of the latter and model name are given the
        parameter *config* is ignored.

        Parameters:
            name (str): model name.
            config (str): model configuration file.
            params (str or list): model parameter file or files for each model 
                source.
            setup (str): model setup file.
            source_names (list): names of the model sources.
            locs (list): location of each model source.
        """
        # Initialize 
        self.config = myConfigParser(interpolation=ExtendedInterpolation())
        self.params = OrderedDict()
        self.setup = None
        self.images = None

        # Load configuration
        if config and not params and not setup:
            self.load_config(config, name)
        elif params:
            assert len(source_names)==len(locs)
            self.config['DEFAULT'] = {'name': name or 'model', 'setup': setup}
            for i,(name,loc) in enumerate(zip(source_name, locs)):
                assert len(loc)==3
                self.config[name] = {'loc': loc}
                pname = os.path.realpath(os.path.expanduser(params[i]))
                if os.path.isfile(pname):
                    self.config[name]['params'] = pname
                elif os.path.isfile(params):
                    self.config[name]['params'] = params
            self.load_params()
            if setup:
                self.load_setup(setup)
            if images:
                raise NotImplementedError

    @abstractmethod
    def fill_grids(self):
        return True

    @abstractmethod
    def load_params(self):
        """Load parameter file for each model source."""
        for section in self.config.sections():
            self.logger.info('Loading parameters for: %s', section)
            self.params[section] = YSO(self.config[section]['params'],
                    loc=self.config.getfloatlist(section, 'loc'))

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
        parser = self._load_configparser(parser, filename)

        return parser

    def load_config(self, config, name=None):
        """Load configuration file.

        Parameters:
            config: configuration file name.
            name: model name.
        """
        # Load file and update name
        self.logger.info('Loading model configuration file')
        self.config = self._load_configparser(self.config, config)
        if name:
            self.name = name
        self.logger.info('Model name: %s', self.name)

        # Load setup and params
        if self.config.get('DEFAULT','setup'):
            self.load_setup(self.config.get('DEFAULT','setup'))
        if self.config.get('DEFAULT','images'):
            self.load_image_setup(self.config.get('DEFAULT','images'))
        self.load_params()

    def load_setup(self, filename):
        """Load setup file.

        Parameters:
            filename: name of the setup file.
        """
        self.logger.info('Loading model setup')
        self.setup = self._load_parser(filename)

    def load_image_setup(self, filename):
        """Load the image setup file.

        Parameters:
            filename: name of the image setup file.
        """
        self.logger.info('Loading image setup')
        self.images = self._load_parser(filename)

    def load_grids(self, rt):
        """Load the grids from file.
        
        Parameters:
            rt: RT transfer config section.
        """
        grids = []
        for grid in self.setup.getlist(rt, 'grids'):
            self.logger.info('Loading grid: %s', os.path.basename(grid))
            fname = os.path.realpath(os.path.expanduser(grid))
            grid = load_struct_array(fname, usecols=None)
            grids += [grid]

        cell_sizes = self.setup.getintlist(rt, 'cell_sizes')

        # Sort grid sizes
        self.logger.info('Sorting grids by cell size (a->Z)')
        ind = np.argsort(cell_sizes)
        grids = map(lambda x: grids[x], ind)
        cell_sizes = map(lambda x: cell_sizes[x], ind)

        return grids, cell_sizes

