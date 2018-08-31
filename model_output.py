from astroSource.container import Container
from myutils.decorators import REGISTERED_FUNCTIONS

class modelOutput(Container):
    """Model output class

    Attributes:
        name (str): model name
        config (myConfigParser): output configuration file
        data (OrderedDict): output data
        logger: logging manager
    """
    logger = get_logger(__name__)

    def __init__(self, name, config):
        """Defines a new model output.

        Parameters:
            name (str): model name
            config (str): configuration file name
        """
        super(modelOutput, self).__init__(name, config)

    def load_data(self, data_name, model_file, PA=0., vlsr=0.):
        assert data_name in self.config.sections()

        pipe = self.config[data_name]['pipe']
        fn = REGISTERED_FUNCTIONS[pipe.lower()]
        self.data[data_name] = fn(self.config[data_name], model_file, PA=PA,
                vlsr=vlsr)

    def load_all(self, model_files, PA=0., vlsr=0.):
        for key in self.config.sections():
            if key not in model_files:
                if 'group' in self.config[key] and \
                        self.config[key]['group'] in model_files:
                    realkey = self.config[key]['group']
                else:
                    raise KeyError('No models with key: %s' % key)
            else:
                realkey = key

            self.load_data(key, model_files[realkey], PA=PA, vlsr=vlsr)
