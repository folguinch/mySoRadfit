from abc import ABCMeta, abstractmethod

class Distribution(object):
    """Base distribution class.

    Attributes:
        __params (myConfigparser): model parameters
    """
    __metaclass__ = ABCMeta

    def __init__(self, param_file):
        """Initialize distribution.

        Parameters:
            param_file (str): distribution parameter file.
        """
        self.__params = self.load_config(param_file)

    def __getitem__(self, key):
        """Get a parameter as a quantity.
        
        Only keys with a section and parameter are allowed."""
        if len(key)!=2:
            raise KeyError('Key must be length 2: %r' % key)
        self.__validate_keys(*keys)
        try:
            # Load a quantity by default
            value = self.__params.getquantity(*key)

            # For debugging
            assert hasattr(value, 'unit')

        except ValueError:
            # If value is a string
            value = self.__params.get(*key)

        return value

    @abstractmethod
    def update(self, section, param, value):
        """Update the value of a parameter"""
        # Check units
        old = self[section, param]
        if not hasattr(value, 'unit') and hasattr(old, 'unit'):
            self.__params[section][param] = value*old.unit
        elif hasattr(value, 'unit') and hasattr(old, 'unit'):
            self.__params[section][param] = value.to(old.unit)
        elif hasattr(value, 'unit') and not hasattr(old, 'unit'):
            raise ValueError('The parameter *%s* in %s does not have unit' % \
                    (param, section))
        else:
            self.__params[section][param] = value

    @property
    def params(self):
        return self.__params

    @property
    def sections(self):
        return self.params.sections()

    #@property
    #def param_list(self):
    #    """Compile a list of all the parameters"""
    #    params = []
    #    for section in self.params.sections():
    #        params += params[section].options()
    #    return params

    @staticmethod
    def load_config(filename):
        """Load the parameter configuration file.

        Parameters:
            filename (str): YSO configuration file name.
        """
        # Verify file
        name = os.path.realpath(os.path.expanduser(filename))
        assert os.path.isfile(name)

        # Load file
        parser = myConfigParser(interpolation=ExtendedInterpolation())
        parser.read(name)

        return parser

    def __validate_keys(self, section, param=None):
        """Validate keys"""
        # Check section
        if section not in self.sections:
            raise KeyError('Section *%s* does not exist' % section)
        elif param is None:
            return True
        else:
            pass

        # Check param
        if param not in self.__params.options(section):
            raise KeyError('Param *%s* does not exist' % param)
        else:
            return True

    def flatten(self, ignore_params=[], get_db_format=False):
        """Return 2 arrays containg parameters and values

        Parameters names in the DEFAULT are renamed *<parameter>_<section>*.
        Parameters in list format are passed as a string

        Parameters:
            ignore_params (list, optional): list of prameters to ignore
            get_db_format (bool, optional): get an array with formats for
                databases
        """
        keys, vals, fmts = [], [], []
        for section in self.sections():
            for param in self.params[section]:
                # Ignores
                if param in ignore_params:
                    continue

                # Keys
                if param in self.params['DEFAULT']:
                    key += ['%s_%s' % (param, section.lower())]
                else:
                    keys += [param]

                # Values
                val = self[section, param]
                if hasattr(val, 'unit'):
                    try:
                        aux = len(val)
                        vals += ['%s' % (list(val.value),)]
                        fmts += ['TEXT']
                    except TypeError:
                        vals += [val.value]
                        fmts += ['REAL']
                else:
                        vals += [val]
                        fmts += ['TEXT']

        if get_db_format:
            return (keys, vals, fmts)
        else:
            return (keys, vals)


