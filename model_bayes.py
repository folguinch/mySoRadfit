from collections import OrderedDict

import numpy as np
from myutils.time_utils import get_timestamp
from myutils.database import load_database

from .model_base import BaseModel

class modelBayes(BaseModel):
    """Manages the models for a bayesian fit

    Attributes:
        config (myConfigParser): model configuration file name.
        params (OrderedDict): model physical parameters.
        setup (myConfigParser): model rt setup.
        images (myConfigParser): model images setup.
        priors (OrderedDict): parameters for the prior calculation.
        db (sqlite3): database.
        logger: logging system.
    """

    def __init__(self, config, database='models_db.sqlite'):
        """Initialize model.

        A timestamp will be assigned as the model name suffix.

        Parameters:
            config (str): configuration file name.
        """
        # Initialize model
        name = get_timestamp()
        super(modelBayes, self).__init__(name=name, config=config)
        self.config = self.config['source']

        # Database
        self.db = None
        self.db = self._load_database(database)

        # Likelihood
        self.priors = self.load_priors()

        # Filter parameters
        #self.fixed = filter(self._filter_params, self.param_list)

    def load_params(self):
        super(modelBayes, self).load_params()

    @property
    def param_list(self):
        pars = []
        for yso in self.params.values():
            pars += yso.param_list
        return pars

    def _filter_params(self, param):
        return param not in self.prior and \
                param not in self.config[db_ignore_params]

    def load_priors(self):
        """Load the prior for the parameters of each source"""
        self.priors = OrderedDict.fromkeys(self.config.sections())
        for key in self.priors.keys():
            self.priors[key] = self._load_parser(
                    self.config[key]['prior'])

    def get_logprior(self, source, param):
        """Get the natural logarithm of the prior value for a parameter.

        Parameters:
            source (str): model source name
            param (str): parameter name
        """
        if self.priors[source]['type']=='uniform':
            section = self.priors[source][param]['section']
            val = self.params[source][section, param]
            vmin = self.priors[source].getquantity(param, 'min')
            vmax = self.priors[source].getquantity(param, 'max')
            if vmin < val < vmax:
                return 0.
            else:
                return np.inf
        else:
            raise NotImplementedError

    def get_logpriors(self, source=None):
        """Get the natural logarithm of the model parameter priors.

        If a source is not specified the priors for each source are given in an
        OrderedDict.

        Parameters:
            source (str, optional): source name
        Returns:
            priors (float): the natural logarithm of all the priors.
        """
        assert source in self.config

        prior = 0.

        for src, pr in self.priors.items():
            if source and src==source:
                prior = sum(self.get_logprior(src, key) \
                        for key in pr.sections())
                break
            elif source:
                continue
            else:
                prior += sum(self.get_logprior(src, key) \
                        for key in pr.sections())

        return prior

    # This only works for a model with one source
    def get_pa(self):
        return self.params['source']['Geometry','pa']

    # This only works for a model with one source
    def get_vlsr(self):
        return self.params['source']['Velocity','vlsr']

    def update_params(self, values):
        """Update the parameters of the model

        The order of the parameters in values follow model source as in the
        config file first and then parameter order as in the prior file.
        Therefore only parameters in the prior can be changed.

        If the model is in the database, the model name will be updated to that
        of the model in the database and return True.

        Parameters:
            values (list): new values of the parameters.

        Returns:
            validate (bool or str): if the model is in the database it returns
                the model ID else True.

        Notes:
            If the values do not have units, the units of the current parameter
            are assigned by the class managing the distribution.
        """
        # Get timestamp
        self.name = get_timestamp()

        # Change model
        for src in self.params.keys():
            for key in self.priors[src].sections():
                section = self.priors[src].get(key, 'section', fallback=None)
                if section is None or section.lower()=='all':
                    val = values.pop(0)
                    for section in self.params[key].sections():
                        self.params[key].update(section, key, val)
                else:
                    self.params[key].update(section, key, values.pop(0))

            # Update database
            if len(self.params.keys())==1:
                validate = self._update_database(table='models')
            else:
                validate = self._update_database(table=src, validate=False)

        if len(values)!=0:
            raise IndexError('There are unspecified parameters')

        return validate

    def _load_database(self, filename, table='models'):
        """Load models database

        Parameters:
            filename (str): database file name
            table (str, optional): table name
        """
        if len(self.params.keys())==1:
            key = self.params.keys()[0]
            keys, vals, fmts = self.params[key].flatten(get_db_format=True)
            keys = ['ID'] + keys
            vals = [self.name] + vals
            fmts = ['TEXT PRIMARY KEY']
            self.db = load_database(filename, table, keys, fmts)
            self._update_database(table=table, keys=keys, vals=vals,
                    validate=False)
        else:
            # Models with more than 1 source are exceptional and probably unique
            for src, yso in self.params.items():
                keys, vals, fmts = yso.flatten()
                if self.db is None:
                    self.db = load_database(filename, src, keys, fmts)
                else:
                    self.db.create_table(src, keys, fmts)

                self._update_database(table=src, keys=keys, vals=vals,
                        validate=False)

    def _update_database(self, table='models', keys=None, vals=None, 
            validate=True):
        """Update the database with the current model parameters
        
        Parameters:
            table (str, optional): table name
            keys (list, optional): table column names
            vals (list, optional): values to store
            validate (bool, optional): whether to validate the parameters
        """

        if vals is None or keys is None:
            keys, val = self.params[source or table].flatten()
            keys = ['ID'] + keys
            vals = [self.name] + vals

        if validate:
            valid = self._validate_db_params(table, keys, vals)
            if valid:
                self.db.update(table, vals)
            return valid
        else:
            self.db.update(table, vals)
            return True

    def _validate_db_params(self, table, keys, vals):
        """Query the database to find suplicates.

        Parameters:
            keys (list): column names. First comlumn is used as ID.
            vals (list): values to test
        """
        cond = ' AND '.join('{0}=:{0}'.format(key) for key in keys)
        self.db.query(table, 'ID', cond, dict(zip(keys, vals)))
        res = self.db.fetchone()

        if res is None:
            return True
        else:
            return res['ID']


