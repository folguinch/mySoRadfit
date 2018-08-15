from myutils.time import get_timestamp

from .model_base import BaseModel

class modelBayes(BaseModel):
    """Manages the models for a bayesian fit

    Attributes:
        config (myConfigParser): model configuration file name.
        params (OrderedDict): model physical parameters.
        setup (myConfigParser): model rt setup.
        images (myConfigParser): model images setup.
        db (sqlite3): database.
        logger: logging system.
    """

    def __init__(self, config, database='models_db.sqlite'):
        """Initialize model.

        A timestamp will be assigned as the model name.

        Parameters:
            config (str): configuration file name.
        """
        # Initialize model
        super(modelBayes, self).__init__(name=get_timestamp(), config=config)

        # Database
        self.db = load_database(database, self.param_list)
