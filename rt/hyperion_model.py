import os

import numpy as np
import astropy.units as u
from hyperion.model import ModelOutput

from ..objects.container import Container
from ..objects.sed import SED
from ..objects.register import REGISTERED_CLASSES

class Hyperion(Container):
    """Creates a container for Hyperion output data.

    Attributes:
        name: name of the model.
        config: output configuration of the model.
        data: data in the model.
    """

    def load_data(self, key, file_name=None, source=None, incl=None, angle=None,
            dtype=None):
        """Load data for key.

        Parameters:
            key: data to load.
            filename: file to open.
            source: source object.
            incl: inclination index.
            angle: inclination angle (model config must be pre-loaded)
            dtype: data type.
        """

        if key=='model':
            assert os.path.isfile(file_name)
            self.data[key] = ModelOutput(file_name)
        elif file_name and dtype:
            assert os.path.isfile(file_name)
            self.data[key] = load_data_by_type(file_name, dtype.lower(), 
                    REGISTERED_CLASSES)
        elif key=='sed':
            assert incl is not None or (angle is not None and \
                    self.config is not None)
            wlg, F = self.data['model'].get_sed(group=0,
                    distance=source.distance.cgs.value, inclination=incl,
                    units='Jy')
            data = np.array(zip(wlg, F[0]), dtype=[('wlg', float),('F', float)])
            self.data[key] = SED(data=data, 
                    units={'wlg':1.*u.micron, 'F':1.*u.Jy})
        else:
            raise NotImplementedError

