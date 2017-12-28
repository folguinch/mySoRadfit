import os
from configparser import ExtendedInterpolation

import numpy as np
from myutils.myconfigparser import myConfigParser

import ulrich, outflow, discs

def hyperion_yso_density(x, y, z, params, loc=(0.,0.,0.)):
    """Claculates the density distribution of a hyperion analytical YSO.

    Hyperion replaces the cavity values only for the envelope.

    Parameters:
        x: x cartesian coordinate.
        y: y cartesian coordinate.
        z: z cartesian coordinate.
        params: model parameters.
        loc: (x,y,z) location of the source
    """
    # Convert coordinates to spherical
    #xy = np.sqrt((x-loc[0])**2 + (y-loc[1])**2)
    r = np.sqrt((x-loc[0])**2 + (y-loc[1])**2 + (z-loc[2])**2)
    th = np.arccos((z-loc[2])/r)

    # Disk
    disc = discs.flared(r, th, params)

    # Envelope and cavity
    envelope = ulrich.density(r, th, params)
    cavity, mask = outflow(r, th, params)
    cavity[cavity>envelope] = envelope[cavity>envelope]
    envelope[~mask] = 0.

    return disc + envelope + cavity
    
class YSO(object):
    """Create a YSO object for the model.

    The YSO object manages the pgysical properties of the yso. When the YSO
    pbject is called with the coordinates of a grid point, the total density is
    returned.

    Attributes:
        params (myConfigParser): physical properties of the YSO.
        loc (iterable): location of the source in the grid.
    """

    def __init__(self, params, loc=(0,0,0)):
        """Initialize the YSO.

        Parameters:
            params (str): YSO configuration file.
            loc (iterable): location of the source.
        """
        assert len(loc)==3
        self.params = self.load_config(params)
        self.loc = loc

    def __call__(self, x, y, z):
        """Density distribution.

        Parameters:
            x, y, z (floats): position where the density is evaluated.
        """
        # Convert coordinates to spherical
        r = np.sqrt((x-self.loc[0])**2 + (y-self.loc[1])**2 + \
                (z-self.loc[2])**2)
        th = np.arccos((z-self.loc[2])/r)

        # Disk
        if 'Disc' in self.parmas:
            disc = discs.flared(r, th, params)
        else:
            disc = 0.

        # Envelope and cavity
        if 'Envelope' in self.params:
            envelope = ulrich.density(r, th, params)
            if 'Cavity' in self.params:
                cavity, mask = outflow(r, th, params)
                cavity[cavity>envelope] = envelope[cavity>envelope]
                envelope[~mask] = 0.
        else:
            envelope = cavity = 0.

        return disc + envelope + cavity

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


