import os
from configparser import ExtendedInterpolation

import numpy as np
import astropy.units as u
from myutils.myconfigparser import myConfigParser
from myutils.coordinates import cart_to_sph, vel_sph_to_cart

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
        r, th, phi = cart_to_sph(x, y, z, pos0=self.loc)
        th = th*u.rad
        #r = np.sqrt((x-self.loc[0])**2 + (y-self.loc[1])**2 + \
        #        (z-self.loc[2])**2)
        ##th = np.arccos((z-self.loc[2])/r)
        #th = np.arctan2(np.sqrt((x-self.loc[0])**2 + (y-self.loc[1])**2),
        #        z-self.loc[2])

        # Disk
        if 'Disc' in self.params:
            disc = discs.flared(r, th, self.params)
        else:
            disc = 0.

        # Envelope and cavity
        if 'Envelope' in self.params:
            envelope = ulrich.density(r, th, self.params)
            if 'Cavity' in self.params:
                cavity, mask = outflow.density(r, th, self.params)
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

    def velocity(self, x, y, z, min_height_to_disc=None):
        """Velocity distribution.

        Parameters:
            x, y, z (floats): position where the density is evaluated.
            min_height_to_disc (float): the velocity of points below this
                height are set to the disc velocity.
        """
        # Convert coordinates to spherical
        r, th, phi = cart_to_sph(x, y, z, pos0=self.loc)

        # Disc radius
        rdisc = self.params.getquantity('Disc', 'rmax')

        # Velocity in Envelope and cavity
        if self.params.get('Velocity', 'envelope', fallback='').lower() == 'ulrich':
            # Envelope
            vr_env, vth_env, vphi_env = ulrich.velocity(r, th, self.params)
        # Cavity
        vr_out, vth_out, vphi_out, mask = outflow.velocity(r, th, self.params)
        vr_env[~mask] = 0.
        vth_env[~mask] = 0.
        vphi_env[~mask] = 0.

        # Disc
        vr_disc, vth_disc, vphi_disc = discs.keplerian_rotation(r, th, 
                self.params)
        #vr_disc[~mask] = 0.
        #vth_disc[~mask] = 0.
        #vphi_disc[~mask] = 0.

        # Combine
        ind = mask & (r.cgs<=rdisc.cgs)
        if min_height_to_disc is not None:
            ind2 = (r.cgs<=rdisc.cgs) & \
                    (np.abs(z) < min_height_to_disc.to(z.unit))
            ind = ind | ind2
        vr = vr_env.cgs + vr_out.cgs
        vth = vth_env.cgs + vth_out.cgs
        vphi = vphi_env.cgs + vphi_out.cgs
        vr[ind] = vr_disc.cgs[ind]
        vth[ind] = vth_disc.cgs[ind]
        vphi[ind] = vphi_disc.cgs[ind]

        return vel_sph_to_cart(vr, vth, vphi, th, phi)


