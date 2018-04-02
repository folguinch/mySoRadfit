import numpy as np
import astropy.units as u
from myutils.math import map_sph_to_cart_axisym

def get_temp_func(params, temp, r, th, phi=None, extrapolate=True, **kwargs):
    """Get a temperature function from a grid.

    The returned function gives the temperature at the new grid points. Only
    spherically symmetric grids are allowed at the moment.

    Parameters:
        params (configparser): model parameters.
        temp (np.array): original temperature array.
        r, th (np.array): coordinates of the original cells.
        phi (np.array, optional): azimuthal spherical coordinates 
            (not implemented)
        extrapolate (bool, default=True): extrapolate inside the dust
            sublimation radius. If False, a temperature of 2.7 K is used inside
            this radius.
        kwargs: parameters for the scipy.map_coordiantes function.
    """
    def temp_func(x, y, z):
        # Convert to spherical
        R = np.sqrt(x**2 + y**2 + z**2)
        TH = np.arctan2(np.sqrt(x**2+y**2), z)

        # Evaluate function
        if phi==None:
            temp = map_sph_to_cart_axisym(temperature[0,:,:], r, th, x, y, z) 
        else:
            raise NotImplementedError
        temp = temp * u.K
        
        # Fill outside grid
        rstar = params.getquantity('Star','r')
        rmin = params.getfloat('Disc','rmin') *\
               params.getquantity('Disc','rsub')
        rout = params.getquantity('Envelope','rout')
        tstar = params.getquantity('Star','t')

        # Inside star
        ind = R.cgs<=rstar.cgs
        temp[ind] = tstar

        # Inside sublimation radius
        ind = (R.cgs>rstar.cgs) & (R.cgs<rmin.cgs)
        if extrapolate:
            temp[ind] = tstar*(rstar.cgs/R.cgs[ind])**(1./2.1)
        else:
            temp[ind] = 2.7*u.K

        # Outside the envelope
        ind = R.cgs > rout.cgs
        temp[ind] = 2.7*u.K

        # Checks
        temp[temp<2.7*u.K] = 2.7*u.K

        return temp

    return temp_func
