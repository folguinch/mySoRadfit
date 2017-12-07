import numpy as np
import astropy.units as u

def density(r, th, params):
    """Cavity density distribution.

    The distribution is assumed to be a power law, so constant density has a
    power exponent equal zero. Code taken from hyperion.

    Parameters:
        r: spherical coordinate radial distance.
        th: polar angle.
        params: model parameters.

    Returns:
        density: the density distribution.
        mask: mask where the density is defined.
    """
    
    # Density
    theta0 = params.getquantity('Cavity','theta0')
    r0 = params.getquantity('Cavity','zref').cgs / \
            np.cos(theta0.to(u.radian))
    density = params.getquantity('Cavity', 'rho0').cgs
    density = density * (r.cgs/r0)**(-1.*params.getfloat('Cavity','rho_exp'))

    # Limit values
    rmin = params.getfloat('Cavity','rmin') *\
            params.getquantity('Cavity','rsub').cgs
    rmax = params.getquantity('Cavity','rout').cgs
    density[r.cgs>rmax] = 0.
    density[r.cgs<rmin] = 0.

    # Mask
    z0 = r0 * np.cos(theta0.to(u.radian))
    R0 = r0 * np.sin(theta0.to(u.radian))
    z = r.cgs * np.cos(theta0.to(u.radian))
    R = r.cgs * np.sin(theta0.to(u.radian))
    zcav = z0 * (R/R0)**params.getfloat('Cavity','power')
    mask = np.abs(z) < zcav
    density[mask] = 0.

    return density, mask
