import numpy as np
from hyperion.util.integrate import integrate_powerlaw

def flared(r, th, params):
    """Calculate the density distribution of a flared disc.

    Function from hyperion.

    Parameters:
        r: spherical coordinate radial distance.
        th: polar angle.
        params: model parameters.
    """

    # Reference density
    rmin = params.getfloat('Disc','rmin') *\
            params.getquantity('Disc','rsub').cgs
    rmax = params.getquantity('Disc','rmax').cgs
    beta = params.getfloat('Disc','beta')
    p = beta - params.getfloat('Disc','alpha')
    r0 = params.getquantity('Disc','r0').cgs
    h0 = params.getquantity('Disc','h0').cgs
    intg = integrate_powerlaw(rmin, rmax, 1.0 + p)
    intg = int1 * r0**-p
    intg = (2.*np.pi)**1.5 * h0 * intg
    rho0 = params.getquantity('Disc','m').cgs/intg

    # Convert coordinates
    R = r.cgs * np.sin(th.to(u.rad))
    z = r.cgs * np.cos(th.to(u.rad))

    # Scale height
    h = h0 * (R/r0)**beta

    # Density
    density = rho0 * (r0/R)**(beta-p) * np.exp(-0.5*(z/h)**2)
    density[r.cgs<rmin] = 0.
    density[r.cgs>rmax] = 0.

    return density
