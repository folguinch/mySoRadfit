import numpy as np
import astropy.units as u
import astropy.constants as ct
from hyperion.util.integrate import integrate_powerlaw

def flared(r, th, params, ignore_rim=False):
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
    intg = intg * r0**-p
    intg = (2.*np.pi)**1.5 * h0 * intg
    rho0 = params.getquantity('Disc','m').cgs/intg

    # Convert coordinates
    R = r.cgs * np.sin(th.to(u.rad))
    z = r.cgs * np.cos(th.to(u.rad))

    # Scale height
    h = h0 * (R/r0)**beta

    # Density
    density = rho0 * (r0/R)**(beta-p) * np.exp(-0.5*(z/h)**2)
    if not ignore_rim:
        density[r.cgs<rmin] = 0.
    else:
    #    rmin = params.getfloat('Velocity','rmin') *\
    #            params.getquantity('Disc','rsub').cgs
        rmin = params.getquantity('Star','r').cgs
        density[r.cgs<rmin] = 0.
    density[r.cgs>rmax] = 0.

    return density

def keplerian_rotation(r, th, params, ignore_rim=False):
    """Calculate the velocity distribution of a Keplerian disc.

    At the moment it only uses the stellar mass and not all the mass inside a
    radius *r*. Infall motion is not included at the moment.

    Parameters:
        r: spherical coordinate radial distance.
        th: polar angle.
        params: model parameters.

    Returns:
        vr, vth, vphi: the velocity for each spherical coordinate direction.
    """
    # Velocity components
    if params.getboolean('Velocity', 'rotation_only'):
        rot_dir = params.getfloat('Velocity','rot_dir')
        r_cyl = r.cgs * np.sin(th)
        vphi = np.sqrt(ct.G.cgs * params.getquantity('Star','m').cgs / r_cyl)
        vphi = rot_dir * vphi

        # Outside the disc
        rmax = params.getquantity('Disc','rmax').cgs
        vphi[r.cgs>rmax] = 0.
        assert vphi.cgs.unit == u.cm/u.s

        # Inner rim
        if ignore_rim:
            rmin = params.getquantity('Star','r').cgs
            vphi[r.cgs<rmin] = 0.
        else:
            rmin = params.getfloat('Velocity','rmin') *\
                    params.getquantity('Disc','rsub').cgs
            vphi[r.cgs<rmin] = 0.

        # Other components
        vr = np.zeros(vphi.shape) * vphi.unit
        vth = np.zeros(vphi.shape) * vphi.unit
    
    else:
        raise NotImplementedError


    return vr, vth, vphi
