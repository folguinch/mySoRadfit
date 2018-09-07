import numpy as np
import astropy.units as u
import astropy.constants as ct
from hyperion.util.integrate import integrate_powerlaw

def flared(r, th, params, component='dust'):
    """Calculate the density distribution of a flared disc.

    Function from hyperion.

    Parameters:
        r: spherical coordinate radial distance.
        th: polar angle.
        params: model parameters.
        component (default=dust): component for the inner rim.

    Note:
        The inner rim is defined by the dust sublimation radius times a factor
        which can depend in the *component* parameter. If a different component
        is included you need to define ```rim_<component>``` in the model
        parameter configuration file.
    """

    # Reference density
    dust_rmin = params.getfloat('Disc','rmin_dust')
    rsub = params.getquantity('Disc','rsub').cgs
    #rmin = params.getfloat('Disc','rmin') *\
    #        params.getquantity('Disc','rsub').cgs
    rmin = dust_rmin * rsub
    rmax = params.getquantity('Disc','rdisc').cgs
    beta = params.getfloat('Disc','beta')
    p = beta - params.getfloat('Disc','alpha')
    r0 = params.getquantity('Disc','r0').cgs
    h0 = params.getquantity('Disc','h0').cgs
    # This integral is to set the dust mass of the disc, thus the inner rim of
    # the dust should be used.
    intg = integrate_powerlaw(rmin, rmax, 1.0 + p)
    intg = intg * r0**-p
    intg = (2.*np.pi)**1.5 * h0 * intg
    rho0 = params.getquantity('Disc','mdisc').cgs/intg

    # Convert coordinates
    R = r.cgs * np.sin(th.to(u.rad))
    z = r.cgs * np.cos(th.to(u.rad))

    # Scale height
    h = h0 * (R/r0)**beta

    # Density evaluated in all points
    density = rho0 * (r0/R)**(beta-p) * np.exp(-0.5*(z/h)**2)

    # Special cases
    if component!='dust':
        comp_rmin = params.getfloat('Disc','rmin_%s' % component)
        density[r.cgs<rsub.cgs*comp_rmin] = 0.
    else:
        density[r.cgs<rmin.cgs] = 0.

    # Inside the star
    rstar = params.getquantity('Star','r').cgs
    density[r.cgs<rstar] = 0.

    # Outside the disc
    density[r.cgs>rmax] = 0.

    return density

def keplerian_rotation(r, th, params, component='dust'):
    """Calculate the velocity distribution of a Keplerian disc.

    At the moment it only uses the stellar mass and not all the mass inside a
    radius *r*. Infall motion is not included at the moment.

    Parameters:
        r: spherical coordinate radial distance.
        th: polar angle.
        params: model parameters.
        component (default=dust): component for the inner rim.

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
        rmax = params.getquantity('Disc','rdisc').cgs
        vphi[r.cgs>rmax] = 0.
        assert vphi.cgs.unit == u.cm/u.s

        # Inner rim
        rsub = params.getquantity('Disc','rsub').cgs
        comp_rmin = params.getfloat('Disc','rmin_%s' % component)
        vphi[r.cgs<rsub*comp_rmin] = 0.

        # Inside the star
        rstar = params.getquantity('Star','r').cgs
        vphi[r.cgs<rstar] = 0.

        # Other components
        vr = np.zeros(vphi.shape) * vphi.unit
        vth = np.zeros(vphi.shape) * vphi.unit
    
    else:
        raise NotImplementedError


    return vr, vth, vphi
