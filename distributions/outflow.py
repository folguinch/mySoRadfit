import numpy as np
import astropy.units as u

def get_mask(r, th, params):
    """Get the points which are not in the cavity.

    Parameters:
        r: spherical coordinate radial distance.
        th: polar angle.
        params: model parameters.
    """
    theta0 = params.getquantity('Cavity','theta0')
    r0 = params.getquantity('Cavity','zref').cgs / \
            np.cos(theta0.to(u.radian))

    # Mask
    z0 = r0 * np.cos(theta0.to(u.radian))
    R0 = r0 * np.sin(theta0.to(u.radian))
    z = r.cgs * np.cos(th.to(u.radian))
    R = r.cgs * np.sin(th.to(u.radian))
    zcav = z0 * (R/R0)**params.getfloat('Cavity','power')
    mask = np.abs(z) < zcav

    return mask

def density(r, th, params, component='dust'):
    """Cavity density distribution.

    The distribution is assumed to be a power law, so constant density has a
    power exponent equal zero. Code taken from hyperion.

    Parameters:
        r: spherical coordinate radial distance.
        th: polar angle.
        params: model parameters.
        component (default=dust): component for the inner rim.

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
    comp_rmin = params.getfloat('Cavity','rmin_%s' % component)
    rsub = params.getquantity('Cavity','rsub').cgs
    rmax = params.getquantity('Cavity','rcav').cgs
    density[r.cgs>rmax] = 0.
    density[r.cgs<rsub*comp_rmin] = 0.

    # Check inside the star
    rstar = params.getquantity('Star','r').cgs
    density[r.cgs<rstar] = 0.

    # Mask
    mask = get_mask(r, th, params)
    density[mask] = 0.

    return density, mask

def velocity(r, th, params, component='dust'):
    """Outflow velocity distribution.

    Parameters:
        r: spherical coordinate radial distance.
        th: polar angle.
        params: model parameters.
        component (default=dust): component for the inner rim.

    Returns:
        vr, vth, vphi: the velocity for each spherical coordinate direction.
        mask: mask where the density is defined.
    """
    # Velocity
    if params.get('Velocity','outflow').lower()=='hubble':
        # Reference radius
        theta0 = params.getquantity('Cavity','theta0')
        r0 = params.getquantity('Cavity','zref').cgs / \
                np.cos(theta0.to(u.radian))

        # Hubble flow
        vr = params.getquantity('Velocity','v0_cav') * (r.cgs/r0.cgs)
        assert vr.cgs.unit == u.cm/u.s
        vth = np.zeros(vr.shape) * vr.unit
        vphi = np.zeros(vr.shape) * vr.unit

        # Mask
        rmax = params.getquantity('Cavity','rcav').cgs
        vr[r.cgs>rmax] = 0.
        mask = get_mask(r, th, params)
        vr[mask] = 0.

        # Inside rim
        comp_rmin = params.getfloat('Cavity','rmin_%s' % component)
        rsub = params.getquantity('Cavity','rsub').cgs
        vr[r.cgs<rsub*comp_rmin] = 0.

        # Inside Star
        rstar = params.getquantity('Star','r').cgs
        vr[r.cgs<rstar] = 0.

        # Velocity cap
        try:
            vcap = params.getquantity('Velocity','v_cap')
            vr[vr.to(vcap.unit)>vcap] = vcap.to(vr.unit)
        except KeyError:
            pass

    else:
        raise NotImplementedError

    return vr.cgs, vth.cgs, vphi.cgs, mask

