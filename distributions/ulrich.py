import numpy as np
import astropy.units as u
import astropy.constants as ct
from hyperion.densities.ulrich_envelope import solve_mu0

def density(r, th, params, ignore_rim=True):
    """Calculate the density in the given position.

    The fixes and warnings are obtained from hyperion.

    Parameters:
        r: spherical coordinate radial distance.
        th: polar angle.
        params: models parameters.
    """
    # Reference density
    rho0 = params.getquantity('Envelope', 'mdot').cgs / \
            (4.*np.pi * np.sqrt(ct.G.cgs*params.getquantity('Star','m').cgs*\
            params.getquantity('Envelope','rc').cgs**3))

    # Angle variables
    mu = np.cos(th.to(u.rad).value)
    rc = params.getquantity('Envelope','rc').cgs
    mu0 = solve_mu0(r.cgs.value/rc.value, mu)

    # Density
    density = (r.cgs/rc)**-1.5 * (1.+mu/mu0)**-0.5 * \
            (mu/mu0 + 2*mu0**2 * rc/r.cgs)**-1
    density = rho0 * density

    # Close to the singularity
    mid1 = (np.abs(mu) < 1.e-10) & (r.cgs.value < rc.value)
    density[mid1] = rho0 / np.sqrt(r[mid1].cgs / rc) \
            / (1. - r[mid1].cgs / rc) / 2.
    mid2 = (np.abs(mu) < 1.e-10) & (r.cgs.value > rc.value)
    density[mid2] = rho0 / np.sqrt(2. * r[mid2].cgs / rc - 1) \
            / (r[mid2].cgs / rc - 1.)
    if np.any((np.abs(mu) < 1.e-10) & (r.cgs == rc)):
        raise OverflowError('Point close to singularity')

    # Outside the envelope
    if ignore_rim:
        #rmin = params.getfloat('Velocity','rmin') *\
        #        params.getquantity('Disc','rsub').cgs
        rmin = params.getquantity('Star','r').cgs
    else:
        rmin = params.getfloat('Envelope','rmin') *\
                params.getquantity('Envelope','rsub').cgs
    mask = (r.cgs > params.getquantity('Envelope','rout').cgs) | \
            (r.cgs < rmin)
    density[mask] = 0.

    return density

def velocity(r, th, params, ignore_rim=False):
    """Calculate the Ulrich velocity distribution in the given points.

    Parameters:
        r: spherical coordinate radial distance.
        th: polar angle.
        params: model parameters.

    Returns:
        vr, vth, vphi: the components of the velocity along each direction.
    """
    # Keplerian
    v0 = np.sqrt(ct.G.cgs * params.getquantity('Star','m').cgs / r.cgs)
    assert v0.cgs.unit == u.cm/u.s

    # Angle variables
    mu = np.cos(th.to(u.rad).value)
    rc = params.getquantity('Envelope','rc').cgs
    mu0 = solve_mu0(r.cgs.value/rc.value, mu)
    theta0 = np.arccos(mu0) * u.rad

    # Velocity components
    rot_dir = params.getfloat('Velocity','rot_dir')
    assert np.abs(rot_dir)==1.
    vr = -1. * np.sqrt(1. + mu/mu0)
    vth = (mu0-mu)/np.sin(th) * np.sqrt(1. + mu/mu0)
    vphi = rot_dir * np.sin(theta0)/np.sin(th) * np.sqrt(1. - mu/mu0)

    # Outside the envelope
    if ignore_rim:
        rmin = params.getquantity('Star','r').cgs
    else:
        rmin = params.getfloat('Velocity','rmin') *\
                params.getquantity('Envelope','rsub').cgs
    mask = (r.cgs > params.getquantity('Envelope','rout').cgs) | \
            (r.cgs < rmin)
    vr[mask] = 0.
    vth[mask] = 0.
    vphi[mask] = 0.

    return v0*vr, v0*vth, v0*vphi
