import numpy as np
import astropy.units as u
import astropy.constants as ct
from hyperion.densities.UlrichEnvelope import solve_mu0

def density(r, th, params):
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
            params.getquantity('rc').cgs**3))

    # Angle variables
    mu = np.cos(th.to(u.rad).value)
    rc = params.getquantity('Envelope','rc').cgs
    mu0 = solve_mu0(r.cgs.value/rc.value, mu)

    # Density
    density = (r.cgs/rc)**-1.5 * (1.+mu/mu0)**-0.5 * \
            (mu/m0 + 2*mu0**2 * rc/r.cgs)**-1
    density = rho0 * density

    # Close to the singularity
    mid1 = (np.abs(mu) < 1.e-10) & (r.cgs.value < rc.value)
    density[mid1] = rho0 / np.sqrt(r[mid1].cgs / rc) \
            / (1. - r[mid1].cgs / rc) / 2.
    mid2 = (np.abs(mu) < 1.e-10) & (r > self.rc)
    density[mid2] = rho0 / np.sqrt(2. * r[mid2].cgs / rc - 1) \
            / (r[mid2].cgs / rc - 1.)
    if np.any((np.abs(mu) < 1.e-10) & (r.cgs == rc)):
        raise OverflowError('Point close to singularity')

    # Outside the envelope
    rmin = params.getfloat('Envelope','rmin') *\
            params.getquantity('Envelope','rsub').cgs
    mask = (r.cgs > params.getquantity('Envelope','rout').cgs) | \
            (r.cgs < rmin)
    density[mask] = 0.

    return density

def velocity(r, th, phi, params):
    pass
