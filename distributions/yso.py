import os
from itertools import product
from configparser import ExtendedInterpolation

import numpy as np
import astropy.units as u
from astropy.io import fits
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

    def __call__(self, x, y, z, ignore_rim=False):
        """Density distribution.

        Parameters:
            x, y, z (floats): position where the density is evaluated.
        """
        disc, envelope, cavity = self.density(x, y, z, ignore_rim=ignore_rim)

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

    def from_fits(self, quantity, section='DEFAULT'):
        """Load a quantity from a FITS file.

        The FITS files have a standard name for each of the quantities.

        Parameters:
            quantity (str): quantity to load.
        """
        dirname = os.path.expanduser(self.params.get(section, 'grids_library'))
        if 'density' in quantity:
            fname = '{0}_rc{1:d}.fits'.format(quantity,
                    self.params.getquantity('Envelope','rc').to(u.au).value)
        fname = os.path.join(dirname, fname)
        assert os.path.isfile(fname)
        self.logger('Loading: %s', os.path.basename(fname))

        return fits.open(fname)[0]

    def to_fits(self, quantity, value, section='DEFAULT'):
        """Save a quantity grid into a fits file.

        The FITS files have a standard name for each of the quantities.

        Parameters:
            quantity (str): quantity to save.
            value (ndarray): values to save.
            section (str): section with the directory name.
        """
        assert value.ndim==3
        if 'density' in quantity:
            rc = self.params.getquantity('Envelope','rc').to(u.au).value
            nz, ny, nx = value.shape
            fname = '{0}_nx{1}_ny{2}_nz{3}_rc{4:d}.fits'.format(quantity, nx, ny,
                    nz, rc)
        dirname = os.path.expanduser(self.params.get(section, 'grids_library'))
        fname = os.path.join(dirname, fname)
        hdu = fits.PrimaryHDU(value)
        hdu.writeto(fname, overwrite=True)


    def density(self, x, y, z, ignore_rim=False, save_components=False,
            from_file=False):
        """Calculates and returns the density distribution by components.

        Parameters:
            x, y, z (floats): position where the density is evaluated.
        """
        # Load from file if needed
        if from_file:
            try:
                dims = 'nx%i_ny%i_nz%i' % x.shape[::-1]
                disc = self.from_fits('density_disc_%s' % dims)
                envelope = self.from_fits('density_envelope_%s' % dims)
                cavity = self.from_fits('density_cavity_%s' % dims)
                return disc, envelope, cavity
            except AssertionError:
                pass
        else:
            pass

        # Convert coordinates to spherical
        r, th, phi = cart_to_sph(x, y, z, pos0=self.loc)

        # Disk
        if 'Disc' in self.params:
            disc = discs.flared(r, th, self.params, ignore_rim=ignore_rim)
        else:
            disc = 0.

        # Envelope and cavity
        if 'Envelope' in self.params:
            envelope = ulrich.density(r, th, self.params, ignore_rim=ignore_rim)
            if 'Cavity' in self.params:
                cavity, mask = outflow.density(r, th, self.params,
                        ignore_rim=ignore_rim)
                cavity[cavity>envelope] = envelope[cavity>envelope]
                envelope[~mask] = 0.
        else:
            envelope = cavity = 0.

        # Overlap
        #if min_height_to_disc is not None:
        #    rdisc = self.params.getquantity('Disc', 'rmax')
        #    ind = (~mask) & (np.abs(z) <= min_height_to_disc.to(z.unit)) & \
        #            (r.cgs<=rdisc.cgs) 
        #    cavity[ind] = (cavity[ind]+disc[ind]+envelope[ind])/2.

        #    ind = mask & (np.abs(z) <= min_height_to_disc.to(z.unit)) & \
        #            (r.cgs<=rdisc.cgs)
        #    disc[ind] = disc[ind]/2
        #    envelope[ind] = envelope[ind]/2
        #envelope[~mask] = 0
        #disc[~mask] = 0
        if save_components:
            self.to_fits('density_disc', disc)
            self.to_fits('density_envelope', envelope)
            self.to_fits('density_cavity', cavity)

        return disc, envelope, cavity

    def velocity(self, x, y, z, ignore_rim=False, disc_dens=None, 
            env_dens=None, cav_dens=None):
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
                self.params, ignore_rim=ignore_rim)
        ind = r.cgs<=rdisc.cgs
        vr_env[ind] = 0.
        vth_env[ind] = 0.
        vphi_env[ind] = 0.
        
        # Combine
        if disc_dens is None and env_dens is None and cav_dens is None:
            disc, envelope, cavity = self.density(x, y, z, ignore_rim=ignore_rim)
        else:
            disc, envelope, cavity = disc_dens, env_dens, cav_dens
        envelope[ind] = 0.
        dens = envelope + cavity
        vr = vr_env + vr_out
        vth = vth_env + vth_out
        vphi = vphi_env + vphi_out
        vr = (vr*dens + vr_disc*disc) / (dens+disc)
        vth = (vth*dens + vth_disc*disc) / (dens+disc)
        vphi = (vphi*dens + vphi_disc*disc) / (dens+disc)
        vr[np.isnan(vr.value)] = 0.
        vth[np.isnan(vth.value)] = 0.
        vphi[np.isnan(vphi.value)] = 0.
        vr[np.isinf(vr.value)] = 0.
        vth[np.isinf(vth.value)] = 0.
        vphi[np.isinf(vphi.value)] = 0.

        #vr_disc[~mask] = 0.
        #vth_disc[~mask] = 0.
        #vphi_disc[~mask] = 0.

        ## Combine
        ## Overlaping section
        #if min_height_to_disc is not None:
        #    #ind2 = (r.cgs<=rdisc.cgs) & \
        #    #        (np.abs(z) < min_height_to_disc.to(z.unit))
        #    #ind = ind | ind2
        #    ind = (~mask) & (np.abs(z) <= min_height_to_disc.to(z.unit)) & \
        #            (r.cgs<=rdisc.cgs) 
        #    # Get the weights
        #    if disc_density is None:
        #        disc_density = discs.flared(r, th, self.params, ignore_rim=True)
        #    if cavity_density is None:
        #        cavity_density, mask2 = outflow.density(r, th, self.params,
        #                ignore_rim=True)

        #    # Weighted average where the cavity and outflow overlaps
        #    vr_out[ind] = (disc_density[ind].cgs*vr_disc[ind].cgs +
        #            cavity_density[ind].cgs*vr_out[ind].cgs) / \
        #                    (disc_density[ind].cgs+cavity_density[ind].cgs)
        #    vth_out[ind] = (disc_density[ind].cgs*vth_disc[ind].cgs +
        #            cavity_density[ind].cgs*vth_out[ind].cgs) / \
        #                    (disc_density[ind].cgs+cavity_density[ind].cgs)
        #    vphi_out[ind] = (disc_density[ind].cgs*vphi_disc[ind].cgs +
        #            cavity_density[ind].cgs*vphi_out[ind].cgs) / \
        #                    (disc_density[ind].cgs+cavity_density[ind].cgs)
        #vr = vr_env.cgs + vr_out.cgs
        #vth = vth_env.cgs + vth_out.cgs
        #vphi = vphi_env.cgs + vphi_out.cgs
        #ind = mask & (r.cgs<=rdisc.cgs)
        #vr[ind] = vr_disc.cgs[ind]
        #vth[ind] = vth_disc.cgs[ind]
        #vphi[ind] = vphi_disc.cgs[ind]

        return vel_sph_to_cart(vr, vth, vphi, th, phi)

    def get_density_velocity(self, x, y, z, temperature=None, ignore_rim=False,
            nquad=1):
        """Optimized function for obtaining the 3-D density and the velocity
        simultaneously.

        This avoid recalculating the density when the velocity function is
        called. The function can also be evaluated in different quadrants to
        avoid slowing the code with larger grids.

        Parameters:
            x, y, z (floats): position where the density is evaluated.
            ignore_rim (bool): if True the inner (dust) rim is ignored and the
                distributions are calculated to the surface of the star for the
                density and the inner velocity rim for the velocity.
            nquad (int, default=1): number of quadrants to devide the x, y, z 
                grids.
        """
        assert x.shape==y.shape==z.shape
        assert x.ndim==3

        # Initialize grids
        disc = np.zeros(x.shape)
        envelope = np.zeros(x.shape)
        cavity = np.zeros(x.shape)
        vx = np.zeros(x.shape)
        vy = np.zeros(x.shape)
        vz = np.zeros(x.shape)
        if temperature is not None:
            temp = np.zeros(x.shape) * u.K

        # Indices of divisions
        xinds = np.arange(-1, x.shape[2]+1, x.shape[2]/nquad)
        xinds[0] = 0
        xinds = zip(xinds[:-1], xinds[1:])
        yinds = np.arange(-1, x.shape[1]+1, x.shape[1]/nquad)
        yinds[0] = 0
        yinds = zip(yinds[:-1], yinds[1:])
        zinds = np.arange(-1, x.shape[0]+1, x.shape[0]/nquad)
        zinds[0] = 0
        zinds = zip(zinds[:-1], zinds[1:])

        for xind, yind, zind in product(xinds, yinds, zinds):
            # Sub axes
            subx = x[zind[0]:zind[1]+1,yind[0]:yind[1]+1,xind[0]:xind[1]+1] 
            suby = y[zind[0]:zind[1]+1,yind[0]:yind[1]+1,xind[0]:xind[1]+1]
            subz = z[zind[0]:zind[1]+1,yind[0]:yind[1]+1,xind[0]:xind[1]+1]

            # Evaluate
            dens = self.density(subx, suby, subz, ignore_rim=ignore_rim)
            vel = self.velocity(subx, suby, subz, ignore_rim=ignore_rim,
                    disc_dens=dens[0], env_dens=dens[1], cav_dens=dens[1])
            if temperature is not None:
                t = temperature(subx, suby, subz)
                temp[zind[0]:zind[1]+1,yind[0]:yind[1]+1,xind[0]:xind[1]+1] = \
                        t

            # Put units
            if not hasattr(disc, 'unit'):
                disc = disc * dens[0].unit
                envelope = envelope * dens[1].unit
                cavity = cavity * dens[2].unit
                vx = vx * vel[0].unit
                vy = vy * vel[1].unit
                vz = vz * vel[2].unit

            # replace
            disc[zind[0]:zind[1]+1,yind[0]:yind[1]+1,xind[0]:xind[1]+1] = \
                    dens[0]
            envelope[zind[0]:zind[1]+1,yind[0]:yind[1]+1,xind[0]:xind[1]+1] = \
                    dens[1]
            cavity[zind[0]:zind[1]+1,yind[0]:yind[1]+1,xind[0]:xind[1]+1] = \
                    dens[2]
            vx[zind[0]:zind[1]+1,yind[0]:yind[1]+1,xind[0]:xind[1]+1] = vel[0]
            vy[zind[0]:zind[1]+1,yind[0]:yind[1]+1,xind[0]:xind[1]+1] = vel[1]
            vz[zind[0]:zind[1]+1,yind[0]:yind[1]+1,xind[0]:xind[1]+1] = vel[2]

        if temperature is not None:
            return disc+envelope+cavity, vx, vy, vz, temp
        else:
            return disc+envelope+cavity, vx, vy, vz

    def abundance(self, temperature):
        """Molecular abundance.

        Parameters:
            temperature: the temperature distribution.
        """
        Ts = self.params.getquantity('Abundance', 't')
        abn = np.array(self.params.getfloatlist('Abundance','abn'))
        abundance = np.ones(temperature.shape) * np.min(abn)

        for i,T in enumerate(Ts):
            if i==0:
                abundance[temperature<T] = abn[i]
            elif i==len(Ts)-1:
                ind = (temperature<T) & (temperature>=Ts[i-1])
                abundance[ind] = abn[i]
                abundance[temperature>=T] = abn[i+1]
            else:
                ind = (temperature<T) & (temperature>=Ts[i-1])
                abundance[ind] = abn[i]

        return abundance

