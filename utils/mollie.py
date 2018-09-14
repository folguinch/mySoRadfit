import os, time
import struct
from string import Template

import numpy as np
from astropy.io import fits
import astropy.units as u
import astropy.constants as ct
from hyperion.model import ModelOutput
from scipy.interpolate import griddata, interp2d
from myutils.logger import get_logger
from myutils.math import rebin_irregular_nd, map_sph_to_cart_axisym, rebin_regular_nd
from myutils.classes.data_3d import Data3D
from myutils.decorators import timed
from Plotter.mesh_plotter import NMeshPlotter

from ..distributions.temperature import get_temp_func

def fill_inner_radius(r, th, temperature, density, yso):
    rstar = yso.params.getquantity('Star','r')
    rmin = yso.params.getfloat('Disc','rmin') *\
            yso.params.getquantity('Disc','rsub')
    tstar = yso.params.getquantity('Star','t').value
    # Inside star
    ind = r.cgs<=rstar.cgs
    temperature[ind] = tstar
    # Inside sublimation radius
    ind = (r.cgs>rstar) & (r.cgs<rmin.cgs)
    temperature[ind] = tstar*((rstar.cgs/r.cgs[ind])**(1./2.1)).value
    # Fill density
    x = r*np.sin(th)
    z = r*np.cos(th)
    density[ind] = yso(x[ind], np.zeros(z[ind].shape), z[ind], ignore_rim=True)

    return temperature, density

def write_fits(dirname, **kwargs):
    fitsnames = []
    for key, val in kwargs.items():
        hdu = fits.PrimaryHDU(val)
        fitsnames += [os.path.join(dirname, key+'.fits')]
        hdu.writeto(fitsnames[-1], overwrite=True)

    return fitsnames

def get_walls(*args):
    walls = []
    for val in args:
        dv = np.min(val[val>0.*u.cm])
        vw = val - dv
        vw = np.append(vw.value, vw[-1].value+2.*dv.value)*vw.unit
        walls += [vw]

    return walls

def rebin_2Dsph_to_cart(val, new_pos, yso, pos=None, rebin=False, interp=False,
        min_height_to_disc=None, logger=get_logger(__name__), **kwargs):
    """Function for rebin or interpolate an input grid into a new grid.

    Parameters:
        val: 2-D spherically symmetric grid of values for 1st and 4th quadrant.
        new_pos: new grid centres.
        yso (YSO): YSO object containg the physical parameters.
        pos: original positions in the grid in spherical coords.
        rebin: rebin original grid
        interp: interpolate original grid
        kwargs: keyword arguments for the rebin or interpolate functions.
    """
    # Checks
    kwargs.setdefault('statistic', 'average')

    # Rebin if needed
    if rebin and pos is not None:
        logger.info('Rebinning grid')
        # Current position mesh
        R, TH = np.meshgrid(pos[0], pos[1])
        YI = R * np.sin(TH)
        ZI = R * np.cos(TH)

        # Extend new positions until maximum radial distance
        maxr = np.nanmax(np.sqrt(new_pos[1].value**2+new_pos[0].value**2))
        delta = np.abs(new_pos[1][0].value-new_pos[1][1].value)
        nymax = int(np.ceil(maxr/delta))
        extra = np.nanmax(new_pos[1].value) + \
                np.arange(1,abs(nymax-int(np.sum(new_pos[1]>0)))+1)*delta
        new_pos0 = np.append(-1.*extra, new_pos[1].value)
        new_pos0 = np.append(new_pos0, extra)
        new_pos0 = np.sort(new_pos0) * new_pos[1].unit

        # Meshes
        #YN, ZN = np.meshgrid(new_pos[1], new_pos[2])
        YN, ZN = np.meshgrid(new_pos0, new_pos[2])

        # Evaluate the function
        if interp:
            #val1 = rebin_2Dsph_to_cart(val, (new_pos[0], np.array([0.])*new_pos[1].unit,
            #    new_pos[2]), yso, pos=pos, interp=True, logger=logger)
            val1 = rebin_2Dsph_to_cart(val, (new_pos0, np.array([0.])*new_pos[1].unit,
                new_pos[2]), yso, pos=pos, interp=True, logger=logger)
            val1 = val1[:,0,:]
        else:
            val1 = yso(YN, np.zeros(YN.shape), ZN)
            assert val1.cgs.unit == u.g/u.cm**3
            val1 = val1.cgs.value

        # Walls
        #walls = get_walls(new_pos[1], new_pos[-1])
        walls = get_walls(new_pos0, new_pos[-1])
        walls = [wall.to(pos[0].unit).value for wall in walls]

        # Rebin 2-D
        val2 = rebin_irregular_nd(val[0,:,:],
                walls[::-1], ZI, YI, statistic=kwargs['statistic'],
                weights=kwargs.get('weights'))

        # Replace nans
        ind = np.isnan(val2)
        val2[ind] = val1[ind]

        # Interpolate
        ZN3, YN3, XN = np.meshgrid(*new_pos[::-1], indexing='ij')
        NEWX = np.sqrt(XN**2+YN3**2)
        val2 = griddata((YN.value.flatten(), ZN.value.flatten()),
                val2.flatten(), (NEWX.to(ZN.unit).value.flatten(), 
                    ZN3.to(ZN.unit).value.flatten()),
                method=kwargs.get('method', 'nearest'))

        return val2.reshape(XN.shape)
    elif interp and pos is not None:
        logger.info('Interpolating grid')
        val1 = map_sph_to_cart_axisym(val[0,:,:], pos[0], pos[1], *new_pos)
        return val1
    else:
        logger.info('Evaluating grid')
        # Meshes
        ZN, YN, XN = np.meshgrid(*new_pos[::-1], indexing='ij')

        # Evaluate the function
        val1 = yso(XN, YN, ZN, ignore_rim=True,
                min_height_to_disc=min_height_to_disc)
        assert val1.cgs.unit == u.g/u.cm**3
        val1 = val1.cgs.value

        return val1

def get_quadrant_ind(x, y, z):
    ind = [(x<=0) & (y<=0) & (z<=0)]
    ind += [(x<=0) & (y<=0) & (z>0)]
    ind += [(x<=0) & (y>0) & (z<=0)]
    ind += [(x<=0) & (y>0) & (z>0)]
    ind += [(x>0) & (y<=0) & (z<=0)]
    ind += [(x>0) & (y<=0) & (z>0)]
    ind += [(x>0) & (y>0) & (z<=0)]
    ind += [(x>0) & (y>0) & (z>0)]

    return ind

def get_bins_quadrant(bins):
    """
    Bins assumed sorted and odd length (even number of cells).
    """
    xmid = (len(bins[0]) - 1)/2
    ymid = (len(bins[1]) - 1)/2
    zmid = (len(bins[2]) - 1)/2

    newbins = [[bins[0][:xmid+1], bins[1][:ymid+1], bins[2][:zmid+1]]]
    newbins += [[bins[0][:xmid+1], bins[1][:ymid+1], bins[2][zmid:]]]
    newbins += [[bins[0][:xmid+1], bins[1][ymid:], bins[2][:zmid+1]]]
    newbins += [[bins[0][:xmid+1], bins[1][ymid:], bins[2][zmid:]]]
    newbins += [[bins[0][xmid:], bins[1][:ymid+1], bins[2][:zmid+1]]]
    newbins += [[bins[0][xmid:], bins[1][:ymid+1], bins[2][zmid:]]]
    newbins += [[bins[0][xmid:], bins[1][ymid:], bins[2][:zmid+1]]]
    newbins += [[bins[0][xmid:], bins[1][ymid:], bins[2][zmid:]]]

    return newbins

def eval_by_quadrant(x, y, z, func, nret=1, **kwargs):
    def replace(val, aux, ind):
        for i in range(nret):
            try:
                val[i][ind] = aux[i]
            except IndexError:
                val[i][ind] = aux
        return val
    val = [np.zeros(x.shape) for i in range(nret)]

    inds = get_quadrant_ind(x, y, z)
    for ind in inds:
        aux = func(x[ind], y[ind], z[ind], **kwargs)
        val = replace(val, aux, ind)
    
    for i in range(nret):
        try:
            val[i] = val[i] * aux[i].unit
        except IndexError:
            val[i] = val[i] * aux.unit

    return val

def rebin_by_quadrant(val, x, y, z, bins):
    inds1 = get_quadrant_ind(x, y, z)
    bins_quadrants = get_bins_quadrant(bins)
    ZN, YN, XN = np.meshgrid(z, y, x, indexing='ij')
    inds3 = get_quadrant_ind(XN, YN, ZN)
    newx = (bins[0][1:] + bins[0][:-1])/2.
    newy = (bins[1][1:] + bins[1][:-1])/2.
    newz = (bins[2][1:] + bins[2][:-1])/2.
    ZN, YN, XN = np.meshgrid(newz, newy, newx, indexing='ij')
    newval = np.zeros(XN.shape)
    inds2 = get_quadrant_ind(XN, YN, ZN)
    for ind1, ind2, ind3, qbins in zip(inds1, inds2, inds3, bins_quadrants):
        print newval[ind2].shape
        newvali[ind2] = rebin_regular_nd(val[ind3], z[ind1], y[ind1], x[ind1], 
                bins=qbins[::-1], statistic='sum')
        print newvali.shape
        newval[ind1] = newvali
    return newval

def rebin_by_chunks(val, x, y, z, bins, chunks=1):
    return 0

@timed
def phys_oversampled_cart(x, y, z, yso, temp_func, oversample=5, 
        logger=get_logger(__name__)):
    """Calculate the density and velocity distributions by oversampling the
    grid.

    This function first creates a new oversampled grid and then rebin this grid
    to the input one by taking weighted averages of the physical quantities.

    Parameters:
        yso (YSO object): the object containing the model parameters.
        xlim (tuple): lower and upper limits of the x-axis.
        ylim (tuple): lower and upper limits of the y-axis.
        zlim (tuple): lower and upper limits of the z-axis.
        cellsize (float): physical size of the cell in the coarse grid.
        oversample (int, default=5): oversampling factor.
        logger (logging): logger manager.
    """
    # Hydrogen mass
    mH = ct.m_p + ct.m_e

    # Special case
    if oversample==1:
        logger.info('Evaluating the grid')
        ZN, YN, XN = np.meshgrid(z, y, x, indexing='ij')
        n, vx, vy, vz, temp = yso.get_all(XN, YN, ZN,
                temperature=temp_func, component='gas')
        n = n / (2.33 * mH)
        temp[temp<2.7*u.K] = 2.7*u.K
        assert n.cgs.unit == 1/u.cm**3
        assert temp.unit == u.K
        assert vx.cgs.unit == u.cm/u.s
        assert vy.cgs.unit == u.cm/u.s
        assert vz.cgs.unit == u.cm/u.s
        return n, (vx, vy, vz), temp

    logger.info('Resampling grid')
    # Create new grid
    dx = np.abs(x[0]-x[1])
    dy = np.abs(y[0]-y[1])
    dz = np.abs(z[0]-z[1])
    xw,xstep = np.linspace(np.min(x)-dx/2., np.max(x)+dx/2., num=len(x)*oversample+1,
            endpoint=True, retstep=True)
    xover = (xw[:-1] + xw[1:]) / 2.
    logger.info('Resampled grid x-step = %s', xstep.to(u.au))
    yw,ystep = np.linspace(np.min(y)-dy/2., np.max(y)+dy/2., num=len(y)*oversample+1,
            endpoint=True, retstep=True)
    yover = (yw[:-1] + yw[1:]) / 2.
    logger.info('Resampled grid y-step = %s', ystep.to(u.au))
    zw,zstep = np.linspace(np.min(z)-dz/2., np.max(z)+dz/2., num=len(z)*oversample+1,
            endpoint=True, retstep=True)
    zover = (zw[:-1] + zw[1:]) / 2.
    logger.info('Resampled grid z-step = %s', zstep.to(u.au))
    ZN, YN, XN = np.meshgrid(zover, yover, xover, indexing='ij')

    # Volumes
    vol = dx*dy*dz
    vol_over = xstep*ystep*zstep

    # Number density
    bins = get_walls(z, y, x)
    n_over, vx_over, vy_over, vz_over, temp = yso.get_all(XN, YN, ZN,
            temperature=temp_func, component='gas', nquad=oversample)
    n_over = n_over / (2.33 * mH)
    assert n_over.cgs.unit == 1/u.cm**3
    N = vol_over * rebin_regular_nd(n_over.cgs.value, zover, yover, xover, bins=bins,
            statistic='sum') * n_over.cgs.unit
    dens = N / vol
    assert dens.cgs.unit == 1/u.cm**3

    # Temperature
    temp = rebin_regular_nd(temp.value*n_over.cgs.value, zover, yover, 
            xover, bins=bins, statistic='sum') * \
                    temp.unit * n_over.cgs.unit
    temp = vol_over * temp / N
    temp[temp<2.7*u.K] = 2.7*u.K
    assert temp.unit == u.K

    # Velocity
    assert vx_over.cgs.unit == u.cm/u.s
    assert vy_over.cgs.unit == u.cm/u.s
    assert vz_over.cgs.unit == u.cm/u.s
    v = []
    for vi in (vx_over, vy_over, vz_over):
        vsum = rebin_regular_nd(vi.cgs.value*n_over.cgs.value, zover, yover, 
                xover, bins=bins, statistic='sum') * \
                        vi.cgs.unit * n_over.cgs.unit
        vsum = vol_over * vsum / N
        vsum[np.isnan(vsum)] = 0.
        assert vsum.cgs.unit == u.cm/u.s
        v += [vsum]

    return dens, v, temp

@timed
def get_physical_props_single(yso, grids, cell_sizes, save_dir, oversample=[3],
        dust_out=None, logger=get_logger(__name__)):
    """Calculate and write the physical properties a model with one source.

    Parameters:
        yso: the model parameters.
        grids: grids where the model will be evaluated
        template: filename of the *define_model.c* file.
        logger: logging system.
    """
    # Models with more than one source should be treated in another function
    # because the oversampling should be different.

    # FITS list
    fitslist = []

    # Validate oversample
    if len(oversample)==1 and len(oversample)!=len(cell_sizes):
        oversample = oversample * len(cell_sizes)
    elif len(oversample)==len(cell_sizes):
        pass
    else:
        raise ValueError('The length of oversample != number of grids')

    # Temperature function
    if yso.params.get('DEFAULT', 'quantities_from'):
        hmodel = yso.params.get('DEFAULT', 'quantities_from')
        logger.info('Loading Hyperion model: %s', os.path.basename(hmodel))
        hmodel = ModelOutput(os.path.expanduser(hmodel))
        q = hmodel.get_quantities()
        temperature = np.sum(q['temperature'].array[1:], axis=0)
        r, th = q.r*u.cm, q.t*u.rad
        temp_func = get_temp_func(yso.params, temperature, r, th)
    elif dust_out is not None:
        logger.info('Loading Hyperion model: %s', os.path.basename(dust_out))
        hmodel = ModelOutput(os.path.expanduser(dust_out))
        q = hmodel.get_quantities()
        temperature = np.sum(q['temperature'].array[1:], axis=0)
        r, th = q.r*u.cm, q.t*u.rad
        temp_func = get_temp_func(yso.params, temperature, r, th)
    else:
        raise NotImplementedError

    # Start from smaller to larger grid
    inv_i = len(grids) - 1
    for i,(grid,cellsz) in enumerate(zip(grids, cell_sizes)):
        # Initialize grid axes
        print '='*80
        logger.info('Working on grid: %i', i)
        logger.info('Oversampling factor: %i', oversample[i])
        logger.info('Grid cell size: %i', cellsz)
        # Multiply by units
        x = grid[0]['x'] * grid[1]['x']
        y = grid[0]['y'] * grid[1]['y']
        z = grid[0]['z'] * grid[1]['z']
        xi = np.unique(x)
        yi = np.unique(y)
        zi = np.unique(z)

        # Density and velocity
        dens, (vx, vy, vz), temp = phys_oversampled_cart(
                xi, yi, zi, yso, temp_func, oversample=oversample[i], 
                logger=logger)

        # Replace out of range values
        dens[dens.cgs<=0./u.cm**3] = 10./u.cm**3
        temp[np.logical_or(np.isnan(temp.value), temp.value<2.7)] = 2.7 * u.K

        # Replace the inner region by rebbining the previous grid
        if i>0:
            # Walls of central cells
            j = cell_sizes.index(cellsz)
            xlen = cell_sizes[j-1] * len(xprev) * u.au
            nxmid = int(xlen.value) / cellsz
            xw = np.linspace(-0.5*xlen.value, 0.5*xlen.value, nxmid+1) * u.au
            ylen = cell_sizes[j-1] * len(yprev) * u.au
            nymid = int(ylen.value) / cellsz
            yw = np.linspace(-0.5*ylen.value, 0.5*ylen.value, nymid+1) * u.au
            zlen = cell_sizes[j-1] * len(zprev) * u.au
            nzmid = int(zlen.value) / cellsz
            zw = np.linspace(-0.5*zlen.value, 0.5*zlen.value, nzmid+1) * u.au
            if nxmid==nymid==nzmid==0:
                logger.warning('The inner grid is smaller than current grid size')
            else:
                logger.info('The inner %ix%ix%i cells will be replaced', nxmid,
                        nymid, nzmid)

                # Rebin previous grid
                # Density
                vol_prev = (cell_sizes[j-1]*u.au)**3
                vol = (cellsz * u.au)**3
                N_cen = vol_prev.cgs * rebin_regular_nd(dens_prev.cgs.value, 
                        zprev.cgs.value, yprev.cgs.value, xprev.cgs.value, 
                        bins=(zw.cgs.value,yw.cgs.value,xw.cgs.value), 
                        statistic='sum') * dens_prev.cgs.unit
                dens_cen = N_cen / vol
                dens_cen = dens_cen.to(dens.unit)
                # Temperature
                T_cen = rebin_regular_nd(temp_prev.value * dens_prev.cgs.value,
                        zprev.cgs.value, yprev.cgs.value, xprev.cgs.value, 
                        bins=(zw.cgs.value,yw.cgs.value, xw.cgs.value), 
                        statistic='sum') * temp_prev.unit * dens_prev.cgs.unit
                T_cen = vol_prev.cgs * T_cen / N_cen.cgs
                T_cen = T_cen.to(temp.unit)

                # Replace
                dens[len(zi)/2-nzmid/2:len(zi)/2+nzmid/2,
                        len(yi)/2-nymid/2:len(yi)/2+nymid/2,
                        len(xi)/2-nxmid/2:len(xi)/2+nxmid/2] = dens_cen
                temp[len(zi)/2-nzmid/2:len(zi)/2+nzmid/2,
                        len(yi)/2-nymid/2:len(yi)/2+nymid/2,
                        len(xi)/2-nxmid/2:len(xi)/2+nxmid/2] = T_cen
        dens_prev = dens
        temp_prev = temp
        xprev = xi
        yprev = yi
        zprev = zi

        # Linewidth and abundance
        linewidth = yso.linewidth(x, y, z, temp).to(u.cm/u.s)

        # Abundance per molecule
        abns = {}
        abn_fmt = 'abn_%s_%i'
        j = 1
        for section in yso.params.sections():
            if not section.lower().startswith('abundance'):
                continue
            mol = yso[section, 'molecule']
            abns[abn_fmt % (mol, inv_i)] = yso.abundance(x, y, z, temp, 
                    index=j, ignore_min=False)
            j = j+1

        # Write FITS
        kw = {'temp%i'%inv_i: temp.value, 'dens%i'%inv_i: dens.cgs.value, 
                'vx%i'%inv_i: vx.cgs.value, 'vy%i'%inv_i: vy.cgs.value, 
                'vz%i'%inv_i: vz.cgs.value, 
                'lwidth%i'%inv_i: linewidth.cgs.value}
        kw.update(abns)
        fitsnames = write_fits(os.path.expanduser(save_dir), 
                **kw)
        fitslist += fitsnames
        inv_i = inv_i - 1

    return fitslist

# Keep for backwards compatibility with script
def set_physical_props(yso, grids, cell_sizes, save_dir, oversample=3,
        dust_out=None, logger=get_logger(__name__)):
    """Calculate and write the physical properties of the model.

    Parameters:
        yso: the model parameters.
        grids: grids where the model will be evaluated
        template: filename of the *define_model.c* file.
        logger: logging system.
    """
    # Load temperature function
    if yso.params.get('DEFAULT', 'quantities_from'):
        hmodel = yso.params.get('DEFAULT', 'quantities_from')
        logger.info('Loading Hyperion model: %s', os.path.basename(hmodel))
        hmodel = ModelOutput(os.path.expanduser(hmodel))
        q = hmodel.get_quantities()
        temperature = np.sum(q['temperature'].array[1:], axis=0)
        r, th = q.r*u.cm, q.t*u.rad
        temp_func = get_temp_func(yso.params, temperature, r, th)
    elif dust_out is not None:
        logger.info('Loading Hyperion model: %s', os.path.basename(dust_out))
        hmodel = ModelOutput(os.path.expanduser(dust_out))
        q = hmodel.get_quantities()
        temperature = np.sum(q['temperature'].array[1:], axis=0)
        r, th = q.r*u.cm, q.t*u.rad
        temp_func = get_temp_func(yso.params, temperature, r, th)
    else:
        raise NotImplementedError

    # Open template
    fitslist = []
    # Start from smaller to larger grid
    for i,grid,cellsz in zip(range(len(grids))[::-1], grids, cell_sizes):
        print '='*80
        logger.info('Working on grid: %i', i)
        logger.info('Grid cell size: %i', cellsz)
        x = grid[0]['x'] * grid[1]['x']
        y = grid[0]['y'] * grid[1]['y']
        z = grid[0]['z'] * grid[1]['z']
        xi = np.unique(x)
        yi = np.unique(y)
        zi = np.unique(z)

        # Density and velocity
        dens, (vx, vy, vz), temp = phys_oversampled_cart(xi, yi, zi, yso,
                temp_func, oversample=oversample if i!=2 else 5, logger=logger)
        dens[dens.cgs<=0./u.cm**3] = 10./u.cm**3
        temp[np.isnan(temp.value)] = 2.7 * u.K

        # Replace the inner region by rebbining the previous grid
        if i<len(grids)-1:
            # Walls of central cells
            j = cell_sizes.index(cellsz)
            xlen = cell_sizes[j-1] * len(xprev) * u.au
            nxmid = int(xlen.value) / cellsz
            xw = np.linspace(-0.5*xlen.value, 0.5*xlen.value, nxmid+1) * u.au
            ylen = cell_sizes[j-1] * len(yprev) * u.au
            nymid = int(ylen.value) / cellsz
            yw = np.linspace(-0.5*ylen.value, 0.5*ylen.value, nymid+1) * u.au
            zlen = cell_sizes[j-1] * len(zprev) * u.au
            nzmid = int(zlen.value) / cellsz
            zw = np.linspace(-0.5*zlen.value, 0.5*zlen.value, nzmid+1) * u.au
            if nxmid==nymid==nzmid==0:
                logger.warning('The inner grid is smaller than current grid size')
            else:
                logger.info('The inner %ix%ix%i cells will be replaced', nxmid,
                        nymid, nzmid)

                # Rebin previous grid
                # Density
                vol_prev = (cell_sizes[j-1]*u.au)**3
                vol = (cellsz * u.au)**3
                N_cen = vol_prev.cgs * rebin_regular_nd(dens_prev.cgs.value, 
                        zprev.cgs.value, yprev.cgs.value, xprev.cgs.value, 
                        bins=(zw.cgs.value,yw.cgs.value,xw.cgs.value), 
                        statistic='sum') * dens_prev.cgs.unit
                dens_cen = N_cen / vol
                dens_cen = dens_cen.to(dens.unit)
                # Temperature
                T_cen = rebin_regular_nd(temp_prev.value * dens_prev.cgs.value,
                        zprev.cgs.value, yprev.cgs.value, xprev.cgs.value, 
                        bins=(zw.cgs.value,yw.cgs.value, xw.cgs.value), 
                        statistic='sum') * temp_prev.unit * dens_prev.cgs.unit
                T_cen = vol_prev.cgs * T_cen / N_cen.cgs
                T_cen = T_cen.to(temp.unit)

                # Replace
                dens[len(zi)/2-nzmid/2:len(zi)/2+nzmid/2,
                        len(yi)/2-nymid/2:len(yi)/2+nymid/2,
                        len(xi)/2-nxmid/2:len(xi)/2+nxmid/2] = dens_cen
                temp[len(zi)/2-nzmid/2:len(zi)/2+nzmid/2,
                        len(yi)/2-nymid/2:len(yi)/2+nymid/2,
                        len(xi)/2-nxmid/2:len(xi)/2+nxmid/2] = T_cen
        dens_prev = dens
        temp_prev = temp
        xprev = xi
        yprev = yi
        zprev = zi

        # Abundance
        abundance = yso.abundance(temp)
        
        # Linewidth
        amu = 1.660531e-24 * u.g
        atoms = yso.params.getfloat('Velocity', 'atoms')
        c_s2 = ct.k_B * temp / (atoms * amu)
        linewidth = np.sqrt(yso.params.getquantity('Velocity', 'linewidth')**2
                + c_s2)

        # Write FITS
        fitsnames = write_fits(os.path.expanduser(save_dir), 
                **{'temp%i'%i: temp.value, 'dens%i'%i: dens.cgs.value, 
                'vx%i'%i: vx.cgs.value, 'vy%i'%i: vy.cgs.value, 'vz%i'%i:
                vz.cgs.value, 'abn%i'%i: abundance,
                'lwidth%i'%i: linewidth.cgs.value})
        fitslist += fitsnames

    return fitslist

def write_setup(section, model, template, rt='mollie', incl=0*u.deg,
        phi=0*u.deg):
    # Open template
    with open(os.path.expanduser(template)) as ftemp:
        temp = Template(ftemp.read())

    # Data to write:
    fmt = '%.8f'
    radius = model.setup.getquantity(rt, 'radius').to(u.pc)
    kwd = {'chwidth': model.images.getquantity(section, 'chwidth').cgs.value, 
            'nchan': model.images.get(section, 'nchan'),
            'minlwidth': model.images.getquantity(section,
                'minlwidth').cgs.value,
            'maxlwidth': model.images.getquantity(section, 
                'maxlwidth').cgs.value,
            'radius': fmt % radius.value, 
            'nphi': model.setup.get(rt, 'nphi'), 
            'nth': model.setup.get(rt, 'nth'), 
            'nlines': model.images.get(section, 'nlines'),
            'line_id': model.images.get(section, 'line_id'),
            'th': incl.to(u.deg).value,
            'phi': phi.to(u.deg).value,
            'npix': model.images.get(section, 'npix'),
            'pixsize': model.images.getquantity(section,
                'pixsize').to(u.pc).value,
            'n_state': model.images.get(section,'n_state'),
            'n_lines': model.images.get(section,'n_lines')}

    # Write file
    dirname = os.path.dirname(os.path.expanduser(template))
    fname = 'setup.c'
    while True:
        try:
            with open(os.path.join(dirname, fname), 'w') as out:
                out.write(temp.substitute(**kwd))
            break
        except IOError:
            print 'Re-trying saving setup.c'
            time.sleep(2)
            continue

def load_model(model, source, filename, logger, old=False, older=False, 
               write='combined_jy', velocity=True, pa=None):
    """Load a Mollie model.

    Written by: K. G. Johnston.
    Modified by: F. Olguin

    Parameters:
        model: model file name.
        source (astroSource): source information.
        logger: logger.
        write: type of image to write (default: data cube in Jy).
        velocity: output cube 3rd axis in velocity or frequency.
        pa: source position angle.
    """

    maxname = 20
    if older: maxname = 16

    logger.info('Opening file: %s', model)
    f = open(model, "rb")
    endian = '>'
    byte = f.read(4)
    nlines = struct.unpack(endian+'l',byte)[0]

    swap_bytes = False
    if nlines > 200:
        swap_bytes = True
        endian = '<'
    logger.info('Swapping bytes? %s', swap_bytes)

    nlines = struct.unpack(endian+'l',byte)[0]
    logger.info('There are %i lines in this data set', nlines)

    nchan = np.zeros(nlines,dtype=int)
    restfreq = np.zeros(nlines)

    nviews = struct.unpack(endian+'l',f.read(4))[0]
    logger.info('There are %i views in this data set', nviews)

    for line in range(nlines):
        nch = struct.unpack(endian+'l',f.read(4))[0]
        nchan[line] = nch
    logger.info('The numbers of channels are %i', nchan)

    nx = struct.unpack(endian+'l',f.read(4))[0]
    ny = struct.unpack(endian+'l',f.read(4))[0]
    cellx = struct.unpack(endian+'f',f.read(4))[0] * u.pc
    celly = struct.unpack(endian+'f',f.read(4))[0] * u.pc
    beamx = struct.unpack(endian+'f',f.read(4))[0] * u.pc
    beamy = struct.unpack(endian+'f',f.read(4))[0] * u.pc

    logger.info('Grid size nx=%i, ny=%i', nx, ny)
    logger.info('Cell size %sx%s', cellx, celly)

    linename = (nlines)*['']

    for i in range(nlines):
        for j in range(maxname):
            bytvar = struct.unpack(endian+'c',f.read(1))[0]
            linename[i] += bytvar.decode('ascii')
        linename[i] = linename[i].strip()

    if (not old) and (not older):
        for i in range(nlines):
            restfreq[i] = struct.unpack(endian+'d',f.read(8))[0]

    logger.info('The lines are:')
    restfreq = restfreq * u.Hz
    for i in range(nlines): 
        logger.info('%s at %s', linename[i], restfreq[i].to(u.GHz))

    maxch = max(nchan)
    chvel = np.zeros((nlines,maxch))

    for i in range(nlines):
        for n in range(nchan[i]):
            chvel[i,n] = struct.unpack(endian+'f',f.read(4))[0]

    chvel = chvel * u.cm/u.s
    for i in range(nlines):
        logger.info('Velocities for line %s:  %s, %s', linename[i], 
                    chvel[i,0].to(u.km/u.s), chvel[i,nchan[i]-1].to(u.km/u.s))

    lng = np.zeros(nviews)
    lat = np.zeros(nviews)

    for i in range(nviews):
        lng[i] = struct.unpack(endian+'f',f.read(4))[0]
        lat[i] = struct.unpack(endian+'f',f.read(4))[0]

    logger.info('Longitudes: %r', lng)

    # fix inclination convention to Hyperion
    lat = 90 - np.array(lat)
    logger.info('Latitudes: %r', lat)

    xc = np.zeros(nx)
    for i in range(nx):
        xc[i] = struct.unpack(endian+'f',f.read(4))[0]

    yc = np.zeros(ny)
    for i in range(ny):
        yc[i] = struct.unpack(endian+'f',f.read(4))[0]

    data = np.ones((nlines,nviews,nx,ny,maxch)) * np.nan
    for l in range(nlines):
        data_bytes = f.read(4 * nviews * nx * ny * nchan[l])
        data[l,:,:,:,:nchan[l]] = np.fromstring(data_bytes, dtype=endian + \
                                                'f4').reshape(nviews, nx, 
                                                              ny, nchan[l])
        logger.info('Max in line %i is %.3e', l, np.nanmax(data[l]))

    f.close()

    logger.info('Min and max brightness in data set: %.3e, %.3e',
                np.nanmin(data), np.nanmax(data))

    if linename[0].strip().lower().startswith("mc"):
        logger.warn('Remember to run CH3CN_CUBE after LOAD_MODEL\n' + \
                    'in order to stack the CH3CN lines into a single spectrum')

    # Set up header common to all files

    rightascension, declination = source.position.ra, source.position.dec
    distance = source.distance.to(u.pc) 

    header_template = fits.Header()

    header_template['OBJECT'] = source.name

    header_template['TELESCOP'] = 'MOLLIE'
    header_template['INSTRUME'] = 'MOLLIE'
    header_template['OBSERVER'] = 'MOLLIE'

    header_template['CTYPE1'] = 'RA---SIN'
    header_template['CTYPE2'] = 'DEC--SIN'

    header_template['CUNIT1'] = 'degree'
    header_template['CUNIT2'] = 'degree'

    header_template['CRPIX1'] = nx / 2.
    header_template['CRPIX2'] = ny / 2. + 1.

    header_template['CDELT1'] = -1.*np.abs(np.degrees((cellx.si/distance.si).value))
    header_template['CDELT2'] = np.degrees((celly.si / distance.si).value)

    header_template['CRVAL1'] = rightascension.to(u.deg).value
    header_template['CRVAL2'] = declination.to(u.deg).value
    header_template['EPOCH'] = 2000

    if pa or source.get_quantity('pa') is not None:
        header_template['CROTA1'] = 0
        if pa is None:
            pa = source.get_quantity('pa')
        header_template['CROTA2'] = (360*u.deg - pa.to(u.deg)).value

    #header_template['BMAJ'] = np.degrees((beamx.si / distance.si).value * 2.35)
    #header_template['BMIN'] = np.degrees((beamx.si / distance.si).value * 2.35)

    header_template['EQUINOX'] = 2000.

    for v in range(nviews):
        for l in range(nlines):

            ### Set up header ###
            header = header_template.copy()
            header['LINE'] = linename[l]

            cdelt3 = (chvel[l,1] - chvel[l,0]).to(u.km/u.s)
            if velocity:
                header['CTYPE3'] = 'VELO-LSR'
                header['CUNIT3'] = 'KM/S'
                crval3 = chvel[l,0].to(u.km/u.s)
                #cdelt3 = (chvel[l,1]/1.e5 - chvel[l,0]/1.e5)
                logger.info("Channel width %s", cdelt3)
            else:
                header['CTYPE3'] = 'FREQ'
                header['CUNIT3'] = 'Hz'
                crval3 = (restfreq[0]*(1.e0 - cdelt3.cgs/ct.c.cgs)).to(u.Hz)
                cdelt3 = (restfreq[0]*cdelt3.cgs/ct.c.cgs).to(u.Hz)
                logger.info("Channel width %s", cdelt3.to(u.MHz))

            header['CRPIX3'] = 1.
            header['CDELT3'] = cdelt3.value
            header['CRVAL3'] = crval3.value

            header['RESTFREQ'] = restfreq[l].to(u.Hz).value

            #logger.info("Channel width %.3e km/s", header['CDELT3'])

            ### Write cube in native units (K) ###

            header_K = header.copy()
            header_K['BUNIT'] = 'K'
            if write == 'indiv_K':
                fits.writeto(filename, data[l,v,:,:,:].transpose(), 
                             header, clobber=True)

            ### Write cube in Jy/pixel ###

            #K_to_Jy_per_beam = (header['RESTFREQ'] / 1e9) ** 2 * \
            #        header['BMAJ']* header['BMIN'] * 3600 ** 2 / 1.224e6

            #header_Jy_per_beam = header.copy()
            #header_Jy_per_beam['BUNIT'] = 'Jy/beam'
            #if write == 'indiv_jy_beam':
            #    fits.writeto(filename, 
            #                 data[l,v,:,:,:].transpose() * K_to_Jy_per_beam, 
            #                 header, clobber=True)

            # Avoid referring to beam since for some cases beam = 0
            pixel_area_deg = np.abs(header['CDELT1']) * header['CDELT2']
            #beam_area_deg = 1.1331 * header['BMAJ'] * header['BMIN']
            #pixels_per_beam = beam_area_deg / pixel_area_deg
            #K_to_Jy_old = K_to_Jy_per_beam / pixels_per_beam
            K_to_Jy = (header['RESTFREQ']*u.Hz).to(u.GHz).value** 2 * 3600 ** 2 / 1.224e6 /\
                    1.1331 * pixel_area_deg

            #logger.info("K to Jy/beam = %.3e", K_to_Jy_per_beam)
            logger.info("K to Jy/pixel = %.3f", K_to_Jy)
            #logger.info("K to Jy/pixel [old]= %.3f", K_to_Jy_old)

            header_Jy = header.copy()
            header_Jy['BUNIT'] = 'Jy'
            if write == 'indiv_jy':
                fits.writeto(filename,
                             data[l,v,:,:,:].transpose() * K_to_Jy, 
                             header, clobber=True)

    ### Now make a stacked cube ###

    logger.info("Writing out stacked cube...")

    for line in range(nlines):
        datamax = np.nanmax(data[line,:,:,:,:])
        logger.info('Line %s, datamax: %.3e', line, datamax)

    minimum_line = nlines - 1
    minimum_velocity = chvel[nlines - 1,0]

    # Now figure out the velocity shift due to line frequencies

    velocity_shift = -1. * ct.c.cgs * (restfreq - restfreq[minimum_line]) / \
            restfreq[minimum_line]

    for line in range(nlines):
        logger.info('Velocity shift: %s', velocity_shift[line].to(u.km/u.s))

    maximum_velocity = chvel[0,nchan[0]-1] + velocity_shift[0]
    logger.info('Min and max velocities: %s, %s',
                minimum_velocity.to(u.km/u.s), maximum_velocity.to(u.km/u.s))

    # Make a new velocity array starting at the minimum velocity
    dv = (chvel[0,1]-chvel[0,0])
    logger.info('Channel width %s', dv.to(u.km/u.s))

    number_channels = int(((maximum_velocity.cgs - minimum_velocity.cgs) /\
            dv.cgs).value) + 1
    logger.info('Number of channels: %i', number_channels)

    velo = np.linspace(0., number_channels-1., number_channels)*dv + \
            minimum_velocity
    logger.info('New velocity range [%s:%s]', velo[0].to(u.km/u.s),
                velo[number_channels-1].to(u.km/u.s))

    # New array to hold 1 line with all the spectra
    # from the 11 lines
    data1 = np.zeros((1,number_channels))

    for i in range(number_channels):
        data1[0,i] = velo[i].cgs.value

    data2 = np.zeros((nviews,nx,ny,number_channels))

    for line in range(nlines):
        for i in range(nchan[line]):
            j = int(((chvel[line,i].cgs + velocity_shift[line].cgs -\
                    velo[0].cgs)/dv.cgs).value)
            data2[:,:,:,j] = (data[line,:,:,:,i] + data2[:,:,:,j])

    nchan[0] = number_channels

    images = []
    for v in range(nviews):

        header = header_template.copy()

        header['LINE'] = 'CH3CN (all)'

        cdelt3 = dv.to(u.km/u.s)
        if velocity:
            header['CTYPE3'] = 'VELO-LSR'
            header['CUNIT3'] = 'KM/S'
            crval3 = minimum_velocity.to(u.km/u.s)
            #cdelt3 = dv/1.e5
            logger.info("Channel width %s", cdelt3)
        else:
            header['CTYPE3'] = 'FREQ'
            header['CUNIT3'] = 'Hz'
            crval3 = (restfreq[0]*(1.e0 - cdelt3.cgs/ct.c.cgs)).to(u.Hz)
            cdelt3 = (restfreq[0]*cdelt3.cgs/ct.c.cgs).to(u.Hz)
            logger.info("Channel width %s", cdelt3.to(u.MHz))

        header['CRPIX3'] = 1.
        header['CDELT3'] = cdelt3.value
        header['CRVAL3'] = crval3.value

        header['RESTFREQ'] = restfreq[l].to(u.Hz).value

        #logger.info("Channel width %.3e km/s", cdelt3)

        header_K['BUNIT'] = 'K'

        if write == 'combined_K':
            fits.writeto(filename, data2[v,:nchan[0]].transpose(), 
                         header, clobber=True)

        #K_to_Jy_per_beam = (header['RESTFREQ']*u.Hz).to(u.GHz).value** 2 * header['BMAJ'] *\
        #        header['BMIN'] * 3600 ** 2 / 1.224e6

        #header_Jy_per_beam = header.copy()
        #header_Jy_per_beam['BUNIT'] = 'Jy/beam'
        #if write == 'combined_jy_beam':
        #    fits.writeto(filename,
        #                 data2[v,:nchan[0]].transpose() * K_to_Jy_per_beam, 
        #                 header, clobber=True)

        # Avoid referring to beam since for some cases beam = 0
        pixel_area_deg = np.abs(header['CDELT1']) * header['CDELT2']
        #beam_area_deg = 1.1331 * header['BMAJ'] * header['BMIN']
        #pixels_per_beam = beam_area_deg / pixel_area_deg
        #K_to_Jy_old = K_to_Jy_per_beam / pixels_per_beam
        K_to_Jy = (header['RESTFREQ']*u.Hz).to(u.GHz).value**2 * 3600**2 / 1.224e6 / \
                1.1331 * pixel_area_deg

        if write == 'combined_jy':
            header['BUNIT'] = 'JY/PIXEL'
            #data2[v] = data[v][:,::-1,:]
            #newdata2 = data2[v,:,:,:nchan[0]].transpose()
            fits.writeto(filename,
                         data2[v,:,:,:nchan[0]].transpose() * K_to_Jy, 
                         header, overwrite=True)
            images += [Data3D(filename)]

    return images

def load_model_cube(model_file, source, filename, pa=0*u.deg, 
        logger=get_logger(__name__,__package__+'.log'), velocity=True, 
        bunit=u.Jy):
    """Load a Mollie model.

    Written by: K. G. Johnston.
    Modified by: F. Olguin
    References to older versions were removed

    Parameters:
        model_file (str): model file name.
        source (astroSource): source information.
        filename (str): output file name.
        pa (astropy.quantity, default=0.): source position angle.
        logger (logging, optional): logger.
        velocity (bool, default=True): output cube 3rd axis in velocity or frequency.
        bunit (astropy.unit, default=Jy): output flux unit.
    """

    # Open file
    logger.info('Opening file: %s', os.path.basename(model_file))
    f = open(model_file, "rb")
    endian = '>'
    byte = f.read(4)
    nlines = struct.unpack(endian+'l', byte)[0]
    
    # Swap bytes
    swap_bytes = False
    if nlines > 200:
        swap_bytes = True
        endian = '<'
    logger.info('Swapping bytes? %s', swap_bytes)

    # Number of lines
    nlines = struct.unpack(endian+'l', byte)[0]
    logger.info('Number of lines: %i', nlines)

    # Number of viewing angles
    nviews = struct.unpack(endian+'l',f.read(4))[0]
    logger.info('Number of viewing angles: %i', nviews)

    # Number of channels
    nchan = np.zeros(nlines,dtype=int)
    for line in range(nlines):
        nch = struct.unpack(endian+'l',f.read(4))[0]
        nchan[line] = nch
    logger.info('Number of channels: %r', nchan)

    # Grid parameters
    nx = struct.unpack(endian+'l',f.read(4))[0]
    ny = struct.unpack(endian+'l',f.read(4))[0]
    cellx = struct.unpack(endian+'f',f.read(4))[0] * u.pc
    celly = struct.unpack(endian+'f',f.read(4))[0] * u.pc
    beamx = struct.unpack(endian+'f',f.read(4))[0] * u.pc
    beamy = struct.unpack(endian+'f',f.read(4))[0] * u.pc
    logger.info('Grid size nx=%i, ny=%i', nx, ny)
    logger.info('Cell size %sx%s', cellx, celly)

    # Line names
    maxname = 20
    linename = (nlines)*['']
    for i in range(nlines):
        for j in range(maxname):
            bytvar = struct.unpack(endian+'c',f.read(1))[0]
            linename[i] += bytvar.decode('ascii')
        linename[i] = linename[i].strip()

    # Rest frequencies
    restfreq = np.zeros(nlines) * u.Hz
    for i in range(nlines):
        restfreq[i] = struct.unpack(endian+'d',f.read(8))[0]*u.Hz
        logger.info('The lines are:')
        logger.info('%s at %s', linename[i], restfreq[i].to(u.GHz))

    # Channel velocity ranges
    maxch = max(nchan)
    chvel = np.zeros((nlines,maxch)) * u.cm/u.s
    for i in range(nlines):
        for n in range(nchan[i]):
            chvel[i,n] = struct.unpack(endian+'f',f.read(4))[0] * u.cm/u.s
        logger.info('Velocity range for line %s:  %s, %s', linename[i],
                    chvel[i,0].to(u.km/u.s), chvel[i,nchan[i]-1].to(u.km/u.s))

    # Viewing angles
    lng = np.zeros(nviews) * u.deg
    lat = np.zeros(nviews) * u.deg
    for i in range(nviews):
        lng[i] = struct.unpack(endian+'f',f.read(4))[0] * u.deg
        lat[i] = struct.unpack(endian+'f',f.read(4))[0] * u.deg
    logger.info('Longitudes: %s', lng)

    # Fix inclination convention to Hyperion
    lat = 90.*u.deg - lat
    logger.info('Latitudes: %s', lat)

    # RA axis
    xc = np.zeros(nx)
    for i in range(nx):
        xc[i] = struct.unpack(endian+'f',f.read(4))[0]

    # Dec axis
    yc = np.zeros(ny)
    for i in range(ny):
        yc[i] = struct.unpack(endian+'f',f.read(4))[0]

    # Data
    data = np.ones((nlines,nviews,nx,ny,maxch)) * np.nan
    for l in range(nlines):
        data_bytes = f.read(4 * nviews * nx * ny * nchan[l])
        data[l,:,:,:,:nchan[l]] = np.fromstring(data_bytes, 
                dtype=endian + 'f4').reshape(nviews, nx, ny, nchan[l])
        logger.info('Max in line %i is %.3e', l, np.nanmax(data[l]))
    f.close()
    logger.info('Min and max brightness in data set: %.3e, %.3e',
                np.nanmin(data), np.nanmax(data))

    # Set up header common to all files
    ra, dec = source.position.ra, source.position.dec
    distance = source.distance.to(u.pc)
    header_template = fits.Header()
    header_template['OBJECT'] = 'MODEL'
    header_template['TELESCOP'] = 'MOLLIE'
    header_template['INSTRUME'] = 'MOLLIE'
    header_template['OBSERVER'] = 'MOLLIE'
    header_template['CTYPE1'] = 'RA---SIN'
    header_template['CTYPE2'] = 'DEC--SIN'
    header_template['CUNIT1'] = 'degree'
    header_template['CUNIT2'] = 'degree'
    header_template['CRPIX1'] = nx / 2.
    header_template['CRPIX2'] = ny / 2. + 1.
    header_template['CDELT1'] = -1.*np.abs(np.degrees((cellx.si/distance.si).value))
    header_template['CDELT2'] = np.degrees((celly.si/distance.si).value)
    header_template['CRVAL1'] = ra.to(u.deg).value
    header_template['CRVAL2'] = dec.to(u.deg).value
    header_template['EPOCH'] = 2000
    header_template['EQUINOX'] = 2000.
    if pa or source.get_quantity('pa') is not None:
        header_template['CROTA1'] = 0
        if pa is None:
            pa = source.get_quantity('pa')
        header_template['CROTA2'] = (360*u.deg - pa.to(u.deg)).value

    # Line indices
    minimum_line = nlines - 1
    minimum_velocity = chvel[nlines - 1,0]

    # Now figure out the velocity shift due to line frequencies
    velocity_shift = -1. * ct.c.cgs * (restfreq - restfreq[minimum_line]) / \
            restfreq[minimum_line]
    for line in range(nlines):
        logger.info('Velocity shift for line %s: %s', linename[line],
                velocity_shift[line].to(u.km/u.s))
    
    # Maximum velocity
    maximum_velocity = chvel[0,nchan[0]-1] + velocity_shift[0]
    logger.info('Min and max velocities: %s, %s',
                minimum_velocity.to(u.km/u.s), maximum_velocity.to(u.km/u.s))

    # Make a new velocity array starting at the minimum velocity
    dv = (chvel[0,1]-chvel[0,0])
    logger.info('Channel width %s', dv.to(u.km/u.s))
    number_channels = int(((maximum_velocity.cgs - minimum_velocity.cgs) /\
            dv.cgs).value) + 1
    logger.info('Number of channels: %i', number_channels)
    velo = np.linspace(0., number_channels-1., number_channels)*dv + \
            minimum_velocity
    logger.info('New velocity range [%s:%s]', velo[0].to(u.km/u.s),
                velo[number_channels-1].to(u.km/u.s))

    # New array to hold 1 line with all the spectra
    cube = np.zeros((nviews,nx,ny,number_channels))
    for line in range(nlines):
        for i in range(nchan[line]):
            j = int(((chvel[line,i].cgs + velocity_shift[line].cgs -\
                    velo[0].cgs)/dv.cgs).value)
            cube[:,:,:,j] = (data[line,:,:,:,i] + cube[:,:,:,j])
    nchan[0] = number_channels

    # Save images per viewing angle
    images = []
    for v in range(nviews):
        # Header
        header = header_template.copy()
        header['LINE'] = '%s (all)' % linename[0].split('(')[0]

        # Save velocity or frequency
        cdelt3 = dv.to(u.km/u.s)
        if velocity:
            header['CTYPE3'] = 'VELO-LSR'
            header['CUNIT3'] = 'KM/S'
            crval3 = minimum_velocity.to(u.km/u.s)
            logger.info("Channel width %s", cdelt3)
        else:
            header['CTYPE3'] = 'FREQ'
            header['CUNIT3'] = 'Hz'
            crval3 = (restfreq[0]*(1. - cdelt3.cgs/ct.c.cgs)).to(u.Hz)
            cdelt3 = (restfreq[0]*cdelt3.cgs/ct.c.cgs).to(u.Hz)
            logger.info("Channel width %s", cdelt3.to(u.MHz))
        header['CRPIX3'] = 1.
        header['CDELT3'] = cdelt3.value
        header['CRVAL3'] = crval3.value
        header['RESTFREQ'] = restfreq[l].to(u.Hz).value

        # Save file
        logger.info('Cube file name: %s', os.path.basename(filename))
        if bunit is u.K:
            header['BUNIT'] = 'K'
            logger.info('Saving cube with units: K')
            fits.writeto(filename, cube[v,:nchan[0]].transpose(), 
                         header, overwrite=True)
            images += [Data3D(filename)]
        elif bunit is u.Jy:
            pixel_area_deg = np.abs(header['CDELT1']) * header['CDELT2']
            K_to_Jy = (header['RESTFREQ']*u.Hz).to(u.GHz).value**2 * \
                    3600**2 / 1.224e6 / 1.1331 * pixel_area_deg
            header['BUNIT'] = 'JY/PIXEL'
            logger.info('Saving cube with units: Jy/pixel')
            fits.writeto(filename,
                    cube[v,:,:,:nchan[0]].transpose() * K_to_Jy, header, 
                    overwrite=True)
            images += [Data3D(filename)]
        else:
            raise NotImplementedError

    return images

