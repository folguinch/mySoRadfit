import os
from string import Template

import numpy as np
import astropy.units as u
import astropy.constants as ct
from hyperion.model import ModelOutput
from scipy.interpolate import griddata
from myutils.logger import get_logger
#from Plotter.mesh_plotter import NMeshPlotter
from myutils.math import rebin_irregular_nd, map_sph_to_cart_axisym

def get_walls(*args):
    walls = []
    for val in args:
        dv = np.min(val[val>0.*u.cm])
        vw = val - dv
        vw = np.append(vw.value, vw[-1].value+2.*dv.value)*vw.unit
        walls += [vw]

    return walls

def rebin_2Dsph_to_cart(val, new_pos, yso, pos=None, rebin=False, interp=False,
        logger=get_logger(__name__), **kwargs):
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

        # Meshes
        YN, ZN = np.meshgrid(new_pos[1], new_pos[2])

        # Evaluate the function
        if interp:
            val1 = rebin_2Dsph_to_cart(val, (new_pos[0], np.array([0.])*new_pos[1].unit,
                new_pos[2]), yso, pos=pos, interp=True, logger=logger)
            val1 = val1[:,0,:]
        else:
            val1 = yso(YN, np.zeros(YN.shape), ZN)
            assert val1.cgs.unit == u.g/u.cm**3
            val1 = val1.cgs.value

        # Walls
        walls = get_walls(new_pos[1], new_pos[-1])
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
        val1 = yso(XN, YN, ZN)
        assert val1.cgs.unit == u.g/u.cm**3
        val1 = val1.cgs.value

        return val1

def set_physical_props(yso, grids, template, logger=get_logger(__name__)):
    """Calculate and write the physical properties of the model.

    Parameters:
        yso: the model parameters.
        grids: grids where the model will be evaluated
        template: filename of the *define_model.c* file.
        logger: logging system.
    """
    # Load quantities
    if yso.params.get('DEFAULT', 'quantities_from'):
        hmodel = yso.params.get('DEFAULT', 'quantities_from')
        logger.info('Loading Hyperion model: %s', os.path.basename(hmodel))
        hmodel = ModelOutput(os.path.expanduser(hmodel))
        q = hmodel.get_quantities()
        density = np.sum(q['density'].array, axis=0)
        temperature = np.sum(q['temperature'].array[1:], axis=0)
        r, th = q.r*u.cm, q.t*u.rad
        rw, tw = q.r_wall*u.cm, q.t_wall*u.rad

        # Volumes
        dr = np.abs(rw[1:] - rw[:-1])
        dt = np.abs(tw[1:] - tw[:-1])
        R, TH = np.meshgrid(r, th)
        DR, DTH = np.meshgrid(dr, dt)
        vol_sph = R**2 * np.sin(TH) * DR * DTH
    else:
        raise NotImplementedError

    ## Density plot
    #styles = ['jet']
    #ylevs = np.array([-1000., 0., 1000., 10000., 20000., 150000]) * u.au
    #zlevs = np.array([-1000., 0., 1000.]) * u.au
    #vmin = {'den':[-25, -21, -20], 'temp':[5, 5, 30]}
    #vmax = {'den':[-12, -14, -12], 'temp':[200, 200, 300]}
    #fig1 = NMeshPlotter(styles=styles, nxcbar=1, rows=3, cols=len(ylevs), 
    #        sharey=True, sharex=True, xsize=4., ysize=4., left=0.8, right=1, 
    #        top=0.6, hspace=0.3)
    #fig2 = NMeshPlotter(styles=styles, nxcbar=1, rows=3, cols=len(ylevs), 
    #        sharey=True, sharex=True, xsize=4., ysize=4., left=0.8, right=1, 
    #        top=0.6, hspace=0.3)
    #fig3 =  NMeshPlotter(styles=styles, nxcbar=1, rows=4, cols=len(ylevs), 
    #        sharey=True, sharex=True, xsize=4., ysize=4., left=0.8, right=1, 
    #        top=0.6, hspace=0.3)
    #fig4 =  NMeshPlotter(styles=styles, nxcbar=1, rows=4, cols=len(zlevs), 
    #        sharey=True, sharex=True, xsize=4., ysize=4., left=0.8, right=1, 
    #        top=0.6, hspace=0.3)
    #axes1,axes2,axes3,axes4 = [], [], [], []

    # Open template
    with open(os.path.expanduser(template)) as ftemp:
        templ = Template(ftemp.read())

    for i, grid in enumerate(grids):
        x = grid[0]['x'] * grid[1]['x']
        y = grid[0]['y'] * grid[1]['y']
        z = grid[0]['z'] * grid[1]['z']
        xi = np.unique(x)
        yi = np.unique(y)
        zi = np.unique(z)

        # Density
        dens = rebin_2Dsph_to_cart(density, (xi, yi, zi), yso, pos=(r, th), 
                rebin=i==2, interp=True, weights=vol_sph, logger=logger)
        dens = dens * u.g/u.cm**3
        mH = ct.m_p + ct.m_e
        dens = dens / (2.33 * mH)
        dens[dens.cgs<=0./u.cm**3] = 10./u.cm**3

        # Temperature
        temp = rebin_2Dsph_to_cart(temperature, (xi, yi, zi), yso, pos=(r, th), 
                rebin=i==2, interp=True, weights=density[0,:,:], logger=logger)
        temp = temp*u.K
        temp[temp<2.7*u.K] = 2.7*u.K

        # Velocity
        if i==2:
            dz = np.abs(zi[0]-zi[1])
        else:
            dz = None
        vx, vy, vz = yso.velocity(x, y, z, min_height_to_disc=dz)
        #vx = vx.to(u.km/u.s)
        #vy = vy.to(u.km/u.s)
        #vz = vz.to(u.km/u.s)

        ## Density
        #den_log = np.log10(dens)
        #den_log[np.isnan(den_log)] = -30
        #axes1 = plot_mesh_slices(fig1, den_log, xi, yi, zi, ylevs.to(yi.unit), 
        #        axes=axes1)
        ## Temperature
        #axes2 = plot_mesh_slices(fig2, temp, xi, yi, zi, ylevs.to(yi.unit), 
        #        axes=axes2, cblabel='Temperature (K)', vmin=vmin['temp'],
        #        vmax=vmax['temp'])
        ## Velocity over density
        #axes3 = plot_field_slices(fig3, vx, vy, vz, x, y, z, den_log, ylevs=ylevs,
        #        axes=axes3)
        #axes4 = plot_field_slices(fig4, vx, vy, vz, x, y, z, den_log, zlevs=zlevs,
        #        axes=axes4)

        # Abundance
        abundance = yso.abundance(temp)
        
        # Linewidth
        amu = 1.660531e-24 * u.g
        atoms = yso.params.getfloat('Velocity', 'atoms')
        c_s = ct.k_B * temp / (atoms * amu)
        linewidth = np.sqrt(yso.params.getquantity('Velocity', 'linewidth')**2
                + c_s)

        # Write template
        sci_fmt = lambda x: '%.8e' % x
        flt_fmt = lambda x: '%.8f' % x
        dens = ','.join(map(sci_fmt, np.ravel(dens.cgs.value,order='F')))
        temp = ','.join(map(flt_fmt, np.ravel(temp.value,order='F')))
        vx = ','.join(map(flt_fmt, np.ravel(vx.cgs.value,order='F')))
        vy = ','.join(map(flt_fmt, np.ravel(vy.cgs.value,order='F')))
        vz = ','.join(map(flt_fmt, np.ravel(vz.cgs.value,order='F')))
        abundance = ','.join(map(sci_fmt, np.ravel(abundance,order='F')))
        linewidth = ','.join(map(flt_fmt, np.ravel(linewidth.cgs.value,order='F')))
        templ = templ.safe_substitute(**{'temp%i'%i: temp, 'dens%i'%i: dens, 
            'vx%i'%i: vx, 'vy%i'%i: vy, 'vz%i'%i: vz, 'abn%i'%i: abundance,
            'linewidth%i'%i: linewidth})
        if i!=len(grids)-1:
            templ = Template(templ)

    # Save file
    dirname = os.path.dirname(os.path.expanduser(template))
    fname = 'define_model.c'
    with open(os.path.join(dirname,fname),'w') as out:
        out.write(templ)

    #fig1.savefig('hyperion_grid_3d_density.png')
    #fig2.savefig('hyperion_grid_3d_temperature.png')
    #fig3.savefig('hyperion_grid_3d_velocity_xz.png')
    #fig4.savefig('hyperion_grid_3d_velocity_xy.png')
