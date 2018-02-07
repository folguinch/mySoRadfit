import os
import struct
from string import Template

import numpy as np
from astropy.io import fits
import astropy.units as u
import astropy.constants as ct
from hyperion.model import ModelOutput
from scipy.interpolate import griddata
from myutils.logger import get_logger
from myutils.math import rebin_irregular_nd, map_sph_to_cart_axisym
from myutils.classes.data_3d import Data3D
from Plotter.mesh_plotter import NMeshPlotter

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
    for key, val in kwargs.items():
        hdu = fits.PrimaryHDU(val)
        hdu.writeto(os.path.join(dirname, key+'.fits'), overwrite=True)

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
        temperature[0,:,:], density[0,:,:] = fill_inner_radius(R, TH,
                temperature[0,:,:], density[0,:,:], yso)
    else:
        raise NotImplementedError

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
        # Number density
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
        vx = np.reshape(vx, dens.shape, order='F')
        vy = np.reshape(vy, dens.shape, order='F')
        vz = np.reshape(vz, dens.shape, order='F')

        # Abundance
        abundance = yso.abundance(temp)
        
        # Linewidth
        amu = 1.660531e-24 * u.g
        atoms = yso.params.getfloat('Velocity', 'atoms')
        c_s = ct.k_B * temp / (atoms * amu)
        linewidth = np.sqrt(yso.params.getquantity('Velocity', 'linewidth')**2
                + c_s)

        # Write FITS
        write_fits(os.path.dirname(os.path.expanduser(template)), 
                **{'temp%i'%i: temp.value, 'dens%i'%i: dens.cgs.value, 
                'vx%i'%i: vx.cgs.value, 'vy%i'%i: vy.cgs.value, 'vz%i'%i:
                vz.cgs.value, 'abn%i'%i: abundance,
                'lwidth%i'%i: linewidth.cgs.value})

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

def write_setup(section, model, template, rt='mollie'):
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
            'nlines_rt': model.setup.get(rt, 'nlines_rt'),
            'line_id': model.images.get(section, 'line_id'),
            'phi': model.images.get(section, 'phi'),
            'th': model.images.get(section, 'th'),
            'npix': model.images.get(section, 'npix'),
            'pixsize': model.images.getquantity(section,
                'pixsize').to(u.pc).value}

    # Write file
    dirname = os.path.dirname(os.path.expanduser(template))
    fname = 'setup.c'
    with open(os.path.join(dirname, fname), 'w') as out:
        out.write(temp.substitute(**kwd))

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

    header_template['BMAJ'] = np.degrees((beamx.si / distance.si).value * 2.35)
    header_template['BMIN'] = np.degrees((beamx.si / distance.si).value * 2.35)

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

            K_to_Jy_per_beam = (header['RESTFREQ'] / 1e9) ** 2 * \
                    header['BMAJ']* header['BMIN'] * 3600 ** 2 / 1.224e6

            header_Jy_per_beam = header.copy()
            header_Jy_per_beam['BUNIT'] = 'Jy/beam'
            if write == 'indiv_jy_beam':
                fits.writeto(filename, 
                             data[l,v,:,:,:].transpose() * K_to_Jy_per_beam, 
                             header, clobber=True)

            # Avoid referring to beam since for some cases beam = 0
            pixel_area_deg = np.abs(header['CDELT1']) * header['CDELT2']
            beam_area_deg = 1.1331 * header['BMAJ'] * header['BMIN']
            pixels_per_beam = beam_area_deg / pixel_area_deg
            K_to_Jy_old = K_to_Jy_per_beam / pixels_per_beam
            K_to_Jy = (header['RESTFREQ']*u.Hz).to(u.GHz).value** 2 * 3600 ** 2 / 1.224e6 /\
                    1.1331 * pixel_area_deg

            logger.info("K to Jy/beam = %.3e", K_to_Jy_per_beam)
            logger.info("K to Jy/pixel = %.3f", K_to_Jy)
            logger.info("K to Jy/pixel [old]= %.3f", K_to_Jy_old)

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

        K_to_Jy_per_beam = (header['RESTFREQ']*u.Hz).to(u.GHz).value** 2 * header['BMAJ'] *\
                header['BMIN'] * 3600 ** 2 / 1.224e6

        header_Jy_per_beam = header.copy()
        header_Jy_per_beam['BUNIT'] = 'Jy/beam'
        if write == 'combined_jy_beam':
            fits.writeto(filename,
                         data2[v,:nchan[0]].transpose() * K_to_Jy_per_beam, 
                         header, clobber=True)

        # Avoid referring to beam since for some cases beam = 0
        pixel_area_deg = np.abs(header['CDELT1']) * header['CDELT2']
        beam_area_deg = 1.1331 * header['BMAJ'] * header['BMIN']
        pixels_per_beam = beam_area_deg / pixel_area_deg
        K_to_Jy_old = K_to_Jy_per_beam / pixels_per_beam
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


