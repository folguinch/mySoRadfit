import os

import astropy.units as u
from myutils.decorators import register_function
from myutils.logger import get_logger
from myutils.montage import rotate

from ..utils.mollie import load_model_cube

def mollie_base_pipe(model_file, source, filename=None, PA=0.*u.deg,
        logger=get_logger(__name__)):
    # Load & save FITS file
    if filename is None:
        filename = model_file + '_raw.fits'
        rotated = model_file + '_rotated.fits'
    else:
        rotated = filename.replace('.fits', '_rotated.fits')
    cube = load_model_cube(model_file, source, filename, 
        logger=logger, pa=PA)

    # Rotate & save
    if PA!=0.*u.deg:
        logger.info('Rotating cube')
        rotate(filename, rotated, z=1.2, cube=True)
        logger.info('Cube written: %s', os.path.basename(rotated))
        return rotated
    else:
        return filename

def casa_base_pipe(config, model_file, source, PA):
    # Run basic pipeline
    filename = mollie_base_pipe(model_file, source, PA=PA)

    # Run CASA

    return concatname

@register_function
def mollie_casa_line_pv(config, model_file, source, PA=0.*u.deg,
        vlsr=0.*u.km/u.s, logger=get_logger(__name__)):
    # Run casa
    filename = casa_base_pipe(config, model_file, source, PA)
    print filename
    exit()

    # Get PV maps

    return data


