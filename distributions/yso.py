import numpy as np

import ulrich, outflow, cavity

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
    
