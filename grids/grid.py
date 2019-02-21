from abc import ABCMeta, abstractmethod

import numpy as np

from .cell import Cell

class GridABC(object):
    """Defines a grid ABC.

    Attibutes:
        grid: grid cells.
        density: grid density.
    """
    __metaclass__ = ABCMeta

    def __init__(self):
        """Initialize a Grid object."""
        self.grid = None
        self.density = None

    def set_density(self, func, params):
        x, y, z = self.coords
        self.density = func(x, y, z, params)

    @property
    def coords(self):
        centers = []
        for cell in grid:
            centers += [cell.center]
        centers = np.array(centers)

        return centers[:,0], centers[:,1], centers[:,2]

    @abstractmethod
    def build(self):
        pass


