from itertools import product

import numpy as np

class Cell(object):
    """Defines a cell object.

    Attributes:
        xwalls: x-coordinate of the cell walls.
        ywalls: y-coordinate of the cell walls.
        zwalls: z-coordinate of the cell walls.
        props: cell physical properties.
    """

    def __init__(self, *walls, **props):
        """Create a new cell object.

        Parameters:
            walls (tuple): the cell walls.
            props (dict): cell physical properties.
        """
        assert len(walls)==6
        self.xwalls = walls[0:2]
        self.ywalls = walls[2:4]
        self.zwalls = walls[4:6]
        self.props = props

    def __getitem__(self, key):
        return self.props[key]

    def __setitem__(self,key, val):
        self.props[key] = val

    @property
    def center(self):
        """Coordinates of the center."""
        center = []
        for walls in self:
            center += [(walls[0]+walls[1])/2]
        return center

    @property
    def volume(self):
        """Cell volume."""
        vol = 1.
        for walls in self:
            step = np.abs(walls[0]-walls[1])
            vol = vol * step
        return vol

    @property
    def mass(self):
        """Cell mass."""
        assert 'density' in self.props
        return self['density']*self.volume

    @property
    def sizes(self):
        """Cell size per axis."""
        return [np.abs(walls[0]-walls[1]) for walls in self]

    def __iter__(self):
        return self

    def next(self):
        for walls in [self.xwalls, self.ywalls, self.zwalls]:
            return walls

    def divide_by(self, n):
        """Divide each cell axis by n.

        The walls of the new cells are returned and are iterated first along x
        then y and finally z as needed by hyperion.

        Parameters:
            n: number of new cells per axis.

        Returns:
            walls: the walls of the new cells.
        """
        # New walls
        wall_pairs = []
        for walls in self:
            step = np.abs(walls[0]-walls[1]) / n
            aux = min(walls)
            wall_pairs += [[(aux+step*i,aux+step*(i+1)) for i in range(int(n))]]
         
        # New cells
        cells = []
        for x,y,z in product(*wall_pairs[::-1]):
            cells += [x+y+z]

        return cells



