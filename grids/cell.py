from itertools import product

import numpy as np

class Cell(object):
    """Defines a cell object.

    Attributes:
        xwalls: x-coordinate of the cell walls.
        ywalls: y-coordinate of the cell walls.
        zwalls: z-coordinate of the cell walls.
    """

    def __init__(self, *walls):
        """Create a new cell object.

        Parameters:
            walls (tuple): the cell walls.
        """
        assert len(walls)==6
        self.xwalls = walls[0:2]
        self.ywalls = walls[2:4]
        self.zwalls = walls[4:6]

    @property
    def center(self):
        """Coordinates of the center."""
        center = []
        for walls in self:
            center += [(walls[0]+walls[1])/2]
        return center

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
            aux = np.min(walls)
            wall_pairs += [[(aux+step*i,aux+step*(i+1)) for i in range(int(n))]]
         
        # New cells
        cells = []
        for x,y,z in product(*wall_pairs[::-1]):
            cells += [x+y+z]

        return cells



