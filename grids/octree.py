import numpy as np
from hyperion.grid import OctreeGrid

from .grid import GridABC

class Octree(GridABC):
    """Defines a Octree grid.

    Attributes:
        center: coordinate center.
        refined: refined cells.
        grid: grid cells.
        density: grid density.
    """

    def __init__(self, xsize, ysize, zsize, center=(0.,0.,0.,)):
        """Initialize a grid object.

        Parameters:
            xsize: size of the x-axis.
            ysize: size of the y-axis.
            zsize: size of the z-axis.
            center: coordinate system center.
        """
        # Check that sizes are in the same units
        ysize = ysize.to(xsize.unit)
        zsize = zsize.to(xsize.unit)

        # Initialize grid
        self.refined = [True]
        self.grid = [Cell(center[0]-xsize/2., center[0]+xsize/2.,
            center[0]-ysize/2., center[0]+ysize/2.,
            center[0]-zsize/2., center[0]+zsize/2., density=0.)]

    def build(self, params, density, criterion, max_depth=10):
        """Build the grid.

        For building the grid, the criterion is used to determine whether the
        grid needs to be refined. The criterion function has only one input
        which is a Cell object.

        The density function input are the coordinates of a point and the
        parameters defining the grid. The density is the total density
        independent from the dust properties.

        Parameters:
            params (myConfigParser): density configuration parameters.
            density (function): density function.
            criterion (function): defines when the recursion should stop.
            max_depth (int, default=10): maximum levels of recursion.
        """
        assert len(self.grid)==1

        self.refined, self.grid = build_grid(self.refined, self.grid, 
                density, criterion, max_depth, 1)

    @staticmethod
    def build_grid(refined, grid, params, density, criterion, max_depth, step):
        grid0 = grid[-1]
        for cell in grid0.divide_by(2.):
            # Evaluate density
            aux = Cell(*cell)
            aux['density'] = density(*aux.center, params)
            grid.append(aux)
            
            # Determine whether grid needs refinement
            refined.append(criterion(aux) and step<=max_depth)
            if refined[-1]:
                refined, grid = build_grid(refined, grid, params, density, 
                        criterion, max_depth, step+1)
        return refined, grid

    def to_hyperion(self, model):
        """Create a hyperion Octree grid.

        The density for each dust component must be set separatedly.
        
        Parameters:
            model (hyperion.Model): hyperion model to fill.
        """
        # Create grid
        x, y, z = self.grid[0].center
        dx, dy, dz = [di/2. for di in self.grid[0].sizes]
        model.set_octree_grid(x, y, z, dx.cgs, dy.cgs, dz.cgs, self.refined)

        return model
