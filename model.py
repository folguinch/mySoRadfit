import argparse

from .model_base import BaseModel
from .grids.octree import Octree
from .distributions.yso import YSO

class Model(BaseModel):
    """Model class.

    Attributes:
        config (myConfigParser): model configuration file name.
        params (OrderedDict): model physical parameters.
        setup (myConfigParser): model rt setup.
        images (myConfigParser): model images setup.
        logger: logging system.
    """

    def fill_grids(self):
        raise NotImplementedError

    def load_params(self):
        """Load parameter file for each model source."""
        for section in self.config.sections():
            self.logger.info('Loading parameters for: %s', section)
            self.params[section] = YSO(self.config[section]['params'],
                    loc=self.config.getfloatlist(section, 'loc'))

    def build_grid(self, criterion, max_depth=10):
        """Build the grid.
        
        Parameters:
            criterion (function): refinement criterion.
            max_depth (int, default=10): maximum level of recursion.
        """
        if self.setup['grid']['type'] == 'octree':
            self.logger.info('Building Octree grid')
            grid = Octree(self.setup.getquantity('grid','xsize'),
                    self.setup.getquantity('grid','ysize'),
                    self.setup.getquantity('grid','zsize'))
            grid.build(self.get_density(), criterion, max_depth=max_depth)

        return grid

    def load_grids(self, rt):
        """Load the grids from file.
        
        Parameters:
            rt: RT transfer config section.
        """
        grids = []
        for grid in self.setup.getlist(rt, 'grids'):
            self.logger.info('Loading grid: %s', os.path.basename(grid))
            fname = os.path.realpath(os.path.expanduser(grid))
            grid = load_struct_array(fname, usecols=None)
            grids += [grid]

        cell_sizes = self.setup.getintlist(rt, 'cell_sizes')

        return grids, cell_sizes

    def get_density(self):
        """Get the density function."""
        def density(x, y, z):
            den = 0.
            for yso in self.sources:
                den = den + yso(x, y, z)
            return den
        return density

class ModelfromConfig(argparse.Action):
    """Load model from configuration file as command line parameter"""

    def __call__(self, parser, namespace, values, option_string=None):
        config = os.path.realpath(os.path.expanduser(values))
        model = Model(config=config)
        return setattr(namespace, self.dest, model)
