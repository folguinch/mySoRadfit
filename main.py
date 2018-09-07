import os, argparse
from configparser import ConfigParser

#from myutils.logger import get_logger
import astropy.units as u
from myutils.argparse_actions import startLogger, LoadConfig
from astroSource.source import LoadSourcefromConfig

#from mySoRadfit import ModelfromConfig
from mySoRadfit.pipes.bayes_fitting import bayes_pipe
#from utils.parsers.args import LoadConfig

#def criterion(cell):
#    if cell.mass >= 1*u.M_sun:
#        return True
#    else:
#        return False

def get_parser():
    parser = argparse.ArgumentParser()
    #config = parser.add_argument('--config', action=LoadConfig, default=config)
    parser.add_argument('--stats', action='store_true', 
            help='Compute stats only')
    parser.add_argument('--bayes', action='store_true', 
            help='Use MCMC+bayes model fitting')
    parser.add_argument('--logger', action=startLogger, 
            help='Main logger file name')
    parser.add_argument('--use', nargs='*', type=list, default=['all'], 
            help='Use these data')
    parser.add_argument('--nwalkers', nargs=1, type=int, default=[2],
            help='Number of walkers for MCMC as a function of parameters')
    parser.add_argument('--steps', nargs=1, type=int, default=[100],
            help='Number of MCMC steps')
    parser.add_argument('--source', action=LoadSourcefromConfig,
            help='Source configuration file')
    parser.add_argument('--outconfig', action=LoadConfig,
            help='Postprocessing configuration')
    #parser.add_argument('model', action=ModelfromConfig, 
    #        help='Model configuration file')
    parser.add_argument('model', type=str, 
            help='Model configuration file')
    return parser

def main():
    """Main function for the RT fitter"""

    # Read command line arguments
    parser = get_parser()
    args = parser.parse_args()

    # Workflow
    if args.bayes:
        args.logger.info('Running Bayesian modelling')
        bayes_pipe(args)
    else:
        raise NotImplementedError
    
if __name__=='__main__':
    main()
