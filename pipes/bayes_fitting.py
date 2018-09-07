import emcee as mc

from myutils.logger import get_logger

from .rt_pipe import rt_pipe
from ..model_bayes import modelBayes
from ..model_output import modelOutput

def log_likelihood(data, model, outconfig, logger=get_logger(__name__)):
    # Setup RTs and run RTs
    files = rt_pipe(model, logger=logger)
    exit()

    # Process output
    modelout = modelOutput(model.name, config=outconfig)
    modelout.load_all(files, PA=model.get_pa(), vlsr=model.get_vlsr())

    # Compare
    for d, mod in zip_str(data, modelout):
        likelihood += -0.5*np.sum((d-mod)**2/e**2 + np.log(2.*np.pi*e**2))

def log_posterior(theta, *args, **kwargs):
    # Update model
    model = kwargs['model']
    validate = model.update_params(theta)

    # Get priors
    logprior = model.get_logpriors()

    # Get likelihood
    loglikelihood = log_likelihood(args[0], model, kwargs['outconfig'],
            logger=kwargs['logger'])
    exit()

    # Posterior
    return logprior + loglikelihood

def bayes_pipe(args):
    """Bayesian fit pipe.

    Parameters:
        args (argparse.parser): argument parser
    """
    # Create initial model
    print('-'*80)
    args.logger.info('Loading Bayes inital model')
    model = modelBayes(args.model)
    ndim = model.get_dimensions()
    nwalkers = args.nwalkers[0] * ndim
    
    # Configure sampler
    args.logger.info('Setting up the sampler:')
    args.logger.info('Number of walkers: %i', nwalkers)
    args.logger.info('Number of parameters: %i', ndim)
    sampler = mc.EnsembleSampler(nwalkers, ndim, log_posterior,
            args=(args.source,),
            kwargs={'model': model, 'outconfig':args.outconfig,
                'logger':args.logger})

    # Initial guesses
    p0 = model.get_random_initial(nwalkers)

    # Run MCMC
    args.logger.info('Running MCMC')
    pos, prob, state = sampler.run_mcmc(p0, args.steps[0])

