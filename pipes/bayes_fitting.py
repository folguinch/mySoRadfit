import emcee as mc

from .rt_pipe import rt_pipe
from ..model_bayes import modelBayes
from ..model_output import modelOutput

def log_likelihood(data, model, outconfig):
    # Setup RTs and run RTs
    files = rt_pipe(model)

    # Process output
    modelout = modelOutput(model.name, outconfig)
    modelout.load_all(files, PA=model.get_pa(), vlsr=model.get_vlsr())

    # Compare
    for d, mod in zip_str(data, modelout):
        likelihood += -0.5*np.sum((d-mod)**2/e**2 + np.log(2.*np.pi*e**2))

def log_posterior(theta, *args, **kwargs):
    # Update model
    model = kwargs['model']
    model.update(theta)

    # Get priors
    logprior = model.get_logpriors()

    # Get likelihood
    loglikelihood = log_likelihood(args[0], model, kwargs['outconfig'])

    # Posterior
    return logprior + loglikelihood

def bayes_pipe(args):
    """Bayesian fit pipe.

    Parameters:
        args (argparse.parser): argument parser
    """
    # Load data

    # Create initial model
    model = modelBayes(args.model)
    
    # Run MCMC
    sampler = mc.EnsembleSampler(nwalkers, ndim, log_posterior, args=(data,),
            kwargs={'model': model, 'outconfig':args.outconfig})
    sampler.run_mcmc()



