import numpy as np
import faps as fp
from pprint import pprint
import sys
from time import time, strftime
from tqdm import tqdm

from amajusmating.dispersal import dispersal_GND

from scipy.stats import beta
from scipy.stats import gamma as gma

def check_parameters(model, sigma):
    """
    Check that variable names are ones that `run_MCMC` will recognise, 
    and that initial parameters are valid. The add placeholder values for
    log prior, likelihood and posterior.

    Parameters
    ----------
    model: dict
        Dictionary of starting model parameters. Keys should be a subset of
        ['missing', 'shape', 'scale', 'mixture', 'assortment'] and values 
        floats giving initial values for those parameters, within appropriate
        boundaries.
    sigma: dict
        Dictionary of standard deviations from which to draw normal
        deviates to apply to each value in current_model.
    """
    if not isinstance(model, dict):
        raise TypeError("`model` should be a dictionary with parameter names as keys and starting values as values.")
    if not isinstance(sigma, dict):
        raise TypeError("`sigma` should be a dictionary with parameter names as keys and standard deviations as values.")
    if model.keys() != sigma.keys():
        raise ValueError("The keys in `model` do not match those in `sigma`")

    valid_parameters = ['missing', 'shape', 'scale', 'mixture', 'assortment']
    if not all([x in valid_parameters for x in model.keys()]):
        raise ValueError("Not all parameter names are valid. Valid names are 'missing', 'shape', 'scale', 'mixture', and 'assortment'.")

    # Check things stay within boundaries
    if model['mixture'] > 1.0 or model['mixture'] < 0:
        raise ValueError('"mixture" parameter should be between 0 and 1.')
    if model['missing'] > 1.0 or model['missing'] < 0:
        raise ValueError('"missing" parameter should be between 0 and 1.')
    if model['shape'] <= 0:
        raise ValueError("'shape' parameter should be positive.")
    if model['scale'] <= 0:
        raise ValueError("'scale' parameter should be positive.")
    if "assortment" in model.keys():
        if model['assortment'] > 1.0 or model['assortment'] < 0:
            raise ValueError('"assortment" parameter should be between 0 and 1.')

    output = dict(model)
    # set initial log probabilities to very small numbers
    output['loglik']        = -10e12
    output['log_prior']     = -10e12
    output['log_posterior'] = -10e12

    return output

def update_parameters(current_model:dict, sigma:dict):
    """
    Apply Gaussian peturbation to a dictionary of model parameters.

    Parameters
    ----------
    current_model: dict
        Dictionary of current parameter values
    sigma: dict
        Dictionary of standard deviations from which to draw normal
        deviates to apply to each value in current_model. At least
        one key must intersect with the keys of current_model.
    """    
    new_model = dict(current_model)
    # random deviates for each parameter
    dev = {k: np.random.normal(0, v) for k, v in sigma.items()}
    # apply deviates to each parameter
    for k in sigma.keys():
        new_model[k] += dev[k]
    
    # Check things stay within boundaries
    if new_model['mixture'] > 1.0: new_model['mixture'] = 1.0
    if new_model['mixture'] <   0: new_model['mixture'] = 0

    if new_model['missing'] > 1.0: new_model['missing'] = 1.0
    if new_model['missing'] <   0: new_model['missing'] = 0

    if new_model['shape'] <   0: new_model['shape'] = 0
    if new_model['scale'] <   0: new_model['scale'] = 0

    if "assortment" in current_model.keys():
        if new_model['assortment'] > 1.0: new_model['assortment'] = 1.0
        if new_model['assortment'] <   0: new_model['assortment'] = 0

    return new_model

def mh_ratio(current, new):
    """
    Decide whether to accept a proposed model over a previous model
    based on their respective (log) probabilities and the
    Metropolis-Hastings ratio.

    Parameters
    ----------
    current: float
        Log probability of the old model.
    new: float
        Log probability of the proposed model.

    Returns
    -------
    Logical. Returns True if the new model is accepted.
    """
    if new > current:
        accept = True
    elif bool(np.random.binomial(1, np.exp(new - current))):
        accept = True
    else:
        accept = False

    return accept

def setup_output(path, model, proposal_sigma, nreps, thin, max_distance = None):
    """
    Set up a blank text file to store model parameters
    Parameters
    ----------
    path: str
        Path and filename where dat should be stored. The suffix
        '.txt' will be automatically appended.
    model: list
        List of parameter names
    Returns
    -------
    A text file is created at the path specified with column headers
    only, showing iteration and time taken at each iteration,
    followed by columns for each parameter in the list model.
    """
    # Create a blank output table on disk.
    out = 'iter\thours'
    for k in sorted(model.keys()): out = out + '\t' + k
    out = out + '\n'
    output = open('{}'.format(path) + '.out', 'w')
    output.write(out)
    output.close()

    # Set up log file
    log_file = open(path + ".log", 'w')
    log_file.write('Metropolis-Hasting analysis mating in A. majus begun {} using FAPS version {}.\n'.format(strftime("%Y-%m-%d %H:%M:%S"), fp.__version__))
    log_file.write('Parameters to initiate the chain:\n')
    pprint(model, log_file)
    if max_distance is not None:
        log_file.write('Maximum dispersal distance allowed: {} metres.\n'.format(max_distance))
    log_file.write('\nGaussian noise is applied to these values at each iteration with standard deviations:\n')
    pprint(proposal_sigma, log_file)
    log_file.write("\nPerforming a total of {} steps, thinning every {} iteration. Output will be saved to:\n{}".format(nreps, thin, path + '.out'))
    log_file.write('\nAnalysis begun {}.\n\n'.format(strftime("%Y-%m-%d %H:%M:%S")))
    log_file.write('MCMC set up. Beginning Markov chain...\n')
    log_file.close()


def write_output(path, i, model, decimals = 3, time0=None):
    """
    Write the output of a model from the current iteration to a
    new line in the output file.
    Parameters
    ----------
    path: str
        Path and filename where the output file is stored.
    i: int
        Iteration index.
    model: dict
        Dictionary of model values with the same parameters
        as the target file.
    decimals: float
        Number of figures to report after the decimal point.
    time0: float
        Optional time from the start of the analysis. This
        will be compared with te current time, and the
        difference reported as the time taken to reach this
        iteration.
    Returns
    -------
    A line is added to the target file with the model values.
    """
    # Format time data
    if time0 is None:
        t = 'NA'
    else:
        t = np.round((time() - time0) / 3600, decimals)

    # Prepare a string of output data to be exported.
    out = str(i) + '\t' + str(t)
    for k in sorted(model.keys()):
            out = out + '\t' + str(round(model[k], decimals))
    # write iteration to disk
    with open('{}'.format(path), 'a') as f:
        f.write(out + '\n')
    f.close()

def run_MCMC(data, initial_parameters, proposal_sigma, priors, nreps, output_dir, chain_name, thin=1, max_distance = np.inf):
        """
        A wrapper function to run Metropolis-Hastings MCMC for paternity and dispersal.

        Parameters
        ----------
        initial_parameters: dict
            Dictionary of starting values for the parameters, with the names as keys and 
            values as floats for the parameter values. Usually this would include the
            'shape' and 'scale' parameters of for the generalised normal distribution,
            the proportion of 'missing' fathers and mixture parameter 'mixture'. Optionally,
            this can include 'assortment' probability p, giving the probability that individuals
            of the same flower colour mate, while individuals of different genotypes mate
            with probabilty 1-p
        proposal_sigma: dict
            Dictionary giving the standard deviation of the normal distribution by which
            to peturb parameters at each iteration. Should have the same keys as 
            `initial_parameters`.
        nreps: int
            Number of iterations to run
        output_dir: str
            Directory to save the output.    
        chain_name: str
            Name for the output file, without a suffix.
        thin: int
            Optional thinning argument. If >1, every one in `thin` samples will be written
            to disk.
        max_distance: float, int
            Maximum distance from the mother a candidate may be. Candidates further than
            this value will have their posterior probability of paternity set to zero.
            This is equivalent to setting a threshold prior on distance.
        
        """
        
        # Check that parameters supplied can be evaluated.
        # For the first iteration, start with an arbitary probabilities
        # It doesn't matter, because for MH-ratio skips the first iteration by default.
        current_model = check_parameters(initial_parameters, proposal_sigma)

        # Set up the datafiles
        out_file = output_dir + chain_name
        setup_output(out_file, current_model, proposal_sigma, nreps, thin)
        t0 = time()

        # RUN THE MCMC.
        for i in tqdm(range(nreps)):
            # UPDATE PARAMETERS
            new_model = update_parameters(current_model, proposal_sigma)
            data.update_covariate_probs(model=new_model, max_distance = max_distance)
            
            # INFER FAMILIES
             # Cluster into families and get likelihoods
            data.sibship_clustering(ndraws=100, use_covariates=True)

            # Probability components for the new model.
            prior_probs = list(priors(new_model).values())
            new_model['log_prior'] = np.log(prior_probs).sum()
            new_model['loglik'] = np.array([fp.alogsumexp(s.lik_partitions) for s in data.sibships.values()]).sum()
            new_model['log_posterior'] = new_model['loglik'] + new_model['log_prior']

            # Decide whether to accept the new model.
            accept = mh_ratio(
                current = current_model['log_posterior'],
                new     = new_model['log_posterior']
            )
            if accept:
                current_model = new_model

            # write iteration to disk
            if(i in np.arange(start = 0, stop = nreps, step=thin)):
                write_output(out_file + '.out', i, current_model, decimals=3, time0=t0)
        
        log_file = open(output_dir + chain_name + ".log", 'a')
        log_file.write('\nMCMC completed {}.\n'.format(strftime("%Y-%m-%d %H:%M:%S")))
        log_file.close()