import numpy as np
import faps as fp
from pprint import pprint
import sys
from time import time, strftime

def update_parameters(current_model, sigma):
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

def setup_output(path, model, proposal_sigma, nreps, thin):
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
    log_file.write('Metropolis-Hasting analysis mating in A. majus begun {}.\n'.format(strftime("%Y-%m-%d %H:%M:%S")))
    log_file.write('Parameters to initiate the chain:\n')
    pprint(model, log_file)
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
    out = str(i) + '\t' + '\t' + str(t)
    for k in sorted(model.keys()):
            out = out + '\t' + str(round(model[k], decimals))
    # write iteration to disk
    with open('{}'.format(path), 'a') as f:
        f.write(out + '\n')
    f.close()


