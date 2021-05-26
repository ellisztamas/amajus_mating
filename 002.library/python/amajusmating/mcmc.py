import numpy as np
import faps as fp
from pprint import pprint
import sys
from time import time, strftime
from tqdm import tqdm

from amajusmating.dispersal import dispersal_GND

from scipy.stats import beta
from scipy.stats import gamma as gma

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

def run_MCMC(faps_data, initial_parameters, proposal_sigma, priors, nreps, output_dir, chain_name, thin=1, max_distance = np.inf):
        """
        A wrapper function to run Metropolis-Hastings MCMC for paternity and dispersal.

        Parameters
        ----------
        initial_parameters: dict
            Dictionary of starting values for the parameters, with the names as keys and 
            values as floats for the parameter values. Usually this would include the
            'shape' and 'scale' parameters of for the generalised normal distribution,
            the proportion of 'missing' fathers and mixture parameter 'lambda'.
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
        current_model = initial_parameters

        # For the first iteration, start with an arbitary probabilities
        # It doesn't matter, because for MH-ratio skips the first iteration by default.
        current_model['log_prior']     = -10e12
        current_model['loglik']        = -10e12
        current_model['log_posterior'] = -10e12
        # Set up the datafiles
        out_file = output_dir + chain_name
        setup_output(out_file, current_model, proposal_sigma, nreps, thin)

        #  IMPORT AND FORMAT DATA.
        t0 = time()

        # RUN THE MCMC.
        for i in tqdm(range(nreps)):
            # UPDATE PARAMETERS
            new_model = update_parameters(current_model, proposal_sigma)
            # Calculate log prior probabilities of the new model.
            prior_probs = list(priors(new_model).values())
            new_model['log_prior'] = np.log(prior_probs).sum()

            # LOG PROBABILITIES OF PATERNITY FOR EACH PARAMETER
            # Update proportion of missing fathers
            faps_data.update_missing_dads(new_model['missing'])
            # Update dispersal
            faps_data.update_dispersal_probs(
                scale = new_model['scale'],
                shape = new_model['shape'],
                mixture = new_model['mixture']
            )
            # Identify candidates who are further than the distance threshold
            # and set their log likelihoods to negative infinity
            ix = faps_data.distances > max_distance
            faps_data.covariates['dispersal'][ix] = -np.inf

            # INFER FAMILIES
            # Incorporate covariate information into paternity_arrays
            cov = sum(faps_data.covariates.values())
            for (p,s) in zip(faps_data.paternity, cov):
                faps_data.paternity[p].add_covariate(s)
            # Cluster into families and get likelihoods
            faps_data.sibship_clustering(ndraws=100, use_covariates=True)

            # Probability components for the new model.
            new_model['loglik'] = np.array([fp.alogsumexp(s.lik_partitions) for s in faps_data.sibships.values()]).sum()
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

# def run_MCMC(paternity_arrays, distance_matrix, initial_parameters, proposal_sigma, nreps, output_dir, chain_name, thin=1, max_distance = np.inf):
#     """
#     A wrapper function to run Metropolis-Hastings MCMC for paternity and dispersal.


#     Parameters
#     ----------
#     paternity_arrays: dict
#         Dictionary of paternityArray objects for each full sibling family.
#     distance_matrix: array
#         Array of distances between the mothers and all candidate fathers.
#     initial_parameters: dict
#         Dictionary of starting values for the parameters, with the names as keys and 
#         values as floats for the parameter values. Usually this would include the
#         'shape' and 'scale' parameters of for the generalised normal distribution,
#         the proportion of 'missing' fathers and mixture parameter 'mixture'.
#     proposal_sigma: dict
#         Dictionary giving the standard deviation of the normal distribution by which
#         to peturb parameters at each iteration. Should have the same keys as 
#         `initial_parameters`.
#     nreps: int
#         Number of iterations to run
#     output_dir: str
#         Directory to save the output.    
#      chain_name: str
#         Name for the output file, without a suffix.
#     thin: int
#         Optional thinning argument. If >1, every one in `thin` samples will be written
#         to disk.
#     max_distance: float, int
#         Maximum distance from the mother a candidate may be. Candidates further than
#         this value will have their posterior probability of paternity set to zero.
#         This is equivalent to setting a threshold prior on distance.
    
#     """
#     current_model = initial_parameters

#     # Set up the datafiles
#     out_file = output_dir + chain_name + ".txt"
#     mcmc.setup_output(out_file, current_model, proposal_sigma, nreps, thin)
#     # setup_output(out_file, current_model.keys())
    
#     # # Set up log file
#     # log_file = open(output_dir + chain_name + ".log", 'w')
#     # log_file.write('Metropolis-Hasting analysis mating in A. majus begun {}.'.format(strftime("%Y-%m-%d %H:%M:%S")))
#     # log_file.write('Initial model:\n')
#     # pprint(current_model, log_file)
#     # log_file.write('\nGaussian noise is applied at each iteration with standard deviations:\n')
#     # pprint(proposal_sigma, log_file)
#     # log_file.write("\nPerforming a total of {} steps, thinning every {} iteration. Output will be saved to:\n{}".format(nreps, thin, out_file))
#     # log_file.write('\nAnalysis begun {}.\n\n'.format(strftime("%Y-%m-%d %H:%M:%S")))

#     #  IMPORT AND FORMAT DATA.
#     t0 = time()

#     # RUN THE MCMC.
#     log_file.write('MCMC set up. Beginning Markov chain...\n')
#     log_file.close()
#     t0 = time()
#     track_accept = 0 # to record mean MH-acceptance rate
#     for i in tqdm(range(nreps)):
#         # UPDATE PARAMETERS
#         new_model = update_parameters(current_model, proposal_sigma)
#         # # Update distance travelled using new parameters
#         # new_model['mean_dist'] = stdev_GND(scale = new_model['scale'],
#         #                                    shape = new_model['shape'])

#         # LOG PROBABILITIES OF PATERNITY FOR EACH PARAMETER
#         # Update proportion of missing fathers
#         for p in paternity_arrays.keys():
#             paternity_arrays[p].missing_parents = new_model['missing']
#         # Update dispersal probabilities
#         if new_model['mixture'] > 1.0: new_model['mixture'] = 1.0
#         # Probability of drawing each male under GND dispersal
#         prob_drawn = dispersal_GND(
#             x     = distance_matrix,
#             scale = new_model['scale'],
#             shape = new_model['shape'],
#             w     = new_model['mixture'])
#         # Identify candidates who are further than the threshold distance
#         # and set their log likelihoods to negative infinity
#         ix = distance_matrix > max_distance
#         prob_drawn[ix] = -np.inf
#         # Incorporate into paternity_arrays
#         for (p,s) in zip(paternity_arrays.keys(), prob_drawn):
#             paternity_arrays[p].add_covariate(s)

#         # INFER FAMILIES
#         # Cluster into families and get likelihoods
#         sc = fp.sibship_clustering(paternity_arrays, ndraws=100, use_covariates=True)
#         new_model['loglik'] = np.array([fp.alogsumexp(s.lik_partitions) for s in sc.values()]).sum()

#         # PRIORS
#         # Incorporate prior probabilities.
#         prior_probs = {
#             'missing'   : beta.pdf(new_model['missing'],a=3,   b=15),
#             'mixture'    : beta.pdf(new_model['mixture'], a=1.1, b=1.1),
#             'shape'     : gma.pdf(new_model['shape'],   a=10, scale = 1/5),
#             # 'mean_dist' : gma.pdf(new_model['mean_dist'],      a=2,   scale= 200),
#             'scale'     : gma.pdf(new_model['shape'],   a=6, scale = 50),
#         }
#         # Log prior probabilities
#         prior_probs = {k: np.log(v) for k,v in prior_probs.items()}
#         log_prior   = np.array(list(prior_probs.values())).sum()
#         # Sum log probabilities and incorporate
#         new_model['loglik'] += log_prior

#         # Decide whether to accept the new model.
#         accept = mh_ratio(
#             current = current_model['loglik'],
#             new     = new_model['loglik']
#         )
#         if accept:
#             current_model = new_model
#             track_accept += 1

#         mean_acceptance = float(track_accept) / (i+1)
#         # write iteration to disk
#         if(i in np.arange(start = 0, stop = nreps, step=thin)):
#             write_output(out_file, i, mean_acceptance, current_model, time0=t0)

#     log_file = open(output_dir + chain_name + ".log", 'a')
#     log_file.write('\nMCMC completed {}.\n'.format(strftime("%Y-%m-%d %H:%M:%S")))
#     log_file.close()
