# from warnings import warn
# from faps.alogsumexp import alogsumexp
import numpy as np
import pandas as pd
import faps as fp
from glob import glob
from tqdm import tqdm
import os

def simulate_mating(data, model, ndraws = 1000, max_distance = np.inf):
    """
    Posterior samples of mating events.

    Draw a sample of mating events from the posterior distributions of sibship stuctures
    and paternity, given observed genotype and dispersal data. Also draw a sample of mating
    events based on dispersal only, to simulate what would be 'expected' if mating was
    random, and based only on distance.

    Parameters
    ==========
    data: `faps_data` class object
        Data about the population.
    model: dict
        Dictionary of starting model parameters. Keys should be a subset of
        ['missing', 'shape', 'scale', 'mixture', 'assortment'] and values 
        floats giving initial values for those parameters, within appropriate
        boundaries.
    ndraws: int, optional
        Number of Monte Carlo draws to perform for sibship clustering. Defaults
        to 1000.
    max_distance: float, int
        Maximum distance from the mother a candidate may be. Candidates further than
        this value will have their posterior probability of paternity set to zero.
        This is equivalent to setting a threshold prior on distance.
    """
    # Update data with parameter values
    data.update_covariate_probs(model = model, max_distance = max_distance)
    # We need to set dispersal probabilities for mothers to zero
    # Otherwise the mother would be drawn many times when drawing from the dispersal
    # kernel only, and cause many invalid partitions.
    jx = [np.where(x == np.array(data.candidates))[0][0] for x in data.mothers]
    for i,j in zip(range(len(data.paternity)), jx):
        data.covariates['dispersal'][i, j] = -np.inf

    # Cluster into sibships, if not already done.
    data.sibship_clustering(ndraws = 1000, use_covariates = True)
    # Draw a posterior data set for observed data, and expected under dispersal only
    obs = fp.posterior_mating(data.sibships, covariates_only=False, ndraws=ndraws)
    exp = fp.posterior_mating(data.sibships, covariates_only=True, ndraws=ndraws)

    return {
        'obs' : obs,
        'exp' : exp
    }

# def observed_sires(data, folder, burnin, ndraws = 1000, max_distance = np.inf):
#     """
#     List possible mating events with their probabilities, distance
#     between parents and flower colours of the parents, integrating out
#     uncertainty in sibship structure.

#     This concatenates MCMC chains in the directory `folder`, then runs 
#     `faps.summarise_sires` on the data set for each draw from the resulting
#     posterior distribution. The output is a concatenation of the siring
#     events data for each set of dispersal paraemeters, which you can 
#     parse out with the column `iter`.

#     Parameters
#     ----------
#     data: `faps_data` class object
#         Data about the population.
#     folder: str
#         Path to a directory containing one or more files containing 
#         output of MCMC runs. Files should have the suffix `.out`.
#     burnin: int
#         Integer number of rows in each MCMC input file to discard.
#     ndraws: int, optional
#         Number of Monte Carlo draws to perform for sibship clustering. Defaults
#         to 1000.
#     max_distance: float, int
#         Maximum distance from the mother a candidate may be. Candidates further than
#         this value will have their posterior probability of paternity set to zero.
#         This is equivalent to setting a threshold prior on distance.

#     Returns
#     -------
#     A data frame listing plausible mating events for each draw from the posterior
#     distribution.
#     * iter: step in the MCMC chain; this indexes draws from the posterior distribution
#     * scale, shape, missing: dispersal parameter values for the posterior draw in
#     question.
#      * mother, father: ID of parents for each putative mating event.
#     * prob: probability that mating between pairs occurred (proportion of valid
#     partitions in which the candidate appears)
#     * offspring: most likely number of offspring sired within the family (weighted
#     average over valid partitions)
#     """
#     # Import MCMC chains
#     chains = glob(folder+"/*out")
#     posterior = {}
#     for chain in chains:
#         k = os.path.basename(chain)
#         posterior[k] = pd.read_csv(chain, sep="\t").loc[lambda x: x.iter >= burnin]
#     posterior = pd.concat(posterior).reset_index()
#     # For this script we are only interested in the generalised-Gaussian part
#     # of the dispersal kernel, so set `mixture` to 1 for all rows.
#     posterior['mixture'] = 1

#     # empty list to store mating events for each iteration
#     siring = {} 
#     # Loop over steps in the MCMC chain and get siring events for each.
#     for i in tqdm(posterior.index):
#         model = posterior.loc[i]
#         # Update data with parameter values
#         data.update_covariate_probs(model = model, max_distance = max_distance)
#         # Cluster into sibships, if not already done.
#         data.sibship_clustering(ndraws = 1000, use_covariates = True)
#         # Call siring events
#         data.sires = fp.summarise_sires(data.sibships)
#         # Add details about the iteration to save for later.
#         data.sires['scale'] = model['scale']
#         data.sires['shape'] = model['shape']
#         data.sires['missing'] = model['missing']
#         data.sires['iter'] = model['iter']    
#         # Send to siring
#         siring[i] = data.sires.filter(items = ['iter', 'scale', 'shape', 'missing', 'mother', 'father', 'prob', 'offspring'])

#     return pd.concat(siring)

# def random_sires(data, model, ndraws = 1000, max_distance = np.inf, threshold = 0.001):
#     """
#     Calculate likelihoods of mating for arbitrary arrays of mating
#     probabilities.

#     `random_sires` returns probilities that maternal individuals mate with each
#     of a set of candidate males based on a matrix of mating probabilties. These
#     are returned in the same format as `summarise_sires()`, and can be
#     processed in the same way. However, while `summarise_sires()` draws finite
#     samples of canidates for each sibship partition separately (because it uses
#     genetic information, and this is also pertinent to the sibships),
#     `random_sires` gets a probability of mating for *every* mother-candidate pair
#     in the population. This is also why the output is much longer than 
#     `summarise_sires`.

#     Parameters
#     ----------
#     Parameters
#     ----------
#     data: `faps_data` class object
#         Data about the population.
#     model: dict
#         Dictionary of starting model parameters. Keys should be a subset of
#         ['missing', 'shape', 'scale', 'mixture', 'assortment'] and values 
#         floats giving initial values for those parameters, within appropriate
#         boundaries.
#     ndraws: int, optional
#         Number of Monte Carlo draws to perform for sibship clustering. Defaults
#         to 1000.
#     max_distance: float, int
#         Maximum distance from the mother a candidate may be. Candidates further than
#         this value will have their posterior probability of paternity set to zero.
#         This is equivalent to setting a threshold prior on distance.

#     Returns
#     -------
#     DataFrame giving mother (taken from the keys of the input dictionary),
#     fathers (inherited from each sibshipCluster object), and probabilties of
#     having sired at least one offspring.
#     """
#     # Update data with parameter values
#     data.update_covariate_probs(model = model, max_distance = max_distance)
#     # Get a single matrix of mating given all covariate probabilty matrices
#     probs = sum(data.covariates.values())
#     # set probabilities for selfed offspring to zero, or this will inflate
#     # apparent assortative mating
#     self_positions = [np.where([x == k for x in data.candidates])[0][0] for k in data.paternity.keys()]
#     for i,j in zip(range(60), self_positions):
#         probs[i, j] = - np.inf
#     probs = probs - fp.alogsumexp(probs, axis=1)[:, np.newaxis]

#     # Cluster into sibships, if not already done.
#     if not hasattr(data, 'sibships'):
#         data.sibship_clustering(ndraws = ndraws, use_covariates = True)

#     # Expected number of mating events for each maternal family.
#     n_sires = [x.mean_nfamilies() for x in data.sibships.values()]
#     n_sires = np.log(n_sires)
#     # Multiply mating probabilities for each candidate by the number of opportunities to mate
#     exp_liks = probs + n_sires[:, np.newaxis]
#     # If there are likelihoods above 1, set threshold to 1
#     # exp_liks[exp_liks > 0] = 0
#     # Make it a dictionary so we can iterate over keys later
#     exp_liks = {k: v for k,v in zip(data.sibships.keys(), exp_liks)}

#     # Create a table in the same format as summarise_sires() would generate.
#     output = []
#     for k,v in data.sibships.items():
#         this_df = pd.DataFrame({
#             'mother'   : k,
#             'father'   : v.candidates,
#             'log_prob' : exp_liks[k],
#             'prob'     : np.exp(exp_liks[k])
#         })
#         output = output + [this_df]
#     output = pd.concat(output)

#     return output
