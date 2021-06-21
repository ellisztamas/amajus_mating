from warnings import warn
from faps.alogsumexp import alogsumexp
import numpy as np
import pandas as pd
import faps as fp

def observed_sires(data, model, ndraws = 1000, max_distance = np.inf):
    """
    List possible mating events with their probabilities, distance
    between parents and flower colours of the parents, integrating out
    uncertainty in sibship structure.

    This runs sibship clustering on each family, including covariate,
    information on each family, then calls `sires()` on each.

    Parameters
    ----------
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

    Returns
    -------
    A data frame listing plausible mating events, giving probability that mating
    between pairs occurred (proportion of valid partitions in which the candidate
    appears), most likely number of offspring sired within the family (weighted 
    average over valid partitions), the distance between mother and father, their
    flower colours, and whether colours match.
    """
    # Update data with parameter values
    data.update_covariate_probs(model = model, max_distance = max_distance)

    # Cluster into sibships, if not already done.
    if not hasattr(data, 'sibships'):
        data.sibship_clustering(ndraws = ndraws, use_covariates = True)
    # Call siring events
    data.sires = fp.summarise_sires(data.sibships)

    # Add data on distance between mother and father and flower colours.
    # Turn distance matrix into a data frame so we can use .loc on it.
    distance_df = pd.DataFrame({
        'mother'   : np.repeat(list(data.mothers), data.n_candidates),
        'father'   : np.tile(data.candidates, len(data.mothers)),
        'distance' : data.distances.flatten()
    })
    # Merge siring events with distances and phenotypes.
    output = (data.sires.merge(distance_df, how="left", on=['mother', 'father']).
    merge(data.flower_colours, how="left", left_on="mother", right_index=True).
    merge(data.flower_colours, how="left", left_on="father", right_index=True, suffixes = ['_mother', "_father"])
    )
    output.assign(match = output['simple_colour_mother'] == output['simple_colour_father'])

    return output

def random_sires(data, model, ndraws = 1000, max_distance = np.inf, threshold = 0.001):
    """
    Calculate likelihoods of mating for arbitrary arrays of mating
    probabilities.

    `random_sires` returns probilities that maternal individuals mate with each
    of a set of candidate males based on a matrix of mating probabilties. These
    are returned in the same format as `summarise_sires()`, and can be
    processed in the same way. However, while `summarise_sires()` draws finite
    samples of canidates for each sibship partition separately (because it uses
    genetic information, and this is also pertinent to the sibships),
    `random_sires` gets a probability of mating for *every* mother-candidate pair
    in the population. This is also why the output is much longer than 
    `summarise_sires`.

    Parameters
    ----------
    Parameters
    ----------
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

    Returns
    -------
    DataFrame giving mother (taken from the keys of the input dictionary),
    fathers (inherited from each sibshipCluster object), and probabilties of
    having sired at least one offspring.
    """
    # Update data with parameter values
    data.update_covariate_probs(model = model, max_distance = max_distance)
    # Get a single matrix of mating given all covariate probabilty matrices
    probs = sum(data.covariates.values())
    # set probabilities for selfed offspring to zero, or this will inflate
    # apparent assortative mating
    self_positions = [np.where([x == k for x in data.candidates])[0][0] for k in data.paternity.keys()]
    for i,j in zip(range(60), self_positions):
        probs[i, j] = - np.inf
    probs = probs - fp.alogsumexp(probs, axis=1)[:, np.newaxis]

    # Cluster into sibships, if not already done.
    if not hasattr(data, 'sibships'):
        data.sibship_clustering(ndraws = ndraws, use_covariates = True)

    # Expected number of mating events for each maternal family.
    n_sires = [x.mean_nfamilies() for x in data.sibships.values()]
    n_sires = np.log(n_sires)
    # Multiply mating probabilities for each candidate by the number of opportunities to mate
    exp_liks = probs + n_sires[:, np.newaxis]
    # If there are likelihoods above 1, set threshold to 1
    # exp_liks[exp_liks > 0] = 0
    # Make it a dictionary so we can iterate over keys later
    exp_liks = {k: v for k,v in zip(data.sibships.keys(), exp_liks)}

    # Create a table in the same format as summarise_sires() would generate.
    output = []
    for k,v in data.sibships.items():
        this_df = pd.DataFrame({
            'mother'   : k,
            'father'   : v.candidates,
            'log_prob' : exp_liks[k],
            'prob'     : np.exp(exp_liks[k])
        })
        output = output + [this_df]
    output = pd.concat(output)
    return output
