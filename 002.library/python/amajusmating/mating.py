from warnings import warn
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
    # Update missing fathers
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

def random_sires(data, probs):
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
    sibships: dictionary
        Dictionary of sibshipCluster objects for multiple maternal families.
        Keys should normally indicate maternal family name.
    probs: array
        Array of (log) probabilites of mating between mothers and each candidate
        father. This should have a row for each mother and a column for each
        candidate. Rows should sum to 1.

    Returns
    -------
    DataFrame giving mother (taken from the keys of the input dictionary),
    fathers (inherited from each sibshipCluster object), and probabilties of
    having sired at least one offspring.
    """
    # Cluster into sibships, if not already done.
    if not hasattr(data, 'sibships'):
        data.sibship_clustering(ndraws = ndraws, use_covariates = True)

    # Check dictionary of sibships.
    if not all([isinstance(x, sibshipCluster) for x in data.sibships.values()]):
        raise TypeError("Not all elements of `data.sibships` are sibshipCluster objects.")
    if isinstance(data.sibships, sibshipCluster):
        raise TypeError("`summarise_sires` is intended to work on a dictionary of sibshipCluster objects, but a single sibshipCluster object was supplied. In this case, call `sires()` directly on the onject, i.e. object.sires().")

    # Check proabbilties are log
    if (probs > 0).any():
        warn("Matrix of mating probabilities contains positive numbers. " \
            "These ought to be log probabilties, which are negative. " \
            "random_sires will run, but results are likely to be garbage."
        )
    # Check nrows match number of data.sibships
    if probs.shape[0] != len(data.sibships.keys()):
        raise ValueError(
            "Matrix of probabilities has {} rows but there are {} sibshipCluster objects".format(probs.shape[0],len(data.sibships.keys()))
        )
    # Check n columns matches number of candidates
    ncandidates = [len(x.candidates) for x in data.sibships.values()][0]
    if probs.shape[1] != ncandidates:
        raise ValueError(
            "Matrix of probabilities has {} columns but there are {} sibshipCluster objects".format(probs.shape[1], ncandidates)
        )

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
