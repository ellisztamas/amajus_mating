from warnings import warn
import numpy as np
import pandas as pd
from faps import sibshipCluster

def mating_probabilities(siring_probabilities, genotypes, geno_col, drop_na=False):
    """
    Calculate likelihoods of mating between each combination of genotypes.

    Parameters
    ----------
    siring_probabilities: DataFrame
        Dataframe with columns 'mother', 'father', 'prob', giving IDs for
        the parents, and a probability of mating having occured between them.
        You can get this by running summarise_sibships() over dictionary of
        sibshipArray clusters.
    genotypes: DataFrame
        Dataframe whose index gives the individual ID, and subsequent giving
        genotypes for one or more loci. The genotype column to use is set by
        geno_col
    geno_col: str
        Column in `genotypes` to use.
    drop_na: boolean, optional
        If True, NA genotypes are ignored.

    Returns
    -------
    Dataframe giving likelihoods of mating for each combination of genotypes.
    'prob' shows the probability that a mother of that genotype receives pollen
    from that genotype, normalised to sum to one.
    """
    # If this function makes it into the FAPS package, add additional check here!

    if geno_col not in genotypes.columns:
        raise ValueError("geno_col not found in the columns of genotypes.")

    # Dataframe of paternal gentoypes in the same order as in sires.
    dads_phens = genotypes.loc[siring_probabilities['father']]
    dads_phens = dads_phens.reset_index(drop=True)
    # Same for the mothers
    mums_phens = genotypes.loc[siring_probabilities['mother']]
    mums_phens = mums_phens.reset_index(drop=True)
    # Unique genotypes
    if drop_na:
        geno = genotypes[geno_col].dropna().unique()
    else:
        geno = genotypes[geno_col].unique()

    # Loop over each combination of genotypes and pull out probabilities of mating.
    output = []
    for maternal in geno:
        for paternal in geno:
            gx  = (mums_phens[geno_col] == maternal) & (dads_phens[geno_col] == paternal)
            gx = np.where(gx)[0]
            # Get mating likleihoods
            df = siring_probabilities.iloc[gx] # Subset of matings for this combination of genotypes
            df = df['prob'].sum()
            output = output + [[maternal, paternal, df]]
    # Join them together
    output = pd.DataFrame(output, columns = ['maternal', 'paternal', 'lik'])
    # Give probabilities of mating for each mother that sum to 1.
    xv = output.groupby("maternal").sum().to_numpy() # Sum likelihoods for maternal genotypes
    output['prob'] = output['lik'] / np.repeat(xv, len(geno)) # divide by sums

    return output



def random_sires(sibships, probs):
    """
    Calculate likelihoods of mating for arbitrary arrays of mating probabilities.

    `random_sires` returns probilities that maternal individuals mate with each
    of a set of candidate males based on a matrix of mating probabilties. These
    are returned in the same format as `summarise_sires()`, and can be processed
    in the same way.

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
    # Check dictionary of sibships.
    if not isinstance(sibships, dict):
        raise TypeError("`sibships` should be a dictionary of sibshipCluster onjects.")
    if not all([isinstance(x, sibshipCluster) for x in sibships.values()]):
        raise TypeError("Not all elements of `sibships` are sibshipCluster objects.")
    if isinstance(sibships, sibshipCluster):
        raise TypeError("`summarise_sires` is intended to work on a dictionary of sibshipCluster objects, but a single sibshipCluster object was supplied. In this case, call `sires()` directly on the onject, i.e. object.sires().")

    # Check proabbilties are log
    if (probs > 0).any():
        warn("Matrix of mating probabilities contains positive numbers. " \
            "These ought to be log probabilties, which are negative. " \
            "random_sires will run, but results are likely to be garbage."
        )
    # Check nrows match number of sibships
    if probs.shape[0] != len(sibships.keys()):
        raise ValueError(
            "Matrix of probabilities has {} rows but there are {} sibshipCluster objects".format(probs.shape[0],len(sibships.keys()))
        )
    # Check n columns matches number of candidates
    ncandidates = [len(x.candidates) for x in sibships.values()][0]
    if probs.shape[1] != ncandidates:
        raise ValueError(
            "Matrix of probabilities has {} columns but there are {} sibshipCluster objects".format(probs.shape[1], ncandidates)
        )

    # Expected number of mating events for each maternal family.
    n_sires = [x.mean_nfamilies() for x in sibships.values()]
    n_sires = np.log(n_sires)
    # Multiply mating probabilities for each candidate by the number of opportunities to mate
    exp_liks = probs + n_sires[:, np.newaxis]
    # Make it a dictionary so we can iterate over keys later
    exp_liks = {k: v for k,v in zip(sibships.keys(), exp_liks)}

    # Create a table in the same format as summarise_sires() would generate.
    output = []
    for k,v in sibships.items():
        this_df = pd.DataFrame({
            'mother'   : k,
            'father'   : v.candidates,
            'log_prob' : exp_liks[k],
            'prob'     : np.exp(exp_liks[k])
        })
        output = output + [this_df]
    output = pd.concat(output)

    return output
