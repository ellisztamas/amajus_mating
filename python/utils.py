from time import time
import numpy as np
from scipy.special import gamma

def generalised_normal_PDF(x, scale, shape, gamma_b=None):
    """
    Calculate the PDF of the generalised normal distribution.

    Parameters
    ----------
    x: vector
        Vector of deviates from the mean.
    scale: float
        Scale parameter.
    shape: float
        Shape parameter
    gamma_b: float, optional
        To speed up calculations, values for Euler's gamma for 1/shape
        can be calculated ahead of time and included as a vector.
    """
    xv = np.copy(x)
    if gamma_b:
        return (shape/(2 * scale * gamma_b ))      * np.exp(-(xv/scale)**shape)
    else:
        return (shape/(2 * scale * gamma(1.0/shape) )) * np.exp(-(xv/scale)**shape)

def dispersal_GND(x, scale, shape, w):
    """
    Calculate a probability that each candidate is a sire
    assuming assuming he is either drawn at random form the
    population, or from a generalised normal function of his
    distance from each mother. The relative contribution of the
    two distributions is controlled by mixture parameter c.

    Parameters
    ----------
    x: vector
        Vector of deviates from the mean.
    scale: float
        Scale parameter.
    shape: float
        Shape parameter
    w: float between 0 and 1.
        The proportion of probability mass assigned to the
        generalised normal function.
    """
    prob_GND = generalised_normal_PDF(x, scale, shape)
    prob_GND = prob_GND / prob_GND.sum(axis=1)[:, np.newaxis]

    prob_drawn = (prob_GND * w) + ((1-w) / x.shape[1])
    prob_drawn = prob_drawn / prob_drawn.sum(1, keepdims=True)
    prob_drawn = np.log(prob_drawn)

    return prob_drawn

def stdev_GND(scale, shape):
    """
    Calculate the 2nd moment of the genealised normal distribution.

    Parameters
    ----------
    scale: float
        Scale parameter.
    shape: float
        Shape parameter

    Returns
    -------
    Float.
    """
    return (scale * gamma(2.0/shape)) / (1.0/shape)


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

    # This is a very hacky way to ensure that values are all positive
    new_model = {k: np.sqrt(v**2) for (k,v) in new_model.items()}

    return new_model

def assortment_probs(focal, candidates, params):
    """
    Calculate the log probabilities of mating between pairs of individuals
    given a vector of assortment probabilities.

    It is assumed that mothers of a certain phenotype mate with individuals
    of the same phenotype with some assortment parameter, a, and with
    individuals of any other phenotype with probability 1-a.

    Parameters
    ----------
    focal: vector
        Vector of phenotypes of the mothers. Cannot contain missing data.
    candidates: vector
        Vector of phenotypes of the candidate fathers. Can contain missing
        data.
    params: dict
        Dictionary containing assortment parameters for each genotype.

    Returns
    -------
    Array of log probabilities of mating with a row for each mother and a
    column for each candidate.
    """
    # unique list of all phenotypes in the focal phenotypes.
    phens = focal.unique().tolist()

    # Identify candidates with missing data and randomly assign another phenotype.
    NAs = np.where(candidates.isna())[0]
    replacements = np.random.choice(phens, size=len(NAs))
    candidates[NAs] = replacements

    # Dictionary of assortment parameters for each phenotype.
    # This is necessary if the input dictionary contains other parameters
    # not related to assortment.
    params = {k:v for k,v in params.items() if k in phens}

    if all([x in params for x in phens]):
        ass_matrix = np.zeros([len(focal), len(candidates)])
        for k in phens:
            # indices of mothers and candidates of phenotype k
            ix = np.where(focal      == k)[0][:, np.newaxis]
            jx = np.where(candidates == k)[0][   np.newaxis]
            # Insert assortment parameters into ass_matrix
            ass_matrix[ix] += (1-params[k]) # set all values to 1-a
            ass_matrix[ix, jx] = params[k] # For candidates that match the mother, set values to a.
        # Normalise and log transform
        ass_matrix = ass_matrix / ass_matrix.sum(1, keepdims=True)
        ass_matrix = np.log(ass_matrix)

        return ass_matrix
    else:
        raise ValueError("Not all phenotypes values present in the focal phenotypes have an assortment parameter.")

# def assortment_probs(phenotypes, a):
#   """
#     Calculate the log probabilities of mating between pairs of individuals
#     based on an assortment parameter between two discrete phenotypes.
#
#     Parameters
#     ----------
#     phenotypes: dict
#         Information about combinations of phenotypes from `setup_phenotypes`.
#     a: float
#         Assortment parameter between 0 and 1. The probability that a genotype
#         mates with the same genotype.
#
#     Returns
#     -------
#     Array of log probabilities of mating with a row for each mother and a
#     column for each candidate.
#     """
# # For candidates with missing data, draw True/False matches at random.
# replacements = np.random.binomial(1, 0.5, size=len(phenotypes['missing'])) == 1
# # Insert these replacements into the array of phenotypesotype matches.
# for i, j in zip(phenotypes['missing'], replacements):
#   phenotypes['matches'][:, i] = j
#
# # Create a matrix of mating probabilities based on assortment parameters and phenotypesotype matches.
# ass_matrix = np.zeros(phenotypes['matches'].shape)
# ass_matrix = ass_matrix + (1-a) # set the default to the value for non-matches.
# ass_matrix[phenotypes['matches']] = a # For pairs whose phenotypesotypes match, add these in.
#
# ass_matrix = ass_matrix / ass_matrix.sum(1, keepdims=True)
# ass_matrix = np.log(ass_matrix)
#
# return ass_matrix

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

def setup_output(path, model):
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
    out = 'iter\tacc\thours'
    for k in sorted(model): out = out + '\t' + k
    out = out + '\n'

    output = open('{}'.format(path), 'w')
    output.write(out)
    output.close()

def write_output(path, i, acceptance, model, decimals = 3, time0=None):
    """
    Write the output of a model from the current iteration to a
    new line in the output file.

    Parameters
    ----------
    path: str
        Path and filename where the output file is stored.
    i: int
        Iteration index.
    acceptance: float
        Mean Metropolis-Hastings acceptance rate.
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
    out = str(i) + '\t' + str(round(acceptance, decimals)) + '\t' + str(t)
    for k in sorted(model.keys()):
            out = out + '\t' + str(round(model[k], decimals))
    # write iteration to disk
    with open('{}'.format(path), 'a') as f:
        f.write(out + '\n')
    f.close()
