import numpy as np

def setup_phenotypes(focal, cand):
    """
    Set up arrays for inferences about mating between phenotype classes.
    
    Parameters
    ----------
    focal: array-like
        Array of phenotypes for the focal individuals, usually mothers
    cand: array-like
        Array of phenotypes for the candidate fathers.
    
    Returns
    -------
    A dictionary containing:
    'unique': Vector of possible phenotype combinations.
    'observed': Array of observed mating combinations, with a row for every
        plant in focal, and a column for every plant in cand.
    'transitions': Matrix of mating probabilities between each pair of
        phenotypes.
        
    """
    # Matrix giving the phenotype combinations for all possible matings
    mating_combs = [i + ' x ' + cand for i in focal]
    mating_combs = np.array(mating_combs)

    # All possible pairs of genotypes
    possible_phenotypes = cand.unique()
    possible_mating = np.array([i + ' x ' + j for i in possible_phenotypes for j in possible_phenotypes])

    # Set up a transition matrix.
    assort_array = np.zeros([len(possible_phenotypes), len(possible_phenotypes)])
    assort_array = assort_array + 1.0/len(possible_phenotypes)
    
    return {
        'unique' : possible_mating,
        'observed' : mating_combs,
        'transitions'  : assort_array
    }

def peturb_transitions(trans_mat, sigma):
    """
    Add gaussian noise to a square transition matrix.
    
    Parameters
    ----------
    trans_mat: array
        Matrix of matrix probabilities from `setup_phenotypes`.
    sigma: float
        Standard deviation of the Gaussian function to peturb transition
        probabilities
        
    Returns
    -------
    An array of the same shape as trans_mat.
    """
    # Create a matrix of standard normal deviates
    peturb = np.random.normal(loc=0, scale=sigma, size= len(trans_mat.flatten()))
    peturb = peturb.reshape(trans_mat.shape)
    # Peturb the transition matrix by those deviates
    trans_mat = trans_mat + peturb
    trans_mat = trans_mat / trans_mat.sum(1, keepdims=True) # normalise rows.
    
    return trans_mat

def assortment_probs(phenotypes, sigma):
    """
    Calculate the log probabilities of mating between pairs of individuals
    based on a matrix of transition probabilities.
    
    Parameters
    ----------
    phenotypes: dict
        Information about combinations of phenotypes from `setup_phenotypes`.
    sigma: float
        Standard deviation of the Gaussian function to peturb transition
        probabilities
    
    Returns
    -------
    Array of log probabilities of mating with a row for each mother and a
    column for each candidate.
    """
    
    phenotypes['transitions'] = peturb_transitions(phenotypes['transitions'], sigma)
    
    prob_drawn = np.zeros(phenotypes['observed'].shape)
    for (i,j) in zip(phenotypes['unique'], phenotypes['transitions'].flatten()):
        ix = phenotypes['observed'] == i
        prob_drawn[ix] = j
    
    prob_drawn = prob_drawn / prob_drawn.sum(1, keepdims=True)
    prob_drawn = np.log(prob_drawn)
    
    return prob_drawn
