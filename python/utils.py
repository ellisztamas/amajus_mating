import numpy as np
from scipy.special import gamma

def generalised_normal_PDF(x, a, b, gamma_b=None):
    """
    Calculate the PDF of the generalised normal distribution.
    
    Parameters
    ----------
    x: vector
        Vector of deviates from the mean.
    a: float
        Scale parameter.
    b: float
        Shape parameter
    gamma_b: float, optional
        To speed up calculations, values for Euler's gamma for 1/b
        can be calculated ahead of time and included as a vector.
    """
    xv = np.copy(x)
    if gamma_b:
        return (b/(2 * a * gamma_b ))      * np.exp(-(xv/a)**b)
    else:
        return (b/(2 * a * gamma(1.0/b) )) * np.exp(-(xv/a)**b)

def dispersal_GND(x, a, b, c):
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
    a: float
        Scale parameter.
    b: float
        Shape parameter
    c: float between 0 and 1.
        The proportion of probability mass assigned to the
        generalised normal function.
    """    
    prob_GND = generalised_normal_PDF(x, a, b)
    prob_GND = prob_GND / prob_GND.sum(axis=1)[:, np.newaxis]

    prob_drawn = (prob_GND * c) + ((1-c) / x.shape[1])
    prob_drawn = np.log(prob_drawn)

    return prob_drawn


# Deprecated code to calculate sampling density.
# If fathers are drawn from the population at random, the probability ought to be proportional to sampling density 
#bins = np.arange(0, 4000, 100)
#bin_labels = np.array([np.digitize(dm, bins) for dm in distance_matrix])

#bin_frequencies = np.array([(bin_labels == x).sum(1) for x in range(1, len(bins))]).T
#sampling_density = np.array([x[y-1] for x,y in zip(bin_frequencies, bin_labels)])
#sampling_density = 1.0 * sampling_density / sampling_density.sum(1)[:, np.newaxis]
#sampling_density = np.log(sampling_density)