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