#!/usr/bin/env python3

"""
Tom Ellis, 26th January 2021

Script to run joint analysis of paternity, sibships and dispersal using
priors that are very skeptical about kurtosis (0.4% of the prior mass is 
below 2).
"""

# Import external packages
import numpy as np
from faps import *
import os

# import local modules.
from run_MCMC import run_MCMC


# FAPS objects and distance matrices are generated in a separate script.
exec(open('scripts/setup_FAPS_GPS.py').read())

# INITIALISE THE MODEL
nreps = 10500 # Total number of iterations to run
thin  = 10 # How often to write samples.
np.random.seed(1246)
max_distance = np.inf

# Dictionary listing starting values.
initial_model = {
    'loglik' : -10e12, # set initial likelihood to a very small number
    'missing' : 0.15, # proportion missing fathers
    # Dispersal parameters
    'shape'  : 1,
    'scale'  : 10,
    'lambda' : 0.8
}

# Proposed values are a Gaussian peturbation away from the previous values.
# This is controlled by the sigma of the gaussian, which is defined for each variable
proposal_sigma = {
    'missing' : 0.025,
    'shape'  : 0.05,
    'scale'  : 2,
    'lambda' : 0.025,
}

# Run the MCMC
run_MCMC(
    patlik,
    distance_matrix,
    initial_model,
    proposal_sigma,
    thin=thin,
    nreps=nreps,
    output_dir= os.path.dirname(os.path.abspath(__file__))+'/',
    chain_name = os.path.splitext(os.path.basename(__file__))[0],
    max_distance = max_distance
)
