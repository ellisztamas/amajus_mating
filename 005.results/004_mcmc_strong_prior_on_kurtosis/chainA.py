#!/usr/bin/env python3

"""
Tom Ellis, 26th January 2021

Script to run joint analysis of paternity, sibships and dispersal using
priors that are very skeptical about kurtosis (0.4% of the prior mass is 
below 2).
"""

# Import external packages
import numpy as np
import os
from scipy.stats import beta
from scipy.stats import gamma as gma

# FAPS objects and distance matrices are generated in a separate script.
exec(open('003.scripts/setup_FAPS_GPS.py').read())

# INITIALISE THE MODEL
nreps = 10500 # Total number of iterations to run
thin  = 10 # How often to write samples.
np.random.seed(1246)
max_distance = np.inf

# Dictionary listing starting values.
initial_model = {
    'missing' : 0.15, # proportion missing fathers
    'shape'  : 1,
    'scale'  : 10,
    'mixture' : 0.8
}

# Proposed values are a Gaussian peturbation away from the previous values.
# This is controlled by the sigma of the gaussian, which is defined for each variable
proposal_sigma = {
    'missing' : 0.025,
    'shape'  : 0.05,
    'scale'  : 2,
    'mixture' : 0.025,
}

# PRIORS
priors = (lambda x : {
    'missing' : beta.pdf(x['missing'], a=3,   b=15),
    'mixture' : beta.pdf(x['mixture'], a=1.1, b=1.1),
    'shape'   : gma.pdf(x['shape'],   a=20,  scale = 1/5),
    'scale'   : gma.pdf(x['scale'],   a=6,   scale = 50)
})

prior_probs = {
            'missing'   : beta.pdf(new_model['missing'],a=3,   b=15),
            'lambda'    : beta.pdf(new_model['lambda'], a=1.1, b=1.1),
            'shape'     : gma.pdf(new_model['shape'],   a=2, scale = 1/2),
            # 'mean_dist' : gma.pdf(new_model['mean_dist'],      a=2,   scale= 200),
            'scale'     : gma.pdf(new_model['shape'],   a=2, scale = 3),
        }

# Run the MCMC
am_data.run_MCMC(
    initial_parameters = initial_model,
    proposal_sigma = proposal_sigma,
    priors = priors,
    thin=thin,
    nreps=nreps,
    output_dir= os.path.dirname(os.path.abspath(__file__))+'/',
    chain_name = os.path.splitext(os.path.basename(__file__))[0],
    max_distance = max_distance
)