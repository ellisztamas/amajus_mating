#!/usr/bin/env python3

"""
Tom Ellis, 27th May 2021

Script to run joint analysis of paternity, sibships and dispersal using
the same priors as result 003, but setting hard limits of 200, 300 or 400m
to dispersal distances.
"""
import numpy as np
import os
from scipy.stats import beta
from scipy.stats import gamma

from amajusmating import mcmc

# FAPS objects and distance matrices are generated in a separate script.
exec(open('003.scripts/setup_FAPS_GPS.py').read())

# INITIALISE THE MODEL
nreps = 4500 # Total number of iterations to run
thin  = 10 # How often to write samples.
# output_dir = "005.results/003_mcmc_restrict_kurtosis/output/"
output_dir = os.path.dirname(os.path.abspath(__file__))+'/output/'
os.makedirs(output_dir, exist_ok=True)

np.random.seed(46)

# PRIORS
priors = (lambda x : {
    'missing' : beta.pdf(x['missing'], a=3,   b=15),
    'mixture' : beta.pdf(x['mixture'], a=1.1, b=1.1),
    'shape'   : gamma.pdf(x['shape'],   a=10,  scale = 1/5),
    'scale'   : gamma.pdf(x['scale'],   a=6,   scale = 50)
})

# Proposed values are a Gaussian peturbation away from the previous values.
# This is controlled by the sigma of the gaussian, which is defined for each variable
proposal_sigma = {
    'missing' : 0.025,
    'shape'  : 0.05,
    'scale'  : 2,
    'mixture' : 0.025,
}

for i in [1,2,3,4]:
    mcmc.run_MCMC(
        data= am_data,
        initial_parameters = {
            'missing' : beta.rvs(a=3, b = 15),
            'mixture' : beta.rvs(a=1.1, b = 1.1),
            'shape'   : gamma.rvs(a=10,  scale = 1/5),
            'scale'   : gamma.rvs(a=6,  scale = 50)
        },
        proposal_sigma = proposal_sigma,
        priors = priors,
        thin=thin,
        nreps=nreps,
        output_dir = output_dir,
        chain_name = 'chain' + str(i) + "_200m",
        max_distance = 200
        )

for i in [1,2,3,4]:
    mcmc.run_MCMC(
        data= am_data,
        initial_parameters = {
            'missing' : beta.rvs(a=3, b = 15),
            'mixture' : beta.rvs(a=1.1, b = 1.1),
            'shape'   : gamma.rvs(a=10,  scale = 1/5),
            'scale'   : gamma.rvs(a=6,  scale = 50)
        },
        proposal_sigma = proposal_sigma,
        priors = priors,
        thin=thin,
        nreps=nreps,
        output_dir = output_dir,
        chain_name = 'chain' + str(i) + "_300m",
        max_distance = 300
        )

for i in [1,2,3,4]:
    mcmc.run_MCMC(
        data= am_data,
        initial_parameters = {
            'missing' : beta.rvs(a=3, b = 15),
            'mixture' : beta.rvs(a=1.1, b = 1.1),
            'shape'   : gamma.rvs(a=10,  scale = 1/5),
            'scale'   : gamma.rvs(a=6,  scale = 50)
        },
        proposal_sigma = proposal_sigma,
        priors = priors,
        thin=thin,
        nreps=nreps,
        output_dir = output_dir,
        chain_name = 'chain' + str(i) + "_400m",
        max_distance = 400
        )

for i in [1,2,3,4]:
    mcmc.run_MCMC(
        data= am_data,
        initial_parameters = {
            'missing' : beta.rvs(a=3, b = 15),
            'mixture' : beta.rvs(a=1.1, b = 1.1),
            'shape'   : gamma.rvs(a=10,  scale = 1/5),
            'scale'   : gamma.rvs(a=6,  scale = 50)
        },
        proposal_sigma = proposal_sigma,
        priors = priors,
        thin=thin,
        nreps=nreps,
        output_dir = output_dir,
        chain_name = 'chain' + str(i) + "_500m",
        max_distance = 500
        )

for i in [1,2,3,4]:
    mcmc.run_MCMC(
        data= am_data,
        initial_parameters = {
            'missing' : beta.rvs(a=3, b = 15),
            'mixture' : beta.rvs(a=1.1, b = 1.1),
            'shape'   : gamma.rvs(a=10,  scale = 1/5),
            'scale'   : gamma.rvs(a=6,  scale = 50)
        },
        proposal_sigma = proposal_sigma,
        priors = priors,
        thin=thin,
        nreps=nreps,
        output_dir = output_dir,
        chain_name = 'chain' + str(i) + "_750m",
        max_distance = 750
        )

for i in [1,2,3,4]:
    mcmc.run_MCMC(
        data= am_data,
        initial_parameters = {
            'missing' : beta.rvs(a=3, b = 15),
            'mixture' : beta.rvs(a=1.1, b = 1.1),
            'shape'   : gamma.rvs(a=10,  scale = 1/5),
            'scale'   : gamma.rvs(a=6,  scale = 50)
        },
        proposal_sigma = proposal_sigma,
        priors = priors,
        thin=thin,
        nreps=nreps,
        output_dir = output_dir,
        chain_name = 'chain' + str(i) + "_1000m",
        max_distance = 1000
        )