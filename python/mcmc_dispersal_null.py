import numpy as np
from faps import *
from utils import *
from time import time, strftime
from sys import stdout
from perform_GND_MCMC import *

filename = 'dispersal_null'
# Set up values to test
np.random.seed(37)
nreps = 33333
thin=100
starting_params = [5.2, 0.37, 0.0]

print('\nMCMC analysis for 3-parameter pollen dispersal in A. majus, starting with zero genetic information.')
print('Analysis begun {}.\n\n'.format(strftime("%Y-%m-%d %H:%M:%S")))
t0 = time()

# FAPS objects and distance matrices are generated in a separate script.
execfile('setup_FAPS_GPS.py')
print('Genotype and GPS data imported.')
# calculate the paternityArray objects
patlik = paternity_array(progeny, mothers, adults, mu = 0.0013, missing_parents= 0.1)
print('paternityArrays objects created after {} minutes.'.format(np.round((time() - t0) / 60, 2)))

# Draw values for a, b and mixture parameter
params = np.array([
    np.random.uniform(0.001,400, size=nreps), # Scale parameter for generalised normal distribution
    np.random.uniform(0.006, 3, size=nreps), # shape parameter for generalised normal
    np.random.uniform(0.001, 1, size=nreps)  # Mixture paramter for gen. normal vs. uniform.
])
params = params.T
params[0] = starting_params

# Run the Markov chain
print('MCMC set up. Beginning Markov chain...\n')
perform_GND_MCMC(distance_matrix, patlik, params, 2, filename)
print('\nMCMC completed {}.\n'.format(strftime("%Y-%m-%d %H:%M:%S", )))