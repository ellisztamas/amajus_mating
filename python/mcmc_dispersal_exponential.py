import numpy as np
from faps import *
import pandas as pd
from utils import *
from time import time, strftime
from sys import stdout

filename = 'dispersal_exponential'
# Set up values to test
np.random.seed(456)
nreps = 33333
thin=100
starting_params = [177, 1.0, 0.8]

logfile = open('../output/log_{}.txt'.format(filename), 'w')
print('\nMCMC analysis for 3-parameter pollen dispersal in A. majus, starting with an exponential distribution.')
logfile.write('MCMC analysis for 3-parameter pollen dispersal in A. majus, starting with an exponential distribution.\n')
print('Analysis begun {}.\n\n'.format(strftime("%Y-%m-%d %H:%M:%S")))
logfile.write('Analysis begun {}.\n\n'.format(strftime("%Y-%m-%d %H:%M:%S")))
t0 = time()

# FAPS objects and distance matrices are generated in a separate script.
execfile('setup_FAPS_GPS.py')
print('Genotype and GPS data imported.')
# calculate the paternityArray objects
patlik = paternity_array(progeny, mothers, adults, mu = 0.0013, missing_parents= 0.1)
print('paternityArrays objects created after {} minutes.'.format(np.round((time() - t0) / 60, 2)))

# If candidates are drawn at random, their probabilities are all uniform
prob_random = 1.0/ adults.size

# Draw values for a, b and mixture parameter
params = np.array([
    np.random.uniform(0.001,400, size=nreps), # Scale parameter for generalised normal distribution
    np.random.uniform(0.006, 3, size=nreps), # shape parameter for generalised normal
    np.random.uniform(0.001, 1, size=nreps)  # Mixture paramter for gen. normal vs. uniform.
])
params = params.T
# initialise the model
current_model = params[0] # the first set of parameters
current_lik   = -10e12 # set initial likelihood to a very small number

save_indices = np.arange(0, params.shape[0]*params.shape[1], thin)
output = np.zeros([len(save_indices), params.shape[1]+1])
counter_all = 0
counter_out = 0

print('MCMC set up. Beginning Markov chain...\n')
print('iter\thours\tscale\tshape\tmix\tloglik')
for i in range(nreps): 
    for p in [0,1,2]:
        new_model = np.copy(current_model)
        new_model[p] = params[i,p]
        
        prob_drawn = dispersal_GND(x=distance_matrix,  a=new_model[0], b=new_model[1], c=new_model[2])

        # add dispersal information to paternityArrays
        df = [patlik[s].add_covariate(prob_drawn[s]) for s in range(len(patlik))]
        del df

        sc = sibship_clustering(patlik, ndraws=100, use_covariates=True)
        new_lik = np.array([alogsumexp(s.lik_partitions) for s in sc]).sum()

        if new_lik > current_lik:
            accept = True
        elif bool(np.random.binomial(1, np.exp(new_lik - current_lik))):
            accept = True
        else:
            accept = False

        if accept:
            current_lik   = new_lik
            current_model = new_model
        if counter_all in save_indices:
            output[counter_out] = np.append(current_model, current_lik)
            counter_out += 1
            
            print "{}\t{}\t{}\t{}\t{}\t{}".format(
                counter_all, np.round((time() - t0) / 3600, 2),
                np.round(current_model, 2)[0],
                np.round(current_model, 2)[1],
                np.round(current_model, 2)[2],
                current_lik)
            stdout.flush()
        counter_all += 1
        
logfile.write('\nMCMC completed {}.\n'.format(strftime("%Y-%m-%d %H:%M:%S", )))
print('\nMCMC completed {}.\n'.format(strftime("%Y-%m-%d %H:%M:%S", )))

out_name = "../output/mcmc_{}.csv".format(filename)
output = pd.DataFrame(output, index = save_indices, columns=['scale', 'shape', 'mix', 'loglik'])
output.round(3).to_csv(out_name)

print('Output saved as {}.'.format(out_name))
logfile.write('Output saved as {}.'.format(out_name))
logfile.close()
