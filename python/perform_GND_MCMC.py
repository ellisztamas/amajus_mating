import numpy as np
from faps import *
from utils import *
from sys import stdout
from time import time
    
def perform_GND_MCMC(distance_matrix, parameters, thin, filename):
    output = open('../output/output_{}.txt'.format(filename), 'w')
    output.write('iter\thours\tscale\tshape\tmix\tloglik\n')
    output.close()
    
    # If candidates are drawn at random, their probabilities are all uniform. This only needs doing once.
    prob_random = 1.0/ adults.size
    # initialise the model
    current_model = parameters[0] # the first set of parameters
    current_lik   = -10e12 # set initial likelihood to a very small number

    save_indices = np.arange(0, parameters.shape[0]*parameters.shape[1], thin)
    output = np.zeros([len(save_indices), parameters.shape[1]+1])
    counter_all = 0
    counter_out = 0
        
    for i in range(parameters.shape[0]):
        for p in [0,1,2]:
            new_model = np.copy(current_model)
            new_model[p] = parameters[i,p]
            
            # Probability of drawing each male under GND dispersal
            prob_drawn = dispersal_GND(x=distance_matrix,  a=new_model[0], b=new_model[1], c=new_model[2])
            # add dispersal information to paternityArrays
            df = [patlik[s].add_covariate(prob_drawn[s]) for s in range(len(patlik))]
            del df
            # Cluster into families and get likelihoods
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
                #output[counter_out] = np.append(current_model, current_lik)
                counter_out += 1
                # write iteration to disk
                with open('../output/output_{}.txt'.format(filename), 'a') as f:
                    f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        counter_all,
                        np.round((time() - t0) / 3600, 2),
                        np.round(current_model, 2)[0],
                        np.round(current_model, 2)[1],
                        np.round(current_model, 2)[2],
                        current_lik)
                            )
                    print("{}\t{}\t{}\t{}\t{}\t{}".format(
                        counter_all, np.round((time() - t0) / 3600, 2),
                        np.round(current_model, 2)[0],
                        np.round(current_model, 2)[1],
                        np.round(current_model, 2)[2],
                        current_lik)
                            )
                f.close()
                
            counter_all += 1