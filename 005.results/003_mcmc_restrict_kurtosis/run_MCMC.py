import numpy as np
import sys
import os

# import local modules.
file_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.dirname(os.path.dirname(file_dir))
# Add library folder to path
sys.path.append(os.path.join(project_dir, 'library'))
from dispersal import *
from mcmc import *

def run_MCMC(paternity_arrays, distance_matrix, initial_parameters, proposal_sigma, nreps, output_dir, chain_name, thin=1, max_distance = np.inf):
    """
    A wrapper function to run Metropolis-Hastings MCMC for paternity and dispersal.


    Parameters
    ----------
    paternity_arrays: dict
        Dictionary of paternityArray objects for each full sibling family.
    distance_matrix: array
        Array of distances between the mothers and all candidate fathers.
    initial_parameters: dict
        Dictionary of starting values for the parameters, with the names as keys and 
        values as floats for the parameter values. Usually this would include the
        'shape' and 'scale' parameters of for the generalised normal distribution,
        the proportion of 'missing' fathers and mixture parameter 'lambda'.
    proposal_sigma: dict
        Dictionary giving the standard deviation of the normal distribution by which
        to peturb parameters at each iteration. Should have the same keys as 
        `initial_parameters`.
    nreps: int
        Number of iterations to run
    output_dir: str
        Directory to save the output.    
     chain_name: str
        Name for the output file, without a suffix.
    thin: int
        Optional thinning argument. If >1, every one in `thin` samples will be written
        to disk.
    max_distance: float, int
        Maximum distance from the mother a candidate may be. Candidates further than
        this value will have their posterior probability of paternity set to zero.
        This is equivalent to setting a threshold prior on distance.
    
    """
    current_model = initial_parameters

    # Set up the datafiles
    out_file = output_dir + chain_name + ".txt"
    setup_output(out_file, current_model.keys())
    # Set up log file
    log_file = open(output_dir + chain_name + ".log", 'w')
    log_file.write('Metropolis-Hasting analysis mating in A. majus begun {}.'.format(strftime("%Y-%m-%d %H:%M:%S")))
    log_file.write('Initial model:\n')
    pprint(current_model, log_file)
    log_file.write('\nGaussian noise is applied at each iteration with standard deviations:\n')
    pprint(proposal_sigma, log_file)
    log_file.write("\nPerforming a total of {} steps, thinning every {} iteration. Output will be saved to:\n{}".format(nreps, thin, out_file))
    log_file.write('\nAnalysis begun {}.\n\n'.format(strftime("%Y-%m-%d %H:%M:%S")))

    #  IMPORT AND FORMAT DATA.
    t0 = time()

    # RUN THE MCMC.
    log_file.write('MCMC set up. Beginning Markov chain...\n')
    log_file.close()
    t0 = time()
    track_accept = 0 # to record mean MH-acceptance rate
    for i in tqdm(range(nreps)):
        # UPDATE PARAMETERS
        new_model = update_parameters(current_model, proposal_sigma)
        # # Update distance travelled using new parameters
        # new_model['mean_dist'] = stdev_GND(scale = new_model['scale'],
        #                                    shape = new_model['shape'])

        # LOG PROBABILITIES OF PATERNITY FOR EACH PARAMETER
        # Update proportion of missing fathers
        for p in paternity_arrays.keys():
            paternity_arrays[p].missing_parents = new_model['missing']
        # Update dispersal probabilities
        if new_model['lambda'] > 1.0: new_model['lambda'] = 1.0
        # Probability of drawing each male under GND dispersal
        prob_drawn = dispersal_GND(
            x     = distance_matrix,
            scale = new_model['scale'],
            shape = new_model['shape'],
            w     = new_model['lambda'])
        # Identify candidates who are further than the threshold distance
        # and set their log likelihoods to negative infinity
        ix = distance_matrix > max_distance
        prob_drawn[ix] = -np.inf
        # Incorporate into paternity_arrays
        for (p,s) in zip(paternity_arrays.keys(), prob_drawn):
            paternity_arrays[p].add_covariate(s)

        # INFER FAMILIES
        # Cluster into families and get likelihoods
        sc = fp.sibship_clustering(paternity_arrays, ndraws=100, use_covariates=True)
        new_model['loglik'] = np.array([fp.alogsumexp(s.lik_partitions) for s in sc.values()]).sum()

        # PRIORS
        # Incorporate prior probabilities.
        prior_probs = {
            'missing'   : beta.pdf(new_model['missing'],a=3,   b=15),
            'lambda'    : beta.pdf(new_model['lambda'], a=1.1, b=1.1),
            'shape'     : gma.pdf(new_model['shape'],   a=10, scale = 1/5),
            # 'mean_dist' : gma.pdf(new_model['mean_dist'],      a=2,   scale= 200),
            'scale'     : gma.pdf(new_model['shape'],   a=6, scale = 50),
        }
        # Log prior probabilities
        prior_probs = {k: np.log(v) for k,v in prior_probs.items()}
        # Sum log probabilities and incorporate
        new_model['loglik'] += np.array(list(prior_probs.values())).sum()

        # Decide whether to accept the new model.
        accept = mh_ratio(
            current = current_model['loglik'],
            new     = new_model['loglik']
        )
        if accept:
            current_model = new_model
            track_accept += 1

        mean_acceptance = float(track_accept) / (i+1)
        # write iteration to disk
        if(i in np.arange(start = 0, stop = nreps, step=thin)):
            write_output(out_file, i, mean_acceptance, current_model, decimals=3, time0=t0)

    log_file = open(output_dir + chain_name + ".log", 'a')
    log_file.write('\nMCMC completed {}.\n'.format(strftime("%Y-%m-%d %H:%M:%S")))
    log_file.close()
