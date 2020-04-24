import numpy as np
import faps as fp
from time import time, strftime
from tqdm import tqdm
from pprint import pprint
from utils_dispersal import *


from scipy.stats import beta
from scipy.stats import gamma as gma
from scipy.special import gamma

def update_parameters(current_model, sigma):
    """
    Apply Gaussian peturbation to a dictionary of model parameters.

    Parameters
    ----------
    current_model: dict
        Dictionary of current parameter values
    sigma: dict
        Dictionary of standard deviations from which to draw normal
        deviates to apply to each value in current_model. At least
        one key must intersect with the keys of current_model.
    """
    new_model = dict(current_model)
    # random deviates for each parameter
    dev = {k: np.random.normal(0, v) for k, v in sigma.items()}
    # apply deviates to each parameter
    for k in sigma.keys():
        new_model[k] += dev[k]

    # This is a very hacky way to ensure that values are all positive
    new_model = {k: np.sqrt(v**2) for (k,v) in new_model.items()}

    return new_model

def assortment_probs(focal, candidates, params):
    """
    Calculate the log probabilities of mating between pairs of individuals
    given a vector of assortment probabilities.

    It is assumed that mothers of a certain phenotype mate with individuals
    of the same phenotype with some assortment parameter, a, and with
    individuals of any other phenotype with probability 1-a.

    Parameters
    ----------
    focal: vector
        Vector of phenotypes of the mothers. Cannot contain missing data.
    candidates: vector
        Vector of phenotypes of the candidate fathers. Can contain missing
        data.
    params: dict
        Dictionary containing assortment parameters for each genotype.

    Returns
    -------
    Array of log probabilities of mating with a row for each mother and a
    column for each candidate.
    """
    # unique list of all phenotypes in the focal phenotypes.
    phens = focal.unique().tolist()

    # Identify candidates with missing data and randomly assign another phenotype.
    NAs = np.where(candidates.isna())[0]
    replacements = np.random.choice(phens, size=len(NAs))
    candidates[NAs] = replacements

    # Dictionary of assortment parameters for each phenotype.
    # This is necessary if the input dictionary contains other parameters
    # not related to assortment.
    params = {k:v for k,v in params.items() if k in phens}

    if all([x in params for x in phens]):
        ass_matrix = np.zeros([len(focal), len(candidates)])
        for k in phens:
            # indices of mothers and candidates of phenotype k
            ix = np.where(focal      == k)[0][:, np.newaxis]
            jx = np.where(candidates == k)[0][   np.newaxis]
            # Insert assortment parameters into ass_matrix
            ass_matrix[ix] += (1-params[k]) # set all values to 1-a
            ass_matrix[ix, jx] = params[k] # For candidates that match the mother, set values to a.
        # Normalise and log transform
        ass_matrix = ass_matrix / ass_matrix.sum(1, keepdims=True)
        ass_matrix = np.log(ass_matrix)

        return ass_matrix
    else:
        raise ValueError("Not all phenotypes values present in the focal phenotypes have an assortment parameter.")

def mh_ratio(current, new):
    """
    Decide whether to accept a proposed model over a previous model
    based on their respective (log) probabilities and the
    Metropolis-Hastings ratio.

    Parameters
    ----------
    current: float
        Log probability of the old model.
    new: float
        Log probability of the proposed model.

    Returns
    -------
    Logical. Returns True if the new model is accepted.
    """
    if new > current:
        accept = True
    elif bool(np.random.binomial(1, np.exp(new - current))):
        accept = True
    else:
        accept = False

    return accept

def setup_output(path, model):
    """
    Set up a blank text file to store model parameters

    Parameters
    ----------
    path: str
        Path and filename where dat should be stored. The suffix
        '.txt' will be automatically appended.
    model: list
        List of parameter names

    Returns
    -------
    A text file is created at the path specified with column headers
    only, showing iteration and time taken at each iteration,
    followed by columns for each parameter in the list model.
    """
    out = 'iter\tacc\thours'
    for k in sorted(model): out = out + '\t' + k
    out = out + '\n'

    output = open('{}'.format(path), 'w')
    output.write(out)
    output.close()

def write_output(path, i, acceptance, model, decimals = 3, time0=None):
    """
    Write the output of a model from the current iteration to a
    new line in the output file.

    Parameters
    ----------
    path: str
        Path and filename where the output file is stored.
    i: int
        Iteration index.
    acceptance: float
        Mean Metropolis-Hastings acceptance rate.
    model: dict
        Dictionary of model values with the same parameters
        as the target file.
    decimals: float
        Number of figures to report after the decimal point.
    time0: float
        Optional time from the start of the analysis. This
        will be compared with te current time, and the
        difference reported as the time taken to reach this
        iteration.

    Returns
    -------
    A line is added to the target file with the model values.
    """
    # Format time data
    if time0 is None:
        t = 'NA'
    else:
        t = np.round((time() - time0) / 3600, decimals)

    # Prepare a string of output data to be exported.
    out = str(i) + '\t' + str(round(acceptance, decimals)) + '\t' + str(t)
    for k in sorted(model.keys()):
            out = out + '\t' + str(round(model[k], decimals))
    # write iteration to disk
    with open('{}'.format(path), 'a') as f:
        f.write(out + '\n')
    f.close()

def run_MCMC(paternity_arrays, distance_matrix, initial_parameters, proposal_sigma, nreps, output_dir, chain_name, thin=1, max_distance = np.inf):
    """
    A wrapper function to run Metropolis-Hastings MCMC for paternity and dispersal.


    Parameters
    ----------
    paternity_arrays: dict
        Dictionary of paternityArray objects for each full sibling family.


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
    # FAPS objects and distance matrices are generated in a separate script.
    exec(open('setup_FAPS_GPS.py').read())
    # # Import flower colour data for the population
    # ros_sulf = pd.read_csv('../data_processed/rosea_sulfurea.csv', index_col='id')
    # ros_sulf = ros_sulf.loc[adults.names] # ensure data are in the same order as for the genotype data
    # # Get a subset of flower colours for the mothers
    # mother_id = np.array([m.names[0] for m in mothers])
    # mothers_colour_genotypes = ros_sulf.loc[mother_id]

    # RUN THE MCMC.
    log_file.write('MCMC set up. Beginning Markov chain...\n')
    log_file.close()
    t0 = time()
    track_accept = 0 # to record mean MH-acceptance rate
    for i in tqdm(range(nreps)):
        # UPDATE PARAMETERS
        new_model = update_parameters(current_model, proposal_sigma)
        # Update distance travelled using new parameters
        new_model['mean_dist'] = stdev_GND(scale = new_model['scale'],
                                           shape = new_model['shape'])

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
            'shape'     : gma.pdf(new_model['shape'],   a=1.1, scale = 4),
            'mean_dist' : gma.pdf(new_model['mean_dist'],      a=2,   scale= 200),
            'scale'     : 1
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
