"""
Tom Ellis, March 2021

Class to hold FAPS objects and variables related to dispersal and assortative mating.
"""

import numpy as np
import pandas as pd
import faps as fp
from time import time, strftime
from tqdm import tqdm

from amajusmating import dispersal
from amajusmating import mcmc

class faps_data(object):
    def __init__(self, paternity, gps, flower_colours, params):
        self.paternity = paternity
        self.gps = gps
        self.flower_colours = flower_colours
        self.params = {}

        # Check the mothers are present in data files.
        self.mothers = self.paternity.keys()
        if any([k not in gps.index for k in self.mothers]):
            ValueError("Not all mothers in the paternity array are present in GPS data.")
        if any([k not in flower_colours.index for k in self.mothers]):
            ValueError("Not all mothers in the paternity array are present in flower-colour data.")

        # Check there are the expected numbers of fathers
        n_candidates = [len(v.candidates) for v in self.paternity.values()]
        assert n_candidates.count(n_candidates[0]) == len(n_candidates)
        n_candidates = n_candidates[0]

        # Distance matrix
        self.distances = dispersal.distance_matrix(gps.loc[self.mothers].to_numpy(), self.gps.to_numpy())
        assert self.distances.shape == (len(self.mothers), n_candidates)

        self.covariates = {}

    def update_missing_dads(self, x:float):
        """
        Update the proportion of missing fathers in each paternity array

        Parameters
        ----------
        x: float between 0 and 1
            Proposal value for the proportion of true fathers not in the sample of candidates.
        
        Returns
        -------
        Nothing; missing_parents attribute for and entry in params are (silently) updated.
        """
        if x > 1 or x < 0:
            raise ValueError("Proportion of missing fathers exceeds zero or one.")
        for p in self.paternity.keys():
            self.paternity[p].missing_parents = x
        self.params['missing'] = x

    def update_dispersal_probs(self, scale:float, shape:float, mixture:float, max_distance:float = np.inf):
        """
        Update probabilities of mating assuming a mixture of a generalised normal 
        (GND) pollen dispersal kernel and a uniform distirbution. The weight
        applied to the GND is set with mixture parameter w.

        Parameters
        ----------
        scale: float
            Scale parameter.
        shape: float
            Shape parameter
        w: float between 0 and 1.
            The proportion of probability mass assigned to the
            generalised normal function.

        Returns
        -------
        An array of log probabilities with a row for each mother and a column
        for each candidate father. Rows sum to one when exponentiated.
        """
        # Probability of drawing each male under GND dispersal
        prob_drawn = dispersal.dispersal_GND(
            x     = self.distances,
            scale = scale,
            shape = shape,
            w     = mixture)
        # Identify candidates who are further than the threshold distance
        # and set their log likelihoods to negative infinity
        ix = self.distances > max_distance
        prob_drawn[ix] = -np.inf
        # Update parameters
        self.covariates['dispersal'] = prob_drawn
        self.params['shape'] = shape
        self.params['scale'] = scale
        self.params['mixture'] = mixture

    def sibship_clustering(self, ndraws:int=100, use_covariates:bool=True):
        """
        Cluster the offspring into sibships given genetic and population data.

        Parameters
        ----------
        ndraws: int
            Number of Monte Carlo simulations to run for each partition.
        use_covariates: logical, optional
            If True, information on prbabilities associated with covariates stored
            in paternityArray objects are incorporated into sibship clustering.

        Returns
        -------
        List of sibshipCluster_objects
        """
        # Incorporate into paternity_arrays
        if self.covariates == {}:
            if use_covariates: raise ValueError("use_covariates is True, but no data about covariates have been supplied.")

        else:
            prob_drawn = sum(self.covariates.values())            
            for (p,s) in zip(self.mothers, prob_drawn):
                self.paternity[p].add_covariate(s)

        # INFER FAMILIES
        # Cluster into families and get likelihoods
        self.sibships = fp.sibship_clustering(self.paternity, ndraws=ndraws, use_covariates=use_covariates)
        self.params['loglik'] = np.array([fp.alogsumexp(s.lik_partitions) for s in self.sibships.values()]).sum()

        return None

    def run_MCMC(self, initial_parameters, proposal_sigma, priors, nreps, output_dir, chain_name, thin=1, max_distance = np.inf):
        """
        A wrapper function to run Metropolis-Hastings MCMC for paternity and dispersal.

        Parameters
        ----------
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

        # For the first iteration, start with an arbitary probabilities
        # It doesn't matter, because for MH-ratio skips the first iteration by default.
        current_model['log_prior']     = -10e12
        current_model['loglik']        = -10e12
        current_model['log_posterior'] = -10e12
        # Set up the datafiles
        out_file = output_dir + chain_name
        mcmc.setup_output(out_file, current_model, proposal_sigma, nreps, thin)

        #  IMPORT AND FORMAT DATA.
        t0 = time()

        # RUN THE MCMC.
        for i in tqdm(range(nreps)):
            # UPDATE PARAMETERS
            new_model = mcmc.update_parameters(current_model, proposal_sigma)
            # Calculate log prior probabilities of the new model.
            prior_probs = list(priors(new_model).values())
            new_model['log_prior'] = np.log(prior_probs).sum()

            # LOG PROBABILITIES OF PATERNITY FOR EACH PARAMETER
            # Update proportion of missing fathers
            self.update_missing_dads(new_model['missing'])
            # Update dispersal
            self.update_dispersal_probs(
                scale = new_model['scale'],
                shape = new_model['shape'],
                mixture = new_model['mixture']
            )
            # Identify candidates who are further than the distance threshold
            # and set their log likelihoods to negative infinity
            ix = self.distances > max_distance
            self.covariates['dispersal'][ix] = -np.inf

            # INFER FAMILIES
            # Incorporate covariate information into paternity_arrays
            cov = sum(self.covariates.values())
            for (p,s) in zip(self.paternity, cov):
                self.paternity[p].add_covariate(s)
            # Cluster into families and get likelihoods
            self.sibship_clustering(ndraws=100, use_covariates=True)

            # Probability components for the new model.
            new_model['loglik'] = np.array([fp.alogsumexp(s.lik_partitions) for s in self.sibships.values()]).sum()
            new_model['log_posterior'] = new_model['loglik'] + new_model['log_prior']

            # Decide whether to accept the new model.
            accept = mcmc.mh_ratio(
                current = current_model['log_posterior'],
                new     = new_model['log_posterior']
            )
            if accept:
                current_model = new_model

            # write iteration to disk
            if(i in np.arange(start = 0, stop = nreps, step=thin)):
                mcmc.write_output(out_file + '.out', i, current_model, decimals=3, time0=t0)
        
        log_file = open(output_dir + chain_name + ".log", 'a')
        log_file.write('\nMCMC completed {}.\n'.format(strftime("%Y-%m-%d %H:%M:%S")))
        log_file.close()