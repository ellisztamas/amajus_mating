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

        # Boolean matrix indicating when mothers (rows) have the same flower
        # colour as candidate fathers (columns).
        maternal_colours  = self.flower_colours.loc[self.mothers].to_numpy()[:, np.newaxis]
        candidate_colours = self.flower_colours.to_numpy()[np.newaxis]
        self.matches = maternal_colours == candidate_colours

        # Empty dictionary to store probabilities for covariates
        self.covariates = {}

    def update_missing_dads(self, x):
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
            w     = mixture
            )
        # Identify candidates who are further than the threshold distance
        # and set their log likelihoods to negative infinity
        ix = self.distances > max_distance
        prob_drawn[ix] = -np.inf
        # Update parameters
        self.covariates['dispersal'] = prob_drawn
        self.params['shape'] = shape
        self.params['scale'] = scale
        self.params['mixture'] = mixture

    def update_assortment_probs(self, p):
        """
        Update assortative mating probabilities.

        Flower colours of mothers and candidate fathers mate with probability p
        if they are of the same phenotype (i.e. if their element in `self. assortment`
        is True) and with probability 1-p if not.

        Parameters
        ----------
        p: float between 0 and 1.
        """
        assert(p >= 0)
        assert(p <= 1)
        # Assign probability p for pairs of individuals of the same phenotype
        # 1-p for pairs that are different
        assortment = np.where(self.matches, p, 1-p)
        # Normalise rows and log.
        assortment = assortment / assortment.sum(1, keepdims=True)
        assortment = np.log(assortment)
        # Update parameters
        self.covariates['assortment'] = assortment
        self.params['assortment'] = p

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
        if self.covariates == {} and use_covariates:
            raise ValueError("use_covariates is True, but no data about covariates have been supplied.")

        else:
            # Multiply covariate probabilities, and normalise rows to sum to 1.
            prob_drawn = sum(self.covariates.values())
            prob_drawn - fp.alogsumexp(prob_drawn)

            for (p,s) in zip(self.mothers, prob_drawn):
                self.paternity[p].add_covariate(s)

        # INFER FAMILIES
        # Cluster into families and get likelihoods
        self.sibships = fp.sibship_clustering(self.paternity, ndraws=ndraws, use_covariates=use_covariates)
        self.params['loglik'] = np.array([fp.alogsumexp(s.lik_partitions) for s in self.sibships.values()]).sum()

        return None