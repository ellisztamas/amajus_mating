import numpy as np
import pandas as pd
from glob import glob
from tqdm import tqdm
import os
import wquantiles as wq

from faps import posterior_mating
from amajusmating.mating import observed_sires

burnin = 1500

# FAPS objects and distance matrices are generated in a separate script.
exec(open('003.scripts/setup_FAPS_GPS.py').read())

folder = "005.results/003_mcmc_restrict_kurtosis/output/"

data = am_data
max_distance = np.inf

# Import MCMC chains
chains = glob(folder+"/*out")
posterior = {}
for chain in chains:
    k = os.path.basename(chain)
    posterior[k] = pd.read_csv(chain, sep="\t").loc[lambda x: x.iter >= burnin]
posterior = pd.concat(posterior).reset_index()
# For this script we are only interested in the generalised-Gaussian part
# of the dispersal kernel, so set `mixture` to 1 for all rows.
posterior['mixture'] = 1

siring = {} 
# Loop over steps in the MCMC chain and get siring events for each.
# for i in tqdm(posterior.index):
i = 4
model = posterior.loc[i]
# Update data with parameter values
data.update_covariate_probs(model = model, max_distance = max_distance)
# Cluster into sibships, if not already done.
data.sibship_clustering(ndraws = 1000, use_covariates = True)

# We need to set mating probabilities for selfing to zero
# Otherwise the mother would be drawn many times, and cause many invalid partitions.
jx = [np.where(x == np.array(data.candidates))[0][0] for x in data.mothers]
for i,j in zip(range(60), jx):
    data.covariates['dispersal'][i, j] = -np.inf

np.random.seed(12)
posterior_mating(am_data.sibships, covariates_only=True, ndraws=10000)

covariate = data.covariates['dispersal'][0]
nfathers   = covariate.shape[0]
prob_array = covariate - fp.alogsumexp(covariate)
prob_array = np.tile(prob_array, nfamilies).reshape([nfamilies, len(covariate)])
prob_array = np.exp(prob_array)


# Need to work out why draw_fathers keeps picking father 2094.
path_samples = np.array([np.random.choice(range(nfathers), 100, replace=True, p = prob_array[i]) for i in range(nfamilies)])
path_samples = path_samples.T