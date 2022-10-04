import numpy as np
import pandas as pd
from glob import glob
from tqdm import tqdm
import os
import seaborn as sns


from faps import posterior_mating
from amajusmating.dispersal import generalised_normal_PDF
from amajusmating import mating

np.random.seed(121)
burnin = 500

# FAPS objects and distance matrices are generated in a separate script.
exec(open('003.scripts/setup_FAPS_GPS.py').read())


folder = "005.results/00*_mcmc_restrict_kurtosis/output/"

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

# Loop over steps in the MCMC chain and get siring events for each.
# for i in tqdm(posterior.index):
i = 4
model = posterior.loc[i]

data.mating = mating.simulate_mating(am_data, model, ndraws=100)


# # Add data on distance between mother and father and flower colours.
# # Turn distance matrix into a data frame so we can use .loc on it.
# distance_df = pd.DataFrame({
#     'mother'   : np.repeat(list(data.mothers), data.n_candidates),
#     'father'   : np.tile(data.candidates, len(data.mothers)),
#     'distance' : data.distances.flatten()
# })

# Merge siring events with distances and phenotypes.
# Observed mating events
# data.mating['obs'] = data.mating['obs'].\
#     merge(data.gps['Easting'], how='left', left_on='mother', right_index=True) .\
#     merge(distance_df, how="left", on=['mother', 'father']).\
#     merge(data.flower_colours, how="left", left_on="mother", right_index=True).\
#     merge(data.flower_colours, how="left", left_on="father", right_index=True, suffixes = ['_mother', "_father"]).keys()
# # Random mating
# data.mating['exp'] = data.mating['exp'].\
#     merge(data.gps['Easting'], how='left', left_on='mother', right_index=True).\
#     merge(distance_df, how="left", on=['mother', 'father']).\
#     merge(data.flower_colours, how="left", left_on="mother", right_index=True).\
#     merge(data.flower_colours, how="left", left_on="father", right_index=True, suffixes = ['_mother', "_father"])

# Summarise mating events and dispersal
{
    'n_mating_events'  : np.sum(data.mating['obs']['father'] != "missing"),
    'missing_dads'     : np.sum(data.mating['obs']['father'] == "missing"),
    'orphans'          : data.mating['obs'].loc[data.mating['obs']['father'] == "missing"]['offspring_sired'].sum(),
    'empty_sibships'   : (data.mating['obs'].loc[data.mating['obs']['father'] == "missing"]['frequency'] == 1).sum(),
    'median_dispersal' : np.median(data.mating['obs']['distance'].dropna()),
    'dispersal>100m'    : np.mean(data.mating['obs']['distance'] > 100),
    'dispersal>500m'    : np.mean(data.mating['obs']['distance'] > 500),
    'dispersal>1000m'   : np.mean(data.mating['obs']['distance'] > 1000),
}

# # Other ways to estimate number of offspring whose father was unsampled.
# # 1. From the last column in posterior_paternity_matrix.
# orphans = np.array([np.exp(v.posterior_paternity_matrix()[:, -1]).sum() for v in data.sibships.values()]).sum()
# orphans / np.array([len(v.offspring) for v in data.sibships.values()]).sum()
# # 2. n offspring for whom the most likely candidate is 'missing'
# np.mean(fp.summarise_paternity(data.sibships)['candidate_1'] == "missing")

# Remove cases where the colour of the father is unkonwn
obs = data.mating['obs'].\
    query('simple_colour_father in ["FR", "hybrid", "Ye"]')
exp = data.mating['exp'].\
    query('simple_colour_father in ["FR", "hybrid", "Ye"]')

for geno in ["FR", "hybrid", "Ye"]:
    w_obs = (obs['simple_colour_father'] == geno).mean()
    w_exp = (exp['simple_colour_father'] == geno).mean()
    print(w_obs/w_exp)
# simple assortment
ass_obs = (obs['simple_colour_father'] == obs['simple_colour_mother']).mean()
ass_exp = (exp['simple_colour_father'] == exp['simple_colour_mother']).mean()
ass_obs / ass_exp

obs['match'] = (obs['simple_colour_father'] == obs['simple_colour_mother'])
obs['bin'] = pd.cut(obs['distance'], np.arange(0,2000,32))

obs.groupby('bin').mean()['match']

import seaborn as sns

# Mating wrt flower colour in spatial bins
boundaries  = [-np.inf, -70, 70, 244, 270, np.inf]

region_obs = pd.cut(obs['Easting'], bins = boundaries, labels = np.arange(1, len(boundaries)))
region_exp = pd.cut(exp['Easting'], bins = boundaries, labels = np.arange(1, len(boundaries)))
for r in range(1, len(boundaries)):
    for geno in ["FR", "hybrid", "Ye"]:
        w_obs = (obs.loc[region_obs == r]['simple_colour_father'] == geno).mean()
        w_exp = (exp.loc[region_exp == r]['simple_colour_father'] == geno).mean()
        print(r, geno, w_obs/w_exp)

for r in range(1, len(boundaries)):
    for geno in ["R/R", "R/r", "r/r"]:
        w_obs = (obs.loc[region_obs == r]['rosea_father'] == geno).mean()
        w_exp = (exp.loc[region_exp == r]['rosea_father'] == geno).mean()
        print(r, geno, w_obs/w_exp)

for r in range(1, len(boundaries)):
    for geno in ["S/+", "s/s"]:
        w_obs = (obs.loc[region_obs == r]['sulfurea_father'] == geno).mean()
        w_exp = (exp.loc[region_exp == r]['sulfurea_father'] == geno).mean()
        print(r, geno, w_obs/w_exp)
