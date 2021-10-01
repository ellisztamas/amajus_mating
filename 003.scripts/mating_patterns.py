import numpy as np
from numpy.core.numeric import outer
import pandas as pd
from tqdm import tqdm

# from faps import posterior_mating
# from amajusmating.dispersal import generalised_normal_PDF
from amajusmating import mating

np.random.seed(121)
input_folder = "005.results/004_mcmc_restrict_kurtosis/output/"
output_folder = input_folder
burnin = 500
max_distance = np.inf
spatial_bins =[-827, -70, 70, 244, 270, 560]
# FAPS objects and distance matrices are generated in a separate script.
exec(open('003.scripts/setup_FAPS_GPS.py').read())

data = am_data





# # Other ways to estimate number of offspring whose father was unsampled.
# # 1. From the last column in posterior_paternity_matrix.
# orphans = np.array([np.exp(v.posterior_paternity_matrix()[:, -1]).sum() for v in data.sibships.values()]).sum()
# orphans / np.array([len(v.offspring) for v in data.sibships.values()]).sum()
# # 2. n offspring for whom the most likely candidate is 'missing'
# np.mean(fp.summarise_paternity(data.sibships)['candidate_1'] == "missing")



