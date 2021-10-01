import numpy as np
import pandas as pd
from glob import glob
from tqdm import tqdm
import os
import wquantiles as wq

from amajusmating.mating import observed_sires

# FAPS objects and distance matrices are generated in a separate script.
exec(open('003.scripts/setup_FAPS_GPS.py').read())

observed_sires(
    data = am_data,
    folder = "005.results/003_mcmc_restrict_kurtosis/output/",
    burnin = 3480,
)


chains = glob("005.results/003_mcmc_restrict_kurtosis/output/*.out")
burnin = 1500
max_distance = np.inf # set an upper bound on dispersal distances.



# am_data.sires['prob'].sum()
# am_data.sires.loc[am_data.sires['prob'] >0.99]['prob'].sum()

# wq.median(am_data.sires['distance'], am_data.sires['prob'])
# fp.posterior_mating(am_data.sibships, covariates_only=True, use_covariates = True, ndraws=1000)
