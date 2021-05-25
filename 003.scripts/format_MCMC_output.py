# Tom Ellis, April 2020
#
# Script to concatenate and thin MCMC output

import pandas as pd
import numpy as np

ix = np.arange(500, 10500, 40)
mcmc = pd.concat([
    pd.read_csv("../004.output/mcmc/mcmc_nolimit_repA.txt", sep="\t").iloc[ix],
    pd.read_csv("../004.output/mcmc/mcmc_nolimit_repB.txt", sep="\t").iloc[ix],
    pd.read_csv("../004.output/mcmc/mcmc_nolimit_repC.txt", sep="\t").iloc[ix],
    pd.read_csv("../004.output/mcmc/mcmc_nolimit_repD.txt", sep="\t").iloc[ix],
])
# Add label for chain
mcmc['chain'] = np.repeat(['A', "B", "C", "D"], len(ix))
mcmc = mcmc.reset_index()

# Save to disk
mcmc.to_csv("../004.output/mcmc_thinned.csv", index=False)
