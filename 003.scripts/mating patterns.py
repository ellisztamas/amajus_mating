# from os import WEXITED
import numpy as np
import wquantiles as wq
import pandas as pd
import faps as fp
from amajusmating import mating
from amajusmating.dispersal import generalised_normal_PDF

# FAPS objects and distance matrices are generated in a separate script.
exec(open('003.scripts/setup_FAPS_GPS.py').read())

model = {
    'missing' : 0.3,
    "mixture" : 1,
    "shape"   : 0.6,
    "scale"   : 30
}

data = am_data
max_distance = np.inf
ndraws=1000
# Update data with parameter values
data.update_covariate_probs(model = model, max_distance = max_distance)

# Cluster into sibships, if not already done.
if not hasattr(data, 'sibships'):
    data.sibship_clustering(ndraws = ndraws, use_covariates = True)
# Call siring events
data.sires = fp.summarise_sires(data.sibships)

# Add data on distance between mother and father and flower colours.
# Turn distance matrix into a data frame so we can use .loc on it.
distance_df = pd.DataFrame({
    'mother'   : np.repeat(list(data.mothers), data.n_candidates),
    'father'   : np.tile(data.candidates, len(data.mothers)),
    'distance' : data.distances.flatten()
})
# Mating probs if dispersal was random
distance_df['dispersal_prob'] = generalised_normal_PDF(
    distance_df['distance'],
    scale = model['scale'],
    shape = model['shape']
)

# Merge siring events with distances and phenotypes.
# Observed mating events
data.sires = (data.sires.merge(distance_df, how="left", on=['mother', 'father']).
merge(data.flower_colours, how="left", left_on="mother", right_index=True).
merge(data.flower_colours, how="left", left_on="father", right_index=True, suffixes = ['_mother', "_father"])
)
data.sires = data.sires.assign(
    match = data.sires['rosea_mother'] == data.sires['rosea_father'],
    cross = lambda x: x.rosea_mother + "/" + x.rosea_father
    )
# Random mating
random_sires = distance_df.\
    merge(data.flower_colours, how="left", left_on="mother", right_index=True).\
        merge(data.flower_colours, how="left", left_on="father", right_index=True, suffixes = ['_mother', "_father"])
random_sires = random_sires.assign(
    match = random_sires['rosea_mother'] == random_sires['rosea_father'],
    cross = lambda x: x.rosea_mother + "/" + x.rosea_father
    )

median_dispersal = wq.median(data.sires['distance'], weights = data.sires['prob'])

# data.sires = data.sires.loc[lambda x: x.distance < median_dispersal]
# random_sires = random_sires.loc[lambda x: x.distance < median_dispersal]

# Strict assortative mating returns no differences
mating_probs = pd.DataFrame({
    'obs'    : data.sires.groupby(by='match').sum()['prob'],
    'random' : random_sires.groupby(by='match').sum()['dispersal_prob']
}).transpose().to_numpy()
mating_probs / mating_probs.sum(1)[:, np.newaxis]

# Try by mating type
mating_probs = pd.DataFrame({
    'obs'    : data.sires.loc[lambda x: x.prob > 0.9].groupby(by='cross').sum()['prob'],
    'random' : random_sires.groupby(by='cross').sum()['dispersal_prob']
})#.to_numpy()
mating_probs['obs'][:3] / mating_probs['obs'][:3].sum()
mating_probs['random'][:3] / mating_probs['random'][:3].sum()
# Hint that yellows mate
mating_probs['obs'][3:6] / mating_probs['obs'][3:6].sum()
mating_probs['random'][3:6] / mating_probs['random'][3:6].sum()

mating_probs['obs'][6:9] / mating_probs['obs'][6:9].sum()
mating_probs['random'][6:9] / mating_probs['random'][6:9].sum()

# Siring success
mating_probs = pd.DataFrame({
    'obs'    : data.sires.loc[lambda x: x.prob > 0.9].groupby(by='rosea_father').sum()['prob'],
    'random' : random_sires.groupby(by='rosea_father').sum()['dispersal_prob']
}).transpose().to_numpy()
m = mating_probs / mating_probs.sum(1)[:, np.newaxis]
m[0] / m[1]

np.quantile(data.sires["Easting"], q=[0, 0.2, 0.4, 0.6, 0.8, 1])


data.sires.groupby(by= "rosea_father").median()['distance']

data.sires = data.sires.merge(data.gps['Easting'], how='left', left_on='mother', right_index=True)
random_sires = random_sires.merge(data.gps, how="left", left_on="mother", right_index=True)

o = data.sires.loc[ (data.sires["Easting"] < -70)].groupby(by='rosea_father').sum()['prob'] 
r = random_sires.loc[(random_sires["Easting"] > -70) ].groupby(by='rosea_father').sum()['dispersal_prob']
m1 = (o/o.sum())/(r/r.sum())

o = data.sires.loc[(data.sires["Easting"] > -70) & (data.sires["Easting"] < 70)].groupby(by='rosea_father').sum()['prob'] 
r = random_sires.loc[(random_sires["Easting"] > -70) & (random_sires["Easting"] < 70)].groupby(by='rosea_father').sum()['dispersal_prob']
m2 = (o/o.sum())/(r/r.sum())

o = data.sires.loc[(data.sires["Easting"] > 70) & (data.sires["Easting"] < 244)].groupby(by='rosea_father').sum()['prob'] 
r = random_sires.loc[(random_sires["Easting"] > 70) & (random_sires["Easting"] < 244)].groupby(by='rosea_father').sum()['dispersal_prob']
m3 = (o/o.sum())/(r/r.sum())

o = data.sires.loc[(data.sires["Easting"] > 244) & (data.sires["Easting"] < 261)].groupby(by='rosea_father').sum()['prob'] 
r = random_sires.loc[(random_sires["Easting"] > 244) & (random_sires["Easting"] < 261)].groupby(by='rosea_father').sum()['dispersal_prob']
m4 = (o/o.sum())/(r/r.sum())

o = data.sires.loc[(data.sires["Easting"] > 261)].groupby(by='rosea_father').sum()['prob'] 
r = random_sires.loc[(random_sires["Easting"] > 261)].groupby(by='rosea_father').sum()['dispersal_prob']
m5 = (o/o.sum())/(r/r.sum())

[m1, m2, m3, m4, m5]

data.sires.loc[(data.sires["Easting"] > 244) & (data.sires["Easting"] < 261) & (data.sires["rosea_father"] == "Ye")]

# No evidence that yellow sires have increased heterozygosity or missing SNP data
pd.DataFrame({
    'names' : adults.names,
    'misisng' : adults.missing_data(by="individual") * adults.nloci,
    'het' : adults.heterozygosity(by="individual"),
    'flower_colour' : ros_sulf['flower_colour']
    }
    ).merge(data.sires, how="right", right_on="father", left_index=True).\
    groupby(by="flower_colour").mean()

data.sires.loc[data.sires['prob'] > 0.1]['father'].value_counts().head(20)

data.sires.loc[(data.sires['father'] == "K1768") & (data.sires['prob'] > 0.01)]

wq.median(data.sires["distance"], data.sires["prob"])