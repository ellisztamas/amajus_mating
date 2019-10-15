import numpy as np
import faps as fp
import pandas as pd
from utils import *

# Import data
progeny = fp.read_genotypes('../data_processed/offspring_2012_genotypes.csv', mothers_col=1, genotype_col=2)
adults  = fp.read_genotypes('../data_processed/parents_2012_genotypes.csv')

# data cleaning. See `SNP checking and cleaning.ipynb` for details
# remove individuals with >7.5% midding data
i = 0.075
progeny = progeny.subset(individuals= progeny.missing_data('individual') < i)
adults  = adults.subset(individuals= adults.missing_data('individual') < i)
# remove loci with more than 10% missing data
r = 0.1
lx = (progeny.missing_data('marker') < r) * (adults.missing_data('marker') < r)
progeny = progeny.subset(loci=lx)
adults  = adults.subset(loci=lx)
# remove loci with heterozygosity less than 15%
h = 0.15
progeny = progeny.subset(loci= (adults.heterozygosity('marker') > h))
adults  = adults.subset( loci= (adults.heterozygosity('marker') > h))
# Split into maternal families.
mothers = adults.subset(individuals=np.unique(progeny.mothers))
mothers = mothers.split(np.unique(progeny.mothers))
progeny = progeny.split(by=progeny.mothers)

mu = 0.0013 # genotype error rate.

# Import GPS data
gps = pd.read_csv("../data_processed/GPS_positions.csv", index_col=0)
gps = gps.loc[adults.names] # reorder to match cleaned SNP data
gps = gps[['Easting','Northing']] # remove the column for altitude. We don't need that here.
gps_mothers = gps.loc[mothers.keys()] # GPS coordinates for the mothers only.
# create a matrix of distances between
distance_matrix = (gps_mothers.values[np.newaxis] - gps.values[:, np.newaxis])**2
distance_matrix = distance_matrix.sum(axis=2)
distance_matrix = np.sqrt(distance_matrix).T

# Import flower colour data for the population
ros_sulf = pd.read_csv('../data_processed/rosea_sulfurea.csv', index_col='id')
ros_sulf = ros_sulf.loc[adults.names] # ensure data are in the same order as for the genotype data
# subset flower colour data for the mothers only.
mothers_colour_genotypes = ros_sulf.loc[mothers.keys()]
