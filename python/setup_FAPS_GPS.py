import numpy as np
from faps import *
import pandas as pd
from dispersal_functions import *

progeny = read_genotypes('../data_processed/offspring_2012_genotypes.csv', mothers_col=1, genotype_col=2)
adults  = read_genotypes('../data_processed/parents_2012_genotypes.csv')

# data cleaning. See `SNP checking and cleaning.ipynb` for details
# remove individuals with >7.5% midding data
i = 0.075
progeny = progeny.subset(    individuals= progeny.missing_data('individual') < i)
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

# pull out genotype information on just the mothers
mothers = adults.subset(progeny.parent_index('m', adults.names))
# split into
ox = progeny.mothers
mothers = mothers.split(ox)
progeny = progeny.split(ox)
mother_id = np.array([np.unique(x.mothers) for x in progeny]).squeeze()
# Sanity check to ensure maternal IDs line up.
all(mother_id == np.array([np.unique(x.names) for x in mothers]).squeeze())

#allele_freqs = adults.allele_freqs() # population allele frequencies
mu = 0.0013 # genotype error rate.
del ox, i, h, lx, r
print('Genotype data imported.')



# Import GPS data
gps         = pd.read_csv("../data_processed/GPS_positions.csv", index_col=0)
gps         = gps.loc[adults.names] # reorder to match cleaned SNP data
gps         = gps[['Easting','Northing']] # remove the column for altitude. We don't need that here.
gps_mothers = gps.loc[mother_id]
# create a matrix of distance_matrixances between
distance_matrix = (gps_mothers.values[np.newaxis] - gps.values[:, np.newaxis])**2
distance_matrix = distance_matrix.sum(axis=2)
distance_matrix = np.sqrt(distance_matrix).T

del gps, gps_mothers # we don't need these anymore.
print('GPS data imported.')