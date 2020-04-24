# Tom Ellis, October 2019.
#
# This script imports and prerpares data on
# 1. Genotype data of parents and offspring
# 2. GPS data
# 3. Flower colour phenotypes
#
# It then creates and saves a paternity array object for genotype data.
# If this file exists when this script is called again, the saved version
# is imported
#
# Finally genotype and paternity arrays are split by maternal family.

import numpy as np
import faps as fp
import pandas as pd
from os.path import isfile
#from utils_dispersal import *

# GENOTYPE DATA
# Genotyping error rate.
mu=0.0013
# Import genpotype data
progeny = fp.read_genotypes('../data_processed/offspring_2012_genotypes.csv', mothers_col=1, genotype_col=2)
adults  = fp.read_genotypes('../data_processed/parents_2012_genotypes.csv')
# SNP Data cleaning.
# remove individuals with >7.5% midding data
md = 0.075
progeny = progeny.subset(individuals= progeny.missing_data('individual') < md)
adults  = adults.subset(individuals= adults.missing_data('individual') < md)
# remove loci with more than 10% missing data
r = 0.1
lx = (progeny.missing_data('marker') < r) * (adults.missing_data('marker') < r)
progeny = progeny.subset(loci=lx)
adults  = adults.subset(loci=lx)
# remove loci with heterozygosity less than 15%
h = 0.15
progeny = progeny.subset(loci= (adults.heterozygosity('marker') > h))
adults  = adults.subset( loci= (adults.heterozygosity('marker') > h))

# PATERNITY ARRAYS
# If this script is run for the first time, the paternityArray is saved for later.
# If the file exists, this is imported, which is much faster.
if isfile("../output/paternity_array.csv"):
    # genotypeArray giving the mother of each maternal family.
    mothers = adults.subset(individuals=np.unique(progeny.mothers))
     # Split into maternal families.
    # Import saved paternityArray.
    patlik = fp.read_paternity_array("../output/paternity_array.csv")
else:
    # A single genotypeArray giving the mother of each of 984 offspring individuals.
    mothers = adults.subset(individuals=progeny.mothers)
    # Create the paternity array and save for later.
    patlik = fp.paternity_array(progeny, mothers, adults, mu = mu, missing_parents= 0.15)
    patlik.write("../output/paternity_array.csv")

# SPLIT BY MATERNAL FAMILYs
# genotypeArray giving the mother of each maternal family.
mothers = adults.subset(individuals=np.unique(progeny.mothers))
# Split objects up by maternal family.
mothers = mothers.split(np.unique(progeny.mothers))
patlik  = patlik.split(by=progeny.mothers)
progeny = progeny.split(by=progeny.mothers)

# GPS DATA
# Import GPS data
gps = pd.read_csv("../data_processed/GPS_positions.csv", index_col=0)
gps = gps.loc[adults.names] # reorder to match cleaned SNP data
gps = gps[['Easting','Northing']] # remove the column for altitude. We don't need that here.
gps_mothers = gps.loc[mothers.keys()] # GPS coordinates for the mothers only.
# create a matrix of distances between
distance_matrix = (gps_mothers.values[np.newaxis] - gps.values[:, np.newaxis])**2
distance_matrix = distance_matrix.sum(axis=2)
distance_matrix = np.sqrt(distance_matrix).T

# FLOWER COLOUR
# Import flower colour data for the population
ros_sulf = pd.read_csv('../data_processed/rosea_sulfurea.csv', index_col='id')
ros_sulf = ros_sulf.loc[adults.names] # ensure data are in the same order as for the genotype data
# subset flower colour data for the mothers only.
mothers_colour_genotypes = ros_sulf.loc[mothers.keys()]

del(md, r, h)
