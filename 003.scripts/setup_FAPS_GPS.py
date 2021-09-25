"""
Tom Ellis, October 2019.

This script imports and prepares data on
1. Genotype data of parents and offspring
2. GPS data
3. Flower colour phenotypes

It then creates and saves a paternity array object for genotype data.
If this file exists when this script is called again, the saved version
is imported
Finally genotype and paternity arrays are split by maternal family.

It returns a faps_data object with all the relevant information.
"""

import numpy as np
import faps as fp
import pandas as pd
from os.path import isfile
# Import local modules
from amajusmating import faps_data

print("Setting up genotype, flower-colour and GPS information using FAPS version {}.\n".format(fp.__version__))

# GENOTYPE DATA
# Genotyping error rate.
mu=0.0013
# Import genpotype data
progeny = fp.read_genotypes('001.data/002.processed/offspring_2012_genotypes.csv', mothers_col=1, genotype_col=2)
adults  = fp.read_genotypes('001.data/002.processed/parents_2012_genotypes.csv')
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
# If the file exists already, this is imported, which is much faster.
# Otherwise a new paternityArray object is created.
if isfile("004.output/paternity_array.csv"):
    # genotypeArray giving the mother of each maternal family.
    mothers = adults.subset(individuals=np.unique(progeny.mothers))
     # Split into maternal families.
    # Import saved paternityArray.
    patlik = fp.read_paternity_array("004.output/paternity_array.csv")
else:
    # A single genotypeArray giving the mother of each of 984 offspring individuals.
    mothers = adults.subset(individuals=progeny.mothers)
    # Create the paternity array and save for later.
    patlik = fp.paternity_array(
        progeny, mothers, adults, mu = mu, missing_parents= 0.32
        )
    patlik.write("004.output/paternity_array.csv")

# SPLIT BY MATERNAL FAMILYs
# genotypeArray giving the mother of each maternal family.
mothers = adults.subset(individuals=np.unique(progeny.mothers))
# Split objects up by maternal family.
mothers = mothers.split(np.unique(progeny.mothers))
patlik  = patlik.split(by=progeny.mothers)
progeny = progeny.split(by=progeny.mothers)

# GPS DATA
# Import GPS data
gps = pd.read_csv("001.data/002.processed/GPS_positions.csv", index_col=0)
gps = gps.loc[adults.names] # reorder to match cleaned SNP data
gps = gps[['Easting','Northing']] # remove the column for altitude. We don't need that here.
# gps_mothers = gps.loc[mothers.keys()] # GPS coordinates for the mothers only.

# FLOWER COLOUR
# Import flower colour data for the population
ros_sulf = pd.read_csv('001.data/002.processed/rosea_sulfurea.csv', index_col='id')
ros_sulf = ros_sulf.loc[adults.names] # ensure data are in the same order as for the genotype data
# Simplify flower colours to yellow, full red or hybrid
ros_sulf['simple_colour'] = 'unkown'
ros_sulf.loc[ros_sulf['flower_colour'].isin(['FR',"Ye"]), 'simple_colour'] = ros_sulf['flower_colour']
ros_sulf.loc[ros_sulf['flower_colour'].isin(['FO',"WR", "WO", "Wh"]), 'simple_colour'] = 'hybrid'

# Create a class object to hold the data.
am_data = faps_data(
    paternity=patlik,
    gps = gps,
    flower_colours = ros_sulf[['rosea','sulfurea', 'simple_colour']],
    params = {}
    )

# Remove variables we don't need anymore.
del(progeny, md, r, h, ros_sulf, patlik, lx, adults, mothers)