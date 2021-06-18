# Assortative mating in a snapdragon hybrid zone

This repository documents analysis associated with the manuscript "Pollinators cause assortative mating for flower colour in a snapdragon hybrid zone" by Tom Ellis, David Field and Nick Barton.

## Table of contents

1. [Experimental set up](#experimental-set-up)
2. [Analysis workflow](#analysis-workflow)
3. [Dependencies](#Dependencies)
4. [Data](#data)
    1. [Parental lines](#parental-lines)
    2. [RIL data](#ril-data)
    3. [Derived data](#derived-data)
    4. [R/QTL objects](#rqtl-objects)

## Experimental set up

## Analysis workflow

1. Format raw SNP, GPS and flower colour data
2. Import genotype and GPS information into Python, and create paternityArray objects.
3. MCMC for the joint analysis of sibships, paternity and dispersal.
4. Thinning MCMC chains.
5. Generate and save sibshipCluster objects and mating events across 1000
posterior draws.
6. Quality check on paternity assignments
7. Quantify male fitness and assortative mating in spatial bins.
8. Simulations