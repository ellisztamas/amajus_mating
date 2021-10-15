# Paternity without using dispersal as a covariate

**Date:** 15th October 2021
**Author:** Tom Ellis

## Background
In principle one should get better estimates of both the pedigree and dispersal
by estimating these together. As a comparison, check paternity results without
including dispersal.

## What did you do?
- `missing_dads.ipynb` estimates the proportion of missing pollen donors and the
    number of offspring with unsampled fathers
- `dispersal.Rmd` plots the dispersal kernel, and cumulative distribution of 
    dispersal


## Main conclusion
- 17.5% of pollen donors are absent
- 32.8% of progeny have an unsampled father.
- Weighted median is 42m (cf 31 when you include a kernel)
- The dispersal kernels look really similar to those using dispersal tbh

## Caveats

## Follow-up
