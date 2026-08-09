# correctEpiAge_MSA2EPICv2

This repository is dedicated to the code availability for the paper: **Cross-Platform Calibration of Epigenetic Age Between the EPICv2 and MSA Arrays**

## Scripts

- `epigenetic_clock_msa_vs_epicv2.R`: R script used to generate the results reported in the paper.
- `calibrate_msa2epic.R`: R script containing the calibration function for mapping MSA-based epigenetic clock values to the EPICv2-referenced scale.

## Usage of the calibration function

First, load the calibration function:

```r
source("calibrate_msa2epic.R")
```

### Horvath clock using Model 1

```r
calibrate_msa_to_epicv2(
  msa_value = 46.8,
  clock = "Horvath",
  model = 1
)
```

### PhenoAge using Model 3

```r
dat$calibrated_phenoage <- calibrate_msa_to_epicv2(
  msa_value = dat$PhenoAge_MSA,
  clock = "PhenoAge",
  model = 3,
  age = dat$Age,
  sex = dat$Sex
)
```

For Model 3, `age` denotes chronological age, and `sex` should be coded as `0` for female and `1` for male.
