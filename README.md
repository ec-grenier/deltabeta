# Run Delta Beta for Mortality

This repository is to run delta beta, a diagnostic tool, across all regions and many model specifications.

## How to run:

Before anything, run a single projection for any model. This is to generate an "allcalcs" file containing the long-run income and climate covariates which jointly determine the shape of a region's curve. Age group specific response function coefficients from the main mortality model are located in `data/csvv/`.

1. Extract covariates from allcalcs files from the single projections generated
2. Run `1_mortality_db_script.R` from command line, R console, or using `run_db.sbatch`. See script for command line arguments to include.
3. (Optional) Combine impacts across age groups. Takes the age group share weighted average.
4. Plot response function (beta) and temperature distributions (delta)

