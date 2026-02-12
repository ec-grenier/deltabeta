#=================================================================================#
# Step 0: Extract covariates
#
# Author: Elliot Grenier (egrenier@uchicago.edu)
# 
# Description:
#'   
#'   Pulls long run income and climate covariates needed for delta beta
#'   calculations. Takes long-run income (13 half bartlett kernel) and long-run 
#'   average temperature (30 year average) from mortality projections
#'   
# How to run:
#'
#'  Run from command line `Rscript 0_extract_covariates.R` or in the IDE of your
#'  chooosing.
#  
#=================================================================================#

#=================================================================================#
# packages and paths ----

packages = c("glue", "dplyr", "readr", "purrr")
invisible(lapply(packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))
rm(packages)

# projection system allcalcs file parent directories
root = "/project/cil/gcp/outputs/mortality/impacts-darwin/single/single-diagnostic-jan2026"
out = "/project/cil/home_dirs/egrenier/cil-comms/deltabeta/data/covars"
allcalcs = "tests.configs.mortality.allmodels-allcalcs-Agespec_interaction_response-oldest.csv" # allcalcs covariates are the same across age groups 

# specifications we need covariates for
specs = expand.grid(rcp = c('rcp45', 'rcp85'),
                    gcm = c('CCSM4'),
                    iam = c('low', 'high'),
                    ssp = c('SSP2', 'SSP3'))

#=================================================================================#
# Extract ----

extract_covars = function(root, allcalcs, rcp, gcm, iam, ssp, out) {
  
  message("[Running: ] ", ssp, ", ", iam, ", ", rcp, ", ", gcm)
  
  # get loggdppc and climtas from mortality
  clim_inc = suppressMessages(read_csv(glue("{root}/{rcp}/{gcm}/{iam}/{ssp}/{allcalcs}"), skip = 21, show_col_types = FALSE))
  clim_inc = clim_inc %>% select(region, 2, loggdppc, climtas) %>% rename(year = `year...2`)
  
  # Save out
  dir.create(out, recursive=T, showWarnings = F)
  
  message(glue("[writing: ] mortality-{rcp}-{gcm}-{iam}-{ssp}-econ_clim.csv"))
  write.csv(clim_inc, glue("{out}/mortality-{rcp}-{gcm}-{iam}-{ssp}-econ_clim.csv"), row.names = F)
  message("---- saved ----\n")
  
}

#=================================================================================#
# Run ----

pwalk(specs, function(rcp, gcm, iam, ssp) {
  extract_covars(root, allcalcs, rcp, gcm, iam, ssp, out)
})
