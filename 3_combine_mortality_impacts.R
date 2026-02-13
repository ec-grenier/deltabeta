#==============================================================================#
#' 
#' 
#' Combine mortality delta beta impacts
#' 
#' 
#==============================================================================#

#==============================================================================#
# 0. Source, set packages and paths -----------

packages = c("glue", "tidyverse", "data.table", "reticulate", "scales", "ncdf4")
invisible(lapply(packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))
rm(packages)

REPO = '/project/cil/home_dirs/egrenier/cil-comms/deltabeta'
source(glue('{REPO}/../../repos/regional-scc/utils/get_econ_vars.R'))
source(glue('{REPO}/../../repos/regional-scc/utils/regionutils.R'))

#==============================================================================#
# 1. Define inputs -----------

# ---- set specs---- # 
# available scenarios:
#  - rcp:      rcp45, rcp85
#  - gcm:      CCSM4
#  - ssp:      SSP2, SSP3
#  - year:     2025, 2050, 2099
specs = expand.grid(rcp = c('rcp45', 'rcp85'), 
                    gcm = c('CCSM4'),
                    iam = c('low', 'high'),
                    ssp = c('SSP2', 'SSP3'),
                    year = c(2025, 2050, 2099)) 

#==============================================================================#
# 1. Combine impacts -----------

combine_impacts = function(rcp, gcm, iam, ssp, year){
  
  message(glue('[running: ] {year}, {rcp}, {gcm}, {iam}, {ssp}'))
  
  pop = get_econvar(units = c('pop', 'pop0to4', 'pop5to64', 'pop65plus'),
                    regions = 'all',
                    iam = iam,
                    ssp = ssp,
                    year_list = c(year))
  
  # get country-level pop share
  pop$iso = substr(pop$region, 1, 3)
  pop$young_share = pop$pop0to4 / pop$pop
  pop$older_share = pop$pop5to64 / pop$pop
  pop$oldest_share = pop$pop65plus / pop$pop
  
  # some IRs have 0 pop yielding NA pop shares. apply iso pop shares to all IRs in a country
  skip_isos = c("ATA", "ATF", "BVT", "CL-", "HMD", "IOT", "SGS", "SP-") # <- unless these exceptions
  pop = pop %>%
    group_by(iso) %>%
    mutate(
      young_share = ifelse(iso %in% skip_isos, young_share, first(na.omit(young_share))),
      older_share = ifelse(iso %in% skip_isos, older_share, first(na.omit(older_share))),
      oldest_share = ifelse(iso %in% skip_isos, oldest_share, first(na.omit(oldest_share)))
    ) %>%
    ungroup()
  
  # Loop through agegroups to obtain population weighted effects, and curves.
  # does a bin-level pop-share weighted average
  result = data.frame()
  ages = c("oldest","older","young")
  for(age in ages){  
    
    db = fread(glue('{REPO}/output/db_out/{rcp}/{gcm}/{iam}/{ssp}/mortality-delta_beta-{year}-{rcp}-{gcm}-{iam}-{ssp}-{age}.csv'))
    
    # get the matching share column
    share = pop[[paste0(age, "_share")]]
    names(share) = pop$region  # so we can match by region
    
    # multiply the relevant columns by the population share
    cols_to_weight = c("beta_fa", "beta_ia", "beta_na", "effect_fa", "effect_ia", "effect_na")
    db = db %>% 
      left_join(pop %>% select(region, !!paste0(age, "_share")), by = "region") %>%
      mutate(across(all_of(cols_to_weight), ~ .x * get(paste0(age, "_share")))) %>%
      select(-paste0(age, "_share"))
    
    if(nrow(result) == 0){
      result = db
    } else {
      result = result %>%
        left_join(db, by = c("region", "bin", glue("T_{year}"), "T_1993", "T_2005", "T_diff"), suffix = c("", paste0("_", age))) %>%
        mutate(across(all_of(cols_to_weight), ~ .x + get(paste0(cur_column(), "_", age)))) %>%
        select(-ends_with(paste0("_", age)))
    }
    
  }
  message('\n==== saving ====\n')
  output = glue('{REPO}/output/db_out/{rcp}/{gcm}/{iam}/{ssp}/mortality-delta_beta-{year}-{rcp}-{gcm}-{iam}-{ssp}-combined.csv')
  message(glue('Output: {output}'))
  write.csv(result, output, row.names=F)
  message('\n==== saved ====\n')
}

#==============================================================================#
# 2. Call func -----------

pwalk(specs, function(year, rcp, gcm, iam, ssp) {

  combine_impacts(rcp, gcm, iam, ssp, year)
  
})