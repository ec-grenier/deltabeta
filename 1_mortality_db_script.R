#==============================================================================#
#' Mortality Delta Beta Script
#'
#' to run: Rscript 1_mortality_db_script.R <rcp> <gcm> <iam> <ssp> <year>   
#==============================================================================#

#==============================================================================#
# 0. Source, set packages and paths -----------

library(glue)

REPO = '/project/cil/home_dirs/egrenier/repos'
source(glue("{REPO}/../cil-comms/deltabeta/yellow_purple_package.R"))
source(glue("{REPO}/../cil-comms/deltabeta/2_mortality_db_wrapper.R"))


# Set paths below
csvv_dir = glue("{REPO}/../cil-comms/deltabeta/data/csvv/") # Should have trailing slash
csvv_name = "Agespec_interaction_response.csvv"
cov_dir = glue("{REPO}/../cil-comms/deltabeta/data/covars")
out_dir = glue("{REPO}/../cil-comms/deltabeta/output/db_out")

# get command line args
args = commandArgs(trailingOnly = TRUE)

#==============================================================================#
# 1. Set globals ----

# to map age group coefs from csvv
age_list = c(3,2,1)
het_list = list(age = c(rep.int(1, 12),rep.int(2, 12),rep.int(3, 12))) 

# read CLI args
rcp = args[1]
gcm = args[2]
iam = args[3]
ssp = args[4]
year = as.numeric(args[5])

# rcp='rcp45'
# gcm='CCSM4'
# year=2099
# ssp = "SSP2"
# iam = "low"

slug=''

cov_dir = glue("{cov_dir}/mortality-{rcp}-{gcm}-{iam}-{ssp}-econ_clim.csv") # this needs to be full .csv path

# Get list of regions to loop over 
region_list = list(unique((fread(cov_dir))$region))

# Testing with subsets or specific regions
# region_list = region_list[[1]][22000:22005] # any range in [1,24378]
# region_list = list('TCD.7.21.77')

#-----set args here-----
args = list(years=year,
            base_year=1993,
            rebase_year=2005,
            csvv.dir=csvv_dir,
            csvv.name=csvv_name,
            het.list=het_list,
            cov.dir=cov_dir,
            covarkey='region',
            list.names=c('climtas','loggdppc'),
            covar.names=c('climtas','loggdppc'),
            func=get.clipped.curve.mortality,
            get.covars=T,
            tas_value="tas",
            ncname="1.6",
            TT_upper_bound=70,
            TT_lower_bound=-70,
            TT_step=1,
            do.clipping=T,
            goodmoney.clipping=T, 
            do.diffclip=T,
            full_db=T, # this should be F
            do_global=F)

#==============================================================================#
# 2. Run Function ----

# create output directory 
out = glue("{out_dir}/{rcp}/{gcm}/{iam}/{ssp}")
dir.create(out, recursive = T, showWarnings = F)
message(glue("[Running delta beta: ] {rcp} {gcm} {iam} {ssp}\n"))

# loop through age groups sequentially
for (age in age_list){
  
  agegroup = ifelse(age==3, "oldest", ifelse(age==2, "older", "young"))
  
  message(glue("[age group: ] {agegroup} \n"))
  
  df = lapply(region_list, get_all_db_tables, age, args)
  df = do.call(rbind, df) %>% as.data.frame() %>% select(hierid, everything()) %>% dplyr::rename(region = hierid)
  
  file_name = glue("mortality-delta_beta-{year}-{rcp}-{gcm}-{iam}-{ssp}-{agegroup}{slug}.csv")
  output = glue("{out}/{file_name}")
  
  message('\n[ saving ]\n')
  message(glue("Saving: {file_name}\nOutput directory: {out}\n"))
  write.csv(df, output, row.names=FALSE)
  message("\n[ saved ]\n")
}

message(glue("[ all age groups completed ]"))
