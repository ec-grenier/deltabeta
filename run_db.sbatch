#!/bin/bash
#SBATCH --job-name=mortality-db
#SBATCH --output=/project/cil/home_dirs/egrenier/cil-comms/deltabeta/output/SLURM/mortality-db-%A-%a.out
#SBATCH --account=cil
#SBATCH --partition=caslake
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=0-23

# Time the whole process
start_time=$(date +%s)

# load in required modules
module load netcdf
module load curl/7.79.1
module load udunits
module load geos/3.9.1
module load sqlite/3.36.0
module load proj/7.2
module load gdal/3.3.3
module load gcc/12.2.0 
module load python/anaconda-2022.05
module load R/4.2.1

# navigate to directory
cd /project/cil/home_dirs/egrenier/cil-comms/deltabeta

# define parameter lists
rcps=( rcp45 rcp85 )
iams=( low high )
gcms=( CCSM4 )
ssps=( SSP2 SSP3 )
years=( 2025 2050 2099 )

# compute sizes
num_rcps=${#rcps[@]}
num_iams=${#iams[@]}
num_gcms=${#gcms[@]}
num_ssps=${#ssps[@]}
num_years=${#years[@]}

# total combinations: 2 * 2 * 1 * 2 * 3 = 24
# so --array=0-23

# flatten SLURM_ARRAY_TASK_ID into parameter indices
idx=${SLURM_ARRAY_TASK_ID}

year_index=$(( idx % num_years ))
idx=$(( idx / num_years ))

ssp_index=$(( idx % num_ssps ))
idx=$(( idx / num_ssps ))

gcm_index=$(( idx % num_gcms ))
idx=$(( idx / num_gcms ))

iam_index=$(( idx % num_iams ))
idx=$(( idx / num_iams ))

rcp_index=$(( idx % num_rcps ))

# extract parameter values
rcp=${rcps[$rcp_index]}
iam=${iams[$iam_index]}
gcm=${gcms[$gcm_index]}
ssp=${ssps[$ssp_index]}
year=${years[$year_index]}

echo "Running task: rcp=$rcp, iam=$iam, gcm=$gcm, ssp=$ssp, year=$year"

# Run projection
echo "Running delta beta..."
Rscript 1_mortality_db_script.R "$rcp" "$gcm" "$iam" "$ssp" "$year" 2>&1

# Capture end time after run finishes
end_time=$(date +%s)
duration=$(( end_time - start_time ))
echo "Job completed at $(date)"
echo "Total elapsed time: ${duration} seconds"