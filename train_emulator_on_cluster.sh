#!/bin/bash
#SBATCH --job-name=hetgp_training          # Job name
#SBATCH --output=JobOut/output_%j.txt      # Output file (%j will be replaced by job ID)
#SBATCH --error=JobOut/error_%j.txt        # Error file (%j will be replaced by job ID)
#SBATCH --ntasks=1                         # Number of tasks (processes)
#SBATCH --cpus-per-task=14                 # Number of CPU cores per task
#SBATCH --mem=60G                          # Memory per job
#SBATCH --time=06:00:00                    # Max execution time (format HH:MM:SS)
#SBATCH --partition=scicore                # Partition to use
#SBATCH --array=1-8                        # Array for 8 datasets

module load R/4.3.0-foss-2021a             # Load R module

# List of datasets
datasets=(
  "dataset1_EIR8_age_05.rds"
  "dataset2_EIR8_age_all.rds"
  "dataset3_EIR22_age_05.rds"
  "dataset4_EIR22_age_all.rds"
  "dataset5_EIR88_age_05.rds"
  "dataset6_EIR88_age_all.rds"
  "dataset7_EIR168_age_05.rds"
  "dataset8_EIR168_age_all.rds"
)

# Select dataset based on SLURM_ARRAY_TASK_ID
dataset_path="/scicore/home/pothin/GROUP/OM_emulator/datasets_22_07_vary_fut_cov/${datasets[$SLURM_ARRAY_TASK_ID-1]}"

# Define X and Y columns
X_cols="futITNcov2023,futITNcov2024,futITNcov2025,futITNcov2026,futITNcov2027,futIRScov2023,futIRScov2024,futIRScov2025,futIRScov2026,futIRScov2027,access2023,access2024,access2025,access2026,access2027"
Y_cols="nTreatments1_2023,nTreatments1_2024,nTreatments1_2025,nTreatments1_2026,nTreatments1_2027,prevalenceRate_2023,prevalenceRate_2024,prevalenceRate_2025,prevalenceRate_2026,prevalenceRate_2027"

# Training fraction
train_fraction=0.2

# Pass parameters to the R script
Rscript good3.R $dataset_path $X_cols $Y_cols $train_fraction
