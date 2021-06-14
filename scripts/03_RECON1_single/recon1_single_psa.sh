#!/bin/bash
#SBATCH --job-name=recon1_emt_a549_psa
#SBATCH --mail-user=scampit@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --licenses=gurobi@slurmdb:8
#SBATCH --output=/home/scampit/Turbo/scampit/Software/emt-cobra/scripts/log/recon1_output.log
#SBATCH --error=/home/scampit/Turbo/scampit/Software/emt-cobra/scripts/err/recon1_error.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=4g
#SBATCH --time=03-00:00:00
#SBATCH --account=lsa1
#SBATCH --partition=standard

module load matlab/R2020a
module load gurobi
matlab -nodisplay -r "addpath(genpath('/home/scampit/Turbo/scampit/Software/emt-cobra/)); 
                      run('/home/scampit/Turbo/scampit/Software/emt-cobra/scripts/03_RECON1_single/recon1_scfba_psa.m'); exit"