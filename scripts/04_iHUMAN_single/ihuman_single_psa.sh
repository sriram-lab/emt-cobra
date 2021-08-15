#!/bin/bash
#SBATCH --job-name=ihuman_emt_a549_psa
#SBATCH --mail-user=scampit@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --licenses=gurobi@slurmdb:1
#SBATCH --output=/home/scampit/Turbo/scampit/Software/emt-cobra/scripts/log/ihuman_output.log
#SBATCH --error=/home/scampit/Turbo/scampit/Software/emt-cobra/scripts/err/ihuman_error.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=4g
#SBATCH --time=03-00:00:00
#SBATCH --account=lsa1
#SBATCH --partition=standard

module load matlab/R2020a
module load gurobi
matlab -nodisplay -r "addpath(genpath('/home/scampit/Turbo/scampit/Software/emt-cobra/')); 
                      run('/home/scampit/Turbo/scampit/Software/emt-cobra/scripts/04_iHUMAN_single/ihuman_scfba_psa.m'); exit"
