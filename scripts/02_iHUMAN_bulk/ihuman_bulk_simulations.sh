#!/bin/bash
#SBATCH --job-name=ihuman_emt_a549_1
#SBATCH --mail-user=scampit@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --licenses=gurobi@slurmdb:8
#SBATCH --output=./ihuman_bulk_output.log
#SBATCH --error=./ihuman_bulk_error.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1g
#SBATCH --cpus-per-task=16
#SBATCH --time=00-05:00:00
#SBATCH --account=lsa1
#SBATCH --partition=standard

module load matlab/R2020a
module load gurobi
matlab -nodisplay -r "run('/home/scampit/Turbo/scampit/Software/emt/srv/ihuman_simulations.m'); exit"
