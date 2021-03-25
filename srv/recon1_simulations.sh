#!/bin/bash
#SBATCH --job-name=recon1_emt_a549_1
#SBATCH --mail-user=scampit@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=./recon1_output.log
#SBATCH --error=./recon1_error.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=4g
#SBATCH --time=03-00:00:00
#SBATCH --account=lsa1
#SBATCH --partition=standard

module load matlab/R2020a
matlab -nodisplay -nodesktop -r "./recon1_scCOBRA ; exit"
