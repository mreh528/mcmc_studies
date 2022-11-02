#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:01:00
#SBATCH --job-name=runMCMC
#SBATCH --output=runMCMC.%j.out
#SBATCH --output=runMCMC.%j.err

module purge
cd /projects/mire0265/mcmc_studies/
source setup_env.sh

./bin/mcmcPerformance -c inputs/configs/config.cfg

echo "All Done!"

