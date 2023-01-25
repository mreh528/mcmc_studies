#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:10:00
#SBATCH --job-name=runMCMC
#SBATCH --output=runMCMC.%j.out
#SBATCH --output=runMCMC.%j.err

module purge
module load singularity

SOFTWARE=/projects/mire0265/software
MCMC=/projects/mire0265/mcmc_studies
CONTAINER=${SOFTWARE}/t2kfitterenv_oa2021_r5_centos_7.sif
JOB_SCRIPT=${MCMC}/job_scripts/sh/runMCMC.sh

cd ${MCMC}

OUTFILE=runMCMC.txt
touch ${OUTFILE}
mv ${OUTFILE} ${SLURM_SCRATCH}

SCRATCH_MREH=/scratch/alpine/mire0265/outputs
singularity run -B ${SCRATCH_MREH} ${CONTAINER} ${JOB_SCRIPT} > ${SLURM_SCRATCH}/${OUTFILE} 2>&1

echo "Done! Moving outputs..."
mv ${SLURM_SCRATCH}/${OUTFILE} ${MCMC}/output/${OUTFILE}
mv ${MCMC}/job_scripts/sh/*.out ${MCMC}/job_scripts/out
mv ${MCMC}/job_scripts/sh/*.err ${MCMC}/job_scripts/err

