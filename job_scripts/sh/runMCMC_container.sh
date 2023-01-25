#!/bin/bash

MCMC=/projects/mire0265/mcmc_studies
cd ${MCMC}
source setup_env.sh

CONFIG=${MCMC}/inputs/configs/config.cfg

echo "Running MCMC..."
./bin/runMCMC -c ${CONFIG}

echo ""
echo "MCMC finished. Calculating posterior covariance..."
./bin/calcPosteriorCovmat -c ${CONFIG}

echo ""
echo "Calculated posterior covariance. Analyzing results..."
./bin/analyzePosteriorCov -c ${CONFIG}


