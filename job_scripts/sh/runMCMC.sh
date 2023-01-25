#!/bin/bash

cd /projects/mire0265/mcmc_studies/
CONFIG=/projects/mire0265/mcmc_studies/inputs/configs/config.cfg
RUN=-1

while [[ $# -gt 0 ]]; do
    case $1 in
        -r|--run)
            RUN=$2
            shift
            shift
            ;;
        -*|--*)
            echo "Unknown argument $1"
            exit 1
            ;;
        *)
            echo "Unknown argument $1"
            exit 1
            ;;
    esac
done

if [ "$RUN" -lt 0 ]; then
    echo "Run number not specified!"
    echo "Usage: ./runMCMC.sh -r [RUN]"
    exit 1
fi

echo "Running MCMC..."
./bin/runMCMC -c ${CONFIG} -r ${RUN}

echo ""
echo "MCMC finished. Calculating posterior covariance..."
./bin/calcPosteriorCovmat -c ${CONFIG} -r ${RUN}

echo ""
echo "Calculated posterior covariance. Analyzing results..."
./bin/analyzePosteriorCov -c ${CONFIG} -r ${RUN}

echo "All Done!"

