#!/bin/bash
#
# Create the dyna_clust_predict conda environment.
# Must be submitted from the project root directory.

ENV_FILE="environment.yml"
ENV_NAME="dyna_clust_predict"

if ! conda info --envs | grep -q "${ENV_NAME}"; then
  echo "Creating conda environment '${ENV_NAME}' at: $(date)"
  mamba env create -f "${ENV_FILE}"
else
  echo "Environment '${ENV_NAME}' already exists."
fi

