#!/bin/bash

# to run use alias: screen -dmS matilda bash -c './matilda.sh'
# Get the directory of the current script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Change to the script's directory
cd "$SCRIPT_DIR"

# Activate the conda environment
source ~/anaconda3/etc/profile.d/conda.sh  # Adjust this path if your conda is installed elsewhere
conda activate matilda

# Run the Python script
python processFlyscans.py
