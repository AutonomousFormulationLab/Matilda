#!/bin/bash
echo $(CONDA_PREFIX:-"/APSshare/miniconda/x86_64")
exit 0


# to run use alias: screen -dmS matilda bash -c './matilda.sh'
# Get the directory of the current script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Change to the script's directory
cd "$SCRIPT_DIR"

#MY_DIR=$(realpath "$(dirname $0)")
#LOG_FILE="${MY_DIR}/logfile.txt"


source ${CONDA_PREFIX}/etc/profile.d/conda.sh
CONDA_ENV=matilda
conda activate "${CONDA_ENV}"

# Run the Python script
python matilda.py
