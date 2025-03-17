#!/bin/bash
#SBATCH --job-name=test_netsurfp
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=5gb
#SBATCH --time=00:10:00
#SBATCH --output=test_netsurfp.log

# Print current working directory
echo "Current working directory:"
pwd

# Load necessary modules (if required)
# module load netsurfp  # If you use modules on your system

# Run NetSurfP test with srun
netsurfp -i test1.fsa -a > test1.myrsa

# Output the result
echo "Test completed, check test1.myrsa for results."
