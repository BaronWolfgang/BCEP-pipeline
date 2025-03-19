#!/bin/bash
#SBATCH --job-name=epigraph_test
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=epigraph_test.log

#run your code
echo "Current working directory:"
pwd

echo $SHELL
python3 src/EpiGraph/inference_customPDB.py --pdb "2BIB" --pdb_path "test_data/" --out_path "test_data/test_output"
