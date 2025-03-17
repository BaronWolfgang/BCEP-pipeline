#!/bin/bash
#SBATCH --job-name=bp3_test
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=bp3_test.log

#run your code
echo "Current working directory:"
pwd

echo $SHELL
python3 src/bp3/bepipred3_CLI.py -i "test_data/pdb_protein_sequences_part_2.fasta" -o "test_data/test_output" -pred "vt_pred" 
