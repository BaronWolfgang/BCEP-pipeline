#!/bin/bash
#SBATCH --job-name=iprscan_test
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=logs/iprscan_test.log

#run your code
echo "Job started at:"
date

echo "Current working directory:"
pwd

echo $SHELL
python3 scripts/iprscan.py --fasta temp/bepipred3/fasta/2BIB.fasta --email "jravilab.msu@gmail.com"

echo "Job finished at:"
date