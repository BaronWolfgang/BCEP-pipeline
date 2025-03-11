# BCEP-pipeline

The pipeline takes as input fasta protein sequences and/or pdb files, which are passed into various pretrained models and outputs per-residue epitope propensity scores as a csv.

## Set up environment

- Create a new environment:
`conda create -n BCEP`

- Activate the environment:
`conda activate BCEP`

- To install necessary libraries:
`conda install --file conda-requirements.txt -c conda-forge -c bioconda`

`pip install -r pip-requirements.txt`
