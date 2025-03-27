# BCEP-pipeline

The pipeline takes as input fasta protein sequences and/or pdb files, which are passed into various pretrained models and outputs per-residue epitope propensity scores as a csv.

## Set up environment

- Create a new environment:
`conda create -n BCEP`

- Activate the environment:
`conda activate BCEP`

### Discotope-3.0
Discotope-3.0 requires cloning their repository into into the src folder:
```
cd src
git clone https://github.com/Magnushhoie/discotope3_web/

# Unzip models to use
unzip models.zip
```

- To install necessary libraries:
```
conda install --file conda-requirements.txt -c conda-forge -c bioconda
pip install -r pip-requirements.txt
```