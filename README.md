# BCEP-pipeline

The pipeline takes as input fasta protein sequences and/or pdb files, which are passed into various pretrained models and outputs per-residue epitope propensity scores as a csv.

## Set up environment

- Create a new environment:
```
conda env create -f environment.yml -n BCEP
```
This will create the environment BCEP and install all dependencies in `environment.yml`

### Discotope-3.0
Discotope needs to be setup independently:
```
cd src/discotope3_web/
pip install .

# Unzip models to use
unzip models.zip
```
