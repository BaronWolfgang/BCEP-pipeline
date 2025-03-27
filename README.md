# BCEP-pipeline

The pipeline takes as input fasta protein sequences and/or pdb files, which are passed into various pretrained models and outputs per-residue epitope propensity scores as a csv.

## Set up environment

- Create a new environment:
`conda create -n BCEP`

- Activate the environment:
`conda activate BCEP`

### BepiPred-3.0
Discotope-3.0 requires cloning their repository into into the src folder:
```
cd src
git clone https://github.com/UberClifford/BepiPred3.0-Predictor
```

### Discotope-3.0
Discotope-3.0 requires cloning their repository into into the src folder:
```
cd src
git clone https://github.com/Magnushhoie/discotope3_web/

# Unzip models to use
unzip models.zip
```

## To install necessary libraries:
```
conda install --file conda-requirements.txt -c conda-forge -c bioconda
pip install -r pip-requirements.txt
```

## Run the pipeline:
Run the script `bcep.py` with the following arguments:
- `--tools` select one or more tools to run from (`discotope3`, `bepipred3`)
Depending on the tool(s) selected different inputs are required:
-- `--pdb` file path to a .pdb file (bepipred3)
-- `--fasta` file path to .fasta file (discotope3) 
- `--out_dir` folder name to store output