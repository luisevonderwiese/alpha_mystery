# Cognate Model
Project for the evaluation of a model for cognate data implemented in RAxML-NG

## Requirements:
1. Setup and activate the conda environment: 
```
conda env create -f environment.yml
conda acitvate alpha-mystery
```
2. Install iqtree
Install iqtree following the instructions on their [website](http://www.iqtree.org/doc/Quickstart).
Place the binary `iqtree2` in `bin/

3. Data
Create lexibench repos
```
lexibench --repos data/lexibench download
lexibench --repos data/lexibench lingpy_wordlists
lexibench --repos data/lexibench character_matrices --format bin.phy bin.catg multi.catg
```

## Execution
```
python remove_invariant_sites.py
python experiment.py
python AIC_analysis.py
```
