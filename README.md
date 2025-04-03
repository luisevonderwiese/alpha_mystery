# Cognate Model
Project for the evaluation of a model for cognate data implemented in RAxML-NG

## Requirements:
1. Setup and activate the conda environment: 
```
conda env create python=3.9 -f environment.yml
conda acitvate alpha-mystery
```
2. Install iqtree
Install iqtree following the instructions on their [website](http://www.iqtree.org/doc/Quickstart).
Place the binary `iqtree2` in `bin/

3. Create Lexibench Data 
Clone the [glottolog repo](https://github.com/glottolog/glottolog) to a directory of your choice, then run:
```
lexibench --repos data/lexibench download
lexibench --repos data/lexibench lingpy_wordlists
lexibench --repos data/lexibench character_matrices --formats bin.phy bin.catg multi.catg --glottolog <your_glottolog_path>
```

## Execution
```
python remove_invariant_sites.py
python experiment.py
python AIC_analysis.py
```
