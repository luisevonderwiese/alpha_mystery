# Cognate Model
Project for the evaluation of a model for cognate data implemented in RAxML-NG

## Requirements:
1. Setup and activate the conda environment: 
```
conda env create -f environment.yml
conda acitvate alpha-mystery
```
2. Data:
Lexibench
Ground Truth Difficulties

## Execution
```
python remove_invariant_sites.py
python experiment.py
python AIC_analysis.py
```
