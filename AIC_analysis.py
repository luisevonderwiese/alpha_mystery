import os
from tabulate import tabulate

import raxmlng

results_dir = os.path.join("results", "raxml")

large_datasets =  ["abvdoceanic", "bowernpny", "iecor"]
datasets = os.listdir("data/msa")
datasets = [d for d in datasets if d not in large_datasets]

low_het_bin = 0
low_het_bing = 0
high_het_bin = 0
high_het_bing = 0

for dataset in datasets:
    aic_bing = raxmlng.aic(os.path.join(results_dir, dataset, "BIN+G"))
    aic_bin = raxmlng.aic(os.path.join(results_dir, dataset, "BIN"))
    alpha = raxmlng.alpha(os.path.join(results_dir, dataset, "BIN+G"))
    if alpha < 90:
        if aic_bin < aic_bing:
            high_het_bin += 1
        else:
            high_het_bing +=1
    else:
        if aic_bin < aic_bing:
            low_het_bin += 1
        else:
            low_het_bing +=1

res = []
res.append(["BIN", high_het_bin, low_het_bin])
res.append(["BIN+G", high_het_bing, low_het_bing])
print(tabulate(res, tablefmt="pipe", headers=["model with lower AIC score", "alpha < 90", "alpha > 90"]))
