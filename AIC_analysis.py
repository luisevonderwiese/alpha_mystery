import os
from tabulate import tabulate
import pandas as pd

import raxmlng
import iqtree

results_dir = os.path.join("data", "raxml")
iqtree_dir = os.path.join("data", "iqtree")


metadata_df = pd.read_csv("data/lexibench/character_matrices/stats.tsv", sep = "\t")
datasets = [row["Name"] for _,row in metadata_df.iterrows()]

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
errors = []
for dataset in datasets:
    alpha = raxmlng.alpha(os.path.join(results_dir, dataset, "BIN+G"))
    if alpha <= 90:
        continue
    llh_bing = raxmlng.final_llh(os.path.join(results_dir, dataset, "BIN+G"))
    llh_bin = raxmlng.final_llh(os.path.join(results_dir, dataset, "BIN"))
    errors.append(abs((llh_bin - llh_bing) / llh_bin))
print(max(errors))


def plausible_tree_anaylsis(datasets, reference_model, high_alpha_only):
    results_dir = os.path.join("data", "raxml")
    iqtree_dir = os.path.join("data", "iqtree")

    plausible_count = 0
    identical_count = 0
    non_plausible_count = 0
    for dataset in datasets:
        alpha = raxmlng.alpha(os.path.join(results_dir, dataset, "BIN+G"))
        if alpha <= 90 and high_alpha_only:
            continue

        msa_path = os.path.join("data/lexibench/character_matrices", dataset, "bin.phy")
        ref_tree_path = raxmlng.best_tree_path(os.path.join(results_dir, dataset, reference_model))
        bing_tree_path = raxmlng.best_tree_path(os.path.join(results_dir, dataset, "BIN+G"))
        prefix = os.path.join(iqtree_dir, dataset, reference_model)
        iqtree.run_statstests(msa_path, ref_tree_path, bing_tree_path, prefix)
        plausible, identical = iqtree.is_plausible(prefix)
        if plausible:
            plausible_count +=1
            if identical:
                identical_count += 1
        else:
            non_plausible_count += 1
    
    print("Reference Model", reference_model)
    print("Not plausible: \t\t", str(non_plausible_count))
    print("Plausible: \t\t", str(plausible_count))
    print("Plausible + identical: \t", str(identical_count))


plausible_tree_anaylsis(datasets, "BIN", True)
plausible_tree_anaylsis(datasets, "prob_MULTIxMK+G", True)
