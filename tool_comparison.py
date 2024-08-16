import os
from tabulate import tabulate
import json
import pandas as pd
import matplotlib.pyplot as plt

import util
import rates
import raxmlng
import iqtree


def scatterplot(dfs, prefix_a, column_a, prefix_b, column_b):
    data_a = dfs[prefix_a][column_a]
    data_b = dfs[prefix_b][column_b]
    plt.scatter(data_a, data_b, s=10)
    plt.xlabel(prefix_a + column_a)
    plt.ylabel(prefix_b + column_b)
    name = prefix_a + column_a + "_vs_" + prefix_b + column_b + ".png"
    plt.savefig(os.path.join("results", "tool_comparison", name))
    plt.clf()
    plt.close()

raxml_results_dir = os.path.join("results", "raxml")
iqtree_results_dir = os.path.join("results", "iqtree")


raxml_gamma_models = ["BIN+G", "BIN+FE+G", "BIN+FC+G", "prob_MULTIxMK+G", "prob_BIN+G"]
iqtree_gamma_models = ["JC2+G", "GTR2+G"]
prefixes = ["", "noinvchar_", "noinvsite_"]
large_datasets =  ["abvdoceanic", "bowernpny", "iecor"]
datasets = os.listdir("data/msa")
datasets = [d for d in datasets if d not in large_datasets]

columns = ["dataset", "num_taxa", "num_sites", "entropy_var", "inv_sites_emp", "zero_freq_emp", "raxml_inv_sites_estimate", "raxml_free_rates_var", "iqtree_inv_sites_estimate", "iqtree_free_rates_var"]
columns += ["raxml_alpha_" + model for model in raxml_gamma_models]
columns += ["iqtree_alpha_" + model for model in iqtree_gamma_models]
dfs = {}
for prefix in prefixes:
    df = pd.DataFrame(columns = columns)
    for i, dataset in enumerate(datasets):
        df.at[i, "dataset"] = dataset
        msa_path = os.path.join("data/msa/", dataset, prefix + "bin.phy")
        align = util.save_msa_read(msa_path)
        df.at[i, "num_taxa"] = util.num_taxa(align)
        df.at[i, "num_sites"] = util.num_sites(align)
        df.at[i, "entropy_var"] = util.entropy_var(align)
        df.at[i, "inv_sites_emp"] = util.inv_sites(align)
        df.at[i, "zero_freq_emp"] = util.zero_freq(align)

        df.at[i, "raxml_inv_sites_estimate"] = raxmlng.inv_estimate(os.path.join(raxml_results_dir, dataset, prefix + "BIN+FO+I")) * 100
        var = rates.var(rates.parse_rates(raxmlng.free_rates(os.path.join(raxml_results_dir, dataset, prefix + "BIN+R4"))))
        if var > 10:
            var = float("nan")
        df.at[i, "raxml_free_rates_var"] = var
        df.at[i, "iqtree_inv_sites_estimate"] = iqtree.inv_estimate(os.path.join(iqtree_results_dir, dataset, prefix + "GTR2+I")) * 100
        df.at[i, "iqtree_free_rates_var"] = rates.var(rates.parse_rates(iqtree.free_rates(os.path.join(iqtree_results_dir, dataset, prefix + "GTR2+R4"))))
        for model in raxml_gamma_models:
            if model.startswith("prob_") and prefix != "":
                continue
            raxml_prefix = os.path.join(raxml_results_dir, dataset, prefix + model)
            df.at[i, "raxml_alpha_" + model] = raxmlng.alpha(raxml_prefix)
        for model in iqtree_gamma_models:
            iqtree_prefix = os.path.join(iqtree_results_dir, dataset, prefix + model)
            df.at[i, "iqtree_alpha_" + model] = min(iqtree.alpha(iqtree_prefix), 100)
    dfs[prefix] = df

if not os.path.isdir(os.path.join("results", "tool_comparison")):
    os.makedirs(os.path.join("results", "tool_comparison"))
scatterplot(dfs, "", "raxml_alpha_BIN+G", "", "iqtree_alpha_GTR2+G")
scatterplot(dfs, "", "raxml_alpha_BIN+FE+G", "", "iqtree_alpha_JC2+G")
scatterplot(dfs, "", "raxml_inv_sites_estimate", "", "iqtree_inv_sites_estimate")
scatterplot(dfs, "", "raxml_free_rates_var", "", "iqtree_free_rates_var")
