import os
from tabulate import tabulate
import json
import pandas as pd
import matplotlib.pyplot as plt

import util
import rates
import iqtree


def scatterplot(dfs, prefix_a, column_a, prefix_b, column_b):
    data_a = dfs[prefix_a][column_a]
    data_b = dfs[prefix_b][column_b]
    plt.scatter(data_a, data_b, s=10)
    plt.xlabel(prefix_a + column_a)
    plt.ylabel(prefix_b + column_b)
    name = prefix_a + column_a + "_vs_" + prefix_b + column_b + ".png"
    plt.savefig(os.path.join("results", "plots_iqtree", name))
    plt.clf()
    plt.close()

results_dir = os.path.join("results", "iqtree")

bin_models = ["JC2+G", "GTR2+G", "GTR2+R4", "GTR2+I"]
gamma_models = ["JC2+G", "GTR2+G"]
prefixes = ["", "noinvchar_", "noinvsite_"]
large_datasets =  ["abvdoceanic", "bowernpny", "iecor"]
datasets = os.listdir("data/msa")
datasets = [d for d in datasets if d not in large_datasets]

for dataset in datasets:
    for prefix in prefixes:
        msa_path = os.path.join("data/msa/", dataset, prefix + "bin.phy")
        for model in bin_models:
            iqtree.run_inference(msa_path, model, os.path.join(results_dir, dataset, prefix + model))


columns = ["dataset", "num_taxa", "num_sites", "entropy_var", "inv_sites_emp", "zero_freq_emp", "inv_sites_estimate", "free_rates_var"]
columns += ["alpha_" + model for model in gamma_models]
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

        df.at[i, "inv_sites_estimate"] = iqtree.inv_estimate(os.path.join(results_dir, dataset, prefix + "GTR2+I")) * 100
        df.at[i, "free_rates_var"] = rates.var(rates.parse_rates(iqtree.free_rates(os.path.join(results_dir, dataset, prefix + "GTR2+R4"))))
        for model in gamma_models:
            iqtree_prefix = os.path.join(results_dir, dataset, prefix + model)
            df.at[i, "alpha_" + model] = iqtree.alpha(iqtree_prefix)
    dfs[prefix] = df

if not os.path.isdir(os.path.join("results", "plots_iqtree")):
    os.makedirs(os.path.join("results", "plots_iqtree"))
scatterplot(dfs, "", "alpha_GTR2+G", "", "alpha_JC2+G")
scatterplot(dfs, "", "alpha_GTR2+G", "", "num_taxa")
scatterplot(dfs, "", "alpha_GTR2+G", "", "num_sites")
scatterplot(dfs, "", "alpha_GTR2+G", "", "entropy_var")
scatterplot(dfs, "", "alpha_GTR2+G", "", "inv_sites_emp")
scatterplot(dfs, "", "alpha_GTR2+G", "", "zero_freq_emp")
scatterplot(dfs, "", "alpha_GTR2+G", "", "inv_sites_estimate")
scatterplot(dfs, "", "alpha_GTR2+G", "", "free_rates_var")

scatterplot(dfs, "noinvsite_", "alpha_GTR2+G", "noinvsite_", "entropy_var")
scatterplot(dfs, "", "alpha_GTR2+G", "noinvsite_", "alpha_GTR2+G")
