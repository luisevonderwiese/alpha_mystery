import os
from tabulate import tabulate
import json
import pandas as pd
import matplotlib.pyplot as plt

import util
import rates
import raxmlng


def scatterplot(dfs, prefix_a, column_a, prefix_b, column_b):
    data_a = dfs[prefix_a][column_a]
    data_b = dfs[prefix_b][column_b]
    plt.scatter(data_a, data_b, s=10)
    plt.xlabel(prefix_a + column_a)
    plt.ylabel(prefix_b + column_b)
    name = prefix_a + column_a + "_vs_" + prefix_b + column_b + ".png"
    plt.savefig(os.path.join("results", "plots", name))
    plt.clf()
    plt.close()

results_dir = os.path.join("results", "raxml")
with open('data/MULTI_X.txt') as f:
    x_values = json.loads(f.read())

bin_models = ["BIN+G", "BIN+FE+G", "BIN+FC+G", "BIN+R4", "BIN+FO+I"]
gamma_models = ["BIN+G", "BIN+FE+G", "BIN+FC+G", "prob_MULTIxMK+G", "prob_BIN+G"]
prefixes = ["", "noinvchar_", "noinvsite_", "equalfreq_"]
large_datasets =  ["abvdoceanic", "bowernpny", "iecor"]
datasets = os.listdir("data/msa")
datasets = [d for d in datasets if d not in large_datasets]

for dataset in datasets:
    for prefix in prefixes:
        msa_path = os.path.join("data/msa/", dataset, prefix + "bin.phy")
        for model in bin_models:
            raxmlng.run_inference(msa_path, model, os.path.join(results_dir, dataset, prefix + model))
    raxmlng.run_inference(os.path.join("data/msa/", dataset, "bin.catg"), "BIN+G", os.path.join(results_dir, dataset, "prob_BIN+G"), "--prob-msa on")
    raxmlng.run_inference(os.path.join("data/msa/", dataset, "multi.catg"), "MULTI" + str(x_values[dataset]) + "_MK+G", os.path.join(results_dir, dataset, "prob_MULTIxMK+G"), "--prob-msa on")

columns = ["dataset", "num_taxa", "num_sites", "entropy_var", "inv_sites_emp", "zero_freq_emp", "inv_sites_estimate", "zero_freq_estimate", "free_rates_var"]
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

        df.at[i, "inv_sites_estimate"] = raxmlng.inv_estimate(os.path.join(results_dir, dataset, prefix + "BIN+FO+I")) * 100
        df.at[i, "zero_freq_estimate"] = raxmlng.zero_freq_estimate(os.path.join(results_dir, dataset, prefix + "BIN+G"))
        var = rates.var(rates.parse_rates(raxmlng.free_rates(os.path.join(results_dir, dataset, prefix + "BIN+R4"))))
        if var > 10:
            var = float("nan")
        df.at[i, "free_rates_var"] = var
        for model in gamma_models:
            if model.startswith("prob_") and prefix != "":
                continue
            raxml_prefix = os.path.join(results_dir, dataset, prefix + model)
            df.at[i, "alpha_" + model] = raxmlng.alpha(raxml_prefix)
    dfs[prefix] = df

if not os.path.isdir(os.path.join("results", "plots")):
    os.makedirs(os.path.join("results", "plots"))
scatterplot(dfs, "", "alpha_BIN+G", "", "alpha_BIN+FE+G")
scatterplot(dfs, "", "alpha_BIN+G", "", "alpha_BIN+FC+G")
scatterplot(dfs, "", "alpha_BIN+FE+G", "", "alpha_BIN+FC+G")
scatterplot(dfs, "", "alpha_BIN+G", "", "alpha_prob_BIN+G")
scatterplot(dfs, "", "alpha_BIN+G", "", "alpha_prob_MULTIxMK+G")
scatterplot(dfs, "", "alpha_BIN+G", "", "num_taxa")
scatterplot(dfs, "", "alpha_BIN+G", "", "num_sites")
scatterplot(dfs, "", "alpha_BIN+G", "", "entropy_var")
scatterplot(dfs, "", "alpha_BIN+G", "", "inv_sites_emp")
scatterplot(dfs, "", "alpha_BIN+G", "", "zero_freq_emp")
scatterplot(dfs, "", "alpha_BIN+G", "", "inv_sites_estimate")
scatterplot(dfs, "", "alpha_BIN+G", "", "free_rates_var")

scatterplot(dfs, "", "zero_freq_emp", "", "zero_freq_estimate")
scatterplot(dfs, "", "alpha_BIN+FE+G", "", "zero_freq_emp")
scatterplot(dfs, "", "alpha_BIN+G", "", "zero_freq_estimate")
scatterplot(dfs, "", "alpha_BIN+FE+G", "", "inv_sites_emp")


scatterplot(dfs, "noinvsite_", "alpha_BIN+G", "noinvsite_", "entropy_var")
scatterplot(dfs, "noinvsite_", "alpha_BIN+FE+G", "noinvsite_", "zero_freq_emp")
scatterplot(dfs, "noinvsite_", "alpha_BIN+G", "noinvsite_", "zero_freq_estimate")
scatterplot(dfs, "noinvsite_", "alpha_BIN+G", "noinvsite_", "zero_freq_emp")
scatterplot(dfs, "", "alpha_BIN+G", "noinvsite_", "alpha_BIN+G")

scatterplot(dfs, "", "alpha_BIN+G", "equalfreq_", "alpha_BIN+G")
scatterplot(dfs, "", "zero_freq_emp", "equalfreq_", "zero_freq_emp")
scatterplot(dfs, "", "zero_freq_estimate", "equalfreq_", "zero_freq_estimate")
scatterplot(dfs, "", "inv_sites_emp", "equalfreq_", "inv_sites_emp")
scatterplot(dfs, "", "inv_sites_estimate", "equalfreq_", "inv_sites_estimate")
