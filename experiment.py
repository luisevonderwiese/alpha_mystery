import os
from tabulate import tabulate
import json
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
from scipy.stats import somersd, rankdata


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

def to_binary(alpha_values):
    threshold = 90
    binary = []
    for alpha in alpha_values:
        if alpha > threshold:
            binary.append(1)
        else:
            binary.append(0)
    return binary

def get_auc_score(dfs, prefix_a, column_a, prefix_b, column_b):
    data_a = dfs[prefix_a][column_a]
    data_b = dfs[prefix_b][column_b]
    binary = to_binary(data_a)
    data_b_filtered = []
    binary_filtered = []
    for i, value in enumerate(data_b):
        if value == value:
            data_b_filtered.append(value)
            binary_filtered.append(binary[i])
    zero_values = [data_b_filtered[i] for i, val in enumerate(binary_filtered) if val == 0]
    one_values = [data_b_filtered[i] for i, val in enumerate(binary_filtered) if val == 1]
    auc = roc_auc_score(binary_filtered, data_b_filtered)
    mult = 1
    if auc < 0.5:
        auc =  roc_auc_score(binary_filtered, [-el for el in  data_b_filtered])
        mult = -1
    return [mult * (2 * auc - 1), sum(zero_values) / len(zero_values), sum(one_values) / len(one_values)] 

#Somers’ D = 2 * AUC – 1
def get_somersd(dfs, prefix_a, column_a, prefix_b, column_b):
    data_a = dfs[prefix_a][column_a]
    data_b = dfs[prefix_b][column_b]
    binary = to_binary(data_a)
    data_b_filtered = []
    binary_filtered = []
    for i, value in enumerate(data_b):
        if value == value:
            data_b_filtered.append(value)
            binary_filtered.append(binary[i])

    s = somersd(binary_filtered, data_b_filtered)
    return (s.statistic, s.pvalue)

def statistical_analysis(prefix, alpha_column, other_columns):
    print(prefix, alpha_column)
    res = []
    for column in other_columns:
        r = [column]
        a = get_auc_score(dfs, prefix, alpha_column, prefix, column)
        r.append(round(a[0], 3))
        r.append(round(a[1], 3))
        r.append(round(a[2], 3))
        #d = get_somersd(dfs, prefix, alpha_column, prefix, column)
        #r.append(d[0])
        #r.append(d[1])
        res.append(r)
    headers = ["column", "AUC", "low alpha mean", "high alpha mean"]
           # "Somers' D", "Somers' D - pvalue"]
    print(tabulate(res, tablefmt="latex", headers=headers))


results_dir = os.path.join("results", "raxml")
metadata_df = pd.read_csv("data/metadata.tsv", sep = "\t")
x_values = {}
for i, row in metadata_df.iterrows():
    x_values[row["dataset"]] = row["max_num_cognate_classes"]

#bin_models = ["BIN", "BIN+G", "BIN+FE+G", "BIN+FC+G", "BIN+R4", "BIN+FO+I"]
#gamma_models = ["BIN+G", "BIN+FE+G", "BIN+FC+G", "prob_MULTIxMK+G", "prob_BIN+G"]
bin_models = ["BIN", "BIN+G", "BIN+R4", "BIN+FO+I"]
gamma_models = ["BIN+G", "prob_MULTIxMK+G", "prob_BIN+G"]

prefixes = [""]
large_datasets =  ["abvdoceanic", "bowernpny", "iecor"]
datasets = os.listdir("data/msa")
datasets = [d for d in datasets if d not in large_datasets]

for dataset in datasets:
    for prefix in prefixes:
        msa_path = os.path.join("data/msa/", dataset, prefix + "bin.phy")
        if not os.path.isfile(msa_path):
            continue
        for model in bin_models:
            raxmlng.run_inference(msa_path, model, os.path.join(results_dir, dataset, prefix + model))
    raxmlng.run_inference(os.path.join("data/msa/", dataset, "bin.catg"), "BIN+G", os.path.join(results_dir, dataset, "prob_BIN+G"), "--prob-msa on")
    raxmlng.run_inference(os.path.join("data/msa/", dataset, "multi.catg"), "MULTI" + str(x_values[dataset]) + "_MK+G", os.path.join(results_dir, dataset, "prob_MULTIxMK+G"), "--prob-msa on")
        
difficult_df = pd.read_csv("data/results_lexibench_difficult.csv")

columns = [
        "dataset", 
        "num_sites", 
        "entropy_var", 
        "inv_sites_emp", 
        "zero_freq_emp", 
        "inv_sites_estimate", 
        "free_rates_var", 
        "inv_sites_error",
        ]
columns += ["alpha_" + model for model in gamma_models]
dfs = {}
for prefix in prefixes:
    df = pd.DataFrame(columns = columns)
    for i, dataset in enumerate(datasets):
        df.at[i, "dataset"] = dataset
        msa_path = os.path.join("data/msa/", dataset, prefix + "bin.phy")
        try:
            align = util.save_msa_read(msa_path)
        except:
            continue
        df.at[i, "num_sites"] = util.num_sites(align)
        df.at[i, "entropy_var"] = util.entropy_var(align)
        df.at[i, "inv_sites_emp"] = util.inv_sites(align)
        df.at[i, "zero_freq_emp"] = util.zero_freq(align)

        df.at[i, "inv_sites_estimate"] = raxmlng.inv_estimate(os.path.join(results_dir, dataset, prefix + "BIN+FO+I")) * 100
        var = rates.var(rates.parse_rates(raxmlng.free_rates(os.path.join(results_dir, dataset, prefix + "BIN+R4"))))
        if var > 10:
            var = float("nan")
        df.at[i, "free_rates_var"] = var
        df.at[i, "inv_sites_error"] = abs(df.at[i, "inv_sites_estimate"] - df.at[i, "inv_sites_emp"])

        for model in gamma_models:
            if model.startswith("prob_") and prefix != "":
                continue
            raxml_prefix = os.path.join(results_dir, dataset, prefix + model)
            df.at[i, "alpha_" + model] = raxmlng.alpha(raxml_prefix)
    df = df.merge(metadata_df, on = "dataset")
    df = df.merge(difficult_df, on = "dataset")
    dfs[prefix] = df

if not os.path.isdir(os.path.join("results", "plots")):
    os.makedirs(os.path.join("results", "plots"))

scatterplot(dfs, "", "alpha_BIN+G", "", "alpha_prob_MULTIxMK+G")







statistical_analysis("", "alpha_BIN+G",
        [
            "entropy_var",
            "num_languages",
            "num_sites",
            "num_concepts",
            "columns_per_concept",
            "difficult",
            "free_rates_var",
            "inv_sites_emp",
            "inv_sites_estimate",
            "inv_sites_error",
            #"zero_freq_emp",
            #"bin_entropy",
            ])


res_low = []
res_high = []
for i, row in dfs[""].iterrows():
    r = [row["dataset"], row["num_languages"], row["columns_per_concept"], row["alpha_BIN+G"]]
    if row["alpha_BIN+G"] > 90:
        res_high.append(r)
    else:
        res_low.append(r)

headers = ["dataset", "num_languages", "columns_per_concept", "alpha"]
print("low alpha")
print(tabulate(res_low, tablefmt="pipe", headers=headers))
print("high alpha")
print(tabulate(res_high, tablefmt="pipe", headers=headers))

print(max(dfs[""]["alpha_prob_MULTIxMK+G"]))
df = dfs[""]
print(max(df[df["alpha_BIN+G"] <= 90]["alpha_BIN+G"]))
