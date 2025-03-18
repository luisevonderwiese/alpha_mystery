import os
from tabulate import tabulate
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score

import util
import rates
import raxmlng
import label_wrapper


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


def statistical_analysis(prefix, alpha_column, other_columns):
    print(prefix, alpha_column)
    res = []
    for column in other_columns:
        r = [column]
        a = get_auc_score(dfs, prefix, alpha_column, prefix, column)
        r.append(round(a[0], 3))
        r.append(round(a[1], 3))
        r.append(round(a[2], 3))
        res.append(r)
    headers = ["column", "AUC", "low alpha mean", "high alpha mean"]
    print(tabulate(res, tablefmt="latex", headers=headers))


results_dir = os.path.join("results", "raxml")
label_dir = os.path.join("results", "difficulty_labels")
if not os.path.isdir(label_dir):
    os.makedirs(label_dir)
metadata_df = pd.read_csv("data/lexibench/character_matrices/stats.tsv", sep = "\t")
wl_df = pd.read_csv("data/lexibench/lingpy_wordlists/stats.tsv", sep = "\t")
for i, row in wl_df.iterrows():
    wl_df.at[i, "dataset"] = row["Dataset"] + "-" + row["Family"].lower()
metadata_df = metadata_df.merge(wl_df, on = "dataset") 
datasets = [row["dataset"] for _,row in metadata_df.iterrows()]
datasets = [d for d in datasets if d != "abvdoceanic-austronesian"]
x_values = {}
for i, row in metadata_df.iterrows():
    x_values[row["wordlist"].split("/")[-1].split(".")[0]] = row["cs_max"]

#bin_models = ["BIN", "BIN+G", "BIN+FE+G", "BIN+FC+G", "BIN+R4", "BIN+FO+I"]
#gamma_models = ["BIN+G", "BIN+FE+G", "BIN+FC+G", "prob_MULTIxMK+G", "prob_BIN+G"]
bin_models = ["BIN", "BIN+G", "BIN+R4", "BIN+FO+I"]
gamma_models = ["BIN+G", "prob_MULTIxMK+G", "prob_BIN+G"]

prefixes = [""]

for dataset in datasets:
    for prefix in prefixes:
        msa_path = os.path.join("data/lexibench/character_matrices/", dataset, prefix + "bin.phy")
        if not os.path.isfile(msa_path):
            continue
        label_wrapper.calculate_label(msa_path, os.path.join(label_dir, prefix + dataset, "label"))
        for model in bin_models:
            raxmlng.run_inference(msa_path, model, os.path.join(results_dir, dataset, prefix + model))
    raxmlng.run_inference(os.path.join("data/lexibench/character_matrices/", dataset, "bin.catg"), "BIN+G", os.path.join(results_dir, dataset, "prob_BIN+G"), "--prob-msa on")
    if x_values[dataset] <= 64:
        raxmlng.run_inference(os.path.join("data/lexibench/character_matrices/", dataset, "multi.catg"), "MULTI" + str(x_values[dataset]) + "_MK+G", os.path.join(results_dir, dataset, "prob_MULTIxMK+G"), "--prob-msa on")
        

columns = [
        "dataset", 
        "num_sites", 
        "entropy_var", 
        "inv_sites_emp", 
        "zero_freq_emp", 
        "inv_sites_estimate", 
        "free_rates_var", 
        "inv_sites_error",
        "difficult"
        ]
columns += ["alpha_" + model for model in gamma_models]
dfs = {}
for prefix in prefixes:
    df = pd.DataFrame(columns = columns)
    for i, dataset in enumerate(datasets):
        print(dataset)
        df.at[i, "dataset"] = dataset
        msa_path = os.path.join("data/lexibench/character_matrices", dataset, prefix + "bin.phy")
        try:
            align = util.save_msa_read(msa_path)
        except:
            continue
        print("reading done")
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
        try: 
            df.at[i, "difficult"] = label_wrapper.get_label(os.path.join(label_dir, prefix + dataset, "label"))
        except Exception as e:
            print(e)
        for model in gamma_models:
            if model.startswith("prob_") and prefix != "":
                continue
            raxml_prefix = os.path.join(results_dir, dataset, prefix + model)
            df.at[i, "alpha_" + model] = raxmlng.alpha(raxml_prefix)
    df = df.merge(metadata_df, on = "dataset")
    dfs[prefix] = df

if not os.path.isdir(os.path.join("results", "plots")):
    os.makedirs(os.path.join("results", "plots"))

scatterplot(dfs, "", "alpha_BIN+G", "", "alpha_prob_MULTIxMK+G")

print(dfs[""])





statistical_analysis("", "alpha_BIN+G",
        [
            "entropy_var",
            "Languages",
            "num_sites",
            "Concepts",
            "cs_mean",
            "difficult",
            "free_rates_var",
            "inv_sites_emp",
            "inv_sites_estimate",
            "inv_sites_error",
            "zero_freq_emp",
            #"bin_entropy",
            ])



res_low = []
res_high = []
for i, row in dfs[""].iterrows():
    r = [row["dataset"], row["Languages"], row["cs_mean"], row["alpha_BIN+G"]]
    if row["alpha_BIN+G"] > 90:
        res_high.append(r)
    else:
        res_low.append(r)

headers = ["dataset", "Languages", "cs_mean", "alpha"]
print("low alpha")
print(tabulate(res_low, tablefmt="pipe", headers=headers))
print("high alpha")
print(tabulate(res_high, tablefmt="pipe", headers=headers))

print(max(dfs[""]["alpha_prob_MULTIxMK+G"]))
df = dfs[""]
print(max(df[df["alpha_BIN+G"] <= 90]["alpha_BIN+G"]))
