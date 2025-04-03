import os
from tabulate import tabulate
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score

import util
import rates
import raxmlng
import label_wrapper


def scatterplot(df, column_a, column_b):
    plt.scatter(df[column_a], df[column_b], s=10)
    plt.xlabel(column_a)
    plt.ylabel(column_b)
    name = column_a + "_vs_" + column_b + ".png"
    plt.savefig(os.path.join("data", "plots", name))
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

def get_auc_score(df, column_a, column_b):
    data_a = df[column_a]
    data_b = df[column_b]
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


def statistical_analysis(df, alpha_column, other_columns):
    print(alpha_column)
    res = []
    for column in other_columns:
        r = [column]
        a = get_auc_score(df, alpha_column, column)
        r.append(round(a[0], 3))
        r.append(round(a[1], 3))
        r.append(round(a[2], 3))
        res.append(r)
    headers = ["column", "AUC", "low alpha mean", "high alpha mean"]
    print(tabulate(res, tablefmt="latex", headers=headers))


results_dir = os.path.join("data", "raxml")
label_dir = os.path.join("data", "difficulty_labels")
if not os.path.isdir(label_dir):
    os.makedirs(label_dir)
metadata_df = pd.read_csv("data/lexibench/character_matrices/stats.tsv", sep = "\t")
wl_df = pd.read_csv("data/lexibench/lingpy_wordlists/stats.tsv", sep = "\t")
metadata_df = metadata_df.merge(wl_df, on = "Name") 

bin_models = ["BIN", "BIN+G", "BIN+R4", "BIN+FO+I"]
gamma_models = ["BIN+G", "prob_MULTIxMK+G", "prob_BIN+G"]

for i, row in metadata_df.iterrows():
    dataset = row["Name"]
    label_wrapper.calculate_label(row["bin.phy"], os.path.join(label_dir, dataset, "label"))
    for model in bin_models:
        raxmlng.run_inference(row["bin.phy"], model, os.path.join(results_dir, dataset, model))
    raxmlng.run_inference(row["bin.catg"], "BIN+G", os.path.join(results_dir, dataset, "prob_BIN+G"), "--prob-msa on")
    x = row["cs_max"]
    if x <= 64:
        raxmlng.run_inference(os.path.join("data/lexibench/character_matrices/", dataset, "multi.catg"), 
                "MULTI" + str(x) + "_MK+G", os.path.join(results_dir, dataset, "prob_MULTIxMK+G"), "--prob-msa on")
        

columns = [
        "Name", 
        "num_sites", 
        "bin_entropy",
        "entropy_var", 
        "inv_sites_emp", 
        "zero_freq_emp", 
        "inv_sites_estimate", 
        "free_rates_var", 
        "inv_sites_error",
        "difficult"
        ]
columns += ["alpha_" + model for model in gamma_models]
df = pd.DataFrame(columns = columns)
for i, row in metadata_df.iterrows():
    dataset = row["Name"]
    print(dataset)
    df.at[i, "Name"] = dataset
    try:
        align = util.safe_msa_read(row["bin.phy"])
    except:
        continue
    df.at[i, "num_sites"] = util.num_sites(align)
    df.at[i, "bin_entropy"] = util.bin_entropy(align)
    df.at[i, "entropy_var"] = util.entropy_var(align)
    df.at[i, "inv_sites_emp"] = util.inv_sites(align)
    df.at[i, "zero_freq_emp"] = util.zero_freq(align)

    df.at[i, "inv_sites_estimate"] = raxmlng.inv_estimate(os.path.join(results_dir, dataset, "BIN+FO+I")) * 100
    var = rates.var(rates.parse_rates(raxmlng.free_rates(os.path.join(results_dir, dataset, "BIN+R4"))))
    if var > 10:
        var = float("nan")
    df.at[i, "free_rates_var"] = var
    df.at[i, "inv_sites_error"] = abs(df.at[i, "inv_sites_estimate"] - df.at[i, "inv_sites_emp"])
    try:
        df.at[i, "difficult"] = label_wrapper.get_label(os.path.join(label_dir, dataset, "label"))
    except ValueError:
        continue
    for model in gamma_models:
        raxml_prefix = os.path.join(results_dir, dataset, model)
        df.at[i, "alpha_" + model] = raxmlng.alpha(raxml_prefix)
    
df = df.merge(metadata_df, on = "Name")

if not os.path.isdir(os.path.join("data", "plots")):
    os.makedirs(os.path.join("data", "plots"))

scatterplot(df, "alpha_BIN+G", "alpha_prob_MULTIxMK+G")


statistical_analysis(df, "alpha_BIN+G",
            [
            #"entropy_var",
            "Languages",
            "num_sites",
            "Concepts",
            "cs_mean",
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
for i, row in df.iterrows():
    r = [row["Name"], row["Languages"], row["cs_mean"], row["alpha_BIN+G"], row["difficult"]]
    if row["alpha_BIN+G"] > 90:
        res_high.append(r)
    else:
        res_low.append(r)

headers = ["Name", "Languages", "cs_mean", "alpha", "difficult"]
print("low alpha")
print(tabulate(res_low, tablefmt="pipe", headers=headers))
print("high alpha")
print(tabulate(res_high, tablefmt="pipe", headers=headers))

print(max(df[df["alpha_prob_MULTIxMK+G"].notnull()]["alpha_prob_MULTIxMK+G"]))
print(min(df[df["alpha_BIN+G"] > 90]["alpha_BIN+G"]))
print(max(df[df["alpha_BIN+G"] <= 90]["alpha_BIN+G"]))
