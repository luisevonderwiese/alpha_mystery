import os
from tabulate import tabulate
import json
import pandas as pd
import matplotlib.pyplot as plt

import util
import rates
import raxmlng

results_dir = os.path.join("results", "raxml_rates")
plots_dir = os.path.join("results", "plots_rates")
if not os.path.isdir(plots_dir):
    os.makedirs(plots_dir)
large_datasets =  ["abvdoceanic", "bowernpny", "iecor"]
datasets = os.listdir("data/msa")
datasets = [d for d in datasets if d not in large_datasets]


for i, dataset in enumerate(datasets):
    msa_path = os.path.join("data/msa/", dataset, "bin.phy")
    prefix = os.path.join(results_dir, dataset, "BIN+G")
    #raxmlng.run_inference_adaptive(msa_path, "BIN+G", prefix)
    srlhs = raxmlng.siterate_lhs(prefix)
    lines = [[srlhs[s][r] for s in range(len(srlhs))] for r in range(4)]

    for line in lines:
        plt.plot(line)
    plt.xlabel("sites")
    plt.ylabel("llh")
    name = dataset + ".png"
    plt.savefig(os.path.join(plots_dir, name))
    plt.clf()
    plt.close()
