import os
from ete3 import Tree

exe_path = ""


def best_tree_path(prefix):
    return prefix + ".raxml.bestTree"

def ml_trees_path(prefix):
    return prefix + ".raxml.mlTrees"

def siterate_lhs_path(prefix):
    return prefix + ".raxml.siterateLH"

def siterate_lhs(prefix):
    path = siterate_lhs_path(prefix)
    if not os.path.isfile(path):
        return []
    with open(path, "r") as f:
        lines = f.readlines()
    siteratelhs = []
    for line in lines:
        siteratelhs.append([float(el) for el in line.split(" ")[:-1]])
    return siteratelhs



def alpha(prefix):
    if not os.path.isfile(prefix + ".raxml.log"):
        return float("nan")
    with open(prefix + ".raxml.log", "r") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("   Rate heterogeneity:"):
            return float(line.split(",  ")[1].split(" ")[1])
    return float('nan')


def inv_estimate(prefix):
    if not os.path.isfile(prefix + ".raxml.log"):
        return float("nan")
    with open(prefix + ".raxml.log", "r") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("   P-inv"):
            return float(line.split(" ")[5][:-1])
    return float('nan')


def final_llh(prefix):
    with open(prefix + ".raxml.log", "r") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("Final LogLikelihood: "):
            return float(line.split(": ")[1])
    return float('nan')




def free_rates(prefix):
    with open(prefix + ".raxml.log", "r") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("   Rate heterogeneity: "):
            return line.split("weights&rates: ")[1]
    return ""

def aic(prefix):
    logpath = prefix + ".raxml.log"
    if not os.path.isfile(logpath):
        return [float('nan'), float('nan'), float('nan')]
    with open(logpath, "r") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("AIC"):
            parts = line.split(" / ")
            scores = []
            for part in parts:
                scores.append(float(part.split(" ")[2]))
            return scores
    return [float('nan'), float('nan'), float('nan')]




def run_inference(msa_path, model, prefix, args = ""):
    if exe_path == "":
        print("Please specify raxmlng.exe_path")
        return
    if not os.path.isfile(msa_path):
        print("MSA " + msa_path + " does not exist")
        return
    prefix_dir = "/".join(prefix.split("/")[:-1])
    if not os.path.isdir(prefix_dir):
        os.makedirs(prefix_dir)
    if not os.path.isfile(best_tree_path(prefix)):
        args = args + " --redo"
    command = exe_path
    command += " --msa " + msa_path
    command += " --model " + model
    command += " --prefix " + prefix
    command += " --threads auto --seed 2 --force model_lh_impr"
    command += " " + args
    os.system(command)
