import os


def best_tree_path(prefix):
    return prefix + ".treefile"


def alpha(prefix):
    if not os.path.isfile(prefix + ".log"):
        return float("nan")
    with open(prefix + ".log", "r") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("Gamma shape alpha:"):
            return float(line.split(": ")[1])
    return float('nan')


def inv_estimate(prefix):
    if not os.path.isfile(prefix + ".log"):
        return float("nan")
    with open(prefix + ".log", "r") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("Proportion of invariable sites:"):
            return float(line.split(": ")[1])
    return float('nan')


def final_llh(prefix):
    with open(prefix + ".log", "r") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("BEST SCORE FOUND"):
            return float(line.split(": ")[1])
    return float('nan')


def free_rates(prefix):
    with open(prefix + ".log", "r") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("Site proportion and rates"):
            return line.split(":  ")[1]
    return ""


def run_inference(msa_path, model, prefix, args = ""):
    if not os.path.isfile(msa_path):
        print("MSA " + msa_path + " does not exist")
        return
    prefix_dir = "/".join(prefix.split("/")[:-1])
    if not os.path.isdir(prefix_dir):
        os.makedirs(prefix_dir)
    command = "./bin/iqtree2"
    command += " -s " + msa_path
    command += " -m " + model
    command += " --prefix " + prefix
    if not os.path.isfile(best_tree_path(prefix)):
        command += " --redo"
    os.system(command)
