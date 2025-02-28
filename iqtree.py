import os

from iqtree_statstest_parser import get_iqtree_results

iqtree_path = "../iqtree-2.0.6-Linux/bin/iqtree2"



def run_statstests(msa_path, ref_tree_path, bing_tree_path, prefix):
    if not os.path.isfile(msa_path):
        print("MSA", msa_path, "does not exist")
        return
    if not os.path.isfile(ref_tree_path):
        print("Tree", ref_tree_path, "does not exist")
        return
    if not os.path.isfile(bing_tree_path):
        print("Tree", bing_tree_path, "does not exist")
        return
    if not os.path.isdir(prefix):
        os.makedirs(prefix)
    with open("temp.trees", "w+", encoding = "utf-8") as temp_tree_file:
        temp_tree_file.write(open(ref_tree_path, "r", encoding = "utf-8").read())
        temp_tree_file.write(open(bing_tree_path, "r", encoding = "utf-8").read())
    command = iqtree_path
    command += " -s " + msa_path 
    command += " -st MORPH" 
    command += " -m MK "
    command += " -pre " + os.path.join(prefix, "significance")
    command += " -z temp.trees"
    command += " -te " + ref_tree_path
    command += " -n 0 -zb 10000 -zw -au -nt 1 -treediff -seed 0" 
    #command += " -redo"
    command += " > " + os.path.join(prefix, "significance.iqtree.log")
    os.system(command)


def is_plausible(prefix):
    iqtree_results = get_iqtree_results(os.path.join(prefix, "significance.iqtree"))
    assert (len(iqtree_results) == 2)
    r = iqtree_results[1]
    return (r["plausible"], r["all_identical"])
