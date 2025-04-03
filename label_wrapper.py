import os
import util

def label_command(msa_path, prefix):
    if msa_path != msa_path:
        print("MSA does not exist")
        return
    command = "label"
    command += " -m " + msa_path
    command += " -i  bin/iqtree2"
    command += " -p " + prefix
    os.system(command)

def get_label(prefix):
    with open(prefix + ".labelGen.log", "r", encoding="utf-8") as out_file:
        line = out_file.readlines()[-3]
    if not line.startswith("Ground Truth Difficulty"):
        raise ValueError("Error during label calculation")
    return float(line.split(" ")[-1])

def calculate_label(msa_path, prefix):
    try:
        get_label(prefix)
    except (FileNotFoundError, ValueError):
        print("calculating difficulty label")
        label_command(msa_path, prefix)
        try:
            get_label(prefix)
        except ValueError:
            print("trying with padded msa")
            util.write_padded_msa(msa_path, "temp.phy")
            label_command("temp.phy", prefix)
            os.remove("temp.phy")
