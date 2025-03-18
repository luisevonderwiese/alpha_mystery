import os

def write_padded_msa(msa_path, outpath):
    with open(msa_path, "r") as msa_file:
        msa_string = msa_file.read()
    parts = msa_string.split("\n\n")
    lines = parts[-1].split("\n")
    block_size = len(lines[1].split(" ")[-1])
    if block_size == 10:
        padding_size = 10
        append_string = " ----------"
    else:
        padding_size = 10 - block_size
        append_string = "-" * padding_size
    if len(parts) != 1:
        msa_string = "\n\n".join(parts[:-1] + ["\n".join([line + append_string for line in lines])])
    else:
        msa_string = "\n".join([lines[0]] + [line + append_string for line in lines[1:]])

    parts = msa_string.split("\n")
    sub_parts = parts[0].split(" ")

    msa_string = "\n".join([" ".join(sub_parts[:-1] + [str(int(sub_parts[-1]) + padding_size)])] + parts[1:])

    with open(outpath, "w+") as new_msa_file:
        new_msa_file.write(msa_string)

def label_command(msa_path, prefix):
    command = "label"
    command += " -m " + msa_path
    command += " -i  bin/iqtree2"
    command += " -p " + prefix
    os.system(command)

def get_label(prefix):
    with open(prefix + ".labelGen.log", "r") as out_file:
        line = out_file.readlines()[-3]
    if not line.startswith("Ground Truth Difficulty"):
        raise ValueError("Error during label calculation")
    return float(line.split(" ")[-1])

def calculate_label(msa_path, prefix):
    try:
        get_label(prefix)
    except (FileNotFoundError, ValueError) as e1:
        print("calculating difficulty label")
        label_command(msa_path, prefix)
        try:
            get_label(prefix)
        except ValueError as e:
            print("trying with padded msa")
            write_padded_msa(msa_path, "temp.phy")
            label_command("temp.phy", prefix)
            os.remove("temp.phy")
