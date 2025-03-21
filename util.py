import os
import math
import numpy as np
from Bio import AlignIO


def write_padded_msa(msa_path, outpath):
    with open(msa_path, "r", encoding="utf-8") as msa_file:
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

    with open(outpath, "w+", encoding="utf-8") as new_msa_file:
        new_msa_file.write(msa_string)


def save_msa_read(path):
    try:
        if path.endswith("abvdoceanic-austronesian/bin.phy"):
            raise ValueError("Somehow reading of the MSA never finishes")
        align = AlignIO.read(path, "phylip-relaxed")
        return align
    except:
        write_padded_msa(path, "temp.phy")
        align = AlignIO.read("temp.phy", "phylip-relaxed")
        os.remove("temp.phy")
        return align


def site_entropy(site):
    site = site.replace("-", "")
    entropy = 0
    for char in ["0", "1"]:
        count = site.count(char)
        if count == 0:
            entropy_x = 0
        else:
            prob = count / len(site)
            entropy_x = prob * math.log2(prob)
        entropy += entropy_x
    return -entropy

def is_invariant(site):
    i = 0
    v0 = None
    while (v0 is None and i < len(site)):
        if site[i] != "-":
            v0 = site[i]
        i+=1
    if v0 is None:
        return True
    for j in range(i, len(site)):
        if site[j] == "-":
            continue
        if site[j] != v0:
            return False
    return True

def entropy_var(align):
    return np.var([site_entropy(align[:, site_index]) for site_index in range(align.get_alignment_length())])

def bin_entropy(align):
    return sum([site_entropy(align[:, site_index]) for site_index in range(align.get_alignment_length())])

def num_sites(align):
    return align.get_alignment_length()

def num_taxa(align):
    return len(align)

def inv_sites(align):
    cnt = 0
    for site_index in range(align.get_alignment_length()):
        if is_invariant(align[:, site_index]):
            cnt += 1
    return (cnt / align.get_alignment_length()) * 100


def zero_freq(align):
    ones = 0
    zeros = 0
    for site_index in range(align.get_alignment_length()):
        site = align[:, site_index]
        ones += site.count("1")
        zeros += site.count("0")
    return zeros / (ones + zeros)
