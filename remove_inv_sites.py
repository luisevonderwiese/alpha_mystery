from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.AlignIO.PhylipIO import RelaxedPhylipWriter

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
        msa_string = "\n\n".join(parts[:-1] + ["\n".join([line + append_string for line in lines[:-1]] + [lines[-1]])])
    else:
        msa_string = "\n".join([lines[0]] + [line + append_string for line in lines[1:-1]] + [lines[-1]])

    parts = msa_string.split("\n")
    sub_parts = parts[0].split(" ")

    msa_string = "\n".join([" ".join(sub_parts[:-1] + [str(int(sub_parts[-1]) + padding_size)])] + parts[1:])

    with open(outpath, "w+") as new_msa_file:
        new_msa_file.write(msa_string)

msa_super_dir = "data/msa"
for dataset in os.listdir(msa_super_dir):
    msa_dir = os.path.join(msa_super_dir, dataset)
    bin_msa_path = os.path.join(msa_dir, "bin.phy")
    out_path = os.path.join(msa_dir, "noinvsite_bin.phy")
    try:
        align = AlignIO.read(bin_msa_path, "phylip-relaxed")
    except:
        write_padded_msa(os.path.join("data/msa/", dataset, "bin.phy"), "temp.phy")
        align = AlignIO.read("temp.phy", "phylip-relaxed")
        os.remove("temp.phy")
    sequences = ["" for _ in range(len(align))]
    for site_index in range(align.get_alignment_length()):
        site = align[:, site_index]
        i = 0
        v0 = None
        while(v0 is None and i < len(site)):
            if site[i] != "-":
                v0 = site[i]
            i+=1
        if v0 is None:
            continue
        invariant = True
        for j in range(i, len(site)):
            if site[j] == "-":
                continue
            if site[j] != v0:
                invariant = False
                break
        if invariant:
            continue
        for j in range(len(align)):
            sequences[j] += site[j]
    new_records = [SeqRecord(sequences[j], id=align[j].id) for j in range(len(align))]
    new_align = MultipleSeqAlignment(new_records, annotations={}, column_annotations={})
    with open(out_path, "w+") as f:
        writer = RelaxedPhylipWriter(f)
        writer.write_alignment(new_align)
