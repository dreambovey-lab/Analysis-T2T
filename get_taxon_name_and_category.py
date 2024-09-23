#!/usr/bin/env python

import subprocess
import fire

clinical_categories = {"2": "Bacteria", "4751": "Fungi", "10239": "Viruses"}


def get_lineage(taxids):
    cmd = ("./.conda/envs/bin/taxonkit lineage"
           " --data-dir ./databases/taxdump -j 4 -t")
    cpt = subprocess.run(cmd,
                         input=taxids,
                         shell=True,
                         executable='/bin/bash',
                         capture_output=True,
                         encoding="utf-8")
    taxname_list = []
    taxid_lineage_string_list = []
    for line in cpt.stdout.strip().split("\n"):
        fields = line.split("\t")
        taxname = fields[1].split(";")[-1]
        if taxname:
            taxname_list.append(taxname)
        else:
            taxname_list.append("Unclassified")
        taxid_lineage_string_list.append(fields[2])
    return taxname_list, taxid_lineage_string_list


def get_clinical_category(taxid_lineage_string_list):
    category_list = []
    for taxid_lineage_string in taxid_lineage_string_list:
        category = "Others"
        for taxid in taxid_lineage_string.split(";"):
            if taxid in clinical_categories:
                category = clinical_categories[taxid]
                break
        category_list.append(category)
    return category_list


def main(taxid_list=None):
    if taxid_list:
        taxid_list = list(map(str, taxid_list))
        taxids = "\n".join(taxid_list)
    else:
        with open("TaxIDs.txt", "r") as in_fp:
            taxids = in_fp.read().replace(",", "\n").strip()
        taxid_list = taxids.split("\n")
    taxname_list, taxid_lineage_string_list = get_lineage(taxids)
    category_list = get_clinical_category(taxid_lineage_string_list)
    with open("TaxID2name2ctg.tsv", "w") as out_fp:
        out_fp.write("TaxID\tName\tCategory\n")
        output_list = zip(taxid_list, taxname_list, category_list)
        for output_tuples in output_list:
            out_fp.write("\t".join(output_tuples))
            out_fp.write("\n")


if __name__ == "__main__":
    fire.Fire(main)
