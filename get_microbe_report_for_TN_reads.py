#!/bin/env python

import re
import fire
import subprocess
import polars as pl
from pathlib import Path
from datetime import datetime
from collections import Counter

project_dir = Path("./projects/T2T/")
data_dir = Path(project_dir, "analysis/per_read/")
kraken_res_dir = Path(project_dir, "analysis/microbe_only_classification/")
sequencers = {"MGI", "ILMN", "ONT"}  # {"ONT", "ILMN", "MGI"}
ref_genomes = {"hg38", "T2T", "YH"}
labs_preset = {}
sample_preset = {}

host_taxid_hierarchy = [7711, 89593, 7742, 7776, 117570, 117571, 8287,
                        1338369, 32523, 32524, 40674, 32525, 9347,
                        1437010, 314146, 9443, 376913, 314293, 9526,
                        314295, 9604, 207598, 9605, 9606, 63221, 741158]


def main(ranks="tsg", dry_run=False):
    res = []
    for sequencer in sequencers:
        labs = [x.name for x in Path(data_dir, sequencer).iterdir() if
                x.is_dir() and x.name.startswith(("Lab", "clinical", "POOLING"))]
        if labs_preset:
            labs = [lab for lab in labs if any(sub in lab for sub in labs_preset)]
        if "clinical" in labs:
            labs.remove("clinical")
            labs = ["clinical"] + labs
        for lab in labs:
            now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print("[" + now + "] Processing " + sequencer + "/" + lab, flush=True)
            data_files = list(Path(data_dir, sequencer, lab).glob("*multi_labels*"))
            if sample_preset:
                data_files = [x for x in data_files if any(sub in x.name for sub in sample_preset)]
            if dry_run:
                print([f.name for f in data_files], flush=True)
            else:
                for f in data_files:
                    res.append(get_microbe_taxid(f))
    for r in ranks:
        taxonomy_ranks = {"k": "Kingdom",
                          "p": "Phylum",
                          "c": "Class",
                          "o": "Order",
                          "f": "Family",
                          "g": "Genus",
                          "s": "Species",
                          "S": "Subspecies",
                          "t": "Taxon"}
        rank = taxonomy_ranks[r]
        res_df = pl.DataFrame(columns=[(rank, pl.Int64)])
        for r in res:
            for identifier, taxid_list in r.items():
                taxon_stats = Counter(unify_taxonomy_rank(taxid_list, rank=rank)).most_common(100)
                df = pl.DataFrame(taxon_stats,
                                  orient="row",
                                  columns=[rank, identifier])
                res_df = res_df.join(df, how="outer", on=rank)
        res_df.write_csv(Path(kraken_res_dir, rank + "_report_for_TN_reads.tsv"), sep="\t")
        res_df.melt(
            id_vars=rank, variable_name="Sample"
        ).pivot(
            values="value", index="Sample", columns=rank
        ).write_csv(Path(kraken_res_dir, rank + "_report_for_TN_reads_transposed.tsv"), sep="\t")
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("[" + now + "] " + "Analysis completed", flush=True)


def get_microbe_taxid(in_file):
    sequencer = in_file.parts[-3]
    lab = in_file.parts[-2]
    prefix = in_file.name.split("_multi")[0]
    nh_kraken_result_file = Path(kraken_res_dir, sequencer, lab, prefix + "_nh_kraken.out")
    list1 = []
    with open(nh_kraken_result_file) as f:
        for line in f:
            list1.append(int(line.split("\t")[2]))
    df = pl.scan_csv(in_file, sep="\t")
    # ReadID genome_algorithm h nh kraken_h kraken_nh
    resolved_df = pl.scan_csv(Path(in_file.parent, prefix + "_resolved.tsv"),
                              sep="\t",
                              dtypes={"ReadID": pl.Utf8, "Resolved": pl.Int32})
    df = df.select(
        pl.all().exclude("^(non_)?(kraken_)?n?h$")
    ).join(
        resolved_df, on="ReadID", how="left"
    ).collect().fill_null(1).lazy()
    ml_kraken_result_file = Path(kraken_res_dir, sequencer, lab, prefix + "_ml_kraken.out")
    df2 = pl.read_csv(ml_kraken_result_file,
                      sep="\t",
                      has_header=False,
                      columns=["column_2", "column_3"],
                      # dtypes={"column_2": pl.Utf8, "column_3": pl.Int32},
                      new_columns=["ReadID", "TaxID"]).lazy()
    df = df.join(df2, on="ReadID")
    res = {}
    for genome in ref_genomes:
        list2 = df.select([
            pl.col("^" + genome + ".*$"),
            pl.col("Resolved"),
            pl.col("TaxID")
        ]).with_column(
            pl.concat_list(
                pl.all().exclude(genome + "__kraken2").exclude("Resolved").exclude("TaxID")
            ).alias(genome + "__non_kraken")
        ).filter(
            (pl.col("Resolved") == 0) &
            pl.col(genome + "__kraken2").is_in(host_taxid_hierarchy).is_not() &
            pl.col(genome + "__non_kraken").arr.eval(
                (pl.element() == 0).all(),
                parallel=True
            ).flatten()
        ).select("TaxID").collect().to_series().to_list()
        res["#".join([sequencer, lab, prefix, genome])] = list1 + list2
    return res


def unify_taxonomy_rank(list_of_taxa, rank):
    if rank == "Taxon":
        list_of_species = list_of_taxa
    else:
        taxonomy_rank_labels = {"Kingdom": "k",
                                "Phylum": "p",
                                "Class": "c",
                                "Order": "o",
                                "Family": "f",
                                "Genus": "g",
                                "Species": "s",
                                "Subspecies": "S"}
        my_string = "\n".join(map(str, list_of_taxa))
        cmd = ("./.conda/envs/bin/taxonkit reformat"
               " --data-dir ./databases/taxdump"
               " -j 4 -I 1 -t --miss-rank-repl '__' -f '{" +
               taxonomy_rank_labels[rank] + "}' | cut -f3")
        cpt = subprocess.run(cmd,
                             input=my_string,
                             shell=True,
                             executable='/bin/bash',
                             capture_output=True,
                             encoding="utf-8")
        list_of_species = map(int, map(lambda x: re.sub("^$", "0", x), cpt.stdout.split("\n")[:-1]))
    return list_of_species


if __name__ == "__main__":
    fire.Fire(main)
