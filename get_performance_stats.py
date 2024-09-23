#!/bin/env python

import os
import filecmp
import traceback
import polars as pl
from math import sqrt
from datetime import datetime
from collections import defaultdict
import concurrent.futures as cf

# 
project_dir = "./projects/T2T/"
parent_qc_dir = os.path.join(project_dir, "analysis/QC/")
parent_dir = os.path.join(project_dir, "analysis/per_read/")
overwrite = False
parallel = True

# taxid below chordata(include)
host_taxid_hierarchy = [7711, 89593, 7742, 7776, 117570, 117571, 8287,
                        1338369, 32523, 32524, 40674, 32525, 9347,
                        1437010, 314146, 9443, 376913, 314293, 9526,
                        314295, 9604, 207598, 9605, 9606, 63221, 741158]

sequencers = {"ILMN"}  # ("ONT", "ILMN", "MGI")
labs_preset = {"Lab49"}

"""
qc:
    clean_reads = host_consensus + multi_labels + non_host_consensus

multi_labels:
    multi_labels = loosened_host_consensus + to_resolve
    for method_x in methods:
        multi_labels = labelled_host_by_method_x_in_multi_labels + labelled_non_host_by_method_x_in_multi_labels

to_resolve:
    to_resolve = multi_labels - loosened_host_consensus
    for method_x in methods:
        loosened_host_consensus = concordant_loosened_host_consensus + discordant_loosened_host_consensus

m8:
    to_resolve = resolved_host + resolved_non_host

resolved:
    for method_x in methods:
        resolved_host = concordant_resolved_host + discordant_resolved_host
        resolved_non_host = concordant_resolved_non_host + discordant_resolved_non_host

summary_stats:
    host_ground_truth = TP + FN = host_consensus + loosened_host_consensus + resolved_host

for method_x in methods:
    labelled_host_by_method_x = TP + FP = host_consensus + concordant_loosened_host_consensus(x) + 
        concordant_resolved_host(x) + discordant_resolved_non_host(x)
    labelled_non_host_by_method_x = TN + FN = non_host_consensus + concordant_resolved_non_host(x) + 
        discordant_loosened_host_consensus(x) + discordant_resolved_host(x)
    TP = host_consensus + concordant_loosened_host_consensus(x) + concordant_resolved_host(x)
    FN = discordant_loosened_host_consensus(x) + discordant_resolved_host(x)
    TN = non_host_consensus + concordant_resolved_non_host(x)
    FP = discordant_resolved_non_host(x)

parse qc -> clean_reads
parse non_host -> non_host_consensus
parse multi_labels -> multi_labels
    -> host_consensus = clean_reads - non_host_consensus - multi_labels
    -> loosened_host_consensus =  multi_labels - to_resolve
        -> concordant_loosened_host_consensus + discordant_loosened_host_consensus
parse blast.m8 -> resolved_host + resolved_non_host -> join to_resolve -> 
    concordant_resolved_host + discordant_resolved_host + concordant_resolved_non_host + discordant_resolved_non_host
"""


def main():
    if parallel:
        with cf.ProcessPoolExecutor(max_workers=3) as executor:
            try:
                for _ in executor.map(get_performance_stats, sequencers):
                    pass
            except Exception:
                traceback.print_exc()
    else:
        for sequencer in sequencers:
            get_performance_stats(sequencer)


def get_performance_stats(sequencer):
    if sequencer in ("ONT",):
        qc_file = parent_qc_dir + "ONT_QC_summary.tsv"
    elif sequencer in ("ILMN", "MGI"):
        qc_file = parent_qc_dir + "SR_QC_summary.tsv"
    else:
        raise Exception("Unknown sequencer!")
    qc = pl.scan_csv(qc_file, sep="\t").filter(
        pl.col("Sequencer") == sequencer
    ).select([
        pl.col("Lab"),
        pl.col("Prefix"),
        pl.col("Clean_Reads")
    ])
    labs = set(qc.select("Lab").collect().to_series())
    if labs_preset:
        labs = [l for l in labs if any(sub in l for sub in labs_preset)]
    result = []
    methods = []
    sample_name_dict = get_sample_name_dict()
    for lab in labs:
        prefixes = qc.filter(pl.col("Lab") == lab).select("Prefix").collect().to_series().to_list()
        for prefix in prefixes:
            now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print("[" + now + "] Processing " + sequencer + "/" + lab + "/" + prefix, flush=True)
            clean = qc.filter((pl.col("Lab") == lab) & (pl.col("Prefix") == prefix)).select("Clean_Reads"). \
                collect().row(0)[0]
            non_host_consensus = count_lines(os.path.join(parent_dir,
                                                          sequencer, lab, prefix + "_non_host_reads.tsv")) - 1
            multi_labels_df = pl.scan_csv(os.path.join(parent_dir, sequencer, lab, prefix + "_multi_labels.tsv"),
                                          sep="\t")
            # Header of multi_labels_df is ["ReadID"] + all methods + ["h", "nh", "kraken_h", "kraken_nh"]
            multi_labels_df = multi_labels_df.select([m for m in multi_labels_df.columns
                                                      if not m in ["h", "nh", "kraken_h", "kraken_nh"]])
            multi_labels = multi_labels_df.collect().height
            host_consensus = clean - non_host_consensus - multi_labels
            # If resolved_file already exists, read it. Otherwise, read to_resolve_file and parse m8 file.
            to_resolve_file = os.path.join(parent_dir, sequencer, lab, prefix +
                                           "_to_resolve.tsv")
            resolved_file = os.path.join(parent_dir, sequencer, lab, prefix + "_resolved.tsv")
            if not os.path.exists(resolved_file):
                to_resolve_df = pl.scan_csv(to_resolve_file, sep="\t").select("ReadID")
                try:
                    blast_out = pl.scan_csv(os.path.join(parent_dir, sequencer, lab, prefix + ".m8"),
                                            sep="\t", has_header=False, null_values="N/A",
                                            dtypes={"column_1": pl.Utf8, "column_2": pl.Utf8, "column_3": pl.Float32,
                                                    "column_4": pl.Int32, "column_5": pl.Int32, "column_6": pl.Int32,
                                                    "column_7": pl.Int32, "column_8": pl.Int32, "column_9": pl.Int32,
                                                    "column_10": pl.Int32, "column_11": pl.Float32,
                                                    "column_12": pl.Float32, "column_13": pl.Int32})
                    """columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", 
                    "qend", "sstart", "send", "evalue", "bitscore", "staxid"] """
                    blast_out = blast_out.groupby("column_1").agg(pl.col("column_13"). \
                        is_in(host_taxid_hierarchy).any().cast(pl.Int8).alias(
                        "Resolved"))
                    resolved_df = to_resolve_df.join(blast_out, left_on="ReadID",
                                                     right_on="column_1", how="left"). \
                        collect().fill_null(0)
                except pl.exceptions.NoDataError:
                    resolved_df = to_resolve_df.with_column(pl.lit(0).alias("Resolved")). \
                        collect()
                resolved_df.write_csv(resolved_file, sep="\t")
                resolved_df = resolved_df.lazy()
                os.remove(to_resolve_file)
                if os.path.exists(os.path.join(parent_dir, sequencer, lab, prefix + ".m8")):
                    os.remove(os.path.join(parent_dir, sequencer, lab, prefix + ".m8"))
            else:
                resolved_df = pl.scan_csv(resolved_file, sep="\t")
            """ The resolved_df now has two columns ["ReadID", "Resolved"]. 
            The "Resolved" column being 1 means it is a host read, and 0 means non-host, as determined by BLAST."""
            loosened_host_consensus_df = multi_labels_df.filter(pl.col("ReadID").is_in(resolved_df.
                                                                                       select(
                "ReadID").collect().to_series().to_list()).is_not())
            loosened_host_consensus = loosened_host_consensus_df.collect().height
            resolved_host = resolved_df.filter(pl.col("Resolved") == 1).collect().height
            resolved_df = resolved_df.join(multi_labels_df, on="ReadID", how="left")
            host_reads_ground_truth = host_consensus + loosened_host_consensus + resolved_host
            host_reads_perc = 100 * host_reads_ground_truth / clean
            res_per_sample = [sequencer, lab, prefix, sample_name_dict[prefix], clean, host_consensus,
                              non_host_consensus, multi_labels, loosened_host_consensus, resolved_host,
                              host_reads_ground_truth, host_reads_perc]
            methods = [m for m in multi_labels_df.columns if not m == "ReadID"]
            for m in methods:
                if loosened_host_consensus == 0:
                    concordant_loosened_host_consensus = discordant_loosened_host_consensus = 0
                else:
                    if "kraken" in m:
                        loosened_host_consensus_df = loosened_host_consensus_df.with_column(pl.col(m). \
                                                                                            is_in(
                            host_taxid_hierarchy).alias(m))
                    concordant_loosened_host_consensus = loosened_host_consensus_df.filter(pl.col(m) == 1). \
                        collect().height
                    discordant_loosened_host_consensus = loosened_host_consensus_df.filter(pl.col(m) == 0). \
                        collect().height
                if resolved_df.collect().height == 0:
                    concordant_resolved_host = discordant_resolved_host = 0
                    concordant_resolved_non_host = discordant_resolved_non_host = 0
                else:
                    if "kraken" in m:
                        resolved_df = resolved_df.with_column(pl.col(m).is_in(host_taxid_hierarchy).alias(m))
                    concordant_resolved_host = resolved_df.filter((pl.col("Resolved") == 1) & (pl.col(m) == 1)). \
                        collect().height
                    discordant_resolved_host = resolved_df.filter((pl.col("Resolved") == 1) & (pl.col(m) == 0)). \
                        collect().height
                    concordant_resolved_non_host = resolved_df.filter((pl.col("Resolved") == 0)
                                                                      & (pl.col(m) == 0)).collect().height
                    discordant_resolved_non_host = resolved_df.filter((pl.col("Resolved") == 0)
                                                                      & (pl.col(m) == 1)).collect().height
                TP = host_consensus + concordant_loosened_host_consensus + concordant_resolved_host
                FN = discordant_loosened_host_consensus + discordant_resolved_host
                TN = non_host_consensus + concordant_resolved_non_host
                FP = discordant_resolved_non_host
                sens = TP / (TP + FN)
                spec = TN / (TN + FP)
                PPV = TP / (TP + FP)
                NPV = TN / (TN + FN)
                Yoden = sens + spec - 1
                MCC = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
                F1 = 2 * TP / (2 * TP + FN + FP)
                G1 = 2 * TN / (2 * TN + FN + FP)
                res_per_sample.extend([TP, FN, TN, FP, sens, spec, PPV, NPV, F1, G1, Yoden, MCC])
            result.append(res_per_sample)
    out_dir = os.path.join(parent_dir, sequencer)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    if labs_preset:
        prefix = [p.split("_")[0] for p in labs_preset]
        out_file = os.path.join(out_dir, sequencer + "_" + "_".join(prefix) + \
                                "_performance_stats.tsv")
    else:
        out_file = os.path.join(out_dir, sequencer + "_performance_stats.tsv")
    new_out_file = out_file + ".new"
    with open(new_out_file, "w") as out_fp:
        header = "Sequencer\tLab\tPrefix\tSample\tClean_reads\tHost_consensus\tNon_host_consensus\t" \
                 "Multi_labels\tLoosened_host_consensus\tBLAST_labelled_host\t" \
                 "Host_reads_ground_truth\tHost_reads_perc\t" + \
                 "\t".join([m + "__" + n for m in methods for n in ["TP", "FN", "TN", "FP", "sens",
                                                                    "spec", "PPV", "NPV", "F1", "G1", "Yoden", "MCC"]])
        out_fp.write(header + "\n")
        for line in result:
            out_fp.write("\t".join(map(str, line)) + "\n")
    if overwrite or not os.path.exists(out_file):
        os.rename(new_out_file, out_file)
    else:
        if filecmp.cmp(out_file, new_out_file):
            print("Identical results. No changes will be made!")
            os.remove(new_out_file)
        else:
            print("Discrepancy found in results. New results in " + new_out_file
                  + ". Exiting without overwriting the old file!")
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("[" + now + "] " + sequencer + " analysis completed", flush=True)


def get_sample_name_dict():
    sample_name_dict = defaultdict(str)
    with open(os.path.join(parent_qc_dir, "samples.tsv"), "r") as in_fp:
        in_fp.readline()
        for line in in_fp:
            data = line.rstrip().split('\t')
            sample_name_dict[data[0]] = data[1]
    return sample_name_dict


def count_lines(in_file):
    n_lines = 0
    size = 1024 * 1024
    with open(in_file, "r") as in_fp:
        buffer = in_fp.read(size)
        while buffer:
            n_lines += buffer.count("\n")
            buffer = in_fp.read(size)
    return n_lines


if __name__ == "__main__":
    main()
