#!/bin/env python

import os
import glob
import traceback
import filecmp
import subprocess
import polars as pl
from io import StringIO
from datetime import datetime
import concurrent.futures as cf

test_mode = False
overwrite = True
parallel = True
num_of_threads = 4
# export POLARS_MAX_THREADS=16

project_dir = "./T2T/"
parent_qc_dir = os.path.join(project_dir, "analysis/QC/")
parent_alignment_dir = os.path.join(project_dir, "analysis/host_reads_removal/")
parent_in_dir = os.path.join(project_dir, "analysis/per_read/")
if test_mode:
    parent_out_dir = os.path.join(project_dir, "test/unique_reads/")
else:
    parent_out_dir = os.path.join(project_dir, "analysis/unique_reads/")
sequencers = ["ONT", "ILMN", "MGI"]  # ("ONT", "ILMN", "MGI")
labs_preset = []
samples_preset = []

ref_genomes = ["T2T", "hg38", "YH"]
alignment_algorithms = ["bwa_mem", "bowtie2_vs", "bowtie2_vf", "minimap2", "winnowmap"]

# TaxIDs below chordata(include)
host_taxid_hierarchy = [7711, 89593, 7742, 7776, 117570, 117571, 8287,
                        1338369, 32523, 32524, 40674, 32525, 9347,
                        1437010, 314146, 9443, 376913, 314293, 9526,
                        314295, 9604, 207598, 9605, 9606]

# T2T genome accession number for each chromosome
T2T_acc2chr = {"CP068277.2": "Chr1",
               "CP068276.2": "Chr2",
               "CP068275.2": "Chr3",
               "CP068274.2": "Chr4",
               "CP068273.2": "Chr5",
               "CP068272.2": "Chr6",
               "CP068271.2": "Chr7",
               "CP068270.2": "Chr8",
               "CP068269.2": "Chr9",
               "CP068268.2": "Chr10",
               "CP068267.2": "Chr11",
               "CP068266.2": "Chr12",
               "CP068265.2": "Chr13",
               "CP068264.2": "Chr14",
               "CP068263.2": "Chr15",
               "CP068262.2": "Chr16",
               "CP068261.2": "Chr17",
               "CP068260.2": "Chr18",
               "CP068259.2": "Chr19",
               "CP068258.2": "Chr20",
               "CP068257.2": "Chr21",
               "CP068256.2": "Chr22",
               "CP068255.2": "ChrX",
               "CP086569.2": "ChrY",
               "CP068254.1": "MT"}


def main():
    for sequencer in sequencers:
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print("[" + now + "] Processing " + sequencer, flush=True)
        file_info_list = get_file_info(sequencer)
        if parallel:
            with cf.ProcessPoolExecutor(max_workers=num_of_threads) as executor:
                try:
                    for _ in executor.map(get_genomic_info, file_info_list):
                        pass
                except Exception:
                    traceback.print_exc()
        else:
            for file_info in file_info_list:
                get_genomic_info(file_info)
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("[" + now + "] " + "Analysis completed", flush=True)


def get_file_info(sequencer: str) -> list:
    my_list = []
    labs = get_subdirectories(os.path.join(parent_in_dir, sequencer))
    if labs_preset:
        labs = [l for l in labs if any(sub in l for sub in labs_preset)]
    for lab in labs:
        in_files = glob.glob(os.path.join(parent_in_dir, sequencer, lab, "*_multi_labels.tsv"))
        prefixes = [os.path.basename(f).split("_multi")[0] for f in in_files]
        if test_mode:
            print(sequencer, lab, "Prefixes:", prefixes, flush=True)
        if samples_preset:
            prefixes = [p for p in prefixes if any(sub in p for sub in samples_preset)]
        for prefix in prefixes:
            my_list.append((sequencer, lab, prefix))
    return my_list


def get_subdirectories(path):
    return [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]


def get_genomic_info(file_info: tuple):
    sequencer, lab, prefix = file_info
    in_dir = os.path.join(parent_in_dir, sequencer, lab)
    out_dir = os.path.join(parent_out_dir, sequencer, lab)
    os.makedirs(out_dir, exist_ok=True)

    # columns in multi_labels.tsv: ["ReadID", "genome_algorithm..genome_algorithm", "h", "nh", "kraken_h", "kraken_nh"]
    multi_labels_df = pl.scan_csv(os.path.join(in_dir, prefix + "_multi_labels.tsv"),
                                  sep="\t")
    resolved_df = pl.scan_csv(os.path.join(in_dir, prefix + "_resolved.tsv"),
                              sep="\t",
                              dtypes={"ReadID": pl.Utf8, "Resolved": pl.Int32})
    multi_labels_df = multi_labels_df.select(
        pl.all().exclude("^(non_)?(kraken_)?n?h$")
    ).join(
        resolved_df, on="ReadID", how="left"
    ).collect().fill_null(1).select([
        # Change all kraken2 columns from *_kraken2 to kraken2_*
        pl.all().exclude("^.*__kraken2.*$"),
        pl.col("^.*__kraken2.*$").map_alias(lambda col_name: "kraken2_" + col_name.replace("__kraken2", ""))
    ]).lazy()

    # If none of the methods gives P and the true lable is P, 
    # mark genome as FN (column "genome").
    # If kraken gives N and the true label is P,
    # mark "genome_kraken2" as FN.
    for genome in ref_genomes:
        multi_labels_df = multi_labels_df.with_column(
            pl.concat_list(pl.col("^" + genome + ".*$")).alias(genome + "_non_kraken")
        ).with_columns([
            (
                (pl.col("Resolved") == 1) &
                pl.col("kraken2_" + genome).is_in(host_taxid_hierarchy).is_not()
            ).cast(pl.Int8).alias(genome + "_kraken2"), # Adding "T2T_kraken2", "hg38_kraken2", "YH_kraken2"
            (
                (pl.col("Resolved") == 1) &
                pl.col("kraken2_" + genome).is_in(host_taxid_hierarchy).is_not() &
                pl.col(genome + "_non_kraken").arr.eval(
                    (pl.element() == 0).all(),
                    parallel=True
                ).flatten()
            ).cast(pl.Int8).alias(genome) # Adding "T2T", "hg38", "YH"
        ])

    # Filtering
    res_df = multi_labels_df.with_columns([
        pl.concat_list(
            pl.col(ref_genomes)
        ).arr.eval(
            (pl.element() == 1).any(),
            parallel=True # Too many columns. Insufficient memory to parallelize
        ).flatten().cast(pl.Boolean).alias("non_kraken"),
        pl.concat_list(
            pl.col([g + "_kraken2" for g in ref_genomes])
        ).arr.eval(
            (pl.element() == 1).any(),
            parallel=True # Too many columns. Insufficient memory to parallelize
        ).flatten().cast(pl.Boolean).alias("kraken")
    ]).filter(
        pl.col("non_kraken") | pl.col("kraken")
    ).select([
        "ReadID",
        pl.col(ref_genomes), 
        pl.col("^kraken2_.*$").map_alias(
            lambda col_name: col_name.replace("kraken2_", "") + "_kraken2_TaxID"
        ),
        pl.col("T2T_non_kraken").arr.eval(
            (pl.element() == 1).any(),
            parallel=True
        ).flatten().cast(pl.Boolean).alias("T2T_any_alignment") # Used to filter reads to only keep those that have alignment records, in order to save time in parsing the alignment files
    ]).collect()

    # Collect Length, GC, Chromosome, and Position for each read in res_df
    res_df = res_df.join(
        get_read_stats(
            file_info,
            res_df
        ), on="ReadID", how="left"
    ).join(
        get_genome_loc(
            file_info,
            res_df.filter(pl.col("T2T_any_alignment")).select("ReadID")
        ), on="ReadID", how="left"
    )

    # Output
    # header = ["T2T", "hg38", "YH", 
    # "T2T_kraken2_TaxID", "hg38_kraken2_TaxID", "YH_kraken2_TaxID", 
    # "Length", "GC", "Chr", "Pos"]
    res_df = res_df.select(pl.all().exclude("T2T_any_alignment"))
    my_cols = ["Length", "GC", "Chr", "Pos", "ReadID"]
    out_file = os.path.join(out_dir, prefix + "_genomic_stats.tsv")
    if (not os.path.exists(out_file)) or overwrite:
        res_df.select([
            pl.all().exclude("^.*kraken2.*$").exclude(my_cols),
            pl.col([g + "_kraken2_TaxID" for g in ref_genomes]),
            pl.col(my_cols)
        ]).write_csv(out_file, sep="\t")
    else:
        res_df.select([
            pl.all().exclude("^.*kraken2.*$").exclude(my_cols),
            pl.col([g + "_kraken2_TaxID" for g in ref_genomes]),
            pl.col(my_cols)
        ]).write_csv(out_file + ".new", sep="\t")
        if filecmp.cmp(out_file, out_file + ".new"):
            os.remove(out_file + ".new")
        else:
            print("Discrepancy found in results\nNew results were written to " +
                  out_file + ".new", flush=True)
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("[" + now + "] " + sequencer + "/" + lab + "/" + prefix + " done", flush=True)


def get_read_stats(file_info: tuple, df: pl.DataFrame) -> pl.DataFrame:
    sequencer, lab, prefix = file_info
    out_dir = os.path.join(parent_out_dir, sequencer, lab)
    reads_file = os.path.join(out_dir, prefix + ".reads")
    df.select("ReadID").write_csv(reads_file)
    fq_file = os.path.join(parent_qc_dir, sequencer, lab, prefix + ".clean.fq.gz")
    cmd = f"seqkit fx2tab <(seqkit grep -f {reads_file} {fq_file}) -l -g -n -i -H"
    cp = subprocess.run(cmd,
                        shell=True,
                        executable='/bin/bash',
                        capture_output=True,
                        encoding="utf-8")
    stats = pl.read_csv(StringIO(cp.stdout),
                        sep="\t",
                        new_columns=["ReadID", "Length", "GC"],
                        dtypes=[pl.Utf8, pl.Int32, pl.Float32])
    os.remove(reads_file)
    return stats


def get_genome_loc(file_info: tuple, df: pl.DataFrame) -> pl.DataFrame:
    sequencer, lab, prefix = file_info
    loc = pl.DataFrame()
    for a in alignment_algorithms:
        if not df.shape[0]:
            break
        file_path = os.path.join(parent_alignment_dir, sequencer, lab, "T2T", a)
        if not os.path.exists(file_path):
            continue
        if sequencer == "ONT":
            alignment_df = pl.scan_csv(
                os.path.join(file_path, prefix + ".PAF"),
                sep="\t",
                has_header=False,
                dtypes={"column_1": pl.Utf8, "column_2": pl.Int32, "column_3": pl.Int32,
                        "column_4": pl.Int32, "column_5": pl.Utf8, "column_6": pl.Utf8,
                        "column_7": pl.Utf8, "column_8": pl.Int32, "column_9": pl.Int32,
                        "column_10": pl.Utf8, "column_11": pl.Int32, "column_12": pl.Int32}
            ).select([
                pl.col("column_1").alias("ReadID"),
                pl.col("column_6").apply(lambda s: T2T_acc2chr[s]).alias("Chr"),
                pl.col("column_8").alias("Pos")
            ]).unique(subset="ReadID").collect()
        else:
            alignment_df = pl.scan_csv(
                os.path.join(file_path, prefix + ".txt"),
                sep="\t",
                has_header=False,
                dtypes={"column_1": pl.Utf8, "column_2": pl.Int32, "column_3": pl.Utf8,
                        "column_4": pl.Int32, "column_5": pl.Int32, "column_6": pl.Utf8,
                        "column_7": pl.Int32}
            ).select([
                pl.col("column_1").alias("ReadID"),
                pl.col("column_3").apply(lambda s: T2T_acc2chr[s]).alias("Chr"),
                pl.col("column_4").alias("Pos")
            ]).unique(subset="ReadID").collect()
        loc = loc.vstack(df.join(alignment_df, on="ReadID", how="inner"))
        df = df.join(alignment_df, on="ReadID", how="anti")
    return loc


if __name__ == "__main__":
    main()
