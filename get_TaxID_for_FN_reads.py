#!/bin/env python

import fire
import traceback
import concurrent.futures as cf
import polars as pl
from pathlib import Path
from collections import Counter
from datetime import datetime

project_dir = Path("./projects/T2T/")
qc_dir = Path(project_dir, "analysis/QC/")
data_dir = Path(project_dir, "analysis/per_read/")
kraken_res_dir = Path(project_dir, "analysis/microbe_only_classification/")
ref_genomes = {"hg38", "T2T", "YH"}
sequencers = {"MGI", "ILMN", "ONT"}  # {"ONT", "ILMN", "MGI"}
labs_preset = {}
samples_preset = {}
parallel = True
max_proc = 4

host_taxid_hierarchy = [7711, 89593, 7742, 7776, 117570, 117571, 8287,
                        1338369, 32523, 32524, 40674, 32525, 9347,
                        1437010, 314146, 9443, 376913, 314293, 9526,
                        314295, 9604, 207598, 9605, 9606, 63221, 741158]


def main(dry_run=False):
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
            if samples_preset:
                data_files = [x for x in data_files if any(sub in x.name for sub in samples_preset)]
            if dry_run:
                print([f.name for f in data_files], flush=True)
            else:
                if parallel:
                    with cf.ProcessPoolExecutor(max_workers=max_proc) as executor:
                        res_iter = executor.map(get_microbe_taxid, data_files)
                        try:
                            for r in res_iter:
                                res.append(r)
                        except Exception:
                            traceback.print_exc()
                else:
                    for f in data_files:
                        res.append(get_microbe_taxid(f))
    res_df = pl.DataFrame(columns=[("Taxon", pl.Int64)])
    for r in res:
        for identifier, taxid_list in r.items():
            print("Aggregating " + identifier, flush=True)
            taxon_stats = Counter(taxid_list).most_common()
            print(taxon_stats, flush=True)
            df = pl.DataFrame(taxon_stats,
                              orient="row",
                              columns=[("Taxon", pl.Int64), (identifier, pl.Int64)])
            res_df = res_df.join(df, how="outer", on="Taxon")
    out_dir = Path(kraken_res_dir, "FN")
    out_dir.mkdir(parents=True, exist_ok=True)
    res_df.write_csv(Path(out_dir, "Microbe_report_for_FN_reads.tsv"), sep="\t")
    res_df.melt(
        id_vars="Taxon", variable_name="Sample"
    ).pivot(
        values="value", index="Sample", columns="Taxon"
    ).write_csv(Path(out_dir, "Microbe_report_for_FN_reads_transposed.tsv"), sep="\t")
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("[" + now + "] " + "Analysis completed", flush=True)


def get_microbe_taxid(in_file):
    sequencer = in_file.parts[-3]
    lab = in_file.parts[-2]
    prefix = in_file.name.split("_multi")[0]
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
        taxon_list = df.select([
            pl.col("^" + genome + ".*$"),
            pl.col("Resolved"),
            pl.col("TaxID")
        ]).with_column(
            pl.concat_list(
                pl.all().exclude(genome + "__kraken2").exclude("Resolved").exclude("TaxID")
            ).alias(genome + "__non_kraken")
        ).filter(
            (pl.col("Resolved") == 1) &
            pl.col(genome + "__kraken2").is_in(host_taxid_hierarchy).is_not() &
            pl.col(genome + "__non_kraken").arr.eval(
                (pl.element() == 0).all(),
                parallel=True
            ).flatten()
        ).select("TaxID").collect().to_series().to_list()
        res["#".join([sequencer, lab, prefix, genome])] = taxon_list
    return res


if __name__ == "__main__":
    fire.Fire(main)
