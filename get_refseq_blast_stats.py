#!/bin/usr/env python

import fire
import traceback
import subprocess
import polars as pl
from collections import OrderedDict
from itertools import islice
from pathlib import Path

project_dir = Path("./projects/T2T")
refseq_dir = project_dir / "ref_DB/genomes/RefSeq"
res_dir = project_dir / "analysis/cmp_gnm/blast/"


class SlicableOrderedDict(OrderedDict):
    def __getitem__(self, k):
        if not isinstance(k, slice):
            return OrderedDict.__getitem__(self, k)
        return SlicableOrderedDict(islice(self.items(), k.start, k.stop))


def main(hseq="Winnowmap", mtype="viral", head=0):
    genome2taxid = get_genome_info_mapping(refseq_dir / mtype)
    if head != 0:
        genome2taxid = genome2taxid[:head]
    all_df = pl.DataFrame()
    for genome, taxid in genome2taxid.items():
        df = pl.DataFrame()
        try:
            df = pl.read_csv(Path(res_dir, hseq, mtype, genome + ".m8"),
                             sep="\t", has_header=False, null_values="N/A",
                             columns=[0, 2, 3, 6, 7, 10, 11],
                             new_columns=["contig", "identity", "length", "qstart", "qend", "evalue", "bitscore"],
                             dtypes={"contig": pl.Utf8, "identity": pl.Float32,
                                     "length": pl.Int32, "qstart": pl.Int32, "qend": pl.Int32,
                                     "evalue": pl.Float32, "bitscore": pl.Float32}
                             )
            """columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", 
            "qend", "sstart", "send", "evalue", "bitscore"] """
            """dtypes={"column_1": pl.Utf8, "column_2": pl.Utf8, "column_3": pl.Float32,
                     "column_4": pl.Int32, "column_5": pl.Int32, "column_6": pl.Int32,
                     "column_7": pl.Int32, "column_8": pl.Int32, "column_9": pl.Int32,
                     "column_10": pl.Int32, "column_11": pl.Float32,
                     "column_12": pl.Float32}"""
        except pl.exceptions.NoDataError:
            df = pl.DataFrame({"contig": [None], "identity": [None],
                               "length": [None], "qstart": [None], "qend": [None], "evalue": [None],
                               "bitscore": [None]},
                              columns={"contig": pl.Utf8, "identity": pl.Float32,
                                       "length": pl.Int32, "qstart": pl.Int32, "qend": pl.Int32,
                                       "evalue": pl.Float32, "bitscore": pl.Float32})
        finally:
            # all_df = all_df.concat(df.select([ "column_1", "column_3", "column_4", "column_7", "column_8"
            # ]).rename({ "column_1": "contig", "column_3": "identity", "column_4": "length", "column_7": "qstart",
            # "column_8": "qend" })
            all_df = all_df.vstack(df.with_columns([
                pl.lit(genome).alias("genome"),
                pl.lit(taxid).alias("taxid")
            ]))
    all_df = all_df.rechunk().groupby(
        ["taxid", "genome", "contig", "qstart", "qend"],
        maintain_order=True
    ).agg([
        pl.col("length").first(),
        pl.col("identity").first(),
        (pl.count() - pl.col("contig").null_count()).alias("hits"),
        pl.col("evalue").first(),
        pl.col("bitscore").first()
    ])
    contig_length_file = refseq_dir / mtype / "contig_length.txt"
    if not contig_length_file.exists():
        subprocess.run("seqkit fx2tab -inlHj 5 *.fna.gz > contig_length.txt", shell=True, cwd=refseq_dir / mtype)
    contig_info_df = pl.read_csv(
        contig_length_file, sep="\t",
        has_header=True, new_columns=["contig", "contig_size"]
    )
    all_df.join(
        contig_info_df, on="contig", how="left"
    ).select([
        "taxid", "genome", "contig", "contig_size", "qstart", "qend", "length", "identity", "hits", "evalue", "bitscore"
    ]).write_csv(
        res_dir / hseq / (mtype + "_blast_stats.tsv"), sep="\t", has_header=True, null_value="NA"
    )


def get_genome_info_mapping(in_dir):
    summary_file = in_dir / "assembly_summary.txt"
    genome2taxid = SlicableOrderedDict()
    with open(summary_file, "r") as fi:
        for line in fi:
            if line.startswith("#"):
                continue
            else:
                fields = line.strip().split("\t")
                genome2taxid[fields[0]] = fields[5]
    return genome2taxid


if __name__ == "__main__":
    fire.Fire(main)
