#!/usr/bin/env python

import os
import glob
import re
import subprocess
import traceback
import argparse
import concurrent.futures as cf
from datetime import datetime
from functools import partial
from itertools import cycle
from simplesam import Reader, Writer

# "/nas/Analyse/wanglei/T2T/"
project_dir = "./projects/T2T/"
parent_qc_dir = os.path.join(project_dir, "analysis/QC/")
parent_out_dir = os.path.join(project_dir, "analysis/host_reads_removal/")

sequencers = {"ILMN"}  # ("ONT", "ILMN", "MGI")
short_read_methods = {
    "bowtie2_vf": "bowtie2 --very-fast -t -p 8 -x ./T2T/ref_DB/indices/bowtie2/{0}/{0} -U {1} --threads 8 -S {2}.sam --no-unal -k 1 2> /dev/null",
    "bowtie2_vs": "bowtie2 --very-sensitive -t -p 8 -x ./T2T/ref_DB/indices/bowtie2/{0}/{0} -U {1} --threads 8 -S {2}.sam --no-unal -k 1 2> /dev/null",
    "bwa_mem": "bwa mem -t 8 ./T2T/ref_DB/indices/bwa/{0}/{0} {1} -o {2}.sam > /dev/null"}
long_read_methods = {
    "minimap2": "minimap2 -x map-ont -t 8 --secondary=no ./T2T/ref_DB/indices/minimap2/{0}/{0}.mmi {1} -o {2}.PAF",
    "winnowmap": "winnowmap -W ./T2T/ref_DB/indices/winnowmap/{0}/{0}_rep_k15.txt -x map-ont -t 8 --secondary=no ./T2T/ref_DB/genomes/{0}.fna {1} -o {2}.PAF"}
ref_genomes = {"hg38", "T2T", "YH"}
labs_preset = {"Lab49"}
max_procs = 4
alignment = False
post_processing = True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dryrun", default=False, action="store_true",
                        help='Prints the commands to submit, without actually running them')
    args = parser.parse_args()
    for sequencer in sequencers:
        labs = get_subdirectories(os.path.join(parent_qc_dir, sequencer))
        if labs_preset:
            labs = [l for l in labs if any(sub in l for sub in labs_preset)]
        for lab in labs:
            scripts = []
            if sequencer in ("ILMN", "MGI"):
                methods = short_read_methods
            elif sequencer in ("ONT",):
                methods = long_read_methods
            else:
                methods = {}
            for genome in ref_genomes:
                for method, cmd in methods.items():
                    in_files = glob.glob(os.path.join(parent_qc_dir, sequencer, lab, "*.clean.fq.gz"))
                    out_dir = os.path.join(parent_out_dir, sequencer, lab, genome, method)
                    subprocess.run("mkdir -p {}".format(out_dir), shell=True, check=True)
                    samples = [str(os.path.basename(in_file).split('.clean')[0]) for in_file in in_files]
                    out_prefixes = [os.path.join(out_dir, sample) for sample in samples]
                    scripts.extend([cmd.format(*args) for args in zip(cycle([genome]), in_files, out_prefixes)])

            if args.dryrun:
                for script in scripts:
                    print(script)
            elif alignment:
                with cf.ProcessPoolExecutor(max_workers=max_procs) as executor:
                    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    print("[" + now + "] Start processing " + sequencer + "/" + lab)
                    res = executor.map(partial(subprocess.run, shell=True, check=True), scripts)
                    try:
                        for _ in res:
                            pass
                    except Exception as exc:
                        traceback.print_exc()
                    else:
                        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                        print("[" + now + "] Alignment complete")

            if post_processing and sequencer in ("ILMN", "MGI"):
                paths = [os.path.join(sequencer, lab, genome, method) for genome in ref_genomes for method in
                         methods.keys()]
                with cf.ProcessPoolExecutor(max_workers=max_procs) as executor:
                    res = executor.map(get_sam_stats, paths)
                    try:
                        for _ in res:
                            pass
                    except Exception as exc:
                        traceback.print_exc()
                    else:
                        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                        print("[" + now + "] Postprocessing complete")


def get_sam_stats(path):
    sams = glob.glob(os.path.join(parent_out_dir, path, "*.sam"))
    # "./T2T/analysis/host_reads_removal/ILMN/clinical/T2T/bwa_mem/D1-1-U83_S1_R1_001.sam"
    for sam_file in sams:
        sample = str(re.split('.sam', os.path.basename(sam_file))[0])
        hits = []
        with open(sam_file, "r") as in_fp:
            in_sam = Reader(in_fp)
            for line in in_sam:
                if not line.mapped:
                    continue
                qname = line.safename
                seqlen = len(line.seq) + sum(map(int, re.findall('(\\d+)H', line.cigar)))
                nm = line.tags["NM"]
                hits.append((qname, line.flag, line.rname, line.pos, seqlen, line.cigar, nm))
        out_file = os.path.join(parent_out_dir, path, sample + ".txt")
        with open(out_file, "w") as out_fp:
            for hit in hits:
                out_fp.write("\t".join(map(str, hit)) + "\n")
        os.remove(sam_file)


def get_subdirectories(path):
    return [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]


if __name__ == "__main__":
    main()
