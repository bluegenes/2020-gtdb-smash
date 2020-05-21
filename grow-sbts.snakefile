"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s grow-sbts.snakefile --configfiles grow_gtdb_sbts.yml
"""

import os
import glob
import pandas as pd
#lineage_csv=config["lineage_csv"]

out_dir = config.get("index_dir", "sbt_index")
logs_dir = os.path.join(out_dir, "logs")


output_files=[]

refInfo=config.get("references")
for refname, info in refInfo.items():
    for alphabet, alphainfo in info["alphabet"].items():
        output_files+=expand(os.path.join(out_dir, "{ref}_{alpha}_k{ksize}_scaled{scaled}.sbt.zip"), ref=refname, alpha=alphabet, ksize=alphainfo["ksizes"], scaled=alphainfo["scaled"])


rule all:
    input: output_files

rule grow_sbt_from_fasta:
    input: lambda w: refInfo[w.db_name]["path"]
    output: os.path.join(out_dir,"{db_name}_{alphabet}_k{ksize}_scaled{scaled}.sbt.zip")
    log: os.path.join(logs_dir, "{db_name}_{alphabet}_{scaled}_{ksize}.grow_sbt.log")
    benchmark: os.path.join(logs_dir, "{db_name}_{alphabet}_{scaled}_{ksize}.grow_sbt.benchmark")
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt *20000,
        runtime = 6000,
    conda: "envs/sourmash3.3.yml"
    shell:
        """
        python scripts/grow-sbtmh.py {input} --alphabet {wildcards.alphabet} --ksize {wildcards.ksize} --scaled {wildcards.scaled} --sbt {output} --track-abundance --input-is-directory
        """

