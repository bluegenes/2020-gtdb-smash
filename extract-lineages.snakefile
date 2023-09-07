"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s gtdb-smash.snakefile --use-conda
"""

import os

lineage_csv=config["lineage_csv"]
gtdb_datadir = config.get("gtdb_datadir", "/home/ctbrown/gtdbtk/release89/fastani/database")
out_dir = config.get("out_dir", "output-gtdb-smash")
data_dir= os.path.join("out_dir", "data")
logs_dir = os.path.join("out_dir", "logs")

rule all:
    input: "final_output.txt"

localrules: build_genome_list, get_genomes

rule build_genome_list:
    input: lineage_csv
    output: os.path.join(out_dir, "lineage-list.txt")
    conda: "envs/sourmash.yml"
    log: os.path.join(logs_dir, "extract-robust-lineages.log")
    shell:
        """
        scripts/extract-robust-lineages.py {input} {output} > {log}
        """
