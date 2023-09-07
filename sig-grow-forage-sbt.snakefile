"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s grow-and-forage-sbt.snakefile --configfile config/grow_gtdb_sbts.yml --use-conda
"""

import os
import glob
import pandas as pd

def try_reading_csv(groups_file):
    # autodetect format
    if '.tsv' in groups_file or '.csv' in groups_file:
        separator = '\t'
        if '.csv' in groups_file:
            separator = ','
        try:
            samples = pd.read_csv(groups_file, dtype=str, sep=separator)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {groups_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    elif '.xls' in groups_file:
        try:
            samples = pd.read_excel(groups_file, dtype=str, sep=separator)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {groups_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    return samples

out_dir = config.get("out_dir", "output-pathdistances")
logs_dir = os.path.join(out_dir, "logs")
index_dir = config.get("index_dir", "gtdb-sbts")
dist_dir = os.path.join(out_dir, "distances")

sampleInfo=config["foraging_samples"]

sbt_targets=[]
query_targets=[]

for sample, info in sampleInfo.items():    
    if info.get("grow_sbt"):
        # build sbt targets
        index_dir= info["grow_sbt"]["sbt_outdir"]
        for alpha, alphainfo in info["alphabet"].items():
            sbt_targets+=expand(os.path.join(index_dir,"{sample}.{alphabet}_scaled{scaled}_k{k}.sbt.zip"), sample=sample, alphabet=alpha, scaled=alphainfo["scaled"], k=alphainfo["ksizes"])
    for alpha, alphainfo in info["alphabet"].items():
        query_targets+=expand(os.path.join(dist_dir, "{sample}.{alphabet}_scaled{scaled}_k{k}.jaccard_from_species.csv"), sample=sample, alphabet=alpha, scaled=alphainfo["scaled"], k=alphainfo["ksizes"])

rule all:
    input: sbt_targets + query_targets

rule grow_sbt:
    output: 
        sbt=os.path.join(index_dir,"{sample}.{alphabet}_scaled{scaled}_k{k}.sbt.zip"),
    threads: 1
    params:
        input_dir= lambda w: sampleInfo[w.sample]['grow_sbt']["input_path"],
        subset_csv=lambda w: sampleInfo[w.sample]['grow_sbt'].get('subset_csv', ''),
        subset_info_colname=lambda w: sampleInfo[w.sample]['grow_sbt'].get('subset_info_colname', 'filename'),
        alpha= lambda w: w.alphabet.rsplit("translate_")[1] if w.alphabet.startswith("translate") else w.alphabet, # remove translate
        translate = lambda w: " --translate " if w.alphabet.startswith("translate") else "",
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=6000,
    log: os.path.join(logs_dir, "grow-sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.grow.log")
    benchmark: os.path.join(logs_dir, "grow-sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.grow.benchmark")
    conda: "envs/forage-env.yml"
    shell:
        # escaped quotes allows for "*" in input dir
        """
        python scripts/grow-sbtmh.py \"{params.input_dir}\" --input-is-directory --sbt {output.sbt} --ksize {wildcards.k} --scaled {wildcards.scaled} --alphabet {params.alpha} --subset-csv {params.subset_csv} --subset-info-colname {params.subset_info_colname} {params.translate} 2> {log}
        """

def find_forage_inputs(w):
    # just require this query to have filenames of interest as a column
    query = sampleInfo[w.sample]["query_csv"]
    if sampleInfo[w.sample].get("grow_sbt", False):
        sbt = rules.grow_sbt.output.sbt,
    else:
        sbt=sampleInfo[w.sample]["use_existing_sbt"]
    return {"query_csv": query, "sbt": sbt}


rule calculate_jaccard_from_common_ancestor:
    input: unpack(find_forage_inputs)
    output: 
        csv=os.path.join(dist_dir, "{sample}.{alphabet}_scaled{scaled}_k{k}.jaccard_from_species.csv"),
        boxplot=os.path.join(dist_dir, "plots", "{sample}.{alphabet}_scaled{scaled}_k{k}.jaccard_from_species.svg"),
    params:
        signature_name_column= lambda w: sampleInfo[w.sample].get('query_signature_name_column_name', 'filename')
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=1200,
    log: os.path.join(logs_dir, "forage", "{sample}.{alphabet}_scaled{scaled}_k{k}.forage.log")
    benchmark: os.path.join(logs_dir, "forage", "{sample}.{alphabet}_scaled{scaled}_k{k}.forage.benchmark")
    conda: "envs/forage-env.yml"
    shell:
        """
        python scripts/forage-sbt.py {input.sbt} --query-csv {input.query_csv} --signature-name-column {params.signature_name_column} \
        --distance-from-species-csv {output.csv} --distance-from-species-plot {output.boxplot} 2>{log}
        """

