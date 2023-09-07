"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s gather.snakefile --use-conda
"""

import os
import glob
import pandas as pd
from collections import defaultdict

# for protein signatures, multipy by 3 if necessary before calculating signature (sourmash v3.x)
ksize_multiplier = {"nucleotide": 1, "protein": 3, "dayhoff": 3, "hp":3, "translate_protein": 3, "translate_dayhoff": 3, "translate_hp": 3}
# "rna" is not sourmash cli friendly
moltype_map = {"nucleotide": "dna", "protein":"protein", "dayhoff": "dayhoff", "hp":"hp", "translate_protein": "protein", "translate_dayhoff": "dayhoff", "translate_hp": "hp"}


def aggregate_sigs(w):
    siglist=[]
    input_type = sampleInfo[w.sample]["input_type"] # protein, dna, rna
    sigfile = os.path.join(compute_dir, input_type, w.alphabet, f"k{w.k}", f"{{acc}}_{w.alphabet}_scaled{w.scaled}_k{w.k}.sig")
    siglist=expand(sigfile, acc=sampleInfo[w.sample]["accessions"])
    return siglist


# is there a way  to traverse directory but only use the chosen sigs?
# (similar to --require-taxonomy for lca)
# maybe if no taxonomy is provided, use this rule? Otherwise the other pair?

#rule index_sbt:
#    input: sigs=aggregate_sigs,
#    output: 
#        sbt=os.path.join(index_dir, "sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.index.sbt.zip"),
#    threads: 1
#    params:
#        alpha= lambda w: w.alphabet.rsplit("translate_")[1] if w.alphabet.startswith("translate") else w.alphabet, # remove translate
#        #alpha_cmd= lambda w: " --" + (w.alphabet.rsplit("translate_")[1] if w.alphabet.startswith("translate") else w.alphabet), # remove translate
#        alpha_cmd = lambda w: "--" + moltype_map[w.alphabet],
#        translate = lambda w: " --translate " if w.alphabet.startswith("translate") else "",
#        input_type = lambda w: sampleInfo[w.sample]["input_type"],
#        ksize = lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
#        sigdir= lambda w: os.path.join(compute_dir, sampleInfo[w.sample]["input_type"], w.alphabet, f"k{w.k}")
#    resources:
#        mem_mb=lambda wildcards, attempt: attempt *5000,
#        runtime=6000,
#    log: os.path.join(logs_dir, "index-sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.index-sbt.log")
#    benchmark: os.path.join(logs_dir, "index-sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.index-sbt.benchmark")
#    conda: "envs/forage-env.yml"
#    shell:
#        """
#        sourmash index --ksize {params.ksize} \
#        --scaled {wildcards.scaled} {params.alpha_cmd} \
#        {output.sbt} \
#        {input.sigs} 2> {log}
#        """

rule index_lca:
    input: 
        sigs=aggregate_sigs,
        taxonomy= lambda w: sampleInfo[w.sample]["info_csv"]
    output:
        os.path.join(index_dir,"lca", "{sample}.{alphabet}_scaled{scaled}_k{k}.index.lca.json.gz"),
    threads: 1
    params:
        alpha= lambda w: (w.alphabet.rsplit("translate_")[1] if w.alphabet.startswith("translate") else w.alphabet), # remove translate
        #alpha_cmd= lambda w: " --" + (w.alphabet.rsplit("translate_")[1] if w.alphabet.startswith("translate") else w.alphabet), # remove translate
        alpha_cmd= lambda w: "--" + moltype_map[w.alphabet],
        translate = lambda w: " --translate " if w.alphabet.startswith("translate") else "",
        input_type = lambda w: sampleInfo[w.sample]["input_type"],
        ksize = lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
        report= lambda w: os.path.join(index_dir,"lca", f"{w.sample}.{w.alphabet}_scaled{w.scaled}_k{w.k}.index.lca.report"),
        #output_prefix = lambda w: os.path.join(index_dir,"{w.sample}.{w.alphabet}_scaled{w.scaled}_k{w.k}.index")
        sigdir= lambda w: os.path.join(compute_dir, sampleInfo[w.sample]["input_type"], w.alphabet, f"k{w.k}")
    resources:
        mem_mb= lambda wildcards, attempt: attempt *100000,
        runtime=600,
    log: os.path.join(logs_dir, "index-lca", "{sample}.{alphabet}_scaled{scaled}_k{k}.index-lca.log")
    benchmark: os.path.join(logs_dir, "index-lca", "{sample}.{alphabet}_scaled{scaled}_k{k}.index-lca.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash lca index \
          --ksize {params.ksize} \
          --scaled {wildcards.scaled} \
          --split-identifiers \
          --require-taxonomy \
          --traverse-directory \
          --report {params.report} \
          {params.alpha_cmd} {input.taxonomy} {output} \
          {params.sigdir} 2> {log} 
        """


# lca takes longer to build, but enables you to select sigs within a larger dir via --require-taxonomy
# then can convert the lca to sbt
rule lca_to_sbt:
    input: os.path.join(index_dir, "lca", "{sample}.{alphabet}_scaled{scaled}_k{k}.index.lca.json.gz")
    output: os.path.join(index_dir, "sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.index.sbt.zip")
    log: os.path.join(logs_dir, "lca-to-sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.lca-to-sbt.log")
    benchmark: os.path.join(logs_dir, "lca-to-sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.lca-to-sbt.benchmark")
    resources:
        mem_mb= lambda wildcards, attempt: attempt *100000,
        runtime=600,
    shell:
        """
        python scripts/convert-lca-to-sbt.py {input} {output} 2> {log}
        """
