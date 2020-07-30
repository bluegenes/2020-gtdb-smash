"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
"""

import os
import glob
import pandas as pd
from collections import defaultdict

out_dir = config["out_dir"]
sigs_dir = config["sigs_dir"]
# need to get input type of files used to generate sigs
input_type = config.get("input_type", "protein")
gather_scaled = config.get("gather_scaled", 100)

gather_dir =  os.path.join(out_dir, "gather")
summary_dir =  os.path.join(out_dir, "gather_tophits")
lca_classify_dir =  os.path.join(out_dir, "lca_classify")
lca_summarize_dir =  os.path.join(out_dir, "lca_summarize")
logs_dir = os.path.join(out_dir, "logs")
envs_dir = "envs"

sampleInfo=config["samples"]
accession2filenames = defaultdict(dict)
accession2signame = defaultdict(dict)

sig_targets=[]

db_names = config.get("db_name", {})
db_names = config.get("db_name", {})
gather_targets, lca_classify_targets, lca_summarize_targets=[],[],[]

# for protein signatures, multipy by 3 if necessary before calculating signature (sourmash v3.x)
ksize_multiplier = {"nucleotide": 1, "dna": 1, "protein": 3, "dayhoff": 3, "hp":3, "translate_protein":3, "translate_dayhoff":3, "translate_hp": 3}
moltype_map = {"nucleotide": "dna", "protein":"protein", "dayhoff": "dayhoff", "hp":"hp", "translate_protein": "protein", "translate_dayhoff": "dayhoff", "translate_hp": "hp"}

# first, build reference dictionary
refInfo = config["references"]
# refInfo.keys = all names
#for ref, info in refInfo.items():
#    ref
# name: {lineages csv, k21: {sbt, lca}, k31: {sbt, lca}, k51: {sbt, lca}}

# using the same configfile that can be used for other stuff.. just grab accessions, signames, filenames; build sigs
for sample, info in sampleInfo.items():
    #input_type = info["input_type"] #protein, dna, rna
    gather_csv = pd.read_csv(info["gather_csv"])
    gather_csv["signame"] = gather_csv["accession"] + " " + gather_csv["species"]
    accession2signame = pd.Series(gather_csv.signame.values,index=gather_csv.accession).to_dict()
    # create "gather_accessions" list
    sampleInfo[sample]["accessions"] = gather_csv.accession.to_list()
    for alpha, alphainfo in info["alphabet"].items():
        ref_alpha=alpha
        if "translate" in alpha:
            ref_alpha = alpha.split("translate_")[1]
        ksizes = alphainfo["ksizes"]
        #gather_targets=expand(os.path.join(gather_dir, "{genome}_x_{db_name}.gather_tophits.csv"), genome=accession2signame.keys(), db_name = refInfo.keys())
        gather_targets+=expand(os.path.join(summary_dir, "{sample}_x_{db_name}.{alphabet}-k{ksize}.gather_tophits.csv"), sample=sample, db_name = refInfo[ref_alpha].keys(), alphabet=alpha, ksize=ksizes)
        #lca_classify_targets= expand(os.path.join(lca_classify_dir,"{sample}_x_{db_name}.{alphabet}-k{ksize}.lca-classify.csv"), sample=sample, db_name = refInfo.keys())
        #lca_summarize_targets= expand(os.path.join(lca_summarize_dir,"{sample}_x_{db_name}.{alphabet}-k{ksize}.lca-summarize.csv"), sample=sample, db_name = refInfo.keys())

#import pdb;pdb.set_trace()
rule all:
    input: gather_targets # lca_classify_targets + lca_summarize_targets

# gather each sig
rule gather_sig:
    input:
        # TARA_IOS_MAG_00042.fa.faa.sig or # TARA_ANE_MAG_00001_protein_scaled100_k11.sig
        query= lambda w: os.path.join(sigs_dir, input_type, "{alphabet}", "k{ksize}", "{genome}_{alphabet}_scaled100_k{ksize}.sig"), 
        #query= os.path.join(sigs_dir, "{genome}.faa.sig")
        db= lambda w: refInfo[moltype_map[w.alphabet]][w.db_name][f"k{w.ksize}"]["sbt"]
    output:
        csv = os.path.join(gather_dir, "{db_name}.{alphabet}-k{ksize}","{genome}_x_{db_name}.{alphabet}-k{ksize}.gather.csv"),
        matches = os.path.join(gather_dir, "{db_name}.{alphabet}-k{ksize}","{genome}_x_{db_name}.{alphabet}-k{ksize}.gather.matches"),
        unassigned = os.path.join(gather_dir, "{db_name}.{alphabet}-k{ksize}", "{genome}_x_{db_name}.{alphabet}-k{ksize}.gather.unassigned")
    params:
        alpha_cmd = lambda w: "--" + moltype_map[w.alphabet],
        ksize = lambda w: int(w.ksize) * ksize_multiplier[w.alphabet],
        scaled = gather_scaled,
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=200,
    log: os.path.join(logs_dir, "gather", input_type, "{alphabet}-k{ksize}",  "{genome}_x_{db_name}.{alphabet}-k{ksize}.gather.log")
    benchmark: os.path.join(logs_dir, "gather", input_type,  "{alphabet}-k{ksize}", "{genome}_x_{db_name}.{alphabet}-k{ksize}.gather.benchmark")
    conda: "envs/sourmash-dev.yml"
    group: "gather"
    shell:
        # --ignore-abundance to turn abund off
        """
        sourmash gather {input.query} {input.db} -o {output.csv} {params.alpha_cmd} \
        --save-matches {output.matches} --threshold-bp=0  \
        --output-unassigned {output.unassigned} \
        --scaled {params.scaled} \
        -k {params.ksize} 2> {log}
        touch {output}
        """
# touch empty output file to enable rna ones to fail (to do: handle failures properly downstream)

rule gather_to_tax:
    input:
        #gather_csv = rules.gather_sig.output.csv,
        gather_csv = os.path.join(gather_dir, "{db_name}.{alphabet}-k{ksize}","{genome}_x_{db_name}.{alphabet}-k{ksize}.gather.csv"),
        #lineages_csv = lambda w: refInfo[w.db_name]["lineages_csv"]
        lineages_csv = lambda w: refInfo[moltype_map[w.alphabet]][w.db_name]["lineages_csv"]
    output:
        gather_tax = os.path.join(gather_dir, "{db_name}.{alphabet}-k{ksize}", "{genome}_x_{db_name}.{alphabet}-k{ksize}.gather_summary.csv"),
        top_matches = os.path.join(gather_dir, "{db_name}.{alphabet}-k{ksize}", "{genome}_x_{db_name}.{alphabet}-k{ksize}.gather_tophits.csv")
    log: os.path.join(logs_dir, "gather_to_tax", "{db_name}.{alphabet}-k{ksize}", input_type, "{genome}_x_{db_name}.{alphabet}-k{ksize}.gather-to-tax.log")
    benchmark: os.path.join(logs_dir, "gather_to_tax", "{db_name}.{alphabet}-k{ksize}","{genome}_x_{db_name}.{alphabet}-k{ksize}.gather-to-tax.benchmark")
    group: "gather"
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=200,
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        python scripts/gather-to-tax.py {input.gather_csv} {input.lineages_csv} --tophits-csv {output.top_matches} > {output.gather_tax} 2> {log}
        """


##### traverse directory is troublesome here --> each set of gather requires it's own folder
rule aggregate_gather_to_tax:
    # make spreadsheet: each proteome:: top lineage hit
    input:
        gather_tophits= lambda w: expand(os.path.join(gather_dir, "{{db_name}}.{{alphabet}}-k{{ksize}}", "{genome}_x_{{db_name}}.{{alphabet}}-k{{ksize}}.gather_tophits.csv"), genome=  sampleInfo[w.sample]["accessions"])
    output:
        summary_csv=os.path.join(summary_dir, "{sample}_x_{db_name}.{alphabet}-k{ksize}.gather_tophits.csv"),
    params:
        gather_dir= os.path.join(gather_dir, "{db_name}.{alphabet}-k{ksize}"),
        true_lineages_cmd = lambda w: "--true-lineages-csv " + sampleInfo[w.sample].get("true_lineages_csv") if sampleInfo[w.sample].get("true_lineages_csv") else "",
    #log: os.path.join(logs_dir, "aggregate_gather", "{sample}_x_{db_name}." + f"{alphabet}-k{ksize}.gather_tophits.log")
    log: os.path.join(logs_dir, "aggregate_gather", "{sample}_x_{db_name}.{alphabet}-k{ksize}.gather_tophits.log")
    benchmark: os.path.join(logs_dir, "aggregate_gather", "{sample}_x_{db_name}.{alphabet}-k{ksize}.gather_tophits.benchmark")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=60,
    conda: "envs/forage-env.yml"
    shell:
        #python scripts/aggregate-gather-to-tax-tophits.py --input-is-directory --output-csv {output} {params.gather_dir} --true-lineages-csv {input.true_lineages} 2> {log}
        """
        python scripts/aggregate-gather-to-tax-tophits.py --input-is-directory --output-csv {output} {params.gather_dir} {params.true_lineages_cmd} 2> {log}
        """


'''
## TRAVERSE DIRECTORY IS AN ISSUE HERE --> maybe cp sigs into temporary directory, traverse, then rm?
## Can't pass all sigs into the command line - need to do all at once. 

rule cp_sigs:
    input:
        lambda w: os.path.join(sigs_dir, sampleInfo[w.sample]["input_type"], "{alphabet}", "k{ksize}", "{genome}_{alphabet}_scaled{scaled}_k{ksize}.sig")
    output:
        os.path.join(lca_classify_dir, "{sample}", "{genome}_{alphabet}_scaled{scaled}_k{ksize}.sig")
    shell:
        """
        cp {input} {output}
        """


# lca classify can do all sigs at once
rule lca_classify_sigs:
    input:
        #sigs= expand(os.path.join(sigs_dir, "{genome}.faa.sig"), genome=accessions), # TARA_IOS_MAG_00042.fa.faa.sig
        sigs = lambda w: expand(os.path.join(sigs_dir, "{input_type}", moltype_map[w.alphabet], "k{ksize}", "{genome}_{alphabet}_scaled{scaled}_k{ksize}.sig"), genome=accessions),
        #db= lambda w: refInfo[w.db_name][f"k{w.ksize}"]["lca"]
        db = lambda w: refInfo[moltype_map[w.alphabet]][w.db_name]["k{w.ksize}"]["lca"]
    output:
        csv = os.path.join(lca_classify_dir, "{sample}_x_{db_name}.{alphabet}-k{ksize}.lca-classify.csv"),
    params:
        sigs_dir = sigs_dir,
        #alpha_cmd = lambda w: "--" + alphabet,
        #ksize = int(ksize) * ksize_multiplier[alphabet],
    resources:
        mem_mb=lambda wildcards, attempt: attempt *20000,
        runtime=200,
    log: os.path.join(logs_dir, "lca-classify", "{sample}_x_{db_name}.{alphabet}-k{ksize}.lca-classify.log")
    benchmark: os.path.join(logs_dir, "lca-classify", "{sample}_x_{db_name}.{alphabet}-k{ksize}.lca-classify.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        # --ignore-abundance to turn abund off
        #sourmash lca classify --query {input.query} --db {input.db} \
        #-o {output.csv} {params.alpha_cmd} \
        #-k {params.ksize} 2> {log}
        """
        sourmash lca classify --query {params.sigs_dir} \
        --traverse-directory --db {input.db} \
        -o {output.csv} --threshold 0  2> {log}
        """

rule lca_summarize_sigs:
    input:
        #sigs= expand(os.path.join(sigs_dir, "{genome}.faa.sig"), genome=accessions), # TARA_IOS_MAG_00042.fa.faa.sig
        sigs = expand(os.path.join(sigs_dir, f"{input_type}", f"{alphabet}", f"k{ksize}", "{genome}_" + f"{alphabet}_scaled{scaled}_k{ksize}.sig"), genome=accessions),
        db= lambda w: db_names[w.db_name]["lca"]
    output:
        csv = os.path.join(lca_summarize_dir, "{sample}_x_{db_name}.{alphabet}-k{ksize}.lca-summarize.csv"),
    params:
        sigs_dir = sigs_dir,
    resources:
        mem_mb=lambda wildcards, attempt: attempt *20000,
        runtime=200,
    log: os.path.join(logs_dir, "lca-classify", "{sample}_x_{db_name}.{alphabet}-k{ksize}.lca-summarize.log")
    benchmark: os.path.join(logs_dir, "lca-classify", "{sample}_x_{db_name}.{alphabet}-k{ksize}.lca-summarize.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash lca summarize  --query {params.sigs_dir} \
        --traverse-directory --singleton \
        --db {input.db} \
        -o {output.csv} --threshold 0  2> {log}
        """
'''
