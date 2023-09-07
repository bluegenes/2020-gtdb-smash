"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s classify-sigs_v2.snakefile --use-conda
"""

## to do: 

# modify build-sigs --> take all the alpha types. use "include to leverage it here"
# modify target generation for new config structure in config/classify_nucl.yml

########3
#ACTUALLY, NUKE THIS AND USE GATHER.SNAKEFILE --> GATHER.RULE
# THEN, CLASSIFY.SNAKEFILE JUST INCLUDES APPROPRIATE RULES
# LCA_SUMMARIZE.RULE
# LCA_CLASSIFY.RULE
# ETC !!!!
#####

import os
import glob

out_dir = config["out_dir"]
sigs_dir = config["sigs_dir"]

logs_dir = os.path.join(out_dir, "logs")
envs_dir = "envs"

basename= config["sample_basename"]
genome_info = config.get("genome_list")
#genome_list = [ x.strip() for x in open('tara-delmont-all-list.txt') ]
#genome_list = [ x.strip() for x in open(genome_info) ]
if genome_info:
    genome_list = [ x.strip().rsplit(".fa")[0] for x in open(genome_info) ]
else:
    genome_csv = pd.read_csv(info["genome_csv"])
    # this assumes they're in same input dir & are same input type as the samples above. ok for now. shrug.
    genome_list = gather_csv.accession.to_list()
    
lca_dbs = config.get("lca_db", {})
sbt_dbs = config.get("sbt_db", {})
gather_targets, lca_classify_targets=[],[]
# just enable single alpha-ksize for now
ksize = config.get("ksize", 19)
alphabet = config.get("alphabet", "dayhoff")
scaled= config.get("scaled", 100)
input_type = config.get("input_type", "protein")

# build gather output dir
gather_dir =  os.path.join(out_dir, "gather", f"{alphabet}-k{ksize}")
summary_dir =  os.path.join(out_dir, "gather_tophits")
lca_classify_dir =  os.path.join(out_dir, "lca_classify")
lca_summarize_dir =  os.path.join(out_dir, "lca_summarize")

if sbt_dbs:
    # build gather targets
    # requires a few steps: 1. gather to db, 2. gather-to-tax, 3. aggregate tax results
    gather_targets=expand(os.path.join(gather_dir, "{genome}_x_{db_name}.gather_tophits.csv"), genome=genome_list, db_name = sbt_dbs.keys())
    gather_targets+=expand(os.path.join(summary_dir, "{sample}_x_{db_name}.gather_tophits.csv"), sample=basename, db_name = sbt_dbs.keys())
if lca_dbs:
   # run sourmash lca classify --> lca db
    lca_classify_targets= expand(os.path.join(lca_classify_dir,"{sample}_x_{db_name}.lca-classify.csv"), sample=basename, db_name = lca_dbs.keys())
    lca_summarize_targets= expand(os.path.join(lca_summarize_dir,"{sample}_x_{db_name}.lca-summarize.csv"), sample=basename, db_name = lca_dbs.keys())

rule all:
    input: gather_targets + lca_classify_targets + lca_summarize_targets

# for protein signatures, multipy by 3 if necessary before calculating signature (sourmash v3.x)
ksize_multiplier = {"nucleotide": 1, "dna": 1, "protein": 3, "dayhoff": 3, "hp":3}
moltype_map = {"nucleotide": "dna", "protein":"protein", "dayhoff": "dayhoff", "hp":"hp", "translate_protein": "protein", "translate_dayhoff": "dayhoff", "translate_hp": "hp"}
# gather each sig
rule gather_sig:
    input:
        # TARA_IOS_MAG_00042.fa.faa.sig or # TARA_ANE_MAG_00001_protein_scaled100_k11.sig
        query= os.path.join(sigs_dir, f"{input_type}", f"{alphabet}", f"k{ksize}", "{genome}_" + f"{alphabet}_scaled{scaled}_k{ksize}.sig"), 
        #query= os.path.join(sigs_dir, "{genome}.faa.sig")
        db= lambda w: sbt_dbs[w.sbt_db]["sbt"]
    output:
        csv = os.path.join(gather_dir, "{genome}_x_{sbt_db}.gather.csv"),
        matches = os.path.join(gather_dir, "{genome}_x_{sbt_db}.gather.matches"),
        unassigned = os.path.join(gather_dir, "{genome}_x_{sbt_db}.gather.unassigned")
    params:
        alpha_cmd = lambda w: "--" + moltype_map[alphabet],
        ksize = int(ksize) * ksize_multiplier[alphabet],
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=60,
    log: os.path.join(logs_dir, "gather", "{genome}_x_{sbt_db}.gather.log")
    benchmark: os.path.join(logs_dir, "gather", "{genome}_x_{sbt_db}.gather.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        # --ignore-abundance to turn abund off
        """
        sourmash gather {input.query} {input.db} -o {output.csv} {params.alpha_cmd} \
        --save-matches {output.matches} --threshold-bp=0  \
        --output-unassigned {output.unassigned} \
        -k {params.ksize} 2> {log}
        touch {output}
        """
# touch empty output file to enable rna ones to fail (to do: handle failures properly downstream)

rule gather_to_tax:
    input:
        gather_csv = rules.gather_sig.output.csv,
        lineages_csv = lambda w: sbt_dbs[w.sbt_db]["lineages_csv"]
    output:
        gather_tax = os.path.join(gather_dir, "{genome}_x_{sbt_db}.gather_summary.csv"),
        top_matches = os.path.join(gather_dir, "{genome}_x_{sbt_db}.gather_tophits.csv")
    log: os.path.join(logs_dir, "gather_to_tax", "{genome}_x_{sbt_db}.gather-to-tax.log")
    benchmark: os.path.join(logs_dir, "gather_to_tax", "{genome}_x_{sbt_db}.gather-to-tax.benchmark")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=60,
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        python scripts/gather-to-tax.py {input.gather_csv} {input.lineages_csv} --tophits-csv {output.top_matches} > {output.gather_tax} 2> {log}
        """

rule aggregate_gather_to_tax:
    # make spreadsheet: each proteome:: top lineage hit
    input:
        gather_tophits= lambda w: expand(os.path.join(gather_dir, "{genome}_x_{{sbt_db}}.gather_tophits.csv"), genome=genome_list)
    output:
        summary_csv=os.path.join(summary_dir, "{sample}_x_{sbt_db}.gather_tophits.csv"),
    params:
        gather_dir= gather_dir 
    #log: os.path.join(logs_dir, "aggregate_gather", "{sample}_x_{sbt_db}." + f"{alphabet}-k{ksize}.gather_tophits.log")
    log: os.path.join(logs_dir, "aggregate_gather", "{sample}_x_{sbt_db}.gather_tophits.log")
    benchmark: os.path.join(logs_dir, "aggregate_gather", "{sample}_x_{sbt_db}.gather_tophits.benchmark")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=60,
    conda: "envs/forage-env.yml"
    shell:
        """
        python scripts/aggregate-gather-to-tax-tophits.py --input-is-directory --output-csv {output} {params.gather_dir} 2> {log}
        """


# lca classify can do all sigs at once
rule lca_classify_sigs:
    input:
        #sigs= expand(os.path.join(sigs_dir, "{genome}.faa.sig"), genome=genome_list), # TARA_IOS_MAG_00042.fa.faa.sig
        sigs = expand(os.path.join(sigs_dir, f"{input_type}", f"{alphabet}", f"k{ksize}", "{genome}_" + f"{alphabet}_scaled{scaled}_k{ksize}.sig"), genome=genome_list),
        db= lambda w: lca_dbs[w.lca_db]["lca"]
    output:
        csv = os.path.join(lca_classify_dir, "{sample}_x_{lca_db}.lca-classify.csv"),
    params:
        sigs_dir = sigs_dir,
        #alpha_cmd = lambda w: "--" + alphabet,
        #ksize = int(ksize) * ksize_multiplier[alphabet],
    resources:
        mem_mb=lambda wildcards, attempt: attempt *20000,
        runtime=200,
    log: os.path.join(logs_dir, "lca-classify", "{sample}_x_{lca_db}.lca-classify.log")
    benchmark: os.path.join(logs_dir, "lca-classify", "{sample}_x_{lca_db}.lca-classify.benchmark")
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
        #sigs= expand(os.path.join(sigs_dir, "{genome}.faa.sig"), genome=genome_list), # TARA_IOS_MAG_00042.fa.faa.sig
        sigs = expand(os.path.join(sigs_dir, f"{input_type}", f"{alphabet}", f"k{ksize}", "{genome}_" + f"{alphabet}_scaled{scaled}_k{ksize}.sig"), genome=genome_list),
        db= lambda w: lca_dbs[w.lca_db]["lca"]
    output:
        csv = os.path.join(lca_summarize_dir, "{sample}_x_{lca_db}.lca-summarize.csv"),
    params:
        sigs_dir = sigs_dir,
    resources:
        mem_mb=lambda wildcards, attempt: attempt *20000,
        runtime=200,
    log: os.path.join(logs_dir, "lca-classify", "{sample}_x_{lca_db}.lca-summarize.log")
    benchmark: os.path.join(logs_dir, "lca-classify", "{sample}_x_{lca_db}.lca-summarize.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash lca summarize  --query {params.sigs_dir} \
        --traverse-directory --singleton \
        --db {input.db} \
        -o {output.csv} --threshold 0  2> {log}
        """
