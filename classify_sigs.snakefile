"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s classify-sigs.snakefile --use-conda
"""
import os
import glob

out_dir = config["out_dir"]
sigs_dir = config["sigs_dir"]

logs_dir = os.path.join(out_dir, "logs")
envs_dir = "envs"

basename= config["sample_basename"]
genome_info = config["genome_info"]
#genome_list = [ x.strip() for x in open('tara-delmont-all-list.txt') ]
genome_list = [ x.strip() for x in open(genome_info) ]

lca_dbs = config.get("lca_db", {})
sbt_dbs = config.get("sbt_db", {})
gather_targets, lca_classify_targets=[],{}
# just enable single alpha-ksize for now
ksize = config.get("ksize", 19)
alphabet = config.get("alphabet", "dayhoff")

# build gather output dir
gather_dir =  os.path.join(out_dir, "gather", f"{alphabet}-k{ksize}")

if sbt_dbs:
    # build gather targets
    # requires a few steps: 1. gather to db, 2. gather-to-tax, 3. aggregate tax results
    gather_targets=expand(os.path.join(gather_dir, "{genome}_x_{db_name}.gather_tophits.csv"), genome=genome_list, db_name = sbt_dbs.keys())
#if lca_dbs:
#    # run sourmash lca classify --> lca db
#    lca_classify_targets= expand(os.path.join(gather_dir,"{sample}_x_{db_name}.lca-classify.csv"), sample=sample, db_name = lca_dbs.keys())

rule all:
    input: gather_targets # + lca_classify_targets

# for protein signatures, multipy by 3 if necessary before calculating signature (sourmash v3.x)
ksize_multiplier = {"dna": 1, "protein": 3, "dayhoff": 3, "hp":3}

# gather each sig
rule gather_sig:
    input:
        query= os.path.join(sigs_dir, "{genome}.faa.sig"), # TARA_IOS_MAG_00042.fa.faa.sig 
        db= lambda w: sbt_dbs[w.sbt_db]["sbt"]
    output:
        csv = os.path.join(gather_dir, "{genome}_x_{sbt_db}.gather.csv"),
        matches = os.path.join(gather_dir, "{genome}_x_{sbt_db}.gather.matches"),
        unassigned = os.path.join(gather_dir, "{genome}_x_{sbt_db}.gather.unassigned")
    params:
        alpha_cmd = lambda w: "--" + alphabet,
        ksize = int(ksize) * ksize_multiplier[alphabet],
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=20,
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

# need to modify bc no true lineages for tara genomes...
#rule aggregate_gather_to_tax:
#    # make spreadsheet: each proteome:: top lineage hit
#    input:
#        gather_tophits= lambda w: expand(os.path.join(gather_dir, "{genome}_x_{{sbt_db}}.gather_tophits.csv"), genome=genome_list)
#        #true_lineages=lambda w: sampleInfo[w.sample]["gather_csv"]
#    output:
#         summary_csv=os.path.join(gather_dir, "gather_tophits", "{sample}", "{sample}.{alphabet}_scaled{scaled}_k{k}.gather_tophits.summary.csv"),
#    params:
#        gather_dir= lambda w: os.path.join(gather_dir, w.sample, w.alphabet, f"k{w.k}")
#    log: os.path.join(logs_dir, "gather_tophits", "{sample}.{alphabet}_scaled{scaled}_k{k}.gather_tophits.log")
#    benchmark: os.path.join(logs_dir, "gather_tophits", "{sample}.{alphabet}_scaled{scaled}_k{k}.gather_tophits.benchmark")
#    conda: "envs/forage-env.yml"
#    shell:
#        """
#        python scripts/aggregate-gather-to-tax-tophits.py --input-is-directory --true-lineages-csv {input.true_lineages} --output-csv {output} {params.gather_dir} 2> {log}
#        """


