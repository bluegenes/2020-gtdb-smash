"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s gtdb-smash.snakefile --use-conda
"""

import os
import glob
import pandas as pd

output_extensions = ["sig"]
output_targets, nucleotide_targets, translate_targets, protein_targets = [],[],[],[]

out_dir = config["out_dir"]
#data_dir= os.path.join(out_dir, "data")
logs_dir = os.path.join(out_dir, "logs")
envs_dir = "envs"
compute_dir = os.path.join(out_dir, "compute")
compare_dir = os.path.join(out_dir, "compare")

index_dir = config.get(out_dir, "representative-index")

index_targets, search_targets=[],[]

sampleInfo=config["representative_indices"]
accession2filenames = {}

for sample, info in sampleInfo.items():    
    if info.get("build_index"):
        info_csv = pd.read_csv(info["build_index"]["info_csv"])
        # accession:: fasta file map
        fasta_dir= info["build_index"]["input_path"]
        info_csv["filename"] = info_csv["filename"].apply(lambda x: os.path.join(fasta_dir, x))
        input_type = info["input_type"] #protein, dna, rna
        accession2filenames[input_type] = pd.Series(info_csv.filename.values,index=info_csv.accession).to_dict()
        # add accessions to sampleInfo dict
        sampleInfo[sample]["accessions"] = info_csv.accession.to_list()
        # build sbt targets
        index_dir= info["build_index"]["index_outdir"]
        index_types= info["build_index"]["index_types"] # lca or sbt
        for alpha, alphainfo in info["alphabet"].items():
            if "lca" in index_types:
                index_targets+=expand(os.path.join(index_dir,"{sample}.{alphabet}_scaled{scaled}_k{k}.index.lca.json"), sample=sample, alphabet=alpha, scaled=alphainfo["scaled"], k=alphainfo["ksizes"])
            if "sbt" in index_types:
                index_targets+=expand(os.path.join(index_dir,"{sample}.{alphabet}_scaled{scaled}_k{k}.index.sbt.zip"), sample=sample, alphabet=alpha, scaled=alphainfo["scaled"], k=alphainfo["ksizes"])

    # build search targets if any
    #if info.get("query_csv"):
        #query_csv = pd.read_csv(info["query_csv"]) 
    #for alpha, alphainfo in info["alphabet"].items():
    #    search_targets+=expand(os.path.join(dist_dir, "{sample}.{alphabet}_scaled{scaled}_k{k}.jaccard_from_species.csv"), sample=sample, alphabet=alpha, scaled=alphainfo["scaled"], k=alphainfo["ksizes"])

rule all:
    input: index_targets + search_targets

# for protein signatures, multipy by 3 if necessary before calculating signature (sourmash v3.x)
ksize_multiplier = {"nucleotide": 1, "protein": 3, "dayhoff": 3, "hp":3, "translate_protein": 3, "translate_dayhoff": 3, "translate_hp": 3}
# "rna" is not sourmash cli friendly
moltype_map = {"nucleotide": "dna", "protein":"protein", "dayhoff": "dayhoff", "hp":"hp", "translate_protein": "protein", "translate_dayhoff": "dayhoff", "translate_hp": "hp"}


# compute sigs to general compute directory, regardles of the sbt sample name (so can reuse sigs for other sbts, if desired)
rule sourmash_compute_dna:
    input: lambda w: accession2filenames["dna"][w.accession]
    output: os.path.join(compute_dir, "dna", "{alphabet}", "k{k}", "{accession}_{alphabet}_scaled{scaled}_k{k}.sig")
    params:
        k= lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
        scaled= lambda w: w.scaled,
        compute_moltypes= lambda w: moltype_map[w.alphabet],
        input_is_protein=False,
        track_abundance=True,
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash", "{accession}_{alphabet}_scaled{scaled}_k{k}.dna.compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{accession}_{alphabet}_scaled{scaled}_k{k}.dna.compute.benchmark")
    conda: "envs/sourmash3.3.yml"
    script: "scripts/sourmash-compute.wrapper.py"

rule sourmash_compute_protein:
    input: lambda w: accession2filenames["protein"][w.accession]
    output: os.path.join(compute_dir, "protein", "{alphabet}", "k{k}", "{accession}_{alphabet}_scaled{scaled}_k{k}.sig")
    params:
        k= lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
        scaled= lambda w: w.scaled,
        compute_moltypes= lambda w: moltype_map[w.alphabet],
        input_is_protein=True,
        track_abundance=True,
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash", "{accession}_{alphabet}_scaled{scaled}_k{k}.protein.compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{accession}_{alphabet}_scaled{scaled}_k{k}.protein.compute.benchmark")
    conda: "envs/sourmash3.3.yml"
    script: "scripts/sourmash-compute.wrapper.py"

rule sourmash_compute_rna:
    input: lambda w: accession2filenames["rna"][w.accession]
    output: os.path.join(compute_dir, "rna", "{alphabet}", "k{k}", "{accession}_{alphabet}_scaled{scaled}_k{k}.sig")
    params:
        k= lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
        scaled= lambda w: w.scaled,
        compute_moltypes= lambda w: moltype_map[w.alphabet],
        input_is_protein=False,
        track_abundance=True,
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash", "{accession}_{alphabet}_scaled{scaled}_k{k}.rna.compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{accession}_{alphabet}_scaled{scaled}_k{k}.rna.compute.benchmark")
    conda: "envs/sourmash3.3.yml"
    script: "scripts/sourmash-compute.wrapper.py"

def aggregate_sigs(w):
    siglist=[]
    input_type = sampleInfo[w.sample]["input_type"] # protein, dna, rna
    sigfile = os.path.join(compute_dir, input_type, w.alphabet, f"k{w.k}", f"{{acc}}_{w.alphabet}_scaled{w.scaled}_k{w.k}.sig")
    siglist=expand(sigfile, acc=sampleInfo[w.sample]["accessions"])
    return siglist


rule index_sbt:
    input: aggregate_sigs,
    output: 
        sbt=os.path.join(index_dir,"{sample}.{alphabet}_scaled{scaled}_k{k}.index.sbt.zip"),
    threads: 1
    params:
        alpha= lambda w: w.alphabet.rsplit("translate_")[1] if w.alphabet.startswith("translate") else w.alphabet, # remove translate
        alpha_cmd= lambda w: " --" + (w.alphabet.rsplit("translate_")[1] if w.alphabet.startswith("translate") else w.alphabet), # remove translate
        translate = lambda w: " --translate " if w.alphabet.startswith("translate") else "",
        input_type = lambda w: sampleInfo[w.sample]["input_type"],
        ksize = lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
        runtime=6000,
    log: os.path.join(logs_dir, "index-sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.index-sbt.log")
    benchmark: os.path.join(logs_dir, "index-sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.index-sbt.benchmark")
    conda: "envs/forage-env.yml"
    shell:
        """
        sourmash index --ksize {params.ksize} --scaled {wildcards.scaled} {params.alpha_cmd}  \
        {output.sbt} {compute_dir}/{params.input_type}/*_{params.alpha}_scaled{wildcards.scaled}_k{wildcards.k}.sig  2> {log}
        """

rule index_lca:
    input: 
        sigs=aggregate_sigs,
        taxonomy= lambda w: sampleInfo[w.sample]["build_index"]["info_csv"]
    output:
        lca=os.path.join(index_dir,"{sample}.{alphabet}_scaled{scaled}_k{k}.index.lca.json"),
        report= os.path.join(index_dir,"{sample}.{alphabet}_scaled{scaled}_k{k}.index.lca.report"),
    threads: 1
    params:
        alpha= lambda w: (w.alphabet.rsplit("translate_")[1] if w.alphabet.startswith("translate") else w.alphabet), # remove translate
        alpha_cmd= lambda w: " --" + (w.alphabet.rsplit("translate_")[1] if w.alphabet.startswith("translate") else w.alphabet), # remove translate
        translate = lambda w: " --translate " if w.alphabet.startswith("translate") else "",
        input_type = lambda w: sampleInfo[w.sample]["input_type"],
        ksize = lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
        #output_prefix = lambda w: os.path.join(index_dir,"{w.sample}.{w.alphabet}_scaled{w.scaled}_k{w.k}.index")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *300000,
        runtime=600000,
    log: os.path.join(logs_dir, "index-lca", "{sample}.{alphabet}_scaled{scaled}_k{k}.index-lca.log")
    benchmark: os.path.join(logs_dir, "index-lca", "{sample}.{alphabet}_scaled{scaled}_k{k}.index-lca.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash lca index --ksize {params.ksize} --scaled {wildcards.scaled} --report {output.report} \
        {params.alpha_cmd} {input.taxonomy} {output.lca} \
        {compute_dir}/{params.input_type}/{wildcards.alphabet}/k{wildcards.k}/*_{params.alpha}_scaled{wildcards.scaled}_k{wildcards.k}.sig  2> {log}
        """
