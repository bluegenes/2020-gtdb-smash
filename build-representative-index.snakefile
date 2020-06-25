"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s gtdb-smash.snakefile --use-conda
"""

import os
import glob
import pandas as pd
from collections import defaultdict

output_extensions = ["sig"]
output_targets, nucleotide_targets, translate_targets, protein_targets = [],[],[],[]

out_dir = config["out_dir"]
compute_dir = config["sigs_dir"]
#data_dir= os.path.join(out_dir, "data")
logs_dir = os.path.join(out_dir, "logs")
envs_dir = "envs"
compare_dir = os.path.join(out_dir, "compare")

index_dir = config.get(out_dir, "representative-index")
gather_dir = os.path.join(out_dir, "gather")

index_targets, gather_targets=[],[]

sampleInfo=config["representative_indices"]
accession2filenames = defaultdict(dict)
accession2signame = defaultdict(dict)

#def make_sigfilename(filename):
# from https://github.com/dib-lab/sourmash_databases/blob/master/Snakefile
#    genomefile = os.path.basename(filename)
#    sigfile = os.path.join(sigs_output_location, genomefile) + '.sig'
#    return sigfile


for sample, info in sampleInfo.items():    
    info_csv = pd.read_csv(info["info_csv"])
    # accession:: fasta file map
    fasta_dir= info["input_path"]
    info_csv["filename"] = info_csv["filename"].apply(lambda x: os.path.join(fasta_dir, x))
    info_csv["signame"] = info_csv["accession"] + " " + info_csv["species"]
    input_type = info["input_type"] #protein, dna, rna
    # if want multiple samples, need to add here instead of overwriting
    sample_acc2file = pd.Series(info_csv.filename.values,index=info_csv.accession).to_dict()
    acc2signame = pd.Series(info_csv.signame.values,index=info_csv.accession).to_dict()
    accession2filenames[input_type].update(sample_acc2file)
    accession2signame.update(acc2signame)
    # add accessions to sampleInfo dict
    sampleInfo[sample]["accessions"] = info_csv.accession.to_list()
    # get index info
    index_dir= info["index_outdir"]
    index_types= info["index_types"] # lca or sbt
    # if we have accessions to gather, handle it.
    gather_accessions = ""
    if info.get("gather_csv"):
        gather_csv = pd.read_csv(info["gather_csv"])
        # this assumes they're in same input dir & are same input type as the samples above. ok for now. shrug.
        gather_csv["filename"] = gather_csv["filename"].apply(lambda x: os.path.join(fasta_dir, x))
        gather_csv["signame"] = gather_csv["accession"] + " " + gather_csv["species"]
        gather_acc2file = pd.Series(gather_csv.filename.values,index=gather_csv.accession).to_dict()
        gather2signame = pd.Series(gather_csv.signame.values,index=gather_csv.accession).to_dict()
        # add accession:: filenames to main acc2filename dict
        accession2filenames[input_type].update(gather_acc2file)
        accession2signame.update(gather2signame)
        # create "gather_accessions" list
        gather_accessions = gather_csv.accession.to_list()
    sampleInfo[sample]["gather_accessions"] = gather_accessions
    # build targets
    for alpha, alphainfo in info["alphabet"].items():
        if "lca" in index_types:
            index_targets+=expand(os.path.join(index_dir,"lca", "{sample}.{alphabet}_scaled{scaled}_k{k}.index.lca.json.gz"), sample=sample, alphabet=alpha, scaled=alphainfo["scaled"], k=alphainfo["ksizes"])
            if gather_accessions:
                gather_targets+=expand(os.path.join(gather_dir, "{sample}", "{alphabet}", "k{k}", "{acc}_x_{sample}.{alphabet}_scaled{scaled}_k{k}.lca.gather.csv"), acc=gather_accessions, sample=sample, alphabet=alpha, scaled=alphainfo["scaled"], k=alphainfo["ksizes"])
                gather_targets+=expand(os.path.join(gather_dir, "{sample}", "{alphabet}", "k{k}", "{acc}_x_{sample}.{alphabet}_scaled{scaled}_k{k}.gather_tophits.csv"), acc=gather_accessions, sample=sample, alphabet=alpha, scaled=alphainfo["scaled"], k=alphainfo["ksizes"])
                gather_targets+=expand(os.path.join(gather_dir, "gather_tophits","{sample}","{sample}.{alphabet}_scaled{scaled}_k{k}.gather_tophits.summary.csv"), sample=sample, alphabet=alpha,scaled=alphainfo["scaled"], k=alphainfo["ksizes"])
        if "sbt" in index_types:
            index_targets+=expand(os.path.join(index_dir,"sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.index.sbt.zip"), sample=sample, alphabet=alpha, scaled=alphainfo["scaled"], k=alphainfo["ksizes"])
            if gather_accessions:
                gather_targets+=expand(os.path.join(gather_dir, "{sample}", "{alphabet}", "k{k}", "{acc}_x_{sample}.{alphabet}_scaled{scaled}_k{k}.sbt.gather.csv"), acc=gather_accessions, sample=sample, alphabet=alpha, scaled=alphainfo["scaled"], k=alphainfo["ksizes"])

rule all:
    input: index_targets + gather_targets

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
        alpha_cmd = lambda w: "--" + moltype_map[w.alphabet],
        signame = lambda w: accession2signame[w.accession],
        abund_cmd = "--track-abundance",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash", "{accession}_{alphabet}_scaled{scaled}_k{k}.dna.compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{accession}_{alphabet}_scaled{scaled}_k{k}.dna.compute.benchmark")
    conda: "envs/sourmash3.3.yml"
    shell:
        """
        sourmash compute -k {params.k} --scaled={wildcards.scaled}  \
        {input} -o {output} {params.alpha_cmd} {params.abund_cmd} --merge={params.signame:q} 2> {log}
        """

rule sourmash_compute_protein:
    input: lambda w: accession2filenames["protein"][w.accession]
    output: os.path.join(compute_dir, "protein", "{alphabet}", "k{k}", "{accession}_{alphabet}_scaled{scaled}_k{k}.sig")
    params:
        k= lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
        alpha_cmd = lambda w: "--" + moltype_map[w.alphabet],
        signame = lambda w: accession2signame[w.accession],
        abund_cmd = "--track-abundance",
        #signame_cmd = lambda w: f"--merge={accession2signame[w.accession]}",
        #scaled= lambda w: w.scaled,
        #compute_moltypes= lambda w: moltype_map[w.alphabet],
        #input_is_protein=True,
        #track_abundance=True,
        #signame = lambda wildcards: lookup_name(wildcards.filename),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash", "{accession}_{alphabet}_scaled{scaled}_k{k}.protein.compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{accession}_{alphabet}_scaled{scaled}_k{k}.protein.compute.benchmark")
    conda: "envs/sourmash3.3.yml"
    shell:
        """
        sourmash compute -k {params.k} --scaled={wildcards.scaled} --input-is-protein \
        {input} -o {output} {params.alpha_cmd} {params.abund_cmd} --merge={params.signame:q} 2> {log}
        """
        #{input} -o {output} {params.alpha_cmd} {params.abund_cmd} {params.signame_cmd} 2> {log}
    #script: "scripts/sourmash-compute.wrapper.py"

rule sourmash_compute_rna:
    input: lambda w: accession2filenames["rna"][w.accession]
    output: os.path.join(compute_dir, "rna", "{alphabet}", "k{k}", "{accession}_{alphabet}_scaled{scaled}_k{k}.sig")
    params:
        k= lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
        alpha_cmd = lambda w: "--" + moltype_map[w.alphabet],
        signame = lambda w: accession2signame[w.accession],
        abund_cmd = "--track-abundance",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash", "{accession}_{alphabet}_scaled{scaled}_k{k}.rna.compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{accession}_{alphabet}_scaled{scaled}_k{k}.rna.compute.benchmark")
    conda: "envs/sourmash3.3.yml"
    shell:
        """
        sourmash compute -k {params.k} --scaled={wildcards.scaled}  \
        {input} -o {output} {params.alpha_cmd} {params.abund_cmd} --merge={params.signame:q} 2> {log}
        """
    #script: "scripts/sourmash-compute.wrapper.py"
        #sourmash compute {moltype_cmd} {abund_cmd} --scaled {scaled} -k {k} {snakemake.input} -o {snakemake.output} -p {snakemake.threads} {extra} {log}

def aggregate_sigs(w):
    siglist=[]
    input_type = sampleInfo[w.sample]["input_type"] # protein, dna, rna
    sigfile = os.path.join(compute_dir, input_type, w.alphabet, f"k{w.k}", f"{{acc}}_{w.alphabet}_scaled{w.scaled}_k{w.k}.sig")
    siglist=expand(sigfile, acc=sampleInfo[w.sample]["accessions"])
    return siglist


rule index_sbt:
    input: sigs=aggregate_sigs,
    output: 
        sbt=os.path.join(index_dir, "sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.index.sbt.zip"),
    threads: 1
    params:
        alpha= lambda w: w.alphabet.rsplit("translate_")[1] if w.alphabet.startswith("translate") else w.alphabet, # remove translate
        #alpha_cmd= lambda w: " --" + (w.alphabet.rsplit("translate_")[1] if w.alphabet.startswith("translate") else w.alphabet), # remove translate
        alpha_cmd = lambda w: "--" + moltype_map[w.alphabet],
        translate = lambda w: " --translate " if w.alphabet.startswith("translate") else "",
        input_type = lambda w: sampleInfo[w.sample]["input_type"],
        ksize = lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
        sigdir= lambda w: os.path.join(compute_dir, sampleInfo[w.sample]["input_type"], w.alphabet, f"k{w.k}")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=6000,
    log: os.path.join(logs_dir, "index-sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.index-sbt.log")
    benchmark: os.path.join(logs_dir, "index-sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.index-sbt.benchmark")
    conda: "envs/forage-env.yml"
    shell:
        """
        sourmash index --ksize {params.ksize} --scaled {wildcards.scaled} {params.alpha_cmd} {output.sbt} \
        {input.sigs} 2> {log}
        """
        #{output.sbt} {compute_dir}/{params.input_type}/{wildcards.alphabet}/k{wildcards.k}/*_{params.alpha}_scaled{wildcards.scaled}_k{wildcards.k}.sig  2> {log}

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
        #touch  {output.report}
        #{compute_dir}/{params.input_type}/{wildcards.alphabet}/k{wildcards.k}/*_{params.alpha}_scaled{wildcards.scaled}_k{wildcards.k}.sig  2> {log}

## gather rules

rule gather_lca:
    input:
        query = lambda w: os.path.join(compute_dir, sampleInfo[w.sample]["input_type"], w.alphabet, f"k{w.k}", f"{w.accession}_{w.alphabet}_scaled{w.scaled}_k{w.k}.sig"),
        db=os.path.join(index_dir, "lca", "{sample}.{alphabet}_scaled{scaled}_k{k}.index.lca.json.gz")
        # maybe later expand this to multiple db's at once? would need to change targets, too.
        #dbs = lambda w: sampleInfo[w.sample]["databases"]["lca"].values()
    output:
        csv = os.path.join(gather_dir, "{sample}", "{alphabet}", "k{k}", "{accession}_x_{sample}.{alphabet}_scaled{scaled}_k{k}.lca.gather.csv"),
        matches = os.path.join(gather_dir, "{sample}", "{alphabet}", "k{k}", "{accession}_x_{sample}.{alphabet}_scaled{scaled}_k{k}.lca.gather.matches"),
        unassigned = os.path.join(gather_dir, "{sample}", "{alphabet}", "k{k}", "{accession}_x_{sample}.{alphabet}_scaled{scaled}_k{k}.lca.gather.unassigned"),
    params:
        alpha= lambda w: (w.alphabet.rsplit("translate_")[1] if w.alphabet.startswith("translate") else w.alphabet), # remove translate
        #alpha_cmd= lambda w: " --" + (w.alphabet.rsplit("translate_")[1] if w.alphabet.startswith("translate") else w.alphabet), # remove translate
        alpha_cmd = lambda w: "--" + moltype_map[w.alphabet],
        translate = lambda w: " --translate " if w.alphabet.startswith("translate") else "",
        input_type = lambda w: sampleInfo[w.sample]["input_type"],
        ksize = lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
        #output_prefix = lambda w: os.path.join(index_dir,"{w.sample}.{w.alphabet}_scaled{w.scaled}_k{w.k}.index")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=600,
    log: os.path.join(logs_dir, "gather", "{accession}_x_{sample}.{alphabet}_scaled{scaled}_k{k}.lca-gather.log")
    benchmark: os.path.join(logs_dir, "gather", "{accession}_x_{sample}.{alphabet}_scaled{scaled}_k{k}.lca-gather.benchmark")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=60,
    conda: "envs/sourmash-dev.yml"
    shell:
        # do we want abundance?? --ignore-abundance to turn off
        """
        sourmash gather {input.query} {input.db} -o {output.csv} {params.alpha_cmd} \
        --save-matches {output.matches} --threshold-bp=0  \
        --output-unassigned {output.unassigned} \
        --scaled {wildcards.scaled} -k {params.ksize} 2> {log}
        touch {output}
        """
# touch empty output file to enable rna ones to fail

rule gather_sbt:
    input:
        query = lambda w: os.path.join(compute_dir, sampleInfo[w.sample]["input_type"], w.alphabet, f"k{w.k}", f"{w.accession}_{w.alphabet}_scaled{w.scaled}_k{w.k}.sig"),
        db=os.path.join(index_dir, "sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.index.sbt.zip")
        # maybe later expand this to multiple db's at once? would need to change targets, too.
    output:
        csv = os.path.join(gather_dir, "{sample}", "{alphabet}", "k{k}", "{accession}_x_{sample}.{alphabet}_scaled{scaled}_k{k}.sbt.gather.csv"),
        matches = os.path.join(gather_dir, "{sample}", "{alphabet}", "k{k}", "{accession}_x_{sample}.{alphabet}_scaled{scaled}_k{k}.sbt.gather.matches"),
        unassigned = os.path.join(gather_dir, "{sample}", "{alphabet}", "k{k}", "{accession}_x_{sample}.{alphabet}_scaled{scaled}_k{k}.sbt.gather.unassigned"),
    params:
        alpha= lambda w: (w.alphabet.rsplit("translate_")[1] if w.alphabet.startswith("translate") else w.alphabet), # remove translate
        alpha_cmd = lambda w: "--" + moltype_map[w.alphabet],
        #alpha_cmd= lambda w: " --" + (w.alphabet.rsplit("translate_")[1] if w.alphabet.startswith("translate") else w.alphabet), # remove translate
        translate = lambda w: " --translate " if w.alphabet.startswith("translate") else "",
        input_type = lambda w: sampleInfo[w.sample]["input_type"],
        ksize = lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
        #output_prefix = lambda w: os.path.join(index_dir,"{w.sample}.{w.alphabet}_scaled{w.scaled}_k{w.k}.index")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=600000,
    log: os.path.join(logs_dir, "gather", "{accession}_x_{sample}.{alphabet}_scaled{scaled}_k{k}.sbt-gather.log")
    benchmark: os.path.join(logs_dir, "gather", "{accession}_x_{sample}.{alphabet}_scaled{scaled}_k{k}.sbt-gather.benchmark")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=60,
    conda: "envs/sourmash-dev.yml"
    shell:
        # do we want abundance?? --ignore-abundance to turn off
        """
        sourmash gather {input.query} {input.db} -o {output.csv} {params.alpha_cmd} \
        --save-matches {output.matches} --threshold-bp=0  \
        --output-unassigned {output.unassigned} \
        --scaled {wildcards.scaled} -k {params.ksize} 2> {log}
        """

rule gather_to_tax:
    input:
        gather_csv = rules.gather_lca.output.csv,
        lineages_csv = lambda w: sampleInfo[w.sample]["info_csv"]
    output:
        gather_tax = os.path.join(gather_dir, "{sample}", "{alphabet}", "k{k}", "{accession}_x_{sample}.{alphabet}_scaled{scaled}_k{k}.gather_summary.csv"),
        top_matches = os.path.join(gather_dir, "{sample}", "{alphabet}", "k{k}", "{accession}_x_{sample}.{alphabet}_scaled{scaled}_k{k}.gather_tophits.csv"),
    log: os.path.join(logs_dir, "gather_to_tax", "{accession}_x_{sample}.{alphabet}_scaled{scaled}_k{k}.gather_to_tax.log")
    benchmark: os.path.join(logs_dir, "gather_to_tax", "{accession}_x_{sample}.{alphabet}_scaled{scaled}_k{k}.gather_to_tax.benchmark")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=600,
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        python scripts/gather-to-tax.py {input.gather_csv} {input.lineages_csv} --tophits-csv {output.top_matches} > {output.gather_tax} 2> {log}
        """
rule aggregate_gather_to_tax:
    # make spreadsheet: each proteome:: top lineage hit
    input:
        gather_tophits= lambda w: expand(os.path.join(gather_dir, "{{sample}}/{{alphabet}}/k{{k}}/{acc}_x_{{sample}}.{{alphabet}}_scaled{{scaled}}_k{{k}}.gather_tophits.csv"), acc=sampleInfo[w.        sample]["gather_accessions"]),
        true_lineages=lambda w: sampleInfo[w.sample]["gather_csv"]
    output:
         summary_csv=os.path.join(gather_dir, "gather_tophits", "{sample}", "{sample}.{alphabet}_scaled{scaled}_k{k}.gather_tophits.summary.csv"),
    params:
        gather_dir= lambda w: os.path.join(gather_dir, w.sample, w.alphabet, f"k{w.k}")
    log: os.path.join(logs_dir, "gather_tophits", "{sample}.{alphabet}_scaled{scaled}_k{k}.gather_tophits.log")
    benchmark: os.path.join(logs_dir, "gather_tophits", "{sample}.{alphabet}_scaled{scaled}_k{k}.gather_tophits.benchmark")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=60,
    conda: "envs/forage-env.yml"
    shell:
        """
        python scripts/aggregate-gather-to-tax-tophits.py --input-is-directory --true-lineages-csv {input.true_lineages} --output-csv {output} {params.gather_dir} 2> {log}
        """
