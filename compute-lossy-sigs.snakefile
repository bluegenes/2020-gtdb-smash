"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Build lossy signatures with abundance information, aggregates into LCA, SBT
"""

import os
import glob
import pandas as pd
from collections import defaultdict

out_dir = config["out_dir"]
compute_dir = config["sigs_dir"]
index_dir = config["index_dir"]
logs_dir = os.path.join(out_dir, "logs")
envs_dir = "envs"

sampleInfo=config["samples"]
accession2filenames = defaultdict(dict)
accession2signame = defaultdict(dict)

sig_targets,index_targets=[],[]

# using the same configfile that can be used for other stuff.. just grab accessions, signames, filenames; build sigs
for sample, info in sampleInfo.items():
    info_csv = pd.read_csv(info["info_csv"])
    # accession:: fasta file map
    fasta_dir= info["input_path"]
    info_csv["filename"] = info_csv["filename"].apply(lambda x: os.path.join(fasta_dir, x))
    info_csv["signame"] = info_csv["accession"] + " " + info_csv["species"]
    input_type = info["input_type"] #protein, dna, rna
    sample_acc2file = pd.Series(info_csv.filename.values,index=info_csv.accession).to_dict()
    acc2signame = pd.Series(info_csv.signame.values,index=info_csv.accession).to_dict()
    accession2filenames[input_type].update(sample_acc2file)
    accession2signame.update(acc2signame)
    sampleInfo[sample]["accessions"] = acc2signame.keys()
    for alpha, alphainfo in info["alphabet"].items():
        sig_targets+=expand(os.path.join(compute_dir, input_type, "{alphabet}", "k{k}", "{accession}_{alphabet}_scaled{scaled}_k{k}.sig"), accession=sample_acc2file.keys(), alphabet=alpha,scaled=alphainfo["scaled"], k=alphainfo["ksizes"])
#        index_targets+=expand(os.path.join(index_dir,"lca", "{sample}.{alphabet}_scaled{scaled}_k{k}.index.lca.json.gz"), sample=sample, alphabet=alpha, scaled=alphainfo["scaled"], k=alphainfo["ksizes"])
#        index_targets+=expand(os.path.join(index_dir,"sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.index.sbt.zip"), sample=sample, alphabet=alpha, scaled=alphainfo["scaled"], k=alphainfo["ksizes"])


rule all:
    input: sig_targets + index_targets


#def find_alphabet(alphabet):
#    alpha=alphabet
#    if "_" in alphabet:
#        alpha = alphabet.rsplit("_")[1]
#    return alpha

# compute sigs to general compute directory, regardless of the sbt sample name (so can reuse sigs for other sbts, if desired)
rule sourmash_compute_dna:
    input: lambda w: accession2filenames["dna"][w.accession]
    output: os.path.join(compute_dir, "dna", "{alphabet}", "k{k}", "{accession}_{alphabet}_scaled{scaled}_k{k}.sig")
    params:
        #k= lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
        alpha = lambda w: w.alphabet.rsplit("_", 1)[1] if "_" in w.alphabet else w.alphabet,
        skipmer_cmd = lambda w: "--skipmer" if "skipmer" in w.alphabet else "",
        signame = lambda w: accession2signame[w.accession],
        #abund_cmd = "--track-abundance", # default
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash", "{accession}_{alphabet}_scaled{scaled}_k{k}.dna.compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{accession}_{alphabet}_scaled{scaled}_k{k}.dna.compute.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        #generate-lossy-signature.py does the ksize multiplication internally
        """
        python scripts/generate-lossy-signature.py {input} --ksize {wildcards.k} \
          --input-alphabet nucleotide --output-alphabet {wildcards.alphabet} \
          --scaled {wildcards.scaled} --signame {params.signame:q} \
          --output_signature {output} {params.skipmer_cmd} 2> {log}
        """
        
        #sourmash compute -k {params.k} --scaled={wildcards.scaled}  \
        #{input} -o {output} {params.alpha_cmd} {params.abund_cmd} --merge={params.signame:q} 2> {log}

rule sourmash_compute_protein:
    input: lambda w: accession2filenames["protein"][w.accession]
    output: os.path.join(compute_dir, "protein", "{alphabet}", "k{k}", "{accession}_{alphabet}_scaled{scaled}_k{k}.sig")
    params:
        #k= lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
        alpha = lambda w: w.alphabet.rsplit("_", 1)[1] if "_" in w.alphabet else w.alphabet,
        #alpha_cmd = lambda w: "--" + moltype_map[w.alphabet],
        signame = lambda w: accession2signame[w.accession],
        #abund_cmd = "--track-abundance",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash", "{accession}_{alphabet}_scaled{scaled}_k{k}.protein.compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{accession}_{alphabet}_scaled{scaled}_k{k}.protein.compute.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        #generate-lossy-signature.py does the ksize multiplication internally
        """
        python scripts/generate-lossy-signature.py {input} --ksize {wildcards.k} \
          --input-alphabet protein --output-alphabet {wildcards.alphabet} \
          --scaled {wildcards.scaled} --signame {params.signame:q} \
          --output_signature {output} 2> {log}
        """
        #sourmash compute -k {params.k} --scaled={wildcards.scaled} --input-is-protein \
        #{input} -o {output} {params.alpha_cmd} {params.abund_cmd} --merge={params.signame:q} 2> {log}

rule sourmash_compute_rna:
    input: lambda w: accession2filenames["rna"][w.accession]
    output: os.path.join(compute_dir, "rna", "{alphabet}", "k{k}", "{accession}_{alphabet}_scaled{scaled}_k{k}.sig")
    params:
        #k= lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
        #alpha_cmd = lambda w: "--" + moltype_map[w.alphabet],
        alpha = lambda w: w.alphabet.rsplit("_", 1)[1] if "_" in w.alphabet else w.alphabet,
        signame = lambda w: accession2signame[w.accession],
        skipmer_cmd = lambda w: "--skipmer" if w.alphabet == "skipmer" else "",
        #abund_cmd = "--track-abundance",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash", "{accession}_{alphabet}_scaled{scaled}_k{k}.rna.compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{accession}_{alphabet}_scaled{scaled}_k{k}.rna.compute.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        #generate-lossy-signature.py does the ksize multiplication internally
        """
        python scripts/generate-lossy-signature.py {input} --ksize {wildcards.k} \
          --input-alphabet nucleotide --output-alphabet {wildcards.alphabet} \
          --scaled {wildcards.scaled} --signame {params.signame:q} \
          --output_signature {output} {params.skipmer_cmd} 2> {log}
        """
        #sourmash compute -k {params.k} --scaled={wildcards.scaled}  \
        #{input} -o {output} {params.alpha_cmd} {params.abund_cmd} --merge={params.signame:q} 2> {log}

def aggregate_sigs(w):
    siglist=[]
    input_type = sampleInfo[w.sample]["input_type"] # protein, dna, rna
    sigfile = os.path.join(compute_dir, input_type, w.alphabet, f"k{w.k}", f"{{acc}}_{w.alphabet}_scaled{w.scaled}_k{w.k}.sig")
    siglist=expand(sigfile, acc=sampleInfo[w.sample]["accessions"])
    return siglist


#### should I just say the alpha is dna or protein for now??
#rule index_lca:
#    input:
#        sigs=aggregate_sigs,
#        taxonomy= lambda w: sampleInfo[w.sample]["info_csv"]
#    output:
#        os.path.join(index_dir,"lca", "{sample}.{alphabet}_scaled{scaled}_k{k}.index.lca.json.gz"),
#    threads: 1
#    params:
#        alpha= lambda w: (w.alphabet.rsplit("translate_")[1] if w.alphabet.startswith("translate") else w.alphabet), # remove translate
#        #alpha_cmd= lambda w: " --" + (w.alphabet.rsplit("translate_")[1] if w.alphabet.startswith("translate") else w.alphabet), # remove translate
#        alpha_cmd= lambda w: "--" + moltype_map[w.alphabet],
#        translate = lambda w: " --translate " if w.alphabet.startswith("translate") else "",
#        input_type = lambda w: sampleInfo[w.sample]["input_type"],
#        ksize = lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
#        report= lambda w: os.path.join(index_dir,"lca", f"{w.sample}.{w.alphabet}_scaled{w.scaled}_k{w.k}.index.lca.report"),
#        #output_prefix = lambda w: os.path.join(index_dir,"{w.sample}.{w.alphabet}_scaled{w.scaled}_k{w.k}.index")
#        sigdir= lambda w: os.path.join(compute_dir, sampleInfo[w.sample]["input_type"], w.alphabet, f"k{w.k}")
#    resources:
#        mem_mb= lambda wildcards, attempt: attempt *100000,
#        runtime=600,
#    log: os.path.join(logs_dir, "index-lca", "{sample}.{alphabet}_scaled{scaled}_k{k}.index-lca.log")
#    benchmark: os.path.join(logs_dir, "index-lca", "{sample}.{alphabet}_scaled{scaled}_k{k}.index-lca.benchmark")
#    conda: "envs/sourmash-dev.yml"
#    shell:
#        """
#        sourmash lca index \
#          --ksize {params.ksize} \
#          --scaled {wildcards.scaled} \
#          --split-identifiers \
#          --require-taxonomy \
#          --traverse-directory \
#          --report {params.report} \
#          {params.alpha_cmd} {input.taxonomy} {output} \
#          {params.sigdir} 2> {log}
#        """
#        #touch  {output.report}
#        #{compute_dir}/{params.input_type}/{wildcards.alphabet}/k{wildcards.k}/*_{params.alpha}_scaled{wildcards.scaled}_k{wildcards.k}.sig  2> {log}


#rule lca_to_sbt:
#    input: os.path.join(index_dir, "lca", "{sample}.{alphabet}_scaled{scaled}_k{k}.index.lca.json.gz")
#    output: os.path.join(index_dir, "sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.index.sbt.zip")
#    log: os.path.join(logs_dir, "lca-to-sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.lca-to-sbt.log")
#    benchmark: os.path.join(logs_dir, "lca-to-sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.lca-to-sbt.benchmark")
#    resources:
#        mem_mb= lambda wildcards, attempt: attempt *100000,
#        runtime=600,
#    conda: "envs/sourmash-dev.yml"
#    shell:
#        """
#        python scripts/convert-lca-to-sbt.py {input} {output} 2> {log}
#        """
