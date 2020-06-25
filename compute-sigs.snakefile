"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s compute-sigs.snakefile --use-conda
"""

import os
import glob
import pandas as pd

compute_dir = config["out_dir"]
logs_dir = os.path.join(compute_dir, "logs")
envs_dir = "envs"

genome_info = config["genomes_info"]
fasta_dir = config["genomes_dir"]
input_type = config["input_type"]
if input_type == "protein":
    filenames = [os.path.join(fasta_dir, x.strip() + ".faa") for x in open(genome_info) ]
    signames = [os.path.basename(x).split(".fa.faa")[0] for x in filenames]
elif input_type == "nucleotide":
    filenames = [os.path.join(fasta_dir, x.strip()) for x in open(genome_info) ]
    signames = [os.path.basename(x).split(".fa")[0] for x in filenames]
else:
    print("input type must be protein or nucleotide")

signame2filenames = dict(zip(signames, filenames))

compute_targets = []

for alpha, info in config["alphabets"].items():
    ksize = info["ksize"]
    scaled = info["scaled"]
    compute_targets+=expand(os.path.join(compute_dir, input_type, "{alphabet}", "k{k}",  "{signame}_{alphabet}_scaled{scaled}_k{k}.sig"), signame = signames, alphabet=alpha, k=ksize, scaled=scaled)

rule all:
    input: compute_targets

# for protein signatures, multipy by 3 if necessary before calculating signature (sourmash v3.x)
ksize_multiplier = {"nucleotide": 1, "protein": 3, "dayhoff": 3, "hp":3, "translate_protein": 3, "translate_dayhoff": 3, "translate_hp": 3}
# "rna" is not sourmash cli friendly
moltype_map = {"nucleotide": "dna", "protein":"protein", "dayhoff": "dayhoff", "hp":"hp", "translate_protein": "protein", "translate_dayhoff": "dayhoff", "translate_hp": "hp"}


# compute sigs to general compute directory, regardles of the sbt sample name (so can reuse sigs for other sbts, if desired)
rule sourmash_compute_nucleotide:
    input: lambda w: signame2filenames[w.signame]
    output: os.path.join(compute_dir, "nucleotide", "{alphabet}", "k{k}", "{signame}_{alphabet}_scaled{scaled}_k{k}.sig")
    params:
        k= lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
        alpha_cmd = lambda w: "--" + moltype_map[w.alphabet],
        abund_cmd = "--track-abundance",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=20,
    log: os.path.join(logs_dir, "sourmash", "{signame}_{alphabet}_scaled{scaled}_k{k}.dna.compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{signame}_{alphabet}_scaled{scaled}_k{k}.dna.compute.benchmark")
    conda: "envs/sourmash3.3.yml"
    shell:
        """
        sourmash compute -k {params.k} --scaled={wildcards.scaled}  \
        {input} -o {output} {params.alpha_cmd} {params.abund_cmd} --merge={wildcards.signame:q} 2> {log}
        """

rule sourmash_compute_protein:
    input: lambda w: signame2filenames[w.signame]
    output: os.path.join(compute_dir, "protein", "{alphabet}", "k{k}", "{signame}_{alphabet}_scaled{scaled}_k{k}.sig")
    params:
        k= lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
        alpha_cmd = lambda w: "--" + moltype_map[w.alphabet],
        abund_cmd = "--track-abundance",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=20,
    log: os.path.join(logs_dir, "sourmash", "{signame}_{alphabet}_scaled{scaled}_k{k}.protein.compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{signame}_{alphabet}_scaled{scaled}_k{k}.protein.compute.benchmark")
    conda: "envs/sourmash3.3.yml"
    shell:
        """
        sourmash compute -k {params.k} --scaled={wildcards.scaled} --input-is-protein \
        {input} -o {output} {params.alpha_cmd} {params.abund_cmd} --merge={wildcards.signame:q} 2> {log}
        """
