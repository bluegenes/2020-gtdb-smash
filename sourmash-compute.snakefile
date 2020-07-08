"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Builds signatures with abundance information.
"""

import os
import glob
import pandas as pd
from collections import defaultdict

out_dir = config["out_dir"]
compute_dir = config["sigs_dir"]
logs_dir = os.path.join(out_dir, "logs")
envs_dir = "envs"

sampleInfo=config["samples"]
accession2filenames = defaultdict(dict)
accession2signame = defaultdict(dict)

sig_targets=[]

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
    for alpha, alphainfo in info["alphabet"].items():
        sig_targets+=expand(os.path.join(compute_dir, input_type, "{alphabet}", "k{k}", "{accession}_{alphabet}_scaled{scaled}_k{k}.sig"), accession=sample_acc2file.keys(), alphabet=alpha,scaled=alphainfo["scaled"], k=alphainfo["ksizes"])


rule all:
    input: sig_targets

### compute rules can be used in any other script via `include`
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
