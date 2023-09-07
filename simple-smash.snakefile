"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s gtdb-smash.snakefile --use-conda
"""

import os
import glob
import pandas as pd
#lineage_csv=config["lineage_csv"]

output_extensions = ["sig"]
output_targets, nucleotide_targets, translate_targets, protein_targets = [],[],[],[]

out_dir = config.get("out_dir", "output-gtdb-smash")
data_dir= os.path.join(out_dir, "data")
logs_dir = os.path.join(out_dir, "logs")
envs_dir = "envs"
compute_dir = os.path.join(out_dir, "compute")
compare_dir = os.path.join(out_dir, "compare")
plots_dir = os.path.join(out_dir, "plots")

# databases
gtdb_datadir = config.get("gtdb_datadir", "/home/ctbrown/gtdbtk/release89/fastani/database")
gtdb_protein_datadir = config.get("gtdb_protein_datadir", "")
gtdb_rna_datadir = config.get("gtdb_rna_datadir", "")

translate = config.get("translate", False)
protein_input = config.get("protein_input", False)
encoding_info = config["sourmash_params"]["encoding"]

genome_extension = "_genomic.fna.gz"

samples_csv = config["samples_csv"]
genome_info = pd.read_csv(samples_csv, sep="\t")
genome_lineages = set(genome_info["accession"].tolist())

path2acc= genome_info.groupby(["path"])['accession'].apply(list).to_dict()

## build sourmash targets
nucleotide_encodings = [] 
if config.get("genome_sigs", True):
    nucleotide_encodings+=["dna"]
if config.get("rna_sigs", False):
    nucleotide_encodings+=["rna"]


for enc in nucleotide_encodings:
    nucleotide_targets += expand(os.path.join(compute_dir, "{encoding}", "{sample}_{encoding}_scaled{scaled}_k{k}.{ext}"), sample = genome_lineages, encoding=enc, scaled=encoding_info[enc]["scaled"], k=encoding_info[enc]["ksize"], ext=output_extensions)
    nucleotide_targets += expand(os.path.join(plots_dir, f"{enc}", "{path}_{encoding}_scaled{scaled}_k{k}_compare.np.matrix.pdf"), encoding=enc, scaled=encoding_info[enc]["scaled"], k=encoding_info[enc]["ksize"], path = path2acc.keys())
    if translate: # just for kicks, translating dna seqs too
        for encoding in ["protein", "dayhoff", "hp"]:
            translate_targets+=expand(os.path.join(compute_dir, f"{enc}","{sample}_{encoding}_scaled{scaled}_k{k}.{ext}"), sample = genome_lineages, encoding=encoding, scaled=encoding_info[encoding]["scaled"], k=encoding_info[encoding]["ksize"], ext=output_extensions)
            translate_targets+= expand(os.path.join(plots_dir, f"{enc}", "{path}_{encoding}_scaled{scaled}_k{k}_compare.np.matrix.pdf"), encoding=encoding, scaled=encoding_info[encoding]["scaled"],k=encoding_info[encoding]["ksize"], path = path2acc.keys())
if protein_input:
    for encoding in ["protein", "dayhoff", "hp"]:
        protein_targets+=expand(os.path.join(compute_dir,"protein","{sample}_{encoding}_scaled{scaled}_k{k}.{ext}"), sample = genome_lineages, encoding=encoding, scaled=encoding_info[encoding]["scaled"],k=encoding_info[encoding]["ksize"], ext=output_extensions)
        protein_targets+= expand(os.path.join(plots_dir, "protein", "{path}_{encoding}_scaled{scaled}_k{k}_compare.np.matrix.pdf"), encoding=encoding, scaled=encoding_info[encoding]["scaled"],k=encoding_info[encoding]["ksize"], path = path2acc.keys())

output_targets = nucleotide_targets + protein_targets + translate_targets

rule all:
    input: output_targets

def find_genome_file(w):
    return glob.glob(os.path.join(gtdb_datadir, "*", f"*{w.sample}*"))

rule sourmash_compute_dna:
    input: find_genome_file
    output: os.path.join(compute_dir, "dna", "{sample}_{encoding}_scaled{scaled}_k{k}.sig")
    params:
        k= lambda w: (int(w.k) * ksize_multiplier[w.encoding]),
        scaled= lambda w: w.scaled,
        compute_moltypes= lambda w: w.encoding,
        input_is_protein=False,
        track_abundance=True,
    threads: 1
    resources:
        #mem_mb=3000,
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=600,
    log: os.path.join(logs_dir, "sourmash", "{sample}_{encoding}_scaled{scaled}_k{k}.dna.compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{sample}_{encoding}_scaled{scaled}_k{k}.dna.compute.benchmark")
    conda: "envs/sourmash3.3.yml"
    script: "scripts/sourmash-compute.wrapper.py"

def find_corresponding_protein_file(w):
    return glob.glob(os.path.join(gtdb_protein_datadir, "*", f"*{w.sample}*"))

rule sourmash_compute_protein:
    input: find_corresponding_protein_file 
    output: os.path.join(compute_dir, "protein", "{sample}_{encoding}_scaled{scaled}_k{k}.sig")
    params:
        k= lambda w: w.k,
        scaled= lambda w: w.scaled,
        compute_moltypes= lambda w: w.encoding,
        input_is_protein=True,
        track_abundance=True,
    threads: 1
    resources:
        #mem_mb=3000,
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash", "{sample}_{encoding}_scaled{scaled}_k{k}.protein.compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{sample}_{encoding}_scaled{scaled}_k{k}.protein.compute.benchmark")
    conda: "envs/sourmash3.3.yml"
    script: "scripts/sourmash-compute.wrapper.py"

def find_corresponding_rna_file(w):
    return glob.glob(os.path.join(gtdb_rna_datadir, "*", f"*{w.sample}*"))

# "rna" is not sourmash-friendly
moltype_map = {"rna": "dna", "dna": "dna", "protein":"protein", "dayhoff": "dayhoff", "hp":"hp"}
# now passing in the protein ksize. multipy by 3 if necessary before calculating signature
ksize_multiplier = {"rna": 1, "dna": 1, "protein": 3, "dayhoff": 3, "hp":3}  

rule sourmash_compute_rna:
    input: find_corresponding_rna_file 
    output: os.path.join(compute_dir, "rna", "{sample}_{encoding}_scaled{scaled}_k{k}.sig")
    params:
        k= lambda w: (int(w.k) * ksize_multiplier[w.encoding]),
        scaled= lambda w: w.scaled,
        compute_moltypes= lambda w: moltype_map[w.encoding],
        input_is_protein=False,
        track_abundance=True,
    threads: 1
    resources:
        #mem_mb=3000,
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash", "{sample}_{encoding}_scaled{scaled}_k{k}.rna.compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{sample}_{encoding}_scaled{scaled}_k{k}.rna.compute.benchmark")
    conda: "envs/sourmash3.3.yml"
    script: "scripts/sourmash-compute.wrapper.py"

def aggregate_sigs(w):
    siglist=[]
    accessions=path2acc[w.path]
    for acc in accessions:
        siglist+=[os.path.join(compute_dir, w.moltype, f"{acc}_{w.encoding}_scaled{w.scaled}_k{w.k}.sig")]
    return siglist

rule sourmash_compare:
    input: sigs=aggregate_sigs
    output:
        np=os.path.join(compare_dir, "{moltype}", "{path}_{encoding}_scaled{scaled}_k{k}_compare.np"),
        csv=os.path.join(compare_dir, "{moltype}", "{path}_{encoding}_scaled{scaled}_k{k}_compare.csv"),
    threads:1
    resources:
        mem_mb=2000,
        runtime=120,
    params:
        include_encodings = lambda w: f"{w.encoding}",
        exclude_encodings = ["nucl", "protein", "dayhoff", "hp"], # this will excude everything except for included encoding
        k = lambda w: f"{w.k}",
    conda: "envs/sourmash3.3.yml"
    script: "scripts/sourmash-compare.wrapper.py"

# sourmash plot each compare matrix numpy output
localrules: sourmash_plot

rule sourmash_plot:
    input: os.path.join(compare_dir, "{moltype}", "{path}_{encoding}_scaled{scaled}_k{k}_compare.np")
    output: os.path.join(plots_dir, "{moltype}", "{path}_{encoding}_scaled{scaled}_k{k}_compare.np.matrix.pdf")
    params:
        plot_dir=lambda w: os.path.join(plots_dir, w.moltype)
    conda: "envs/sourmash3.3.yml"
    shell:
        """
        sourmash plot --output-dir {params.plot_dir} --labels --pdf {input}
        """
