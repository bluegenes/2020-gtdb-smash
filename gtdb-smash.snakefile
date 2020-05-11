"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s gtdb-smash.snakefile --use-conda
"""

import os
import glob
#lineage_csv=config["lineage_csv"]

output_extensions = ["sig"]
output_targets, nucleotide_targets, translate_targets, protein_targets = [],[],[],[]

out_dir = config.get("out_dir", "output-gtdb-smash")
data_dir= os.path.join(out_dir, "data")
logs_dir = os.path.join(out_dir, "logs")
compute_dir = os.path.join(out_dir, "compute")

# databases
gtdb_datadir = config.get("gtdb_datadir", "/home/ctbrown/gtdbtk/release89/fastani/database")
gtdb_protein_datadir = config.get("gtdb_protein_datadir", "")
gtdb_rna_datadir = config.get("gtdb_rna_datadir", "")

translate = config.get("translate", False)
protein_input = config.get("protein_input", False)
encoding_info = config["sourmash_params"]["encoding"]

genome_lineage_list =config["genome_lineage_list"] # generated by extract-robust-lineages
genome_extension = "_genomic.fna.gz"
genome_lineages = [line.strip().split(genome_extension)[0] for line in open(genome_lineage_list)][:10]

nucleotide_encodings = [] 
if config.get("genome_sigs", True):
    nucleotide_encodings+=["dna"]
elif config.get("rna_sigs", False):
    nucleotide_encodings+=["rna"]

for enc in nucleotide_encodings:
    nucleotide_targets += expand(os.path.join(compute_dir, "{encoding}", "{sample}_{encoding}_scaled{scaled}_k{k}.{ext}"), sample = genome_lineages, encoding=enc, scaled=encoding_info[enc]["scaled"], k=encoding_info[enc]["ksize"], ext=output_extensions)
    if translate: # just for kicks, translating dna seqs too
        for encoding in ["protein", "dayhoff", "hp"]:
            # ksizes must be multiplied by 3 (give nucleotide ksizes to sourmash)
            translate_ksizes = [int(k)*3 for k in encoding_info[encoding]["ksize"]]
            translate_targets+=expand(os.path.join(compute_dir, {enc},"{sample}_{encoding}_scaled{scaled}_k{k}.{ext}"), sample = genome_lineages, encoding=encoding, scaled=encoding_info[encoding]["scaled"], k=translate_ksizes, ext=output_extensions)

if protein_input:
    for encoding in ["protein", "dayhoff", "hp"]:
        protein_targets+=expand(os.path.join(compute_dir,"protein","{sample}_{encoding}_scaled{scaled}_pk{k}.{ext}"), sample = genome_lineages, encoding=encoding, scaled=encoding_info[encoding]["scaled"],k=encoding_info[encoding]["ksize"], ext=output_extensions)

output_targets = nucleotide_targets + protein_targets + translate_targets

rule all:
    input: output_targets

#localrules: get_genomes, get_proteins

#rule get_genomes:
#    input: os.path.join(gtdb_datadir, "{sample}_genomic.fna.gz")
#    output: os.path.join(data_dir, "{sample}_genomic.fna.gz") 
#    shell: 
#        """    
#        cp {input} {output}
#        """

# too many - don't want to copy all over!
#rule get_proteins:
    #input: find_corresponding_protein_file 
    #output: os.path.join(data_dir, "{sample}_protein.faa.gz")
    #shell:
    #    """
    #    cp {input} {output}
    #    """

def find_genome_file(w):
    return glob.glob(os.path.join(gtdb_datadir, "*", f"*{w.sample}*"))

rule sourmash_compute_dna:
    #input: rules.get_genomes.output
    input: find_genome_file
    output: os.path.join(compute_dir, "dna", "{sample}_{encoding}_scaled{scaled}_k{k}.sig")
    params:
        k= lambda w: w.k,
        scaled= lambda w: w.scaled,
        compute_moltypes= lambda w: w.encoding,
        input_is_protein=False,
        track_abundance=True,
    threads: 1
    resources:
        #mem_mb=3000,
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=200,
    log: os.path.join(logs_dir, "sourmash", "{sample}_{encoding}_scaled{scaled}_k{k}.dna.compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{sample}_{encoding}_scaled{scaled}_k{k}.dna.compute.benchmark")
    conda: "envs/sourmash3.3.yml"
    script: "scripts/sourmash-compute.wrapper.py"

def find_corresponding_protein_file(w):
    return glob.glob(os.path.join(gtdb_protein_datadir, "*", f"*{w.sample}*"))

rule sourmash_compute_protein:
    input: find_corresponding_protein_file 
    output: os.path.join(compute_dir, "protein", "{sample}_{encoding}_scaled{scaled}_pk{k}.sig")
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
        runtime=200,
    log: os.path.join(logs_dir, "sourmash", "{sample}_{encoding}_scaled{scaled}_k{k}.protein.compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{sample}_{encoding}_scaled{scaled}_k{k}.protein.compute.benchmark")
    conda: "envs/sourmash3.3.yml"
    script: "scripts/sourmash-compute.wrapper.py"

def find_corresponding_rna_file(w):
    return glob.glob(os.path.join(gtdb_rna_datadir, "*", f"*{w.sample}*"))

# "rna" is not sourmash-friendly
moltype_map = {"rna": "dna", "dna": "dna", "protein":"protein", "dayhoff": "dayhoff", "hp":"hp"}

rule sourmash_compute_rna:
    input: find_corresponding_rna_file 
    output: os.path.join(compute_dir, "rna", "{sample}_{encoding}_scaled{scaled}_k{k}.sig")
    params:
        k= lambda w: w.k,
        scaled= lambda w: w.scaled,
        compute_moltypes= lambda w: moltype_map[w.encoding],
        input_is_protein=True,
        track_abundance=True,
    threads: 1
    resources:
        #mem_mb=3000,
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=200,
    log: os.path.join(logs_dir, "sourmash", "{sample}_{encoding}_scaled{scaled}_k{k}.rna.compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{sample}_{encoding}_scaled{scaled}_k{k}.rna.compute.benchmark")
    conda: "envs/sourmash3.3.yml"
    script: "scripts/sourmash-compute.wrapper.py"
