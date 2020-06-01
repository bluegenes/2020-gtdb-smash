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

index_dir = config.get(out_dir, "gtdb-sbts")
dist_dir = os.path.join(out_dir, "distances")


sbt_targets, query_targets=[],[]


#for
#for enc in nucleotide_encodings:
    #nucleotide_targets += expand(os.path.join(compute_dir, "{encoding}", "{sample}_{encoding}_scaled{scaled}_k{k}.{ext}"), sample = genome_lineages, encoding=enc, scaled=encoding_info[enc]["scaled"], k=encoding_info[enc]["ksize"], ext=output_extensions)
    #nucleotide_targets += expand(os.path.join(plots_dir, f"{enc}", "{path}_{encoding}_scaled{scaled}_k{k}_compare.np.matrix.pdf"), encoding=enc, scaled=encoding_info[enc]["scaled"], k=encoding_info[enc]["ksize"], path = path2acc.keys())
#    if translate: # just for kicks, translating dna seqs too
#        for encoding in ["protein", "dayhoff", "hp"]:
#            translate_targets+=expand(os.path.join(compute_dir, f"{enc}","{sample}_{encoding}_scaled{scaled}_k{k}.{ext}"), sample = genome_lineages, encoding=encoding, scaled=encoding_info[encoding]["scaled"], k=encoding_info[encoding]["ksize"], ext=output_extensions)
#            translate_targets+= expand(os.path.join(plots_dir, f"{enc}", "{path}_{encoding}_scaled{scaled}_k{k}_compare.np.matrix.pdf"), encoding=encoding, scaled=encoding_info[encoding]["scaled"],k=encoding_info[encoding]["ksize"], path = path2acc.keys())
#if protein_input:
#    for encoding in ["protein", "dayhoff", "hp"]:
#        protein_targets+=expand(os.path.join(compute_dir,"protein","{sample}_{encoding}_scaled{scaled}_k{k}.{ext}"), sample = genome_lineages, encoding=encoding, scaled=encoding_info[encoding]["scaled"],k=encoding_info[encoding]["ksize"], ext=output_extensions)
#        protein_targets+= expand(os.path.join(plots_dir, "protein", "{path}_{encoding}_scaled{scaled}_k{k}_compare.np.matrix.pdf"), encoding=encoding, scaled=encoding_info[encoding]["scaled"],k=encoding_info[encoding]["ksize"], path = path2acc.keys())

#output_targets = nucleotide_targets + protein_targets + translate_targets

#def build_target(folder, samplelist, encoding_set, encoding_info, out_extensions):
#    targets = []
#    for enc in encodings:
#        targets = expand(os.path.join(folder, enc, "{sample}_{encoding}_scaled{scaled}_k{k}.{ext}"), sample = samplelist, encoding=enc, scaled=encoding_info[enc]["scaled"],              k=encoding_info[enc]["ksize"], ext=out_extensions)
#    return targets

sampleInfo=config["foraging_samples"]
accession2filenames = {}

for sample, info in sampleInfo.items():    
    if info.get("grow_sbt"):
        # map genome acession: fasta files find fasta files 
        info_csv = pd.read_csv(info["query_csv"])
        # accession:: fasta file map
        fasta_dir= info["grow_sbt"]["input_path"]
        info_csv["filename"] = info_csv["filename"].apply(lambda x: os.path.join(fasta_dir, x))       #x.replace(os.path.dirname(x), dst))
        input_type = info["input_type"] #protein, dna, rna
        accession2filenames[input_type] = pd.Series(info_csv.filename.values,index=info_csv.accession).to_dict()
        # add accessions to sampleInfo dict
        sampleInfo[sample]["accessions"] = info_csv.accession.to_list()
        # build sbt targets
        index_dir= info["grow_sbt"]["sbt_outdir"]
        for alpha, alphainfo in info["alphabet"].items():
            sbt_targets+=expand(os.path.join(index_dir,"{sample}.{alphabet}_scaled{scaled}_k{k}.sbt.zip"), sample=sample, alphabet=alpha, scaled=alphainfo["scaled"], k=alphainfo["ksizes"])
    for alpha, alphainfo in info["alphabet"].items():
        query_targets+=expand(os.path.join(dist_dir, "{sample}.{alphabet}_scaled{scaled}_k{k}.jaccard_from_species.csv"), sample=sample, alphabet=alpha, scaled=alphainfo["scaled"], k=alphainfo["ksizes"])

rule all:
    input: sbt_targets + query_targets

# for protein signatures, multipy by 3 if necessary before calculating signature (sourmash v3.x)
ksize_multiplier = {"rna": 1, "dna": 1, "protein": 3, "dayhoff": 3, "hp":3, "translate_protein": 3, "translate_dayhoff": 3, "translate_hp": 3}
# "rna" is not sourmash cli friendly
moltype_map = {"rna": "dna", "dna": "dna", "protein":"protein", "dayhoff": "dayhoff", "hp":"hp", "translate_protein": "protein", "translate_dayhoff": "dayhoff", "translate_hp": "hp"}


# compute sigs to general compute directory, regardles of the sbt sample name (so can reuse sigs for other sbts, if desired)
rule sourmash_compute_dna:
    input: lambda w: accession2filenames["dna"][w.accession]
    output: os.path.join(compute_dir, "dna", "{accession}_{alphabet}_scaled{scaled}_k{k}.sig")
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
    output: os.path.join(compute_dir, "protein", "{accession}_{alphabet}_scaled{scaled}_k{k}.sig")
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
    output: os.path.join(compute_dir, "rna", "{accession}_{alphabet}_scaled{scaled}_k{k}.sig")
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
    sigfile = os.path.join(compute_dir, input_type, f"{{acc}}_{w.alphabet}_scaled{w.scaled}_k{w.k}.sig")
    siglist=expand(sigfile, acc=sampleInfo[w.sample]["accessions"])
    #accessions=path2acc[w.path]
    #for acc in accessions:
    #    siglist+=[os.path.join(compute_dir, w.moltype, f"{acc}_{w.encoding}_scaled{w.scaled}_k{w.k}.sig")]
    return siglist

rule grow_sbt:
    input: aggregate_sigs
    output: 
        sbt=os.path.join(index_dir,"{sample}.{alphabet}_scaled{scaled}_k{k}.sbt.zip"),
    threads: 1
    params:
        #input_dir= compute_dir  #lambda w: sampleInfo[w.sample]['grow_sbt']["input_path"],
        #subset_info_colname=lambda w: sampleInfo[w.sample]['grow_sbt'].get('subset_info_colname', 'filename'),
        alpha= lambda w: w.alphabet.rsplit("translate_")[1] if w.alphabet.startswith("translate") else w.alphabet, # remove translate
        translate = lambda w: " --translate " if w.alphabet.startswith("translate") else "",
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
        runtime=6000,
    log: os.path.join(logs_dir, "grow-sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.grow.log")
    benchmark: os.path.join(logs_dir, "grow-sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.grow.benchmark")
    conda: "envs/forage-env.yml"
    shell:
        """
        python scripts/grow-sbtmh.py {input} --sbt {output.sbt} --ksize {wildcards.k} --scaled {wildcards.scaled} --alphabet {params.alpha} {params.translate} 2> {log}
        """

def find_forage_inputs(w):
    # just require this query to have filenames of interest as a column
    query = sampleInfo[w.sample]["query_csv"]
    if sampleInfo[w.sample].get("grow_sbt", False):
        sbt = rules.grow_sbt.output.sbt,
    else:
        sbt=sampleInfo[w.sample]["use_existing_sbt"]
    return {"query_csv": query, "sbt": sbt}


rule calculate_jaccard_from_common_ancestor:
    input: unpack(find_forage_inputs)
    output: 
        csv=os.path.join(dist_dir, "{sample}.{alphabet}_scaled{scaled}_k{k}.jaccard_from_species.csv"),
        boxplot=os.path.join(dist_dir, "plots", "{sample}.{alphabet}_scaled{scaled}_k{k}.jaccard_from_species.svg"),
    params:
        signature_name_column= lambda w: sampleInfo[w.sample].get('query_signature_name_column_name', 'filename')
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=1200,
    log: os.path.join(logs_dir, "forage", "{sample}.{alphabet}_scaled{scaled}_k{k}.forage.log")
    benchmark: os.path.join(logs_dir, "forage", "{sample}.{alphabet}_scaled{scaled}_k{k}.forage.benchmark")
    conda: "envs/forage-env.yml"
    shell:
        """
        python scripts/forage-sbt.py {input.sbt} --query-csv {input.query_csv} --signature-name-column {params.signature_name_column} \
        --distance-from-species-csv {output.csv} --distance-from-species-plot {output.boxplot} 2>{log}
        """

