"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s grow-and-forage-sbt.snakefile --configfile config/grow_gtdb_sbts.yml --use-conda
"""

import os
import glob
import pandas as pd

def csv_reader(groups_file):
    # autodetect format
    if '.tsv' in groups_file or '.csv' in groups_file:
        separator = '\t'
        if '.csv' in groups_file:
            separator = ','
        try:
            samples = pd.read_csv(groups_file, dtype=str, sep=separator)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {groups_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    elif '.xls' in samples_file:
        try:
            samples = pd.read_excel(groups_file, dtype=str, sep=separator)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {groups_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    return samples

out_dir = config.get("out_dir", "output-pathdistances")
logs_dir = os.path.join(out_dir, "logs")
index_dir = config.get("index_dir", "gtdb-sbts")
dist_dir = os.path.join(out_dir, "distances")

sampleInfo=config["foraging_samples"]

sbt_targets=[]
query_targets=[]

for sample, info in sampleInfo.items():    
    if info.get("grow_sbt"):
        # build sbt targets
        sbt_extensions=["sbt.zip", "signame2filenames.csv", "md5_duplicates.csv"]
        index_dir= info["grow_sbt"]["sbt_outdir"]
        for alpha, alphainfo in info["alphabet"].items():
            sbt_targets+=expand(os.path.join(index_dir,"{sample}_{alphabet}_scaled{scaled}_k{k}.{ext}"), sample=sample, alphabet=alpha, scaled=alphainfo["scaled"], k=alphainfo["ksizes"], ext=sbt_extensions)
    for alpha, alphainfo in info["alphabet"].items():
        query_targets+=expand(os.path.join(dist_dir, "{sample}_{alphabet}_scaled{scaled}_k{k}.jaccard_from_species.csv"), sample=sample, alphabet=alpha, scaled=alphainfo["scaled"], k=alphainfo["ksizes"])

rule all:
    input: sbt_targets + query_targets

# ooooof, snakemake finds all these input files just to calculate the dag!
# do this differently:: just pass in subset csv, handle within the pyscript!
# == ALWAYS INDEX FROM DIR
#def find_sbt_input(w):
#    input_files=[]
#    if sampleInfo[w.sample]["grow_sbt"].get("subset_csv"):
#        accs = csv_reader(sampleInfo[w.sample]["grow_sbt"]["subset_csv"])["accession"].tolist()
#        for acc in accs:
#            input_files+=glob.glob(os.path.join(sampleInfo[w.sample]["grow_sbt"].get("input_path", ""), "*", f"*{acc}*"))
#    else:
#        input_files = sampleInfo[w.sample]["grow_sbt"].get("input_path", "")
#    return input_files

rule grow_sbt:
 #   input: find_sbt_input
    output: 
        sbt=os.path.join(index_dir,"{sample}_{alphabet}_scaled{scaled}_k{k}.sbt.zip"),
    threads: 1
    params:
        input_dir= lambda w: sampleInfo[w.sample]['grow_sbt']["input_path"],
        subset_csv=lambda w: sampleInfo[w.sample]['grow_sbt'].get('subset_csv', ''),
        subset_info_colname=lambda w: sampleInfo[w.sample]['grow_sbt'].get('subset_info_colname', 'accession')
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=6000,
    log: os.path.join(logs_dir, "grow-sbt", "{sample}_{alphabet}_scaled{scaled}_k{k}.grow.log")
    benchmark: os.path.join(logs_dir, "grow-sbt", "{sample}_{alphabet}_scaled{scaled}_k{k}.grow.benchmark")
    conda: "envs/forage-env.yml"
    shell:
        """
        python scripts/grow-sbtmh.py {params.input_dir} --input-is-directory --sbt {output.sbt} --ksize {wildcards.k} --scaled {wildcards.scaled} --alphabet {wildcards.alphabet} --subset-csv {params.subset_csv} --subset-info-colname {params.subset_info_colname} 2> {log}
        """

def find_forage_inputs(w):
    # just require this query to have filenames of interest as a column
    query = sampleInfo[w.sample]["query_csv"]
    if sampleInfo[w.sample].get("grow_sbt", False):
        sbt = rules.grow_sbt.output.sbt,
    else:
        sbt=sampleInfo[w.sample]["use_existing_sbt"]["sbt"]
    return {"query_csv": query, "sbt": sbt}


rule calculate_jaccard_from_common_ancestor:
    input: unpack(find_forage_inputs)
    output: 
        csv=os.path.join(dist_dir, "{sample}_{alphabet}_scaled{scaled}_k{k}.jaccard_from_species.csv"),
        boxplot=os.path.join(dist_dir, "plots", "{sample}_{alphabet}_scaled{scaled}_k{k}.jaccard_from_species.svg"),
    params:
        signature_name_column: lambda w: sampleInfo[w.sample].get('query_signature_name_column_name', 'filename')
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=1200,
    log: os.path.join(logs_dir, "forage", "{sample}_{alphabet}_scaled{scaled}_k{k}.forage.log")
    benchmark: os.path.join(logs_dir, "forage", "{sample}_{alphabet}_scaled{scaled}_k{k}.forage.benchmark")
    conda: "envs/forage-env.yml"
    shell:
        ## to do: separate out the query csv processing --> other script. assume this query has a column with appropriate filenames or signature names
        #python scripts/forage-sbt.py {input.sbt} --database_csv {input.database_csv}  --query_csv {input.query_csv} --duplicates_csv {input.duplicates_csv} --distance_from_species_csv {output.csv} --distance_from_species_plot {output.boxplot} 2>{log}
       
       # this is how it _should_end up: 
        """
        python scripts/forage-sbt.py {input.sbt} --query_csv {input.query_csv} --query-signature-name-column {params.signature_name_column}--distance_from_species_csv {output.csv} --distance_from_species_plot {output.boxplot} 2>{log}
        """

