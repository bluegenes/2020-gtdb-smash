"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s sourmash-compare.snakefile --use-conda
"""

import os
import glob
#lineage_csv=config["lineage_csv"]

out_dir = config.get("out_dir", "output-gtdb-smash")
logs_dir = os.path.join(out_dir, "logs")
envs_dir = "envs"
compute_dir = os.path.join(out_dir, "compute")
compare_dir = os.path.join(out_dir, "compare")

# names are getting crazy. let's just use `--traverse-directory` to build compare matrix! Plot won't work at this scale, but could do tsne, as with 2020-pep hashclust


rule all:
    #expand(os.path.join(out_dir, "plass", "plots", "plass_cdhit100_scaled{scaled}_{moltype}_k{k}_above{min_count}_{ctype}_compare.np.matrix.pdf"), scaled=[1], min_count=[2], moltype="dayhoff",k=ksizes, ctype=["cosine", "jaccard"]) 
    expand(os.path.join(out_dir, "compare", "{type}"))




# build all compare matrices: np and csv output
rule sourmash_compare_cosine: # use abundances
    input: sigs=expand(os.path.join(out_dir, "hashclust", "sigs", "{sample}_plass100_scaled{{scaled}}_{{moltype}}_k{{ksize}}_above{{min_count}}_renamed.sig"), sample=SAMPLES)
    output:
        np=os.path.join(out_dir, "plass", "compare", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_cosine_compare.np"),
        csv=os.path.join(out_dir, "plass", "compare", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_cosine_compare.csv"),
    params:
        include_encodings = lambda w: w.moltype,
        exclude_encodings = ["nucl", "protein", "dayhoff", "hp"], # this will exclude everything except for included encoding
        k = lambda w: w.ksize,
    log: os.path.join(logs_dir, "sourmash", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_cosine_compare.log")
    benchmark: os.path.join(logs_dir, "sourmash", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_cosine_compare.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "sourmash-compare.wrapper.py")

rule sourmash_compare_jaccard: #ignore abundances
    input: sigs=expand(os.path.join(out_dir, "hashclust", "sigs", "{sample}_plass100_scaled{{scaled}}_{{moltype}}_k{{ksize}}_above{{min_count}}_renamed.sig"), sample=SAMPLES)
    output:
        np=os.path.join(out_dir, "plass", "compare", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_jaccard_compare.np"),
        csv=os.path.join(out_dir, "plass", "compare", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_jaccard_compare.csv"),
    params:
        ignore_abundance = True,
        include_encodings = lambda w: w.moltype,
        exclude_encodings = ["nucl", "protein", "dayhoff", "hp"], # this will exclude everything except for included encoding
        k = lambda w: w.ksize,
    log: os.path.join(logs_dir, "sourmash", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_jaccard_compare.log")
    benchmark: os.path.join(logs_dir, "sourmash", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_jaccard_compare.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "sourmash-compare.wrapper.py")

rule sourmash_plot_cosine:
    input: rules.sourmash_compare_cosine.output.np
    output: os.path.join(out_dir, "plass", "plots", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_cosine_compare.np.matrix.pdf"),
    params:
        plot_dir=os.path.join(out_dir, "plass", "plots")
    log: os.path.join(logs_dir, "sourmash", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_cosine_plot.log")
    benchmark: os.path.join(logs_dir, "sourmash", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_cosine_plot.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    shell:
        """
        sourmash plot --output-dir {params.plot_dir} --labels --pdf {input} 2> {log}
        """

rule sourmash_plot_jaccard:
    input: rules.sourmash_compare_jaccard.output.np
    output: os.path.join(out_dir, "plass", "plots", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_jaccard_compare.np.matrix.pdf"),
    params:
        plot_dir=os.path.join(out_dir, "plass", "plots")
    log: os.path.join(logs_dir, "sourmash", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_jaccard_plot.log")
    benchmark: os.path.join(logs_dir, "sourmash", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_jaccard_plot.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    shell:
        """
        sourmash plot --output-dir {params.plot_dir} --labels --pdf {input} 2> {log}
        """
