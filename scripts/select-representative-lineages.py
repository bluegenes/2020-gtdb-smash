import os
import sys
import glob
import argparse
import pandas as pd


def split_lineages_to_columns(row, taxo_col):
    #rank_d = {"d__": "superkingdom", "p__": "phylum", "c__": "class", "o__": "order", "f__": "family", "g__": "genus", "s__": "species"}
    # split taxonomy column on ";"
    lineages = row[taxo_col].split(";")
    if lineages[0].startswith("d__"):
        row["superkingdom"] = lineages[0].split("d__")[1]
        row["phylum"] = lineages[1].split("p__")[1]
        row["class"] = lineages[2].split("c__")[1]
        row["order"] = lineages[3].split("o__")[1]
        row["family"] = lineages[4].split("f__")[1]
        row["genus"] = lineages[5].split("g__")[1]
        row["species"] = lineages[6].split("s__")[1]
        row["strain"] = ""
    else:
        row["superkingdom"] = lineages[0]
        row["phylum"] = lineages[1]
        row["class"] = lineages[2]
        row["order"] = lineages[3]
        row["family"] = lineages[4]
        row["genus"] = lineages[5]
        row["species"] = lineages[6]
        row["strain"] = ""
    return row


def main(args):
    #read in csv
    queryDF= pd.read_csv(args.query_csv)
    # if this is an evol paths dataset, drop anchor species (confound analysis bc all genomes w/in a path are picked from them)
    final_column_order = ["accession", "superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain", "signame", "filename"]
    cols_to_add=[]
    if "rank" in queryDF.columns:
        queryDF = queryDF[queryDF["rank"] != "species"]
        cols_to_add = ["path","rank"]

    # if we have a single "lineage" column we need to split:
    if "lineage" in queryDF.columns:
        queryDF = queryDF.apply(split_lineages_to_columns, axis=1, args=(["lineage"]))
        cols_to_add+=["lineage"]
    elif "strain" not in queryDF.columns:
        queryDF["strain"] = ""
    # now select representatives at the chosen level
    rep_col = args.representative_rank.lower()
    if rep_col not in queryDF.columns:
        sys.stderr.write(f"representive rank {rep_col} is not present in the columns of your input query csv\n\n")
        sys.exit(-1)

    # easy way: use pandas groupby to select first unique member of each group
    #re: groupby:: passing as_index=False, will achieve a filtration, which returns the grouped row.
    # take nth row, here row 0 = first row. If expanding to later rows, watch out for Nans
    # do queryDF.groupby([rep_col].cumcount() to see order
    # https://pandas.pydata.org/pandas-docs/stable/user_guide/groupby.html
    select_n = args.nth_to_select
    repDF = queryDF.copy().groupby([rep_col], as_index=False).nth(select_n)
    repDF.dropna(inplace=True)

    # add full filepath to filename
    #repDF["filepath"] = repDF["filename"].apply(lambda x: os.path.join(args.fasta_path, x)) #, axis=1)
    repDF["signame"] = repDF["accession"] + " " + repDF["species"]

    # reorder columns so it works for lca index
    repDF = repDF[final_column_order+ cols_to_add]
    repDF.to_csv(args.output_csv, index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("query_csv")
    p.add_argument("--output-csv")
    p.add_argument("--fasta-path", default="/home/ntpierce/2020-gtdb-smash/gtdb_r89_rep_genomes_faa", help="directory to find fasta files")
    p.add_argument("--representative-rank", default="family", help="rank for which to select representative genomes")
    p.add_argument("--nth-to-select", type=int, default=0, help="select nth genome as the 'representative': default=0 (select first)")
    p.add_argument("--keep-all", action="store_true", help="keep all genomes")
    args = p.parse_args()
    if not args.output_csv:
        if args.keep_all:
            args.output_csv=(args.query_csv).rsplit(".", 1)[0] + f".with-signames.csv"
        else:
            n = str(args.nth_to_select)
            args.output_csv=(args.query_csv).rsplit(".", 1)[0] + f".n{n}th-representative-at-{args.representative_rank}.csv"
    sys.exit(main(args))
