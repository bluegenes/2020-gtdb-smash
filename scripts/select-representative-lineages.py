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
    else:
        row["superkingdom"] = lineages[0]
        row["phylum"] = lineages[1]
        row["class"] = lineages[2]
        row["order"] = lineages[3]
        row["family"] = lineages[4]
        row["genus"] = lineages[5]
        row["species"] = lineages[6]
    return row


def main(args):
    #read in csv
    queryDF= pd.read_csv(args.query_csv)

    # if this is an evol paths dataset, drop anchor species (confound analysis bc all genomes w/in a path are picked from them)
    final_column_order = ["accession", "superkingdom", "phylum", "class", "order", "family", "genus", "species","filename"]
    cols_to_add=[]
    if "rank" in queryDF.columns:
        queryDF = queryDF[queryDF["rank"] != "species"]
        cols_to_add = ["path","rank"]

    # if we have a single "lineage" column we need to split:
    if "lineage" in queryDF.columns:
        queryDF = queryDF.apply(split_lineages_to_columns, axis=1, args=(["lineage"]))
        cols_to_add+=["lineage"]

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
    repDF = queryDF.groupby([rep_col], as_index=False).nth(0)

    # reorder columns so it works for lca index
    repDF = repDF[final_column_order+ cols_to_add]
    repDF.to_csv(args.output_csv, index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("query_csv")
    p.add_argument("--output-csv")
    p.add_argument("--representative-rank", default="family", help="rank for which to select representative genomes")
    args = p.parse_args()
    if not args.output_csv:
        args.output_csv=(args.query_csv).rsplit(".", 1)[0] + f".representative-at-{args.representative_rank}.csv"
    sys.exit(main(args))
