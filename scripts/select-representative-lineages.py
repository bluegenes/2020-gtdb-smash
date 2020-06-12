import os
import sys
import glob
import argparse
import pandas as pd


#def select_representative_lineages(csv):
    # loop through csv
    # select set of unique genomes for each rank
    # print csv for just those genomes
#    return newDF


def split_lineages_to_columns(row, taxo_col):
    rank_dir = {"d__": "superkingdom", "p__": "phylum", "c__": "class", "o__": "order", "f__": "family", "g__": "genus", "s__": "species"}
    # split taxonomy column on ";"
    lineages = row[taxo_col].split(";")
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
    # glob for filenames
    #queryDF = queryDF.apply(find_filename, axis=1, args=(args.sigfile_directory,args.identifier_column, args.full_paths))

    taxo_col = args.taxonomy_column
    if taxo_col not in queryDF.columns:
        sys.stderr.write(f"representive rank {rep_col} is not present in the columns of your input query csv\n\n")
        sys.exit(-1)
    # drop anchor species (confound analysis bc all genomes w/in a path are picked from them)
    taxoDF = queryDF[queryDF["rank"] != "species"]
    taxoDF = taxoDF.apply(split_lineages_to_columns, axis=1, args=([taxo_col]))

    rep_col = args.representative_rank
    if rep_col not in taxoDF.columns:
        sys.stderr.write(f"representive rank {rep_col} is not present in the columns of your input query csv\n\n")
        sys.exit(-1)
    # fast way: use pandas groupby to select first unique member of each group
    #re: groupby:: passing as_index=False, will achieve a filtration, which returns the grouped row.
    # take nth row, here row 0 = first row. If expanding to later rows, watch out for Nans
    # do queryDF.groupby([rep_col].cumcount() to see order
    # https://pandas.pydata.org/pandas-docs/stable/user_guide/groupby.html
    repDF = taxoDF.groupby([rep_col], as_index=False).nth(0)

    # write taxoDF, lineage, filename
    repDF[["accession", "lineage", "filename"]].to_csv(args.output_csv, index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("query_csv")
    p.add_argument("--output-csv")
    p.add_argument("--taxonomy-column", default="lineage", help="column in query_csv that can be used to identify corresponding filename")
    p.add_argument("--representative-rank", default="family", help="rank for which to select representative genomes")
    args = p.parse_args()
    if not args.output_csv:
        args.output_csv=(args.query_csv).rsplit(".", 1)[0] + f".representative-at-{args.representative_rank}.csv"
    sys.exit(main(args))
