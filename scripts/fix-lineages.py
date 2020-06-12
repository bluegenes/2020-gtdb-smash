import os
import sys
import glob
import argparse
import pandas as pd

def try_reading_csv(groups_file):
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
    elif '.xls' in groups_file:
        try:
            samples = pd.read_excel(groups_file, dtype=str, sep=separator)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {groups_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    return samples


def main(args):
    #read in csv
    queryDF= try_reading_csv(args.query_csv)
    linDF= try_reading_csv(args.lineages_csv)
    # join superkingdom--> species columns to "lineage" column
    linDF["lineage"] = linDF["superkingdom"] + ";" + linDF["phylum"]  + ";" + linDF["class"]  + ";" + linDF["order"]  + ";" + linDF["family"] + ";" + linDF["genus"]  + ";" + linDF["species"]

    # drop incorrect lineages
    queryDF.drop(columns = ["lineage"], inplace=True)
    # merge to get correct lineages
    queryDF = queryDF.merge(linDF[["accession", "lineage"]], on="accession")
    rank_order= ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]
    queryDF['rank'] = pd.Categorical(queryDF['rank'], rank_order)
    queryDF.sort_values(by=["path", "rank"], inplace=True)
    # write full csv
    queryDF.to_csv(args.output_csv, index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("query_csv")
    p.add_argument("--lineages-csv")
    p.add_argument("--output-csv")
    args = p.parse_args()
    if not args.output_csv:
        args.output_csv=(args.query_csv).rsplit(".", 1)[0] + ".correct-lineages.csv"
    sys.exit(main(args))
