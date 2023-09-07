import os
import sys
import argparse
import pandas as pd


def main(args):
    #read in csv input files
    genome2gathertophits = pd.concat([pd.read_csv(csv).assign(accession=os.path.basename(csv).rsplit("_x_")[0]) for csv in args.input_files])
    column_order = ["accession","superkingdom","phylum","class","order","family","genus","species"]
    genome2gathertophits = genome2gathertophits[column_order]

    tophits = pd.read_csv(args.gather_tophits_csv)
    lineages= pd.read_csv(args.lineages_csv)

    # use lineages file to add true lineages to tophits
    # for each genome,

    tophits_assessment.to_csv(args.output_csv, index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--gather-tophits-csv", required=True)
    p.add_argument("--lineages-csv", required=True)
    p.add_argument("--output-csv", required=True)
    args = p.parse_args()
    sys.exit(main(args))
