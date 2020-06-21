import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd

#from sourmash.lca.command_index import load_taxonomy_assignments

def split_col_w_potential_nan(row):
    try:
        info = row["match_at_rank"].split(" ")
        row["percent_match"] = info[0]
        row["lineage"] = info[1]
    except:
        row["percent_match"] = 0
        row["lineage"] = "no match"
    return row

rank_to_lineage_map = {"superkingdom": 0,"phylum": 1,"class": 2,"order": 3,"family": 4,"genus": 5,"species":6}

def assess_gather_lineage_at_rank(row):
    match_rank = row["rank"]
    match_rank_index = rank_to_lineage_map[match_rank]
    true_lineage = row["true_lineage"]
    true_lineage_at_match_rank = true_lineage[match_rank_index]
    gather_lineage_at_match_rank = row["lineage"]
    if gather_lineage_at_match_rank == "no match":
        row["lineage_match"] = "no gather match"
    elif true_lineage_at_match_rank in gather_lineage_at_match_rank:
        row["lineage_match"] = "correct lineage"
    else:
        row["lineage_match"] = "incorrect lineage"
    return row


def collect_input_files_from_dir(input_dir, subset_csv=None, subset_info_colname="accession"):
    """
    Scan a directory for input fasta or signature files.
    If subset_csv is provided, glob the directory for the pattern provided in the
    subset_info_colname column.
    """
    input_files=[]
    if subset_csv:
        try:
            match_strings = set(pd.read_csv(subset_csv)[subset_info_colname].tolist())
            for m in match_strings:
                input_files+=glob.glob(os.path.join(input_dir, f"*{m}*.gather_tophits.csv"))
        except:
            sys.stderr.write(f"can't collect input files from dir {input_dir} using subset_csv {subset_csv} column {subset_info_colname}")
            sys.exit()
    else:
        input_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if ".gather_tophits.csv" in f]
    if not input_files:
        sys.stderr.write(f"can't collect input files from dir {input_dir}")
        sys.exit()
    return input_files



def main(args):
    #read in csv input files
    input_files = args.input_files
    if args.input_is_directory:
        input_files = collect_input_files_from_dir(input_files[0], args.true_lineages_csv)

    genome2gathertophits = pd.concat([pd.read_csv(csv).assign(accession=os.path.basename(csv).rsplit("_x_")[0]) for csv in input_files])
    column_order = ["accession","superkingdom","phylum","class","order","family","genus","species"]
    genome2gathertophits = genome2gathertophits[column_order]
    # ok, we should melt + split this dataframe
    # want: accession, rank, percent match, lineage

    # first melt
    gatherMelt = pd.melt(genome2gathertophits, id_vars=["accession"], value_vars=["superkingdom","phylum","class","order","family","genus","species"], var_name='rank', value_name = "match_at_rank")

    # then split each rank on " ". sigh, watch out for nans. Maybe just drop?
    gatherMelt = gatherMelt.apply(split_col_w_potential_nan, axis=1)
    gatherMelt.drop("match_at_rank", axis=1, inplace=True)

    if args.true_lineages_csv:
        # cols will be: accession,superkingdom,phylum,class,order,family,genus,species,strain,signame,filename
        #tax_assign, _ = load_taxonomy_assignments(args.true_lineages_csv,start_column=args.start_column)
        #print(f'loaded {len(tax_assign)} tax assignments.')
        true_lineages= pd.read_csv(args.true_lineages_csv)
        true_lineages["true_lineage"] = true_lineages[["superkingdom","phylum","class","order","family","genus","species"]].values.tolist()
        # first, add true lineage to gatherMelt
        gatherMelt = gatherMelt.merge(true_lineages[["accession", "true_lineage"]], on = "accession")
        gatherMelt.dropna(inplace=True)
        # now assess gather == true lineage?
        gatherMelt = gatherMelt.apply(assess_gather_lineage_at_rank, axis=1)
        # count true/false in lineage_match column
        matches = gatherMelt.groupby("rank")["lineage_match"].value_counts().reset_index(name="num_matches")
        matches_pivot = matches.pivot(index='rank', columns='lineage_match')
        matches_pivot.columns = matches_pivot.columns.droplevel()
        #matches = gatherMelt.groupby("rank")["lineage_match"].agg() #.reset_index()
        #gatherMelt["match_freq"] = gatherMelt.groupby("rank")["lineage_match"].transform("count")

        gatherMelt = gatherMelt.merge(matches_pivot, left_on = "rank", right_index=True)
        gatherMelt = gatherMelt.rename(columns={"correct lineage": "total_correct_lineages", "incorrect lineage": "total incorrect lineages", "no gather match": "total no match"})


        # print result
        matches_pivot.to_csv(args.output_match_info, index=False)

        outDF = gatherMelt
    else:
        outDF = genome2gathertophits
    outDF.to_csv(args.output_csv, index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("input_files", nargs='+')
    p.add_argument("--output-csv", required=True)
    p.add_argument("--output-match-info")
    p.add_argument("--true-lineages-csv")
    p.add_argument("--start-column", type=int, default=2)
    p.add_argument("--input-is-directory", action="store_true")
    args = p.parse_args()
    if not args.output_match_info:
        args.output_match_info = args.output_csv.rsplit(".csv")[0] + ".matchinfo.csv"
    sys.exit(main(args))
