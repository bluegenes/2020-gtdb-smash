import os
import sys
import argparse
import glob
import pprint
import matplotlib.pyplot as plt
import pandas as pd

import screed
import sourmash


def load_sbt_file(tree_file):
    if os.path.exists(tree_file):
        try:
            sbt = sourmash.load_sbt_index(tree_file)
            return sbt
        except:
            sys.stderr.write(f"cannot load sbt file at {tree_file}\n")


def forage_by_name(sbt, query_accessions, threshold=0.1):
    acc2sig = {}
    num_accs= len(query_accessions)
    num_sigs=0
    for sig in sbt.signatures(): # generator with sig info
        #filename = sig.filename
        if sig.name() in query_accessions:
            #print(sig.minhash.get_mins())
            acc2sig[sig.name()] = sig
            num_sigs+=1
            if num_sigs == num_accs:
                break
    return acc2sig


def assess_group_distance(groupD, acc2sig): # provide options for choosing to use abund or not? (compute either one OR both, as is done now)
    for group, acclist in groupD.items():
        siglist = [acc2sig[acc] for acc in acclist] # get signatures
        compare_and_plot(siglist, ignore_abundance=False, filename= f"{group}_cosine") # need ksize, scaled, etc. can we pass in a base name, then add group name?
        compare_and_plot(siglist, ignore_abundance=True, filename= f"{group}_jaccard")
        # do we want to save distances as a csv for later? kinda depends on how long it takes to regenerate this...
        # do we want to compute all-by-all at the end? or just care about the groups?


def compare_and_plot(siglist, ignore_abundance=False, filename=None):
    dist_matrix = sourmash.compare.compare_all_pairs(siglist, ignore_abundance=ignore_abundance)
    plot_mat = sourmash.fig.load_matrix_and_labels(dist_matrix)
    dist_fig = sourmash.fig.plot_composite_matrix(plot_mat)
    if filename:
        plt.savefig(filename) #, format="svg")


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
            samples = pd.read_excel(groupss_file, dtype=str, sep=separator)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {groups_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    return samples



def build_sbt_accession_from_filename(groupInfo, db_infofile, filename=None):
    import pdb; pdb.set_trace()
    # if dataframe doesn't have filename
    if "filename" not in groupInfo.columns:
        dbInfo = csv_reader(db_infofile)
        # add filenames to groupInfo
        groupInfo=groupInfo.merge(dbInfo[["accession", "filename"]], on=["accession"])
    #generate sbt_accession from filenames
    samples["sbt_accession"] = samples["filename"].str.rsplit("_", 1, expand=True)[0]
    #sbt_accession = os.path.basename(fasta_file.rsplit("_", 1)[0])
    if filename:
        if filename.endswith(".csv"):
            groupInfo.to_csv(filename)
        elif filename.endswith(".tsv"):
            groupInfo.to_csv(filename, sep="\t")
        else:
            sys.stdout.write(f"Can't write updated groupinfo file. Filename {filename} must end with '.csv' or '.tsv'")
    return groupInfo


def build_group2accession(group_infofile, db_infofile="gtdb-lineages.csv"):
    groupDF = csv_reader(group_infofile)

    if "sbt_accession" not in groupDF.columns:
        updated_groupInfofile = group_infofile.rsplit(".", 1)[0] + ".sbtinfo.csv"
        groupDF = find_sbt_accessions(groupDF, db_infofile, updated_groupInfofile)

    # build group: accession dictionary
    group2acc = groupDf.groupby('path').agg({'sbt_accession':lambda x: list(x)})
    # build set of unique query accessions
    acc_set = set(groupDF["sbt_accession"].unique())

    return group2acc, acc_set


def main(args):
    sbt = load_sbt_file(args.sbt)
    query_accessionD, all_acc = build_group_to_accession(args.query_csv, args.database_csv)
    accession2sig = forage_by_name(sbt, all_acc, args.threshold)
    assess_group_distance(query_accessionD, accession2sig)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("sbt")
    p.add_argument("--query_csv", required=True)
    p.add_argument("--database_csv", required=True)
    p.add_argument("--threshold", type=int, default=0.1)
    args = p.parse_args()
    sys.exit(main(args))
