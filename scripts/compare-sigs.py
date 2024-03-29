import os
import sys
import argparse
import glob
import pprint
import matplotlib.pyplot as plt
#import pandas as pd
import seaborn as sns

import screed
import sourmash
from sourmash import sourmash_args #import load_file_as_signatures
import sourmash.fig

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

def get_anchor_containment(acclist, siglist):
    # return accession, anchor containment as list of tuples
    contain_tuples = []
    anchor_sig = siglist[0]
    for acc, sig in zip (acclist, siglist):
        containment = anchor_sig.contained_by(sig)
        contain_tuples.append((acc, containment))
    return contain_tuples

def get_anchor_jaccard(acclist, siglist):
    # return accession, anchor containment as list of tuples
    jaccard_tuples = []
    anchor_sig = siglist[0]
    for acc, sig in zip (acclist, siglist):
        jaccard_dist = anchor_sig.jaccard(sig)
        jaccard_tuples.append((acc, jaccard_dist))
    return jaccard_tuples

def assess_jaccard_and_containment(groupD, acc2sig, abund=False): # provide options for choosing to use abund or not? (compute either one OR both, as is done now)
    anchorJaccard,anchorContainment = {},{}
    for group, accInfo in groupD.items():
        steps_to_common_ancestor = {"species": 0, "genus": 1, "family": 2, "order": 3, "class": 4, "phylum": 5, "superkingdom": 6}
        siglist = [""]*7
        acclist = [""]*7
        for rank, acc in accInfo.items():
            idx = steps_to_common_ancestor[rank]
            siglist[idx] = acc2sig[acc]
            acclist[idx] = acc
        # assess jaccard and containment from anchor
        anchorContainment[group] = get_anchor_containment(acclist,siglist)
        anchorJaccard[group] = get_anchor_jaccard(acclist,siglist)

    return anchorContainment, anchorJaccard

## TODO: modify plotting for new format

def plot_all_distances(speciesDist, dist_csv, dist_plot=None):
    # given all the distances generated from each path, boxplot the distances!
    # groupDF['steps_to_common_ancestor'] = groupDF["rank"].map(steps_to_common_ancestor)
    # for now, just turn into pandas DF and box plot
    # to do: set common rank_order / steps to common_ancestor once, earlier, so don't have discordance btwn assess_group_distance and here
    # to do: CLEAN UP!

    # first, build distance dataframe as before (to keep plotting the same for now)
    distance_only = {}
    for key, val in speciesDist.items():
        distance_only[key] = [v[1] for v in val]
    #dist_only = {k,v[:][0] for k, v in speciesDist.items()}
    rank_order= ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]
    # this imports with path as the rownames.
    distDF = pd.DataFrame.from_dict(distance_only, orient="index", columns= rank_order)

    if dist_plot:
        sns.set_style("white")
        plt.figure(figsize=(11,7))
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=12)
        sns.boxplot(data=distDF, palette="GnBu_d")
        sns.stripplot(data=distDF,color='black',alpha=0.7) #, size=8)
        plt.xlabel("Common Ancestor", size=18, labelpad=20)
        plt.ylabel("Jaccard", size=20, labelpad=15)
        plt.savefig(dist_plot)

    ## now, build long form data
    distanceTest = pd.DataFrame.from_dict(speciesDist, orient="index", columns= rank_order)
    distanceTest = distanceTest.reset_index() # make index (path name) a column
    distanceTest = distanceTest.rename(columns={"index": "path"})
    distanceMelt = pd.melt(distanceTest, id_vars=["path"], value_vars=["species", "genus", "family", "order", "class", "phylum", "superkingdom"], var_name='rank')
    distanceMelt[['accession', 'jaccard']] = pd.DataFrame(distanceMelt['value'].tolist(), index=distanceMelt.index)
    distanceMelt.drop("value", axis=1, inplace=True)
    distanceMelt['rank'] = pd.Categorical(distanceMelt['rank'], rank_order)
    distanceMelt.sort_values(by=["path","rank"], inplace=True)
    with open(dist_csv, "w") as out:
        distanceMelt.to_csv(dist_csv, index=False)
    sys.stderr.write(f"wrote distance csv to {dist_csv}\n")


def main(args):
    # from query csv, build dictionary of group:: filenames (signature names)
    groupDF = try_reading_csv(args.paths_csv)
    signame_col = args.signature_name_column
    if signame_col not in groupDF.columns:
        sys.stderr.write(f"your paths csv does not have the {signame_col} column. " \
        "we need this column to find corresponding signatures.\n")
        sys.exit()

    group2signame = (groupDF.groupby('path').apply(lambda x: dict(zip(x['rank'],x[signame_col]))).to_dict())
    signames_to_find = set(groupDF[signame_col].unique())

    # find signatures in sbt
    #signame2sig = forage_by_name(sbt, signames_to_find)#, args.threshold)
    # INSTEAD, load signames from file
    all_sigs = sourmash.load_file_as_signatures(args.from_file)

    # assess and plot distances
    dist_from_species_level = assess_group_distance(group2signame, signame2sig, args.full_compare_csv, args.full_compare_plot)
    #plot_all_distances(dist_from_species_level, args.distance_from_species_csv, args.distance_from_species_plot)
    plot_all_distances(dist_from_species_level, args.distance_from_species_csv, args.distance_from_species_plot)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("paths_csv")
    p.add_argument("--from-file", required=True)
    p.add_argument("--signature-name-column", default="filename", help="column with signature names in the sbt. By default, this should be the fasta file basename")
    p.add_argument("--from-anchor-csv", default="path_species_jaccard_dists.csv")
    p.add_argument("--from-anchor-plot", default="path_species.boxplot.svg")
    args = p.parse_args()
    sys.exit(main(args))
