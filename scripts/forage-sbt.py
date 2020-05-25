import os
import sys
import argparse
import glob
import pprint
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import screed
import sourmash
# why fig doesn't work without this? best practice?
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


def load_sbt_file(tree_file):
    if os.path.exists(tree_file):
        try:
            sbt = sourmash.load_sbt_index(tree_file)
            sys.stderr.write(f"loaded sbt file at {tree_file}\n")
            return sbt
        except:
            sys.stderr.write(f"cannot load sbt file at {tree_file}\n")
            sys.exit()


def forage_by_name(sbt, query_signames, threshold=0.1):
    acc2sig = {}
    num_sigs=0
    num_accs= len(query_signames)
    for sig in sbt.signatures(): # generator with sig info
        if sig.name() in query_signames:
            #print(sig.minhash.get_mins())
            acc2sig[sig.name()] = sig
            num_sigs+=1
            if num_sigs == num_accs:
                break
    return acc2sig


def sourmash_compare_and_plot(siglist, ignore_abundance=False, plot_and_save=False, filename=None):
    dist_matrix = sourmash.compare.compare_all_pairs(siglist, ignore_abundance=ignore_abundance)
    labels = [item.name() for item in siglist]
    #optionally plot with sourmash plot
    if plot_and_save:
        dist_fig = sourmash.fig.plot_composite_matrix(dist_matrix, labels)
        if filename:
            plt.savefig(filename) #, format="svg")
    return dist_matrix, labels


def assess_group_distance(groupD, acc2sig, abund=False): # provide options for choosing to use abund or not? (compute either one OR both, as is done now)
    speciesDist = {}
    for group, accInfo in groupD.items():
        #siglist = [acc2sig[acc] for acc in accInfo.values()] # get signatures
        # take the time or order sigs properly in siglist
        steps_to_common_ancestor = {"species": 0, "genus": 1, "family": 2, "order": 3, "class": 4, "phylum": 5, "superkingdom": 6}
        siglist = [""]*7
        for rank, acc in accInfo.items():
            idx = steps_to_common_ancestor[rank]
            siglist[idx] = acc2sig[acc]
        jaccard_dist, labels=sourmash_compare_and_plot(siglist, ignore_abundance=True) #, filename= f"{group}_jaccard.svg") # need ksize, scaled, etc. can we pass in a base name, then add group name?
        if abund:
            cosine_dist, labels=sourmash_compare_and_plot(siglist, ignore_abundance=False)  #, filename= f"{group}_cosine.svg")
        # build distance:: steps_to_common_ancestor!
        # all we care about is the distance from species-level? so every pairwise with that genome.
        # find the accession that is the species entry
        # grab column (or row) from np matrix / dataframe

        #given proper ordering above, we can assume jaccard_dist is ordered species::superkingdom
        # to do: riple check this.
        jaccard_from_species = jaccard_dist[0, :]
        speciesDist[group] = jaccard_from_species

    return speciesDist


def plot_all_distances(speciesDist, dist_csv, dist_plot=None):
    # given all the distances generated from each path, boxplot the distances!
    # groupDF['steps_to_common_ancestor'] = groupDF["rank"].map(steps_to_common_ancestor)
    # for now, just turn into pandas DF and box plot
    # to do: set common rank_order / steps to common_ancestor once, earlier, so don't have discordance btwn assess_group_distance and here

    rank_order= ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]
    # this imports with path as the rownames.
    distDF = pd.DataFrame.from_dict(speciesDist, orient="index", columns= rank_order)
    # species info is always 100% and therefore useless/compressing plot info. try dropping!
    distDF.drop(['species',], axis=1, inplace=True)
    with open(dist_csv, "w") as out:
        distDF.to_csv(dist_csv,index=False)
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


def main(args):
    # from query csv, build dictionary of group:: filenames (signature names)
    groupDF = try_reading_csv(args.query_csv)
    signame_col = args.signature_name_column
    if signame_col not in groupDF.columns:
        sys.stderr.write(f"your query csv does not have the {signame_col} column. " \
        "we need this column to find corresponding signatures in the sbt.\n")
        sys.exit()

    group2signame = (groupDF.groupby('path').apply(lambda x: dict(zip(x['rank'],x[signame_col]))).to_dict())
    signames_to_find = set(groupDF[signame_col].unique())

    sbt = load_sbt_file(args.sbt)
    # find signatures in sbt
    signame2sig = forage_by_name(sbt, signames_to_find)#, args.threshold)
    # assess and plot distances
    dist_from_species_level = assess_group_distance(group2signame, signame2sig)
    plot_all_distances(dist_from_species_level, args.distance_from_species_csv, args.distance_from_species_plot)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("sbt")
    p.add_argument("--query-csv", required=True)
    p.add_argument("--signature-name-column", default="filename", help="column with signature names in the sbt. By default, this should be the fasta file basename")
    #p.add_argument("--threshold", type=int, default=0.1)
    p.add_argument("--distance-from-species-csv", default="path_species_jaccard_dists.csv")
    p.add_argument("--distance-from-species-plot", default="testpaths.boxplot.svg")
    args = p.parse_args()
    sys.exit(main(args))
