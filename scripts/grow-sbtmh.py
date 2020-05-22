import os
import sys
import argparse
import glob
import pprint

import screed
import sourmash
from sourmash.sbtmh import SigLeaf
import pandas as pd

def csv_reader(groups_file):
    # autodetect format
    samples=""
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


def maybe_read_fasta_file(fasta_file):
    records=[]
    accession=""
    try:
        records=screed.open(fasta_file)
         # build a short name for this signature using the filename
        accession = os.path.basename(fasta_file.rsplit("_", 1)[0])
    except:
        # log the error somehow
        sys.stderr.write(f"\ncannot read fasta file {fasta_file}\n")
    return records, accession


def maybe_load_sbt_file(tree_file):
    sbt=""
    if os.path.exists(tree_file):
        try:
            sbt = sourmash.load_sbt_index(tree_file)
        except:
            sys.stderr.write(f"\ncannot load sbt file at {tree_file}\n")
    else:
        sbt = sourmash.create_sbt_index()
    return sbt


def determine_appropriate_fresh_minhash(ksize, scaled_val, abund, alphabet):
    if alphabet == "dna":
        mh = sourmash.MinHash(ksize=ksize, n=0, scaled=scaled_val, track_abundance=abund, is_protein=False)
    elif alphabet == "protein":
        mh = sourmash.MinHash(ksize=ksize, n=0, scaled=scaled_val, track_abundance=abund, is_protein=True, dayhoff=False, hp=False)
    elif alphabet == "dayhoff":
        mh = sourmash.MinHash(ksize=ksize, n=0, scaled=scaled_val, track_abundance=abund, is_protein=True, dayhoff=True, hp=False)
    elif alphabet == "hp":
        mh = sourmash.MinHash(ksize=ksize, n=0, scaled=scaled_val, track_abundance=abund, is_protein=True, dayhoff=False, hp=True)
    return mh


def load_or_generate_sig(input_file, ksize, scaled, alphabet, abundance):
    sig=""
    if input_file.endswith(".sig"):
        sig = sourmash.load_one_signature(input_file, ksize=ksize)
    else:
        # read file and add sigs
        records, accession = maybe_read_fasta_file(input_file)
        # start with fresh minhash
        mh = determine_appropriate_fresh_minhash(ksize, scaled, abundance, alphabet)
        if records:
            for record in records:
                if alphabet == "dna":
                    mh.add_sequence(record.sequence)
                else:
                    #if translate: ... need to do anything differently?
                    mh.add_protein(record.sequence)
            # minhash --> signature
            sig = sourmash.SourmashSignature(mh, name=accession)
    return sig


def grow_sbt(input_files, sbt_file, ksize, scaled, alpha, abund, input_is_directory, dup_sigs, sig2file, csv_subset):
    sig2filename={}
    duplicated_md5_sigs=set()
    md5sum_set=set()
    if input_is_directory:
        refdir = input_files[0]
        if csv_subset:
            accs = (csv_reader(csv_subset))["accession"].tolist()
            input_files=[]
            for acc in accs:
                input_files+=glob.glob(os.path.join(refdir, f"*{acc}*"))
        else:
            input_files = [os.path.join(refdir, f) for f in os.listdir(refdir)]
    # create or load sbt
    sbt = maybe_load_sbt_file(sbt_file)
    # iterate through input files; add to sbt
    for n, filename in enumerate(input_files):
        # swipe some handy progress reporting code from titus:
        if n % 100 == 0:
            sys.stderr.write(f"... loading {filename} file {n} of {len(input_files)}")
        sig = load_or_generate_sig(filename, ksize, scaled, alpha, abund)
        if sig: # is this necessary?
            if sig.minhash:
                sig2filename[sig.name()]=filename
                md5 = sig.md5sum()
                if md5 not in md5sum_set:
                    leaf = SigLeaf(md5, sig)
                    sbt.add_node(leaf)
                    md5sum_set.add(md5)
                else:
                    # this only records signatures that do not get added.
                    # (original sig with that md5 is in the sbt)
                    duplicated_md5_sigs.add(sig.name())  # don't really need correspondence, do we? just keep set.
                    sys.stderr.write(f"duplicated md5sum {sig.name()}!\n")
                    # this is kinda dumb. don't really want to track this if we can help it.
                    sig2filename[sig.name()]=filename

    # keep namemap, becuase the filenames from gtdb-lineage-csv don't seem to be what we need
    # name and filename are kinda redundant here. decide what we want to keep.
    with open(sig2file, "w") as out:
        out.write(",".join(["sbt_accession", "filename"]) + "\n")
        for name, filename in sig2filename.items():
            out.write(",".join([name, filename]) + "\n")
    with open(dup_sigs, "w") as out:
        for name in duplicated_md5_sigs:
            out.write(",".join([name, sig2filename[name]]) + "\n")
        #out.write("\n".join(duplicated_md5_sigs))
    # save the tree
    # hmm.. overwriting seems complicated but managed? see sourmash/sbt_storage.py
    sbt.save(sbt_file)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("input_files", nargs='+')
    p.add_argument("--sbt")
    p.add_argument("--ksize", type=int, default=31)
    p.add_argument("--scaled", type=int, default=1000)
    p.add_argument("--alphabet", default="dna")
    p.add_argument("--input-is-directory", action="store_true")
    p.add_argument("--subset-csv", default=None)
    p.add_argument("--track-abundance", action="store_true")
    args = p.parse_args()
    if not args.sbt.endswith(".sbt.zip"):
        sys.stderr.write("sbt file must end with .sbt.zip")
        sys.exit()
    else:
        dupes= args.sbt.rsplit(".sbt.zip")[0] + ".md5_duplicates.csv"
        print(dupes)
        signame2filename= args.sbt.rsplit(".sbt.zip")[0] + ".signame2filenames.csv"
        print(signame2filename)

    sys.exit(grow_sbt(args.input_files, args.sbt, args.ksize, args.scaled, args.alphabet, args.track_abundance, args.input_is_directory, dupes, signame2filename, args.subset_csv))
