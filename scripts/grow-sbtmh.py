import os
import sys
import argparse
import glob
import pprint

import screed
import sourmash
from sourmash.sbtmh import SigLeaf


def maybe_read_fasta_file(fasta_file):
    records=[]
    accession=""
    try:
        records=screed.open(fasta_file)
         # build a short name for this signature using the filename
        accession = os.path.basename(fasta_file.rsplit("_", 1)[0])
    except:
        # log the error somehow
        sys.stderr.write(f"cannot read fasta file {fasta_file}\n")
    return records, accession


def maybe_load_sbt_file(tree_file):
    sbt=""
    if os.path.exists(tree_file):
        try:
            sbt = sourmash.load_sbt_index(tree_file)
        except:
            sys.stderr.write(f"cannot load sbt file at {tree_file}\n")
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


def grow_sbt(input_files, sbt_file, ksize, scaled, alpha, abund, input_is_directory):
    if input_is_directory:
        refdir = input_files[0]
        input_files = [os.path.join(refdir, f) for f in os.listdir(refdir)]
    # create or load sbt
    sbt = maybe_load_sbt_file(sbt_file)
    # iterate through input files; add to sbt
    md5sum_set=set()
    for inp in input_files:
        #sys.stderr.write(f"{inp}\n")
        sig = load_or_generate_sig(inp, ksize, scaled, alpha, abund)
        if sig:
            # protein and hp were exiting with ValueError: This will insert duplicated entries .. maybe checking for empty mins will solve?
            if len(sig.minhash.get_mins()) > 0:
                md5 = sig.md5sum()
                if md5 not in md5sum_set:
                    leaf = SigLeaf(md5, sig)
                    sbt.add_node(leaf)
                    md5sum_set.add(md5)
                else:
                    sys.stderr.write(f"duplicated md5sum {md5}!")
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
    p.add_argument("--track-abundance", action="store_true")
    args = p.parse_args()
    sys.exit(grow_sbt(args.input_files, args.sbt, args.ksize, args.scaled, args.alphabet, args.track_abundance, args.input_is_directory))
