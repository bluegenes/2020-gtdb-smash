import os
import sys
import argparse
import glob
import pprint

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
    path_siglist=[]
    num_accs= len(query_accessions)
    num_sigs=0
    for sig in sbt.signatures(): # generator with sig info
        #filename = sig.filename
        if sig.name() in query_accessions:
            #print(sig.minhash.get_mins())
            path_siglist.append(sig)
            num_sigs+=1
            if num_sigs == num_accs:
                break
    return path_siglist

def main(args):
    # create or load sbt
    sbt = load_sbt_file(args.sbt)
    accession_info = forage_by_name(sbt, args.query_accessions, args.threshold)
    print(accession_info)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("sbt")
    p.add_argument("--query_accessions", nargs='+', required=True)
    p.add_argument("--threshold", type=int, default=0.1)
    args = p.parse_args()
    sys.exit(main(args))
