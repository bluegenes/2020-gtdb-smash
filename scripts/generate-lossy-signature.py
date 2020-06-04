import os
import sys
import argparse
import glob
import pprint

import screed
import sourmash
from sourmash._minhash import hash_murmur
from sourmash.logging import notify

import sequence_encodings

def try_reading_fasta_file(fasta_file):
    records=[]
    try:
        records=screed.open(fasta_file)
    except ValueError:
        sys.stderr.write(f"File {fasta_file} doesn't seem to be a sequence file, skipping. \n")
    return records

def translate_sequence(sequence, input_type, molecule):
    if input_type == "nucleotide":
        translator = NUCLEOTIDE_ENCODINGS[molecule]
    elif input_type == "protein":
        translator = PEPTIDE_ENCODINGS[molecule]
    return sequence.translate(translator)

def reencode_sequence(sequence, input_type, molecule):
    # modified from encode_peptide in sencha code
    assert input_type in ["nucleotide", "protein"]
    if input_type == "nucleotide":
        valid_encodings = NUCLEOTIDE_ENCODINGS.keys()
    elif input_type == "protein":
        valid_encodings = PEPTIDE_ENCODINGS.keys()
    if molecule in valid_encodings:
        return translate_sequence(sequence, input_type, molecule)
    elif molecule == input_type:
        return sequence
    else:
        raise ValueError(
            f"{molecule} is not a valid {input_type} encoding, "
            "only "
            f"{', '.join(valid_encodings)} can be used"
        )


def kmers(seqstr, ksize, skip=None):
    # modified from compute-dna-mh-another-way.py
    if skip:
        assert len(skip) == 2 # skip should be a tuple
        # skip = (2, 3) # keep skip_m out of every skip_n
        skip_m, skip_n = skip
        span = int(skip_n * (ksize / skip_m - 1) + skip_m) # rounds down
        # kmerize
        for start in range(len(seqstr) - span + 1):
            # get kmer of length span, which is skip_n * (ksize / skip_m - 1) + skip_m
            skip_idx = 0
            while skip_idx < span:
                skipmer+=seqstr[start+skip_idx : start+skip_idx+skip_m]
                skip_idx+=skip_n
            yield skipmer
    else:
        for start in range(len(seq) - ksize + 1):
            yield seq[start:start + k]


def hash_sequence(seqstr, input_type, ksize, skipinfo=None): #, ignore_abundance=True):
    hashes=[]
    # modify sequence if needed based on alphabet (e.g. protein --> dayhoff) # NOT nucl-> protein translation
    reencoded_seq = reencode_sequence(sequence, input_type, alphabet)
    # check that we can kmerize?
    #if len(seqstr) < ksize:
    #    return []
    for fwd_kmer in kmers(seqstr, ksize, skipinfo):
        if input_type == "nucleotide":
            rev_kmer = reverse(complement(fwd_kmer))
            if fwd_kmer < rev_kmer:
                kmer = fwd_kmer
            else:
                kmer = rev_kmer
            hash = hash_murmur(kmer)
        else:
            hash = hash_murmur(fwd_kmer)
        if hash < 0:
            hash += 2**64
        hashes+=[hash] # use list for now
        #yield hash
    return hashes

    ### how to keep track of abundance with this method? eh, ignore for now
    #mh = sourmash.Minhash(ksize=k)
    #if ignore_abundance:
    #    for kmer in kmers:
    #        mh.add_hash(hash_murmur(kmer))
    #else:
    #    for kmer, abund in kmers:
            #mins (default None) - list of hashvals, or (hashval, abund) pairs
    #        mh.add_hash((hash_murmur(kmer), abund))


def generate_sigs_from_file(input_file, input_type, alphabet, ksize, scaled, singleton=False, skipinfo=None): #, ignore_abundance=False): #, translate=False):
    sigs = []
    # start with fresh minhash
    mh = sourmash.Minhash(ksize=ksize)

    # read file and add sigs
    records = try_reading_fasta_file(input_file)
    if records:
        for record in records:
            mh.add_many(hash_sequence(record.sequence, input_type, ksize, skipmer=skipinfo))
            if singleton:
                # how to set additional info (alphabet, etc?)
                sig+=[sourmash.SourmashSignature(mh, name=record.name)]
                # reset mh
                mh = sourmash.Minhash(ksize=ksize)
        if not singleton:
            # minhash --> signature, using filename as signature name ..i think this happens automatically if don't provide name?
            sigs+=[sourmash.SourmashSignature(mh, name=os.path.basename(input_file))]

    # does scaled happen after generating all hashes? how to implement scaled?
    return sigs


def generate_lossy_sigs(args):
    # define args
    input_fasta = args.input_fasta
    assert os.path.exists(input_fasta)
    ksize = args.ksize
    input_type = args.input_alphabet
    alphabet = args.output_alphabet
    scaled=args.scaled
    singleton = args.singleton
    if args.skipmer:
        #just set default for now
        skip=(3,2)

    # build signature(s) from file
    sigs = generate_sigs_from_file(input_fasta, input_type, alphabet, ksize, scaled, singleton, skip) #, args.ignore_abundance, args.translate)

    # build output file if necessary
    out_sig = args.sinature_file
    if not out_sig:
        out_sig = input_fasta.rsplit(".fa")[0] +f".{alphabet}_k{ksize}_scaled{scaled}.sig"
    # write sigs to file
    with open out_sig as out:
        sourmash.signature.save_signatures(sigs, fp=out)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("input_fasta")
    p.add_argument("--signature-file")
    p.add_argument("--ksize", type=int, default=31)
    p.add_argument("--scaled", type=int, default=1000)
    p.add_argument("--input-alphabet", default="protein", help="input type (nucleotide or protein)")
    p.add_argument("--output-alphabet", default="dayhoff", help="desired alphabet for sourmash signatures. options: nucleotide, protein, dayhoff, or hp")
    #p.add_argument("--ignore-abundance", action="store_true", help="do not store abundance information for signatures")
    p.add_argument("--skipmer", action="store_true", help="compute skipmer kmers with (m=2,n=3)")
    p.add_argument("--singleton", action="store_true", help="with fasta inputs, add one signature per fasta entry, rather than one per file")
    args = p.parse_args()
    sys.exit(generate_lossy_sigs(args))
