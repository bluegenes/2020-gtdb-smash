import os
import sys
import argparse
import glob
import pprint

import screed
import sourmash
from sourmash._minhash import hash_murmur
from sourmash.logging import notify

import sequence_encodings as enc

def try_reading_fasta_file(fasta_file):
    records=[]
    try:
        records=screed.open(fasta_file)
    except ValueError:
        sys.stderr.write(f"File {fasta_file} doesn't seem to be a sequence file, skipping. \n")
    return records

def translate_sequence(sequence, input_type, molecule):
    if input_type == "nucleotide":
        translator = enc.NUCLEOTIDE_ENCODINGS[molecule]
    elif input_type == "protein":
        translator = enc.PEPTIDE_ENCODINGS[molecule]
    return sequence.translate(translator)

def reencode_sequence(sequence, input_type, molecule):
    # modified from encode_peptide in sencha code
    assert input_type in ["nucleotide", "protein"]
    if input_type == "nucleotide":
        valid_encodings = enc.NUCLEOTIDE_ENCODINGS.keys()
    elif input_type == "protein":
        valid_encodings = enc.PEPTIDE_ENCODINGS.keys()
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

def kmers(seq, ksize, skip=None):
    # modified compute-dna-mh-another-way.py
    if skip:
        # see https://github.com/dib-lab/sourmash/pull/309/files
        assert len(skip) == 2 # skip should be a tuple
        # skip = (3, 2) # keep skip_m out of every skip_n
        skip_n, skip_m = skip
        # instead of asserting that it must be divisible, overshoot and then cut down?
        #assert ksize % skip_m == 0
        span = int(skip_n * (ksize / skip_m - 1) + skip_m) # rounds down
        # kmerize
        for start in range(len(seq) - span + 1):
            # get kmer of length span, which is skip_n * (ksize / skip_m - 1) + skip_m
            skip_idx = 0
            skipmer=""
            #while skip_idx < span:
            while skip_idx <= span: # = here allows you to overshoot, then cut down to desired ksize
                skipmer+=seq[start+skip_idx : start+skip_idx+skip_m]
                skip_idx+=skip_n
            # overshoot, then cut down. Any actual benefit to doing this?
            skipmer= skipmer[:ksize]
            #print(len(skipmer))
            #if you need to overshoot, then you will have with an improper kmer on the last one --- break
            if len(skipmer) < ksize:
                break
            yield skipmer
    else:
        # yield regular kmer
        for start in range(len(seq) - ksize + 1):
            yield seq[start:start + ksize]


def hash_sequence(seqstr, input_type, ksize, alphabet, skipinfo=None):
    hashes=[]
    # modify sequence if needed based on alphabet (e.g. protein --> dayhoff) # NOT nucl-> protein translation
    # hmm.. do this by kmer so we can revcomp in nucleotide space, then translate. Otherwise complement doesn't make sense
    #reencoded_seq = reencode_sequence(seqstr, input_type, alphabet)
    # check that we can kmerize?
    if len(seqstr) < ksize:
        return hashes
    for fwd_kmer in kmers(seqstr, ksize, skipinfo):
        if input_type == "nucleotide":
            # for nucleotide input, get reverse-complement, select smaller kmer
            rev_kmer = enc.reverse(enc.complement(fwd_kmer))
            if fwd_kmer < rev_kmer: # just a consistent way to choose a kmer, right?
                kmer = fwd_kmer
            else:
                kmer = rev_kmer
        else:
            # protein input, no need to revcomp
            kmer=fwd_kmer

        # tranlate, then hash
        translated_kmer = reencode_sequence(kmer, input_type, alphabet)
        #print(f"orig: {kmer}")
        #print(f"trans: {translated_kmer}")
        hash = hash_murmur(translated_kmer)
        if hash < 0:
            hash += 2**64
        hashes+=[hash]
        #yield hash
    return hashes


def hash_sequence_list(seqlist, input_type, ksize, alphabet, skipinfo=None):
    seqlist_hashes=[]
    for seq in seqlist:
        seqlist_hashes+=hash_sequence(seq, input_type, ksize, alphabet, skipinfo)
    return seqlist_hashes



def generate_sigs_from_file(input_file, input_type, alphabet, ksize, scaled, singleton=False, skipinfo=None, ignore_abundance=False): #, translate=False):
    sigs = []
    abund = not ignore_abundance

    # start with fresh minhash
    if input_type == "nucleotide":
        mh = sourmash.MinHash(ksize=ksize, n=0, scaled=scaled, track_abundance=abund) #, is_protein=False) # need is_protein??
    elif input_type == "protein":
        k=ksize*3 ## need to multiply by 3 to get same ksize, bc add_protein method does k/3
        mh = sourmash.MinHash(ksize=k, n=0, scaled=scaled, track_abundance=abund) #, is_protein=True) # need is_protein??

    # read file and add sigs
    records = try_reading_fasta_file(input_file)
    file_chunk = []
    if records:
        for n, record in enumerate(records):
            if singleton:
                hashes = hash_sequence(record.sequence, input_type, ksize, alphabet, skipinfo=skipinfo)
                # add_many takes care of scaled, abundance tracking, etc
                mh.add_many(hashes)
                sigs+=[sourmash.SourmashSignature(mh, name=record.name)]
                # reset mh
                mh = sourmash.MinHash(ksize=ksize, n=0, scaled=scaled, track_abundance=abund) #, is_protein=False) # need is_protein??
            else:
                if file_chunk and n % 10000 == 0:
                    hashes = hash_sequence_list(file_chunk, input_type, ksize, alphabet, skipinfo=skipinfo)
                    mh.add_many(hashes)
                    file_chunk=[]
                    #sys.stderr.write(f"... loading {filename} file {n} of {len(input_files)}\n")
                else:
                    file_chunk.append(record.sequence)
        if not singleton:
            if file_chunk:
                hashes = hash_sequence_list(file_chunk, input_type, ksize, alphabet, skipinfo=skipinfo)
                mh.add_many(hashes)
            # minhash --> signature, using filename as signature name ..i think this happens automatically if don't provide name?
            sigs+=[sourmash.SourmashSignature(mh, name=os.path.basename(input_file))]
    # return many sigs (singleton) or one sig per file
    return sigs


def generate_lossy_sigs(args):
    # define args

    # print avail output alphabets
    if args.print_output_alphabets:
        sys.stdout.write("\nNucleotide Options:\n\t")
        sys.stdout.write("\n\t".join(enc.NUCLEOTIDE_ENCODINGS.keys()) + "\n\n")
        sys.stdout.write("Protein Options\n\t")
        sys.stdout.write("\n\t".join(enc.PEPTIDE_ENCODINGS.keys()) + "\n\n")
        sys.exit(0)

    input_fasta = args.input_fasta
    assert os.path.exists(input_fasta)
    ksize = args.ksize
    input_type = args.input_alphabet
    alphabet = args.output_alphabet
    scaled=args.scaled
    singleton = args.singleton
    skip=None
    if args.skipmer:
        if input_type == "protein":
            sys.stdout.write("\nError: Please specify --input-alphabet as nucleotide or remove the ---skipmer option \n\n")
            sys.exit(0)
        #just set default for now
        skip=(3,2)

    # build signature(s) from file
    sigs = generate_sigs_from_file(input_fasta, input_type, alphabet, ksize, scaled, singleton, skip) #, args.ignore_abundance, args.translate)

    # build output file if necessary
    out_sig = args.signature_file
    if not out_sig:
        if singleton:
            out_sig = input_fasta.rsplit(".fa")[0] + f".{alphabet}_k{ksize}_scaled{scaled}_singleton.sig"
        else:
            out_sig = input_fasta.rsplit(".fa")[0] + f".{alphabet}_k{ksize}_scaled{scaled}.sig"
    # write sigs to file
    with open(out_sig, "w") as out:
        sourmash.signature.save_signatures(sigs, fp=out)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("input_fasta")
    p.add_argument("--signature-file")
    p.add_argument("--ksize", type=int, default=11)
    p.add_argument("--scaled", type=int, default=1)
    p.add_argument("--input-alphabet", default="protein", help="input type (nucleotide or protein)")
    p.add_argument("--output-alphabet", default="dayhoff", help="desired alphabet for sourmash signatures")
    p.add_argument("--print-output-alphabets", action="store_true")
    #p.add_argument("--ignore-abundance", action="store_true", help="do not store abundance information for signatures")
    p.add_argument("--skipmer", action="store_true", help="compute skipmer kmers with (m=2,n=3)")
    p.add_argument("--singleton", action="store_true", help="with fasta inputs, add one signature per fasta entry, rather than one per file")
    args = p.parse_args()
    sys.exit(generate_lossy_sigs(args))
