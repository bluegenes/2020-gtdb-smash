# first, gather all '.fna.gz'
import os
import random

configfile: "config.yml"

prefix = config.get('prefix', '')
genomes_location = config['genomes_location']
sigs_output_location = config['sigs_output_location']

signature_scaled = int(config['signature_scaled'])
signature_ksizes = [ int(x) for x in config['signature_ksizes'] ]
lca_scaled = [ int(x) for x in config['lca_scaled'] ]
lca_ksizes = [ int(x) for x in config['lca_ksizes'] ]

if min(lca_scaled) < signature_scaled:
    raise ValueError("smallest lca_scaled needs to be >= signature_scaled")

for ksize in lca_ksizes:
    if ksize not in signature_ksizes:
        raise ValueError("signature_ksizes must include lca ksize {}".format(ksize))

###

all_files = []
for root, dirs, files in os.walk(genomes_location, topdown=False):
    for name in files:
        if name.endswith('_protein.faa.gz'):
            filename = os.path.join(root, name)
            if filename.startswith('./'): filename = filename[2:]
            all_files.append(filename)

print('found {} files under {}'.format(len(all_files), genomes_location))
random.shuffle(all_files)
print('examples:', all_files[:5])

def make_sigfilename(filename):
    genomefile = os.path.basename(filename)
    sigfile = os.path.join(sigs_output_location, genomefile) + '.sig'

    return sigfile

# rule all: make sigs
rule all:
    input:
        [ make_sigfilename(genomefile) for genomefile in all_files ],
        expand(prefix+"gtdb-release89-k{k}-{scaled}.lca.json.gz",
               k=lca_ksizes, scaled=lca_scaled),
        expand(prefix+"gtdb-release89-k{k}.sbt.zip", k=31)
#        prefix+"gtdb-release89-k21-lowrank.lca.json.gz",
#        prefix+"gtdb-release89-k31-lowrank.lca.json.gz",
#        prefix+"gtdb-release89-k51-lowrank.lca.json.gz",
#        prefix+"gtdb-release89-k21.sbt.json",
#        prefix+"gtdb-release89-k31.sbt.json",
#        prefix+"gtdb-release89-k51.sbt.json"

rule make_lineage_csv:
    input:
        config['gtdb_taxonomy_tsv']
    output:
        'gtdb-lineages.csv'
    conda: 'env-sourmash.yml'
    shell:
        '../update-gtdb-taxonomy.py {input} -o {output}'

rule make_lca_db:
    input:
        "gtdb-lineages.csv",
        [ make_sigfilename(genomefile) for genomefile in all_files ]
    output:
        prefix+"gtdb-release89-k{ksize,\d+}-{scaled,\d+}.lca.json.gz"
    params:
        ksize="{ksize}",
        prefix=prefix,
        scaled="{scaled}",
    conda: 'env-sourmash.yml'
    shell:
        """sourmash lca index -k {params.ksize} \
                           --scaled {params.scaled} --require-taxonomy \
                           --split-identifiers \
                           --dayhoff \
                           --report report-gtdb-k31.txt \
                           --traverse-directory -C 3 \
                           gtdb-lineages.csv \
                           {output} \
                           {sigs_output_location}
        """

rule make_sbt_db:
    input:
        [ make_sigfilename(genomefile) for genomefile in all_files ]
    output:
        prefix+"gtdb-release89-k{ksize}.sbt.zip"
    params:
        ksize="{ksize}",
        prefix=prefix
    conda: 'env-sourmash.yml'
    shell:
        """sourmash index -k {params.ksize} \
                           --traverse-directory \
                           {output} {sigs_output_location}
        """

filename_to_names = None
def lookup_name(filename):
    global filename_to_names
    print('XXX', filename)

    if filename_to_names is None:
        import csv
        d = {}
        with open('gtdb-lineages.csv', 'rt') as fp:
            r = csv.DictReader(fp)
            for row in r:
                fn = row['filename']
                name = "{} {}".format(row['accession'], row['species'])
                d[fn] = name
        filename_to_names = d

    key = os.path.basename(filename)
    print('looking up', key, filename_to_names[key])
    return filename_to_names[key]

# generic rule: compute signature
rule compute_sig:
    input:
        os.path.join(genomes_location, "{filename}"),
    output:
        os.path.join(sigs_output_location + "{filename}.sig")
    params:
        signame = lambda wildcards: lookup_name(wildcards.filename),
        ksizes = ",".join([ str(x) for x in signature_ksizes ]),
        scaled = signature_scaled,
    conda: 'env-sourmash.yml'
    shell: """
        sourmash compute -k {params.ksizes} --scaled={params.scaled} \
             {input} -o {output} --dayhoff --no-dna --merge={params.signame:q}
    """

# extract high-rank hash values for a particular ksize
rule extract_hash_values:
    input:
        prefix+"gtdb-release89-k{ksize}.lca.json.gz"
    output:
        prefix+"high-rank-hashes.k{ksize}.hashvals"
    conda: 'env-sourmash.yml'
    shell:
        "~/2019-sourmash-gtdb/extract-high-rank-hashes.py {input} --lowest-rank=phylum -o {output}"

# scrub LCA database of high-rank hash values
rule scrub_lca_db:
    input:
        prefix+"gtdb-release89-k{ksize}.lca.json.gz",
        prefix+"high-rank-hashes.k{ksize}.hashvals"
    output:
        prefix+"gtdb-release89-k{ksize}-lowrank.lca.json.gz"
    conda: 'env-sourmash.yml'
    shell:
        "~/2019-sourmash-gtdb/scrub-lca-db.py {input[0]} {input[1]} -o {output}"
