


python generate-lossy-signature.py ../testdata/GCA_003153015.1_protein.one.fna --ksize 31 --input-alphabet nucleotide --output-alphabet amino-keto --scaled  1 --skipmer
python generate-lossy-signature.py ../testdata/GCA_003153015.1_protein.one.fna --ksize 31 --input-alphabet nucleotide --output-alphabet purine-pyrimidine --scaled  1


python generate-lossy-signature.py ../testdata/GB_GCA_002727195.1_protein.one.faa --ksize 19 --input-alphabet protein --output-alphabet dayhoff --scaled  1
python generate-lossy-signature.py ../testdata/GB_GCA_002727195.1_protein.head200.faa --ksize 19 --input-alphabet protein --output-alphabet dayhoff --scaled 1
python generate-lossy-signature.py ../testdata/GB_GCA_002727195.1_protein.head200.faa --ksize 19 --input-alphabet protein --output-alphabet aa9 --scaled 1


python generate-lossy-signature.py ../testdata/GB_GCA_002727195.1_protein.head200.faa --ksize 19 --input-alphabet protein --output-alphabet dayhoff --scaled 1 --singleton

