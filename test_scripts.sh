



#python scripts/grow-sbtmh.py test_grow_sbt --sbt test_new_growsbt/sbt-7paths.dayhoff_k13_scaled1.sbt.zip --ksize 13 --input-is-directory --scaled 1 --alphabet dayhoff --subset-info-colname accession --subset-csv gtdb_evol_groups_first_seven.tsv
#python scripts/grow-sbtmh.py test_grow_sbt --sbt test_new_growsbt/sbt-7paths.protein_k7_scaled1.sbt.zip --ksize 7 --input-is-directory --scaled 1 --alphabet protein --subset-info-colname accession --subset-csv gtdb_evol_groups_first_seven.tsv

#python scripts/match-accession-to-filename.py gtdb_evol_groups_first_seven.tsv --sigfile-directory gtdb_faa_pathsonly
#python scripts/match-accession-to-filename.py gtdb_evol_groups.tsv --sigfile-directory gtdb_faa_pathsonly

#python scripts/forage-sbt.py test_new_growsbt/sbt-7paths.dayhoff_k13_scaled1.sbt.zip --query-csv gtdb_evol_groups_first_seven.with_filenames.csv --distance-from-species-csv test_new_growsbt/sbt-7paths_dayhoff_k13.pathwise-dist-from-species.csv --distance-from-species-plot test_new_growsbt/sbt-7paths_dayhoff_k13.pathwise-dist-from-species.svg

#python scripts/grow-sbtmh.py test_grow_sbt --sbt test_new_growsbt/sbt-7paths.dayhoff_k13_scaled1_singleton.sbt.zip --ksize 13 --input-is-directory --scaled 1 --alphabet dayhoff --subset-info-colname accession --subset-csv gtdb_evol_groups_first_seven.tsv --singleton

# find dna filenames
#python scripts/match-accession-to-filename.py gtdb_evol_groups_first_seven.tsv --sigfile-directory /home/ctbrown/gtdbtk/release89/fastani/database --output-csv /home/ntpierce/2020-gtdb-smash/gtdb_evol_groups_first_seven.with_dna_filenames.csv
#find rna filenames
#python scripts/match-accession-to-filename.py gtdb_evol_groups_first_seven.tsv --sigfile-directory "gtdb_r89_rep_genomes_protein_fna/*" --output-csv /home/ntpierce/2020-gtdb-smash/gtdb_evol_groups_first_seven.with_rna_filenames.csv

# find dna filenames
#python scripts/match-accession-to-filename.py gtdb_evol_groups.tsv --sigfile-directory /home/ctbrown/gtdbtk/release89/fastani/database --output-csv /home/ntpierce/2020-gtdb-smash/gtdb_evol_groups.with_dna_filenames.csv
#find rna filenames
#python scripts/match-accession-to-filename.py gtdb_evol_groups.tsv --sigfile-directory "gtdb_r89_rep_genomes_protein_fna/*" --output-csv /home/ntpierce/2020-gtdb-smash/gtdb_evol_groups.with_rna_filenames.csv --full-paths



#python scripts/grow-sbtmh.py test_grow_sbt --sbt test_new_growsbt/sbt-7paths.dayhoff_k13_scaled1_singleton.sbt.zip --ksize 13 --input-is-directory --scaled 1 --alphabet dayhoff --subset-info-colname accession --subset-csv gtdb_evol_groups_first_seven.tsv --singleton

#python scripts/grow-sbtmh.py evol_distances/compute/protein/GCA_002336895_dayhoff_scaled200_k19.sig evol_distances/compute/protein/GCA_003506995_dayhoff_scaled200_k19.sig evol_distances/compute/protein/GCA_002424095_dayhoff_scaled200_k19.sig evol_distances/compute/protein/GCA_002429065_dayhoff_scaled200_k19.sig evol_distances/compute/protein/GCA_002437405_dayhoff_scaled200_k19.sig evol_distances/compute/protein/GCA_003527555_dayhoff_scaled200_k19.sig evol_distances/compute/protein/GCA_003506995_dayhoff_scaled200_k19.sig evol_distances/compute/protein/GCF_000739655_dayhoff_scaled200_k19.sig evol_distances/compute/protein/GCF_000172155_dayhoff_scaled200_k19.sig --sbt evol_distances/sbts/evol_paths_pep.dayhoff_scaled200_k19.sbt.zip --ksize 19 --scaled 200 --alphabet dayhoff

#python scripts/forage-sbt.py evol_distances/sbts/evol_paths_pep.dayhoff_scaled200_k19.sbt.zip --query-csv /home/ntpierce/2020-gtdb-smash/gtdb_evol_groups.with_filenames.csv --signature-name-column filename --distance-from-species-csv evol_distances/distances/evol_paths_pep.dayhoff_scaled200_k19.jaccard_from_species.csv --distance-from-species-plot evol_distances/distances/plots/evol_paths_pep.dayhoff_scaled200_k19.jaccard_from_species.svg #2>evol_distances/logs/forage/evol_paths_pep.dayhoff_scaled200_k19.forage.log

#python scripts/grow-sbtmh.py evol_distances/compute/dna/GCA_002437405_nucleotide_scaled1000_k31.sig evol_distances/compute/dna/GCA_003527555_nucleotide_scaled1000_k31.sig evol_distances/compute/dna/GCF_000739635_nucleotide_scaled1000_k31.sig evol_distances/compute/dna/GCA_002336895_nucleotide_scaled1000_k31.sig evol_distances/compute/dna/GCA_003506995_nucleotide_scaled1000_k31.sig evol_distances/compute/dna/GCA_002424095_nucleotide_scaled1000_k31.sig evol_distances/compute/dna/GCA_002429065_nucleotide_scaled1000_k31.sig evol_distances/compute/dna/GCA_002437405_nucleotide_scaled1000_k31.sig evol_distances/compute/dna/GCA_003527555_nucleotide_scaled1000_k31.sig evol_distances/compute/dna/GCA_003506995_nucleotide_scaled1000_k31.sig evol_distances/compute/dna/GCF_000739655_nucleotide_scaled1000_k31.sig evol_distances/compute/dna/GCF_000172155_nucleotide_scaled1000_k31.sig --sbt testdna.sbt.zip --ksize 31 --scaled 1000 --alphabet nucleotide


#python scripts/grow-sbtmh.py --sbt evol_distances/sbts/testpep.dayhoff_scaled200_k21.sbt.zip --ksize 21 --scaled 200 --alphabet dayhoff --input-files evol_distances/compute/protein/GCA_002424095_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002429065_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002437405_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_003506995_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_003527555_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002386905_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002387565_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002424095_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002429065_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002437405_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_003506995_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002307065_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002292185_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_003527555_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002424095_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002429065_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002437405_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_003506995_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002712205_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002862135_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002323995_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002424095_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002429065_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002437405_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_003506995_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002323995_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002696885_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002712205_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002424095_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002429065_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002437405_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_003527555_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCF_000739635_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002336895_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_003506995_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002424095_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002429065_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_002437405_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_003527555_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCA_003506995_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCF_000739655_dayhoff_scaled200_k21.sig evol_distances/compute/protein/GCF_000172155_dayhoff_scaled200_k21.sig

#python scripts/match-accession-to-filename.py gtdb_evol_groups2.tsv --sigfile-directory gtdb_r89_rep_genomes_faa
#python scripts/select-representative-lineages.py gtdb_evol_groups2.with_filenames.csv --representative-rank family 
#python scripts/select-representative-lineages.py gtdb_evol_groups2.with_filenames.csv --representative-rank order

# fix lineages for original gtdb_evol_groups

#python scripts/fix-lineages.py gtdb_evol_groups.tsv --lineages-csv gtdb-lineages.csv --output-csv gtdb_evol_groups.correct-lineages.csv
#python scripts/match-accession-to-filename.py gtdb_evol_groups.correct-lineages.csv --sigfile-directory gtdb_r89_rep_genomes_faa
#python scripts/select-representative-lineages.py gtdb_evol_groups.correct-lineages.with_filenames.csv --representative-rank family 
#python scripts/select-representative-lineages.py gtdb_evol_groups.correct-lineages.with_filenames.csv --representative-rank order

#python scripts/match-accession-to-filename.py gtdb-lineages.csv --sigfile-directory gtdb_r89_rep_genomes_faa --output-csv gtdb-lineages.protein-filenames.csv


#python scripts/select-representative-lineages.py gtdb-lineages.protein-filenames.csv --representative-rank order
#python scripts/select-representative-lineages.py gtdb-lineages.protein-filenames.csv --representative-rank family
#python scripts/select-representative-lineages.py gtdb-lineages.protein-filenames.csv --representative-rank genus

python scripts/select-representative-lineages.py gtdb-lineages.protein-filenames.csv --representative-rank family --nth-to-select 1

