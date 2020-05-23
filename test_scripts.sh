



#python scripts/grow-sbtmh.py test_grow_sbt --sbt test_new_growsbt/sbt-7paths.dayhoff_k13_scaled1.sbt.zip --ksize 13 --input-is-directory --scaled 1 --alphabet dayhoff --subset-info-colname accession --subset-csv gtdb_evol_groups_first_seven.tsv

#python scripts/match-accession-to-filename.py gtdb_evol_groups_first_seven.tsv --sigfile-directory gtdb_faa_pathsonly
#python scripts/match-accession-to-filename.py gtdb_evol_groups.tsv --sigfile-directory gtdb_faa_pathsonly

#python scripts/forage-sbt.py test_new_growsbt/sbt-7paths.dayhoff_k13_scaled1.sbt.zip --query-csv gtdb_evol_groups_first_seven.with_filenames.csv --distance-from-species-csv test_new_growsbt/sbt-7paths_dayhoff_k13.pathwise-dist-from-species.csv --distance-from-species-plot test_new_growsbt/sbt-7paths_dayhoff_k13.pathwise-dist-from-species.svg

python scripts/grow-sbtmh.py test_grow_sbt --sbt test_new_growsbt/sbt-7paths.dayhoff_k13_scaled1_singleton.sbt.zip --ksize 13 --input-is-directory --scaled 1 --alphabet dayhoff --subset-info-colname accession --subset-csv gtdb_evol_groups_first_seven.tsv --singleton

