


#python scripts/grow-sbtmh.py gtdb_r89_rep_genomes_faa/GB_GCA_003149575.1_protein.faa.gz gtdb_r89_rep_genomes_faa/RS_GCF_001072465.1_protein.faa.gz --sbt test_out.sbt.zip --ksize 7 --alphabet protein --track-abundance --scaled 1

#python scripts/forage-sbtmh.py GB_GCA_003149575.1 --sbt test_out.sbt.zip

#python scripts/grow-sbtmh.py test_grow_sbt --sbt test_grow_sbt/sbt-7paths_dayhoff_k13.sbt.zip --ksize 13 --input-is-directory --scaled 1 --alphabet dayhoff --track-abundance


python scripts/forage-sbt.py test_grow_sbt/sbt-7paths_dayhoff_k13.sbt.zip --query_csv gtdb_evol_groups_first_seven.sbtinfo.csv  
