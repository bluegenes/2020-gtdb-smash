


python scripts/grow-sbtmh.py gtdb_r89_rep_genomes_faa/GB_GCA_003149575.1_protein.faa.gz gtdb_r89_rep_genomes_faa/RS_GCF_001072465.1_protein.faa.gz --sbt test_out.sbt.zip --ksize 7 --alphabet protein --track-abundance --scaled 1

python scripts/forage-sbtmh.py GB_GCA_003149575.1 --sbt test_out.sbt.zip
