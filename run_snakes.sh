 #!/bin/bash

 #snakemake -s gtdb-smash.snakefile --configfiles config/config.yml --profile default --nolock -n

#snakemake -s gtdb-smash.snakefile --configfiles config/protein-only.yml --profile default --nolock -n
#snakemake -s gtdb-smash.snakefile --configfiles config/protein-only.yml --profile farm --cluster-config cluster_config.yml --jobs 15

#snakemake -s simple-smash.snakefile --configfiles config/genomes_with_path_ids.yml --profile default --nolock -n

#snakemake -s simple-smash.snakefile --configfiles config/genomes_with_path_ids.yml --profile farm --cluster-config cluster_config.yml --jobs 30


#snakemake -s grow-sbts.snakefile --configfiles config/grow_gtdb_sbts.yml --profile farm --cluster-config cluster_config.yml --jobs 10 --restart-times 0

#snakemake -s grow-and-forage-sbt.snakefile --configfile config/grow_first_seven_sbts.yml --profile default --jobs 17 --nolock

#snakemake -s grow-and-forage-sbt.snakefile --configfile config/grow_gtdb_sbts.yml --profile default  --jobs 1 --restart-times 0


#snakemake -s grow-and-forage-sbt.snakefile --configfile config/grow_gtdb_sbts.yml --profile farm --cluster-config config/grow-sbt-clusterconfig.yml --jobs 25 --nolock

#snakemake -s sig-sbt-forage.snakefile --configfile config/grow_gtdb_sbts.yml --profile farm --cluster-config config/grow-sbt-clusterconfig.yml --jobs 15
#snakemake -s sig-sbt-forage.snakefile --configfile config/grow_gtdb_sbts.yml --profile default --nolock --jobs 1 #-n


#snakemake -s sig-sbt-forage.snakefile --configfile config/grow_gtdb_sbts_rna_dna.yml --profile farm --cluster-config config/bml_clusterconfig.yml --jobs 20 --nolock
#snakemake -s sig-sbt-forage.snakefile --configfile config/grow_gtdb_sbts.yml --profile farm --cluster-config config/bml_clusterconfig.yml --jobs 50

#snakemake -s representative-index.snakefile --configfile config/representative-index.yml --profile default --jobs 1 #--nolock --restart-times 0 -n

#snakemake -s representative-index.snakefile --configfile config/representative-index.yml --profile farm --cluster-config config/grow-sbt-clusterconfig.yml --jobs 1 -n #--nolock --restart-times 0 -n

# set mem limit here so no more than one lca job runs at once
snakemake -s build-representative-index.snakefile --configfile config/build-representative-index.yml --profile farm --cluster-config config/grow-sbt-clusterconfig.yml --jobs 15 --resources mem_mb=300000 #-n #--nolock --restart-times 0 -n

