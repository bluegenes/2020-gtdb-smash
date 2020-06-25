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
#snakemake -s sig-sbt-forage.snakefile --configfile config/grow_gtdb_sbts_rna_dna.yml --profile default  --jobs 4 --nolock
#snakemake -s sig-sbt-forage.snakefile --configfile config/grow_gtdb_sbts.yml --profile farm --cluster-config config/bml_clusterconfig.yml --jobs 50

#snakemake -s representative-index.snakefile --configfile config/representative-index.yml --profile default --jobs 1 #--nolock --restart-times 0 -n

#snakemake -s representative-index.snakefile --configfile config/representative-index.yml --profile farm --cluster-config config/grow-sbt-clusterconfig.yml --jobs 1 -n #--nolock --restart-times 0 -n

# set mem limit here so no more than one lca job runs at once

#snakemake -s build-representative-index.snakefile --configfile config/build-representative-index.yml --profile farm --cluster-config config/grow-sbt-clusterconfig.yml --jobs 40 #--until index_lca # --restart-times 0 --until index_lca --latency-wait 120 #40 #--resources mem_mb=340000 # -n #--nolock --restart-times 0 -n

#snakemake -s build-representative-index.snakefile --configfile config/build-full-gtdb.yml --profile farm --cluster-config config/grow-sbt-clusterconfig.yml --jobs 40 --nolock #--restart-times 0 -n

#snakemake -s build-representative-index.snakefile --configfile config/build-representative-index.yml --profile farm --cluster-config cluster_config.yml --jobs 40 #--until sourmash_compute_protein #--until index_lca #--resources mem_mb=150000 #--latency-wait 120 #--restart-times 0  #--latency-wait 120 #-n

#snakemake -s build-representative-index.snakefile --configfile config/build-representative-index.yml --profile farm --cluster-config config/grow-sbt-clusterconfig.yml --jobs 40 # --restart-times 0 #40 #--resources mem_mb=340000 # -n #--nolock --restart-times 0 -n

#snakemake -s build-representative-index.snakefile --configfile config/build-representative-index.yml --profile default --jobs 1  #-R aggregate_gather_to_tax --restart-times 0

#snakemake -s build-representative-index.snakefile --configfile config/build-representative-index-dev.yml --profile default --jobs --until sourmash_compute_protein #  -R aggregate_gather_to_tax # -R gather_to_tax

#snakemake -s classify_sigs.snakefile --configfile tara_delmont/dayhoff_config.yml --profile default --jobs 1
#snakemake -s classify_sigs.snakefile --configfile tara_delmont/dayhoff_config.yml --profile farm --cluster-config tara_delmont/cluster_config.yml --jobs 20

#snakemake -s classify_sigs.snakefile --configfile tara_delmont/protein_config.yml --profile farm --cluster-config tara_delmont/cluster_config.yml --jobs 10
#snakemake -s classify_sigs.snakefile --configfile tara_delmont/protein_config.yml --profile farm --cluster-config tara_delmont/cluster_config.yml --jobs 40
#snakemake -s classify_sigs.snakefile --configfile tara_delmont/protein_config.yml --profile default --jobs 30 

#snakemake -s compute-sigs.snakefile --configfile tara_delmont/compute-protein-sigs.yml --profile farm --cluster-config tara_delmont/cluster_config.yml --jobs 40
snakemake -s compute-sigs.snakefile --configfile tara_delmont/compute-nucleotide-sigs.yml --profile farm --cluster-config tara_delmont/cluster_config.yml --jobs 40


