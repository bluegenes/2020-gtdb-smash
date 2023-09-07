


#for i in $(seq 3 100)
#  do 
#snakemake -s compute-lossy-sigs.snakefile --configfile config/compute_gtdb_lossy.yml --profile default --conda-create-envs-only #--jobs 20 #--batch all=$i/100 --jobs 20
#  done

snakemake -s compute-lossy-sigs.snakefile --configfile config/compute_gtdb_lossy.yml --profile farm --cluster-config cluster_config.yml --jobs 25 --resources mem_mb=300000

