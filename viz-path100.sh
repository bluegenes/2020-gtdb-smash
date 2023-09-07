
# path 100
# accession,path,rank,lineage,filename
#GCA_003166025,path100,species,d__Bacteria;p__Actinobacteriota;c__Acidimicrobiia;o__Acidimicrobiales;f__Palsa-688;g__Fen-671;s__Fen-671 sp003166025,GB_GCA_003166025.1_protein.faa.gz
#GCA_003158835,path100,genus,d__Bacteria;p__Actinobacteriota;c__Acidimicrobiia;o__Acidimicrobiales;f__Palsa-688;g__Fen-671;s__Fen-671 sp003158835,GB_GCA_003158835.1_protein.faa.gz
#GCA_003156795,path100,family,d__Bacteria;p__Actinobacteriota;c__Acidimicrobiia;o__Acidimicrobiales;f__Palsa-688;g__PALSA-688;s__PALSA-688 sp003156795,GB_GCA_003156795.1_protein.faa.gz
#GCF_900128965,path100,order,d__Bacteria;p__Actinobacteriota;c__Acidimicrobiia;o__Acidimicrobiales;f__Acidimicrobiaceae;g__Ferrithrix;s__Ferrithrix thermotolerans,RS_GCF_900128965.1_protein.faa.gz
#GCA_002727895,path100,class,d__Bacteria;p__Actinobacteriota;c__Acidimicrobiia;o__Microtrichales;f__MedAcidi-G1;g__UBA9410;s__UBA9410 sp002727895,GB_GCA_002727895.1_protein.faa.gz
#GCF_000308455,path100,phylum,d__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Mycobacteriales;f__Mycobacteriaceae;g__Nocardia;s__Nocardia abscessus,RS_GCF_000308455.1_protein.faa.gz
#GCA_002424095,path100,superkingdom,d__Bacteria;p__Firmicutes_A;c__Clostridia;o__4C28d-15;f__UBA1242;g__UBA5489;s__UBA5489 sp002424095,GB_GCA_002424095.1_protein.faa.gz
#
#sourmash compare --dayhoff -k 19 --scaled 100
#sourmash compare --hp -k 33 --scaled 100 


conda activate forage-env

#alphabet=dayhoff
#alphabet=protein
alphabet=hp
#ksize=19
#ksize=11
ksize=35
#smash_k=$(expr $ksize /* 3)
smash_k=$((ksize*3))
sigdir="/home/ntpierce/2020-gtdb-smash/index-gtdb/sigs/protein/${alphabet}/k${ksize}"

sourmash compare --${alphabet} -k ${smash_k} ${sigdir}/GCA_003166025_*.sig ${sigdir}/GCA_003158835_*.sig ${sigdir}/GCA_003156795_*.sig ${sigdir}/GCF_900128965_*.sig ${sigdir}/GCA_002727895_*.sig ${sigdir}/GCF_000308455_*.sig ${sigdir}/GCA_002424095_*.sig --csv path100.${alphabet}-k${ksize}.compare.csv -o path100.${alphabet}-k${ksize}.compare.np

sourmash plot path100.${alphabet}-k${ksize}.compare.np --labels


