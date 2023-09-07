
cut -f 11 -d "," /home/ntpierce/2020-gtdb-smash/gtdb-lineages.protein-filenames.n1th-representative-at-genus.csv | tail -n +2 > /home/ntpierce/thumper/gtdb-genus-n1th/gtdb-pep.rep-genus-n1th.protein-filenames.txt
cut -f 11 -d "," /home/ntpierce/2020-gtdb-smash/gtdb-lineages.rna-filenames.n1th-representative-at-genus.csv | tail -n +2 > /home/ntpierce/thumper/gtdb-genus-n1th/gtdb-pep.rep-genus-n1th.rna-filenames.txt
cut -f 11 -d "," /home/ntpierce/2020-gtdb-smash/gtdb-lineages.dna-filenames.n1th-representative-at-genus.csv | tail -n +2 > /home/ntpierce/thumper/gtdb-genus-n1th/gtdb-pep.rep-genus-n1th.dna-filenames.txt

