#!/bin/bash -l

#$ -l h_rt=24:00:00
#$ -N nif_meta_1
#$ -j y
#$ -m ea
#$ -M jdamashe@hamilton.edu

#$ -pe omp 28
#$ -l mem_per_core=18G

# Run on the Boston University server (SCC) by Julian Damashek (jdamashe@hamilton.edu; juliandamashek@gmail.com)

ml metabolic/4.0

# make list of genomes and decompress fasta files
# input to METABOLIC was .faa files ofr CDS from genomes

#cd metabolic_input_faas
#ls *.faa.gz > files.txt
#while read -r i; do gunzip $i; done < files.txt
#cd ..

# run METABOLIC
METABOLIC-G.pl -in metabolic_input_faas -t 28 -o diazo_gtdb_metabolic_out

# organize
cd diazo_gtdb_metabolic_out && tar -zcf intermediate_files.tar.gz intermediate_files && rm -rf intermediate_files && cd ..
cd metabolic_input_faas && ls *.faa > files.txt && ls *.fasta >> files.txt && ls *.gene >> files.txt
for i in `cat files.txt`; do gzip $i; done && rm -f files.txt && cd ..

# done!
