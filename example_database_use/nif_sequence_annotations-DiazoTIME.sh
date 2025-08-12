#!/bin/bash -l

#$ -l h_rt=4:00:00
#$ -N nif_IMG_dia
#$ -j y
#$ -m ea

#$ -pe omp 28
#$ -l mem_per_core=4G

## Use DIAMOND to search nifHDK amino acid sequences against database of GTDB sequences in the DiazoTIME db ##
## Will output a table for each gene that's queried, giving the top hit in the GTDB database ##
ml diamond/2.0.11

## User files - Only one is required, but the appropriate DIAMOND db needs to be used ##
NIFH_INPUT="your_nifH_sequences.faa.gz";
NIFK_INPUT="your_nifK_sequences.faa.gz";                                      
NIFD_INPUT="your_nifD_sequences.faa.gz";

## DIAMOND DB files ##
DB=~/Nif_DIAMOND_dbs;

## Classifying nifH sequences ##
echo 'Starting DIAMOND search of nifH hits!'
diamond blastp -q $NIFH_INPUT -d $DB/gtdb_r214_NifH.dmnd -o nifH_DiazoTIME_diamond.txt -f 6 qseqid sseqid evalue bitscore -e 0.01 -k 1 --ultra-sensitive --quiet 
-p 28

## Classifying nifK sequences ##
echo 'Starting DIAMOND search of nifK hits!'
diamond blastp -q $NIFK_INPUT -d $DB/gtdb_r214_NifK.dmnd -o nifK_DiazoTIME_diamond.txt -f 6 qseqid sseqid evalue bitscore -e 0.01 -k 1 --ultra-sensitive --quiet 
-p 28

## Classifying nifD sequences ##
echo 'Starting DIAMOND search of nifD hits!'
diamond blastp -q $NIFD_INPUT -d $DB/gtdb_r214_NifD.dmnd -o nifD_DiazoTIME_diamond.txt -f 6 qseqid sseqid evalue bitscore -e 0.01 -k 1 --ultra-sensitive --quiet 
-p 28
