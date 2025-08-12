
####################################### LOAD PACKAGES ####################################### 

library(tidyverse)
library(dplyr)

#################################### PREPARE ENVIRONMENT ####################################

rm(list = ls())
Sys.setenv(LANG = "en")
options(stringsAsFactors = FALSE)


###################################### IMPORT DATA #########################################

## DiazoTIME database - includes GTDB genome accession numbers, taxonomy, and metabolic classification ##
## Genomes downloaded from GTDB r214 ##
diazotime <- read.csv("example_database_use/DiazoTIME_GTDBr214_taxonomy_and_METABOLIC.csv") %>% mutate(accession=gsub('.[^.]*$', '', accession))
 
## Metadata for user sequences ##
## Should include sequence IDs associated with the .faa sequences that were queried against DiazoTIME db ##
meta <- read_csv("input_your_metadata_here.csv")
meta$Sequence_ID <- as.character(meta$Sequence_ID)

## nifH DIAMOND output ##
diamond.nifH <- read_tsv("nifH_DiazoTIME_diamond.txt", col_names=F)
colnames(diamond.nifH) <- c("qseqid","accession","evalue","bitscore")
diamond.nifH <- diamond.nifH %>% mutate(gene="NifH", Sequence_ID=as.character(gsub('_.*', '', diamond.nifH$qseqid)))

## nifK DIAMOND output ##
diamond.nifK <- read_tsv("nifK_DiazoTIME_diamond.txt", col_names=F)
colnames(diamond.nifK) <- c("qseqid","accession","evalue","bitscore")
diamond.nifK <- diamond.nifK %>% mutate(gene="NifK", Sequence_ID=as.character(gsub('_.*', '', diamond.nifK$qseqid)))

## nifD DIAMOND output ##
diamond.nifD <- read_tsv("nifK_DiazoTIME_diamond.txt", col_names=F)
colnames(diamond.nifD) <- c("qseqid","accession","evalue","bitscore")
diamond.nifD <- diamond.nifD %>% mutate(gene="NifD", Sequence_ID=as.character(gsub('_.*', '', diamond.nifD$qseqid)))

## Combine all 3 nif genes together ##
diamond.all.nif <- rbind(diamond.nifH, diamond.nifK, diamond.nifD) 


################################## MERGE WITH DIAZOTIME DB ##################################

## Join with DiazoTIME - taxonomy and metabolic predictions ##
diamond.all.nif <- diamond.all.nif %>% left_join(diazotime, by="accession") %>% 
  left_join(meta, by="Sequence_ID", relationship="many-to-many") ## add in metadata
 