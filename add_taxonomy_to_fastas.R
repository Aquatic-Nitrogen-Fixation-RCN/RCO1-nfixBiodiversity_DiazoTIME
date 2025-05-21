library(tidyverse)
library(phylotools)
library(R.utils)

setwd("~/Documents/data/RCN_N2fix/RCN1/omics/GTDB/")

# add taxonomy to sequence names

# #load genome taxonomy data - updated to GTDB r214
# input file is the genome metadata file downloaded from GTDB r214, parsed to just have the 'accession' and 'taxonomy' columns
genome.tax <- read_tsv("combined_gtdb_meta_taxonomy_only_r214.tsv") %>%
  mutate(Genome=gsub('.[^.]*$', '', accession), .before=2, Genome=str_sub(Genome,4,-1)) %>% select(-accession)

# rename nifH fna file
# read in fasta file
nifH.gtdb.fnas <- read.fasta("gtdb_r214_nifH_renamed.fna.gz")

# make table to rename
nifH.seq.names <- nifH.gtdb.fnas %>% select(seq.name) %>% mutate(gtdb.name=str_sub(seq.name,4,-1), Genome=str_extract(gtdb.name, pattern="^\\w+")) %>% 
  left_join(genome.tax, by="Genome") %>% select(c(seq.name,gtdb_taxonomy)) %>% 
  mutate(new.name=paste(seq.name,gtdb_taxonomy,sep="-")) %>% select(-gtdb_taxonomy)

# rename
nifH.gtdb.fnas.with.tax <- rename.fasta(infile="gtdb_r214_nifH_renamed.fna.gz", ref_table=nifH.seq.names, outfile="gtdb_r214_nifH_with_tax.fna")
gzip(filename="gtdb_r214_nifH_with_tax.fna", destname="gtdb_r214_nifH_with_tax.fna.gz", overwrite=T)



# rename nifD fna file
# read in fasta file
nifD.gtdb.fnas <- read.fasta("gtdb_r214_nifD_renamed.fna.gz")

# make table to rename
nifD.seq.names <- nifD.gtdb.fnas %>% select(seq.name) %>% mutate(gtdb.name=str_sub(seq.name,4,-1), Genome=str_extract(gtdb.name, pattern="^\\w+")) %>% 
  left_join(genome.tax, by="Genome") %>% select(c(seq.name,gtdb_taxonomy)) %>% 
  mutate(new.name=paste(seq.name,gtdb_taxonomy,sep="-")) %>% select(-gtdb_taxonomy)

# rename
nifD.gtdb.fnas.with.tax <- rename.fasta(infile="gtdb_r214_nifD_renamed.fna.gz", ref_table=nifD.seq.names, outfile="gtdb_r214_nifD_with_tax.fna")
gzip(filename="gtdb_r214_nifD_with_tax.fna", destname="gtdb_r214_nifD_with_tax.fna.gz", overwrite=T)



# rename nifK fna file
# read in fasta file
nifK.gtdb.fnas <- read.fasta("gtdb_r214_nifK_renamed.fna.gz")

# make table to rename
nifK.seq.names <- nifK.gtdb.fnas %>% select(seq.name) %>% mutate(gtdb.name=str_sub(seq.name,4,-1), Genome=str_extract(gtdb.name, pattern="^\\w+")) %>% 
  left_join(genome.tax, by="Genome") %>% select(c(seq.name,gtdb_taxonomy)) %>% 
  mutate(new.name=paste(seq.name,gtdb_taxonomy,sep="-")) %>% select(-gtdb_taxonomy)

# rename
nifK.gtdb.fnas.with.tax <- rename.fasta(infile="gtdb_r214_nifK_renamed.fna.gz", ref_table=nifK.seq.names, outfile="gtdb_r214_nifK_with_tax.fna")
gzip(filename="gtdb_r214_nifK_with_tax.fna", destname="gtdb_r214_nifK_with_tax.fna.gz", overwrite=T)



# rename nifH faa file (original names formatted slightly differently)
# read in fasta file
nifH.gtdb.faas <- read.fasta("gtdb_r214_diazo_nifH_genes_renamed.faa.gz")

# make table to rename
nifH.faa.names <- nifH.gtdb.faas %>% select(seq.name) %>% mutate(new.name=str_replace(seq.name, pattern=" \\(", replacement="-"), new.name=str_replace(new.name, pattern="\\)", replacement="")) %>% 
  mutate(Genome=str_extract(new.name, pattern="^\\w+")) %>% left_join(genome.tax, by="Genome") %>% 
  mutate(new.name=paste(new.name,gtdb_taxonomy,sep="-")) %>% select(-c(gtdb_taxonomy,Genome))

# rename
nifH.gtdb.faas.with.tax <- rename.fasta(infile="gtdb_r214_diazo_nifH_genes_renamed.faa.gz", ref_table=nifH.faa.names, outfile="gtdb_r214_nifH_with_tax.faa")
gzip(filename="gtdb_r214_nifH_with_tax.faa", destname="gtdb_r214_nifH_with_tax.faa.gz", overwrite=T)



# rename nifD faa file 
# read in fasta file
nifD.gtdb.faas <- read.fasta("gtdb_r214_diazo_nifD_genes_renamed.faa.gz")

# make table to rename
nifD.faa.names <- nifD.gtdb.faas %>% select(seq.name) %>% mutate(new.name=str_replace(seq.name, pattern=" \\(", replacement="-"), new.name=str_replace(new.name, pattern="\\)", replacement="")) %>% 
  mutate(Genome=str_extract(new.name, pattern="^\\w+")) %>% left_join(genome.tax, by="Genome") %>% 
  mutate(new.name=paste(new.name,gtdb_taxonomy,sep="-")) %>% select(-c(gtdb_taxonomy,Genome))

# rename
nifD.gtdb.faas.with.tax <- rename.fasta(infile="gtdb_r214_diazo_nifD_genes_renamed.faa.gz", ref_table=nifD.faa.names, outfile="gtdb_r214_nifD_with_tax.faa")
gzip(filename="gtdb_r214_nifD_with_tax.faa", destname="gtdb_r214_nifD_with_tax.faa.gz", overwrite=T)



# rename nifK faa file 
# read in fasta file
nifK.gtdb.faas <- read.fasta("gtdb_r214_diazo_nifK_genes_renamed.faa.gz")

# make table to rename
nifK.faa.names <- nifK.gtdb.faas %>% select(seq.name) %>% mutate(new.name=str_replace(seq.name, pattern=" \\(", replacement="-"), new.name=str_replace(new.name, pattern="\\)", replacement="")) %>% 
  mutate(Genome=str_extract(new.name, pattern="^\\w+")) %>% left_join(genome.tax, by="Genome") %>% 
  mutate(new.name=paste(new.name,gtdb_taxonomy,sep="-")) %>% select(-c(gtdb_taxonomy,Genome))

# rename
nifK.gtdb.faas.with.tax <- rename.fasta(infile="gtdb_r214_diazo_nifK_genes_renamed.faa.gz", ref_table=nifK.faa.names, outfile="gtdb_r214_nifK_with_tax.faa")
gzip(filename="gtdb_r214_nifK_with_tax.faa", destname="gtdb_r214_nifK_with_tax.faa.gz", overwrite=T)

