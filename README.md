# nfixBiodiversity_DiazoTIME

## Files:
- metabolic_diazo.sh - script used to run METABOLIC on diazotroph genomes. Formatted for the Boston University server (SCC). 

- GTDB_diazotroph_genomes_r214_DiazoTIME.R - R script for general analyses and graphs of taxonomy and metabolic data.

- change_HMMER_output_names_to_combine_GTDB_and_NCBI.py - takes amino acid fasta files from GTDB (or AnnoTree), which have sequence names that use NCBI contig names, and renames the sequences so they have both GTDB and NCBI accessions. Uses the dictionary “combined_gtdb_r214_genome_contigs_dict.txt”.

- add_taxonomy_to_fastas.R - Added GTDB taxonomy with R. Used taxonomy file downloaded from GTDB, parsed to just have “accession” and “taxonomy” columns: combined_gtdb_meta_taxonomy_only_r214.tsv. 

- Had to manually retrieve .fna files for genes (AnnoTree only have .faa files). This used the scripts filterbyname_Nif?_fna.sh (one for each gene) and a list of contigs names (GTDB_nif?_gene_accessions.txt; one for each gene, containing only the NCBI gene accession numbers from the .faa files originally from AnnoTtree) to retrieve sequences from the CDS nucleic acid sequences file downloaded from the GTDB database (gtdb_proteins_reps_r214_genes.fna).