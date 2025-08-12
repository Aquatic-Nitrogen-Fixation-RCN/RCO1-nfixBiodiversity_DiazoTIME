# nfixBiodiversity_DiazoTIME

The Diazotroph Taxonomic Identity and MEtabolism (DiazoTIME) database contains annotated taxonomy and metabolic predictions for nifH-, nifD-, and nifK- containing genomes (2798 genomes) in the Genome Taxonomy Database (GTDB; r214; Parks et al. 2022). All files for this database are found on Zenodo (https://doi.org/10.5281/zenodo.15311395). This is a supplemental repository that contains shell, python, and R scripts that were used in database generation. Example scripts for using the database are also included in this repository. 


## Contents overview

	|- README                # the top level description of content (this doc)
	|
	|- database_generation/  # code and supplemental reference files used to generate DiazoTIME
	| |- add_taxomony_to_fastas.R    
	| |- change_HMMER_output_names_to_combine_GTDB_and_NCBI.py   
	| |- filterbyname_NifD_fna.sh    
	| |- filterbyname_NifH_fna.sh
	| |- filterbyname_NifK_fna.sh
	| |- metabolic_diazo.sh
	| |+ reference/     
	|
	|- example_database_use/  # code and example files for using DiazoTIME 
	| |- Nif_DIAMOND_dbs/
	| |- DiazoTIME_GTDBr214_taxonomy_and_METABOLIC.csv
	| |- DiazoTIME_GTDBr214_taxonomy_and_METABOLIC.xlsx
	| |- nif_DiazoTIME_data_processing.R
	+ +- nif_sequence_annotations-DiazoTIME.sh     
	
	

## File descriptions

- **database_generation/**: Scripts used in the construction of the DiazoTIME database.
    
    - **metabolic_diazo.sh**: Script used to run METABOLIC on diazotroph genomes. Formatted for the Boston University server (SCC). 

    - **change_HMMER_output_names_to_combine_GTDB_and_NCBI.py**: Takes amino acid fasta files from GTDB (or AnnoTree), which have sequence names that use NCBI contig names, and renames the sequences so they have both GTDB and NCBI accessions. Uses the dictionary “combined_gtdb_r214_genome_contigs_dict.txt”.

    - **add_taxonomy_to_fastas.R**: Added GTDB taxonomy with R. Used taxonomy file downloaded from GTDB, parsed to just have “accession” and “taxonomy” columns: combined_gtdb_meta_taxonomy_only_r214.tsv. 

    - **reference/**: Reference files used in R script.

    - Had to manually retrieve .fna files for genes (AnnoTree only have .faa files). This used the scripts filterbyname_Nif?_fna.sh (one for each gene) and a list of contigs names (GTDB_nif?_gene_accessions.txt; one for each gene, containing only the NCBI gene accession numbers from the .faa files originally from AnnoTtree) to retrieve sequences from the CDS nucleic acid sequences file downloaded from the GTDB database (gtdb_proteins_reps_r214_genes.fna).
    
    
- **example_database_use/**: Example use case of DiazoTIME database. Users can place their own nifH/nifD/nifK sequence files in amino acid format (.faa) inside the directory, run "nif_sequence_annotations-DiazoTIME.sh", and process the output using "nif_DiazoTIME_data_processing.R". Sample metadata files can also be integrated into the output of the R script if provided. 

    - **nif_sequence_annotations-DiazoTIME.sh**: A shell script file that can be used to match nif sequences with their closest representative in the DiazoTIME database of GTDB genomes. 
    
    - **Nif_DIAMOND_dbs/**: Three diamond databases that are called by the "nif_sequence_annotations-DiazoTIME.sh" script.
    
    - **DiazoTIME_GTDBr214_taxonomy_and_METABOLIC.csv** and **DiazoTIME_GTDBr214_taxonomy_and_METABOLIC.csv**: The actual DiazoTIME database downloaded from the Zenodo database repository (https://doi.org/10.5281/zenodo.15311395), which includes GTDB genome accessions, taxonomy, and metabolic predictions. 
    
    - **nif_DiazoTIME_data_processing.R**: R code that takes the output from "nif_sequence_annotations-DiazoTIME.sh" and merges it with the DiazoTIME database ("DiazoTIME_GTDBr214_taxonomy_and_METABOLIC.csv")