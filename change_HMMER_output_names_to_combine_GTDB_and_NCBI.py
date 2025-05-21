#!/usr/bin/env python

# Takes .faa fasta file coming from GTDB files or AnnoTree - MUST end with '.faa' extension
#		and uses a dictionary that has both the GTDB and NCBI names for the GTDB genomes. (format: 'NCBI_name':GTDB_name)
#		Currently the .faa file will only have the NCBI names (from the original contigs).
#		This script saves the name for each sequence as: >GTDB_ACCESSION (NCBI GENE ACCESSION)

# Usage: change_HMMER_output_names_to_combine_GTDB_and_NCBI.py [dictionary] [.faa file to change] 

import sys
import csv
import ast
import collections
import re

#save dictionary
with open(str(sys.argv[1]), 'r') as dict_file:
	dict_annot = ast.literal_eval(dict_file.read())

with open(str(sys.argv[2]), 'r') as seq_file:
	with open(str(sys.argv[2].split('.faa')[0]+'_notes_from_script.txt'), 'w') as notes_file:
		notes_file.write('trailing_characters_removed,NCBI_accession'+'\n')
		with open(str(sys.argv[2].split('.faa')[0]+'_renamed.faa'), 'w') as renamed_file:
			with open(str(sys.argv[2].split('.faa')[0]+'_tax_table.csv'), 'w') as tax_table: #table to store GTDB and NCBI names, to use with left_join in R to attach GTDB metadata (taxonomy, etc.)
				tax_table.write('GTDB_accession,NCBI_accession'+'\n')
											
				seq_counter = 0
				bad_counter = 0
				
				for line in seq_file:
					if line.strip('\r\n').startswith('>'):
						NCBI_accession = line.strip('\r\n').split(".")[0].split('>')[1] #save NCBI annotation (first item in sequence name) and remove initial >
						#print(NCBI_accession)
						try:
							renamed_file.write('>'+dict_annot[NCBI_accession]+' ('+line.strip('\r\n').split(' #')[0].split('>')[1]+')'+'\n') #what to do if the key IS in the dictionary
							tax_table.write(dict_annot[NCBI_accession]+','+line.strip('\r\n').split(' #')[0].split('>')[1]+'\n')
							seq_counter += 1
							#print(NCBI_accession)
						except KeyError:
							notes_file.write('Being difficult: '+str(NCBI_accession)+'\n')
							bad_counter += 1
					else:
						renamed_file.write(line) #write the sequence data			

print('I am done! Go grow yer trees!')
print(('%i sequences were renamed.'%(seq_counter)))
print(('%i sequences were difficult....'%(bad_counter)))
			
