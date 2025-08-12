#PBS -N filt_NifD
#PBS -l procs=1
#PBS -l mem=15gb
#PBS -q una-new
#PBS -M jdamashe@hamilton.edu
#PBS -m ea
#PBS -j oe
#PBS -l walltime=2:00:00
#PBS -r n

cd $PBS_O_WORKDIR

module load python/anaconda-damashek

filterbyname.sh in=/usr/local/damashek/amoeba/databases/gtdb/protein_fna_reps/gtdb_proteins_reps_r214_genes.fna out=GTDB_r214_NifD.fna names=GTDB_nifD_gene_accessions.txt include=t --ignorejunk
