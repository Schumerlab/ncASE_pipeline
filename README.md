This is a pipeline for determining allele specific expression using two colinear reference genomes and non-competitive mapping (ncASE). 

Further updates will extend to non-colinear cases and other mappers.

Dependencies:

samtools v~1.9
bcftools v~1.9
bwa
R
bedtools

All dependencies can be globally installed or set in the configuration file

All parameters are set in the configuration file, see cfg examples and below

################
###PARAMETERS
#################

#full path to the parental genomes
#Note: make sure to unzip them!
genome1=
genome2=

#genome index or not, 0 to index, 1 to skip
skip_genome_index=

#mapping program to use, options are bwa or star
mapping_program=

#path to ncASE_pipeline scripts
path_to_ncASE=./

#path to mapper if not globally installed
path_to_mapper=

#path to samtools if not globally installed
path_to_samtools=

#path to bedtools if not globally installed
#Note: give only the path to the parent bedtools folder, e.g. /usrbin/bedtools2 not /usrbin/bedtools2/bin/intersectBed
path_to_bedtools=

#path and name of gtf file, leave empty if you are using a transcriptome 
gtf_file=

#allow sites where the alternate allele count is 0 (100% ASE). 0 for no, 1 for yes
allow_zero_counts=0

#allelic bias tolerance proportion. For example 0.05 would allow sites where the inferred counts for parent 1 in the two genome references is up to 5% different
allelic_bias_threshold=

#read type, options are SE or PE
read_type=

#read length
read_length=

#list containing full path to the reads, with PE reads tab separated on the same line
read_list=

#number of individuals to run in each job
number_indiv_per_job=

#Aims file listing the chromosome, position, parent 1 base, and parent 2 base,  see swordtail examples for format
provide_AIMs=

#slurm commands for each job
slurm_command_map=
slurm_command_variant_call=
slurm_command_ncASE=
