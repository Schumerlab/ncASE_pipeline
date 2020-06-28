# ncASE 

ncASE is a pipeline for determining allele specific expression using two colinear reference genomes and non-competitive mapping (ncASE). 

Further updates will extend to non-colinear cases and other mappers.

## Dependencies

samtools v1.9
bcftools v1.9
bwa
R
bedtools

All dependencies can be globally installed or set in the configuration file.

This program runs batch jobs through Slurm workload manager, found on many HPC clusters.

## Usage

```bash 
perl ncASE_parallel.pl configuration.cfg
```

## Parameters

All parameters are set in the configuration file, see .cfg examples and below.

The following parameters should appear in your configuration file, one parameter per line:

Full path to the parental genomes (make sure to unzip them!)
```bash
genome1=
genome2=
```

Genome index or not: 0 to index, 1 to skip
```bash
skip_genome_index=
```

Mapping program to use: bwa or star
```bash
mapping_program=
```
Path to ncASE_pipeline scripts
```bash
path_to_ncASE=./
```

Path to mapper (if not globally installed)
```bash
path_to_mapper=
```

Path to samtools (if not globally installed)
```bash
path_to_samtools=
```

Path to bedtools (if not globally installed)
Note: give only the path to the parent bedtools folder, e.g. /usrbin/bedtools2 not /usrbin/bedtools2/bin/intersectBed
```bash
path_to_bedtools=
```

Path and name of gtf annotation file, leave empty if you are using a transcriptome 
```bash
gtf_file=
```

Allow sites where the alternate allele count is 0 (100% ASE)? 0 for no, 1 for yes
```bash
allow_zero_counts=0
```

Allelic bias tolerance proportion. For example, 0.05 would allow sites where the inferred counts for parent 1 in the two genome references is up to 5% different
```bash
allelic_bias_threshold=
```

Read type: SE (single) or PE (paired)
```bash
read_type=
```

Read length
```bash
read_length=
```

List containing full path to the reads. PE reads should be tab-separated on the same line
```bash
read_list=
```

Number of individuals to run in each job sent to job scheduler
```bash
number_indiv_per_job=
```

AIMs file listing the chromosome, position, parent 1 base, and parent 2 base (see swordtail examples for format)
```bash
provide_AIMs=
```

Slurm commands for each job. If you are using another job scheduler, edit these commands so that they are compatible with your scheduler's syntax.
```bash
slurm_command_map=
slurm_command_variant_call=
slurm_command_ncASE=
```
