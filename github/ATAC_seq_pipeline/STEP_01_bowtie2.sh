#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=devjanee.swain.lenz@duke.edu
#SBATCH --array=4
#SBATCH --mem=10G

### define locations of read 1 fasta file (S1), read 2 fast file (S2), output file (S3), and genome (S4)

S1=../*S"$SLURM_ARRAY_TASK_ID"_L006_R1_001.fastq
S2=../*S"$SLURM_ARRAY_TASK_ID"_L006_R2_001.fastq 
S3=../S"$SLURM_ARRAY_TASK_ID"_mapped.sam
S4=../../LENZ_3856_170120B1/hg19

#### working on cluster, so need to load bowtie2 before running

module load bowtie2 

bowtie2 -p 8 -t --phred33 --local -x $S4 -1 $S1 -2 $S2 -S $S3
