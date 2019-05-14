#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=devjanee.swain.lenz@duke.edu
#SBATCH --array=14
#SBATCH --mem=10G

# define locations of input files (S2) and output files (S3)

S2=../S"$SLURM_ARRAY_TASK_ID"_sorted.sam
S3=../S"$SLURM_ARRAY_TASK_ID"_rmdup.bam

module load samtools

samtools rmdup --output-fmt BAM -S $S2 $S3

