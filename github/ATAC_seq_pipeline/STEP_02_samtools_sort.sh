#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=devjanee.swain.lenz@duke.edu
#SBATCH --array=1-14
#SBATCH --mem=10G

# define locations of input files (S1), output files (S2) and temporary files (S3)

S1=../S"$SLURM_ARRAY_TASK_ID"_mapped.sam
S2=../S"$SLURM_ARRAY_TASK_ID"_sorted.sam
S3=../S"$SLURM_ARRAY_TASK_ID"_temp

module load samtools

samtools sort -T $S3 -n -O SAM $S1 -o $S2

