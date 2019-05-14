#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=devjanee.swain.lenz@duke.edu
#SBATCH --array=1-14
#SBATCH --mem=10G

# define location of input file (S1)

S1=../final/S"$SLURM_ARRAY_TASK_ID"_sorted.bam

module load samtools

samtools index $S1

