#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=devjanee.swain.lenz@duke.edu
#SBATCH --array=1-14
#SBATCH --mem=10G

# define location of input files (S1) and output files (S2)
S1=../S"$SLURM_ARRAY_TASK_ID"_rmdup.bam
S2=../S"$SLURM_ARRAY_TASK_ID".bed

module load bedtools2

bedtools bamtobed -i $S1 > $S2

