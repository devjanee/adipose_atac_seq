#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=devjanee.swain.lenz@duke.edu
#SBATCH --array=4,9-14
#SBATCH --mem=10G

# define location of input file (S1), output file (S2), unmapped file (S3), chain (S4)

S1=../S"$SLURM_ARRAY_TASK_ID"pT4.bed
S2=../S"$SLURM_ARRAY_TASK_ID"_hg19_final.bed
S3=../S"$SLURM_ARRAY_TASK_ID"_unmapped_final.bed
S4=../panTro4ToHg19.over.chain

~/bin/liftOver $S1 $S4 $S2 $S3
