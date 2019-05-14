#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=devjanee.swain.lenz@duke.edu
#SBATCH --array=1-14
#SBATCH --mem=10G

S1=../final/S"$SLURM_ARRAY_TASK_ID"_hg19_final.bed
S2=../final/S"$SLURM_ARRAY_TASK_ID".bam

module load bedtools2

bedtools bedtobam -i $S1 -g ~/reference_files/hg19.chrom.sizes > $S2
