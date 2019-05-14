#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=devjanee.swain.lenz@duke.edu
#SBATCH --mem=10G

# define locations of sorted and indexed bam files (S1), bed file (S2) and output file (S3)

S1= ../*sorted.bam
S2= ../all_peaks_merged.bed
S3= ../2018_02_12_all_peaks_count.bed 

module load bedtools2

bedtools multicov -bams $S1 -bed $S2  > $S3 
