#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=devjanee.swain.lenz@duke.edu
#SBATCH --mem=20G

cat *narrowPeak > all_peaks_unsorted.bed

sort -k1,1 -k2,2n all_peaks_unsorted.bed > all_peaks_sorted.bed

module load bedtools2

mergeBed -c 4 -o collapse -i all_peaks_sorted.bed > all_peaks_merged.bed