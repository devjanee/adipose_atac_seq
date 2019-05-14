#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=devjanee.swain.lenz@duke.edu
#SBATCH --mem=10G

# define locations of .bed files (S1) and output suffix (S2) 

S1=../*.bed
S2=../human_all

module load macs2

macs2 callpeak -B -f BED -q 0.01 --nomodel --shift -100 --ext 200 -t $S1 -n $S2 
