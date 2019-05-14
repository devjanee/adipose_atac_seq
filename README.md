# adipose_atac_seq

ATAC_seq_pipeline - This folder contains the shell scripts used to process raw sequencing data into a counts table and is to be run in order of steps 

- Step 1. Mapping with bowtie2 (note, use the proper genome)
- Steps 2 - 3. Removing duplicate PCR reads 
- Step 4. Creating bed files for LiftOver
- Step 5. LiftOver command (note, use proper chain [e.g. going from chimpanzee to human genome will need panTro4ToHg19.over.chain] and if reciprocal liftOver is required, then script will be ran multiple times)
- Step 6. Peak calling script using macs2
- Step 7. Concatenating and merging peaks together from all species
- Step 8 - 9. Creating and indexing bam files for count tables
- Step 10. Creating count table (rows are peaks and columns are raw counts for each sample)
- Step 11. R script to perform DESeq2 analyses

###########################################

gkmSVM_null_v_randomGenome.R - R commands to determine if a set of sequences is similar to 100 sets of 1100 random genomic sequences using the R package "gkmSVM" (http://www.beerlab.org/gkmsvm/)

gkmSVM_OCRs_v_closestNull.R - R commands to classify OCRs to the closest common sites using R package "gkmSVM" (http://www.beerlab.org/gkmsvm/)

TFBSTools_motif_analyses.R - R commands to score motifs in a given set of sequences using  R package "TFBSTools" (https://bioconductor.org/packages/release/bioc/vignettes/TFBSTools/inst/doc/TFBSTools.html) 
