# load all the libraries
library(Rcpp)
library(kernlab)
library(seqinr)
library(ROCR)
library(BSgenome)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
library(GenomeInfoDb)
library(IRanges)
library(BiocGenerics)
library(S4Vectors)
library(gkmSVM)

# remember to load all the dependencies or won't work

# define positive and negative sets

posfn='~/Documents/ATAC_seq/gkmSVM/test/2019_01_30_hu_high.fasta'
negfn='~/Documents/ATAC_seq/gkmSVM/test/hu_high_null.fasta'

# define kernel output and prefix for output

kernelfn='~/Documents/ATAC_seq/gkmSVM/test/hu_high_kernel.test'
svmfnprfx='~/Documents/ATAC_seq/gkmSVM/test/hu_high_test'
outfn='~/Documents/ATAC_seq/gkmSVM/test/hu_high_weights.txt'

# calculate kernel matrix

gkmsvm_kernel(posfn, negfn, kernelfn)

# train SVM and save pdf

# can use the kernel matrix to train other data sets; I am using the same sets to train because I want to find out what the motifs are that separate posfn from negfn

gkmsvm_trainCV(kernelfn, posfn, negfn, svmfnprfx, outputCVpredfn='~/Documents/ATAC_seq/gkmSVM/test/hu_high_test_cvpred.out',outputROCfn='~/Documents/ATAC_seq/gkmSVM/test/hu_high_test_roc.out',nCV=10)

# use non-redundant 10mer to see which kmers are predictive in the classifier
seqfile = '~/nr6mers.fa'
gkmsvm_classify(seqfile, svmfnprfx, outfn,L=6)