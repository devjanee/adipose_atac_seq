###### Script to test if null gkm sets are different than rest of genome ######

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

for (x in 1:100){
  	negfn=paste('~/Documents/ATAC_seq/gkmSVM/chimpanzee/random/random_', x, ".fasta", sep="")
	kernelfn=paste('~/Documents/ATAC_seq/gkmSVM/chimpanzee/ch_high_v_random/ch_high_v_random_',x,'_kernel',sep="")
	svmfnprfx=paste('~/Documents/ATAC_seq/gkmSVM/chimpanzee/ch_high_v_random/ch_high_v_random_',x,sep="")
	outputPDFfn=paste('~/Documents/ATAC_seq/gkmSVM/chimpanzee/ch_high_v_random/ch_high_v_random_AUC_',x,".pdf", sep="")
	gkmsvm_kernel(posfn, negfn, kernelfn)
	auc <- gkmsvm_trainCV(kernelfn, posfn, negfn, svmfnprfx,outputPDFfn=outputPDFfn)
	auc.df[x,1]<-auc$aucss[2]
	auc.df[x,2]<-auc$aucss[3]
	outfn=paste('~/Documents/ATAC_seq/gkmSVM/chimpanzee/ch_high_v_random/ch_high_v_random_',x,".txt",sep="")
}

for (x in 1:100){
	negfn=paste('~/Documents/ATAC_seq/gkmSVM/chimpanzee/random/random_', x, ".fasta", sep="")
	kernelfn=paste('~/Documents/ATAC_seq/gkmSVM/chimpanzee/ch_high_v_random/ch_high_v_random_',x,'_kernel',sep="")
	svmfnprfx=paste('~/Documents/ATAC_seq/gkmSVM/chimpanzee/ch_high_v_random/ch_high_v_random_',x,sep="")
	outputPDFfn=paste('~/Documents/ATAC_seq/gkmSVM/chimpanzee/ch_high_v_random/ch_high_v_random_AUC_',x,".pdf", sep="")

auc <- gkmsvm_trainCV(kernelfn, posfn, negfn, svmfnprfx,outputPDFfn=outputPDFfn)
	auc.df[x,1]<-auc$aucss[2]
	auc.df[x,2]<-auc$aucss[3]
}