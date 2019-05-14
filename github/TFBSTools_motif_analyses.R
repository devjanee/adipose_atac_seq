# TFBSTools analysis on NFIA and PPARG motifs

# load libraries

library(TFBSTools)
library(JASPAR2016)
library("Biostrings")
library(data.table)

# read in fasta files, each ortholog's fasta header is a number [e.g. >1] so that pairwise comparisons are easy

common_hg19<-readDNAStringSet("~/Documents/ATAC_seq/2018_08_01_common_matched_size_numbered_hg19.fasta")

HH<-readDNAStringSet("~/Downloads/2018_09_18_hu_high_seq.fasta")
HL<-readDNAStringSet("~/Downloads/2018_09_18_hu_low_seq.fasta")
CL<-readDNAStringSet("~/Downloads/2018_09_18_ch_low_seq.fasta")
CH<-readDNAStringSet("~/Downloads/2018_09_18_ch_high_seq.fasta")

#create position weight matrix class for NFIA

NFIA <- PWMatrix(ID="M5660_1.02", name="NFIA", 
                matrixClass="SMAD", strand="+",
                bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                tags=list(family="NFI", species="10090",
                          tax_group="vertebrates",medline="7592839", 
                          type="SELEX",ACC="P53762", pazar_tf_id="TF0000003",
                          TFBSshape_ID="11", TFencyclopedia_ID="580"),
                profileMatrix=matrix(c(0.014729848,	0,	0.000134409,	0.000168424,	0.022730452,	0.698634129,	0.22996995,	0.246710884,	0.231629107,	0.132193146,	0.000358715,	0.000363544,	0.000288871,	0.99308857,	0.530965169,
0.03970791,	0.000511934,	0,	3.37E-05,	0.976779969,	0.062830086,	0.363455342,	0.250407742,	0.152878793,	0.097794638,	0.000260884,	0.999127495,	0.999494475,	0.006312173,	0.002400716,
0.004011533,	0.005759263,	0.999697581,	0.999629467,	6.99E-05,	0.099508286,	0.175663444,	0.259069987,	0.387159185,	0.066160869,	0.972085439,	0.000181772,	7.22E-05,	0.000399505,	0.456705729,
0.941550708,	0.993728803,	0.000168011,	0.000168424,	0.000419639,	0.1390275,	0.230911263,	0.243811388,	0.228332915,	0.703851347,	0.027294962,	0.000327189,	0.000144436,	0.000199752,	0.009928385),
                                     byrow=TRUE, nrow=4,
                                     dimnames=list(c("A", "C", "G", "T"))
                                     )
                )

siteset_NFIA_common<-searchSeq(NFIA, common_hg19, min.score="80%",strand="*")
siteset_NFIA_common.df<-as.data.frame(writeGFF3(siteset_NFIA_common,scoreType="relative"))
siteset_NFIA_common.DT<-data.table(siteset_NFIA_common.df)
siteset_NFIA_common_max<-siteset_NFIA_common.DT[,.SD[which.max(score)],by=list(seqname)]
siteset_NFIA_common_sum<-siteset_NFIA_common.DT[,sum(score),by=list(seqname)]

siteset_NFIA_HH<-searchSeq(NFIA, HH, min.score="80%",strand="*")
siteset_NFIA_HH.df<-as.data.frame(writeGFF3(siteset_NFIA_HH,scoreType="relative"))
siteset_NFIA_HH.DT<-data.table(siteset_NFIA_HH.df)
siteset_NFIA_HH_max<-siteset_NFIA_HH.DT[,.SD[which.max(score)],by=list(seqname)]
siteset_NFIA_HH_sum<-siteset_NFIA_HH.DT[,sum(score),by=list(seqname)]

siteset_NFIA_HL<-searchSeq(NFIA, HL, min.score="80%",strand="*")
siteset_NFIA_HL.df<-as.data.frame(writeGFF3(siteset_NFIA_HL,scoreType="relative"))
siteset_NFIA_HL.DT<-data.table(siteset_NFIA_HL.df)
siteset_NFIA_HL_max<-siteset_NFIA_HL.DT[,.SD[which.max(score)],by=list(seqname)]
siteset_NFIA_HL_sum<-siteset_NFIA_HL.DT[,sum(score),by=list(seqname)]

siteset_NFIA_CH<-searchSeq(NFIA, CH, min.score="80%",strand="*")
siteset_NFIA_CH.df<-as.data.frame(writeGFF3(siteset_NFIA_CH,scoreType="relative"))
siteset_NFIA_CH.DT<-data.table(siteset_NFIA_CH.df)
siteset_NFIA_CH_max<-siteset_NFIA_CH.DT[,.SD[which.max(score)],by=list(seqname)]
siteset_NFIA_CH_sum<-siteset_NFIA_CH.DT[,sum(score),by=list(seqname)]

siteset_NFIA_CL<-searchSeq(NFIA, CL, min.score="80%",strand="*")
siteset_NFIA_CL.df<-as.data.frame(writeGFF3(siteset_NFIA_CL,scoreType="relative"))
siteset_NFIA_CL.DT<-data.table(siteset_NFIA_CL.df)
siteset_NFIA_CL_max<-siteset_NFIA_CL.DT[,.SD[which.max(score)],by=list(seqname)]
siteset_NFIA_CL_sum<-siteset_NFIA_CL.DT[,sum(score),by=list(seqname)]

PPARG <- PWMatrix(ID="M1894_1.02", name="PPARG", 
                matrixClass="", strand="+",
                bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                tags=list(family="", species="10090",
                          tax_group="vertebrates",medline="7592839", 
                          type="SELEX",ACC="P53762", pazar_tf_id="TF0000003",
                          TFBSshape_ID="11", TFencyclopedia_ID="580"),
                profileMatrix=matrix(c(0.032064549,	0.031726753,	0.761137105,	0.049094166,	5.50E-03,	0.002189421,	0.380008503,	0.023775278,	0.138228695,
0.017912184,	0.13184225,	0.161989512,	0.939534714,	0.965457597,	0.144833831,	0.139504921,	0.264332904,	0.100492359,
0.015392447,	0.815621705,	0.052607388,	0.001642066,	0,	0.000547355,	0.039856882,	0.114977093,	0.412709934,
0.93463082,	0.020809292,	0.024265994,	0.009729053,	0.029044236,	0.852429393,	0.440629694,	0.596914725,	0.348569013),
                                     byrow=TRUE, nrow=4,
                                     dimnames=list(c("A", "C", "G", "T"))
                                     )
                )

siteset_PPARG_common<-searchSeq(PPARG, common_hg19, min.score="80%",strand="*")
siteset_PPARG_common.df<-as.data.frame(writeGFF3(siteset_PPARG_common,scoreType="relative"))
siteset_PPARG_common.DT<-data.table(siteset_PPARG_common.df)
siteset_PPARG_common_max<-siteset_PPARG_common.DT[,.SD[which.max(score)],by=list(seqname)]
siteset_PPARG_common_sum<-siteset_PPARG_common.DT[,sum(score),by=list(seqname)]

siteset_PPARG_HH<-searchSeq(PPARG, HH, min.score="80%",strand="*")
siteset_PPARG_HH.df<-as.data.frame(writeGFF3(siteset_PPARG_HH,scoreType="relative"))
siteset_PPARG_HH.DT<-data.table(siteset_PPARG_HH.df)
siteset_PPARG_HH_max<-siteset_PPARG_HH.DT[,.SD[which.max(score)],by=list(seqname)]
siteset_PPARG_HH_sum<-siteset_PPARG_HH.DT[,sum(score),by=list(seqname)]

siteset_PPARG_HL<-searchSeq(PPARG, HL, min.score="80%",strand="*")
siteset_PPARG_HL.df<-as.data.frame(writeGFF3(siteset_PPARG_HL,scoreType="relative"))
siteset_PPARG_HL.DT<-data.table(siteset_PPARG_HL.df)
siteset_PPARG_HL_max<-siteset_PPARG_HL.DT[,.SD[which.max(score)],by=list(seqname)]
siteset_PPARG_HL_sum<-siteset_PPARG_HL.DT[,sum(score),by=list(seqname)]

siteset_PPARG_CH<-searchSeq(PPARG, CH, min.score="80%",strand="*")
siteset_PPARG_CH.df<-as.data.frame(writeGFF3(siteset_PPARG_CH,scoreType="relative"))
siteset_PPARG_CH.DT<-data.table(siteset_PPARG_CH.df)
siteset_PPARG_CH_max<-siteset_PPARG_CH.DT[,.SD[which.max(score)],by=list(seqname)]
siteset_PPARG_CH_sum<-siteset_PPARG_CH.DT[,sum(score),by=list(seqname)]

siteset_PPARG_CL<-searchSeq(PPARG, CL, min.score="80%",strand="*")
siteset_PPARG_CL.df<-as.data.frame(writeGFF3(siteset_PPARG_CL,scoreType="relative"))
siteset_PPARG_CL.DT<-data.table(siteset_PPARG_CL.df)
siteset_PPARG_CL_max<-siteset_PPARG_CL.DT[,.SD[which.max(score)],by=list(seqname)]
siteset_PPARG_CL_sum<-siteset_PPARG_CL.DT[,sum(score),by=list(seqname)]


