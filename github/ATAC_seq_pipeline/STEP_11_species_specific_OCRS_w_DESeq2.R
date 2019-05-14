### May 12, 2018
### Devjanee Swain Lenz, PhD
### Laboratory of Greg Wray, Duke University

#########################################################
#														#
# Using DESeq2 to call differential peaks b/t species	#
#														#
#########################################################


### filter out rows with zero counts #####

counts <- read.table("~/2018_02_12_all_peaks_count.bed ")
cts <- counts[,4:15]
row_sub = apply(cts, 1, function(row) all(row !=0 ))
counts.zero.filter<- cts[row_sub,]
write.table(counts.zero.filter,"~/mydata.text",sep="\t",quotes=F,colnames=F,rownames=F)

### load DESeq2 library ###

library("DESeq2")

# define counts table ("cts") ; use peak names as rownames because filtering makes it easier to spot check than numbers of rows

# read in bed file with chr(column 1), start(column 2), end(column 3), peak_name(column 4), raw counts for samples(columns 5-24, in this case)

bed <- read.table("~/mydata.txt")

# make count table with peak names

cts <- bed[,5:24]
rownames(cts)<-bed$name

# read in variables
# coldata explained in DESeq2 vignettes - in this case, rownames are sample names and the one column is the variable "species"

coldata<-read.table("~/coldata.txt")

# run DESeq for the species model

dds <- DESeqDataSetFromMatrix(countData = counts.zero.filter, colData = coldata.1, design = ~ species)
dds <- DEseq(dds)

##### PCA of dds

r.dds<-rlog(dds)
plotPCA(r.dds,intgroup="species")

##### make a heat map

library(pheatmap)
pheatmap(cor(assay(dds.all.3)), border_color=NA)

# compare to a null model

null.v.species <- nbinomWaldTest(dds)

# only take regions that are significant for a species effect for both DEseq and nbinomWaldTest results

res.null.v.species <- res(null.v.species)
null.v.sp.filter <- (which(res.null.v.species$padj<0.05))

res<- results(dds)
res.sig<-(which(res.sig$padj<0.05))

length(intersect(rownames(res.sig),rownames(null.v.sp.filter)))

# same - continue using DEseq(dds) results

# take the contrasts between each species pair

res.2 <- results(dds, contrast=c("species","hu","ch"))
res.3 <- results(dds, contrast=c("species","hu","mac"))
res.4 <- results(dds, contrast=c("species","ch","mac"))

# make vectors of significantly differenential regions between pairs of species

Hu_Ch_diff <- which(res.2$padj<0.05)
Hu_Mac_diff <-which(res.3$padj<0.05)
Ch_Mac_diff <-which(res.4$padj<0.05)

# get regions specific to one species

Hu_spec <- cts.norm[intersect(Hu_Ch_diff,Hu_Mac_diff),]
Ch_spec <- cts.norm[intersect(Hu_Ch_diff,Ch_Mac_diff),]

cts.norm <- counts(dds, normalized=TRUE)
 
# human and chimpanzee different from each other and different from macaque 
# make data frame, and take averages of H,C, and M and figure out increases, decreases

CHLCA_ch_dif_names<-intersect(rownames(Hu_spec),rownames(Ch_spec))
CHLCA_ch_dif_df<- cts.norm[CHLCA_ch_dif_names,]

CHLCA_ch_dif_df<- as.data.frame(cts.norm[CHLCA_ch_dif_names,])

# means for each species
CHLCA_ch_dif_df$H<-rowMeans(subset(CHLCA_ch_dif_df, select = c(H1.1,H2.1,H2.2,H3.1,H1.2,H1.3)), na.rm = TRUE)
CHLCA_ch_dif_df$C<-rowMeans(subset(CHLCA_ch_dif_df, select = c(C1.1,C2.1,C2.2,C1.2,C2.3)), na.rm = TRUE)
CHLCA_ch_dif_df$M<-rowMeans(subset(CHLCA_ch_dif_df, select = c(M1.1,M1.2)), na.rm = TRUE)

# logical statements, which species is higher or lower
CHLCA_ch_dif_df$H_C<-CHLCA_ch_dif_df$H<CHLCA_ch_dif_df$C
CHLCA_ch_dif_df$H_M<-CHLCA_ch_dif_df$H<CHLCA_ch_dif_df$M
CHLCA_ch_dif_df$C_M<-CHLCA_ch_dif_df$C<CHLCA_ch_dif_df$M

# categorize based on logical statements
# H<C<M = H2
# H>C>M = H4
# C<H<M = C2
# C>H>M = C4
# if macaque is in the middle = "x" because we can't figure out what's going on with these regions

CHLCA_ch_dif_df$CHLCA<-apply(CHLCA_ch_dif_df,1,function(x) if(x[17]!=F&x[18]!=F&x[19]!=F){x[20]="H2"}else if(x[17]==F&x[18]==F&x[19]==F){x[20]="H4"}else if(x[17]==F&x[18]!=F&x[19]!=F){x[20]="C2"}else if(x[17]!=F&x[18]==F&x[19]==F){x[20]="C4"}else{x[20]="x"})

# filter out H2 & C4 from hu_spec_low
# filter out H4 & C2 from hu spec hi
# filter out C2 $ H4 from ch spec low
# filter out H4 & C2 from ch spec hi 
# double check to get rid of mac in middle

H2_names<-rownames(CHLCA_ch_dif_df[which(CHLCA_ch_dif_df$CHLCA=="H2"),])
C4_names<-rownames(CHLCA_ch_dif_df[which(CHLCA_ch_dif_df$CHLCA=="C4"),])
Hu_low_2<-Hu_low[!(rownames(Hu_low) %in% H2_names)&!(rownames(Hu_low) %in% C4_names),]

Hu_low_2$H<-rowMeans(subset(Hu_low_2, select = c(H1.1,H2.1,H2.2,H3.1,H1.2,H1.3)), na.rm = TRUE)
Hu_low_2$C<-rowMeans(subset(Hu_low_2, select = c(C1.1,C2.1,C2.2,C1.2,C2.3)), na.rm = TRUE)
Hu_low_2$M<-rowMeans(subset(Hu_low_2, select = c(M1.1,M1.2)), na.rm = TRUE)

Hu_low_2$H_C<-Hu_low_2$H<Hu_low_2$C
Hu_low_2$H_M<-Hu_low_2$H<Hu_low_2$M
Hu_low_2$C_M<-Hu_low_2$C<Hu_low_2$M

Hu_low_2$CHLCA<-apply(Hu_low_2,1,function(x) if(x[17]!=F&x[18]!=F&x[19]!=F){x[20]="H1"}else if(x[17]==F&x[18]==F&x[19]==F){x[20]="H4"}else if(x[17]==F&x[18]!=F&x[19]!=F){x[20]="C3"}else if(x[17]!=F&x[18]==F&x[19]==F){x[20]="C4"}else{x[20]="x"})

# length of hu low where ch and mac are indistinguishable
length(which(Hu_low_2$CHLCA=="H1"))
# [1] 680
# oh this number seems close to the TPM number

# 190 are mac in the middle 

Hu_high_2<-Hu_high[!(rownames(Hu_high) %in% H4_names)&!(rownames(Hu_high) %in% C2_names),]
Hu_high_2<-as.data.frame(Hu_high_2)

Hu_high_2$H<-rowMeans(subset(Hu_high_2, select = c(H1.1,H2.1,H2.2,H3.1,H1.2,H1.3)), na.rm = TRUE)
Hu_high_2$C<-rowMeans(subset(Hu_high_2, select = c(C1.1,C2.1,C2.2,C1.2,C2.3)), na.rm = TRUE)
Hu_high_2$M<-rowMeans(subset(Hu_high_2, select = c(M1.1,M1.2)), na.rm = TRUE)

Hu_high_2$H_C<-Hu_high_2$H<Hu_high_2$C
Hu_high_2$H_M<-Hu_high_2$H<Hu_high_2$M
Hu_high_2$C_M<-Hu_high_2$C<Hu_high_2$M

Hu_high_2$CHLCA<-apply(Hu_high_2,1,function(x) if(x[17]!=F&x[18]!=F&x[19]!=F){x[20]="H2"}else if(x[17]==F&x[18]==F&x[19]==F){x[20]="H3"}else if(x[17]==F&x[18]!=F&x[19]!=F){x[20]="C2"}else if(x[17]!=F&x[18]==F&x[19]==F){x[20]="C4"}else{x[20]="x"})

length(which(Hu_high_2$CHLCA=="H3"))
# [1] 599
# 599 that are Hu high while ch and mac are indistinguishable
# 151 mac in the middle

Ch_high_2<-Ch_high[!(rownames(Ch_high) %in% H2_names)&!(rownames(Ch_high) %in% C4_names),]

Ch_high_2<-as.data.frame(Ch_high_2)

Ch_high_2$H<-rowMeans(subset(Ch_high_2, select = c(H1.1,H2.1,H2.2,H3.1,H1.2,H1.3)), na.rm = TRUE)
Ch_high_2$C<-rowMeans(subset(Ch_high_2, select = c(C1.1,C2.1,C2.2,C1.2,C2.3)), na.rm = TRUE)
Ch_high_2$M<-rowMeans(subset(Ch_high_2, select = c(M1.1,M1.2)), na.rm = TRUE)

Ch_high_2$H_C<-Ch_high_2$H<Ch_high_2$C
Ch_high_2$H_M<-Ch_high_2$H<Ch_high_2$M
Ch_high_2$C_M<-Ch_high_2$C<Ch_high_2$M
Ch_high_2$CHLCA<-""

Ch_high_2$CHLCA<-apply(Ch_high_2,1,function(x) if(x[17]!=F&x[18]!=F&x[19]!=F){x[20]="H2"}else if(x[17]==F&x[18]==F&x[19]==F){x[20]="H3"}else if(x[17]==F&x[18]!=F&x[19]!=F){x[20]="C2"}else if(x[17]!=F&x[18]==F&x[19]==F){x[20]="C3"}else{x[20]="x"})

length(which(Ch_high_2=="C3"))
#[1] 863

# 176 are mac in the middle

Ch_low_2<-Ch_low[!(rownames(Ch_low) %in% H4_names)&!(rownames(Ch_low) %in% C2_names),]

Ch_low_2<-as.data.frame(Ch_low_2)

Ch_low_2$H<-rowMeans(subset(Ch_low_2, select = c(H1.1,H2.1,H2.2,H3.1,H1.2,H1.3)), na.rm = TRUE)
Ch_low_2$C<-rowMeans(subset(Ch_low_2, select = c(C1.1,C2.1,C2.2,C1.2,C2.3)), na.rm = TRUE)
Ch_low_2$M<-rowMeans(subset(Ch_low_2, select = c(M1.1,M1.2)), na.rm = TRUE)

Ch_low_2$H_C<-Ch_low_2$H<Ch_low_2$C
Ch_low_2$H_M<-Ch_low_2$H<Ch_low_2$M
Ch_low_2$C_M<-Ch_low_2$C<Ch_low_2$M

Ch_low_2$CHLCA<-""

Ch_low_2$CHLCA<-apply(Ch_low_2,1,function(x) if(x[17]!=F&x[18]!=F&x[19]!=F){x[20]="H2"}else if(x[17]==F&x[18]==F&x[19]==F){x[20]="H3"}else if(x[17]==F&x[18]!=F&x[19]!=F){x[20]="C1"}else if(x[17]!=F&x[18]==F&x[19]==F){x[20]="C3"}else{x[20]="x"})

# 422 are chimpanzee low where human and macaque are same
# 82 mac in the middle

HLCH_intersect<-(intersect(rownames(Hu_low_2),rownames(Ch_high_2)))
# 2 are different between H, C and M where M is in the middle; otherwise M is in the middle but not significantly different than C

HHCL_intersect<-(intersect(rownames(Hu_high_2),rownames(Ch_low_2)))
# 5 are different between H, C and M where M is in the middle; otherwise M is in the middle but not significantly differentt than H

# Human decreased and Chimpanzee increased have overlaps so need to filter out the things that are unique

# Human increase and chimpanzee decrease use HHCL_intersect
# Human decrease and chimpanzee increase use HLCH_intersect

HuDecreased<-Hu_high_2[!(rownames(Hu_low_2) %in% HLCH_intersect),]
HuIncreased<-Hu_high_2[!(rownames(Hu_high_2) %in% HHCL_intersect),]
ChDecreased<-Hu_high_2[!(rownames(Ch_low_2) %in% HHCL_intersect),]
ChIncreased<-Hu_high_2[!(rownames(Ch_low_2) %in% HLCH_intersect),]


