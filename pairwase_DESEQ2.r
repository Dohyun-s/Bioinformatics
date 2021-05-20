##learn DESEQ 2 for pairwise comparison

library("apeglm")
library("DESeq2")
library("vsn")
library("ggplot2")
library(dplyr)
library("pheatmap")
library("RColorBrewer")
setwd("C:/Users/user/Documents/2021-1/bioinformatic/project")
library(EnhancedVolcano)
GSM465<-read.table("GSM1169465.txt", sep="\t")
colnames(GSM465)=c("ENSEMBL","Gene_name","Count")

GSM466<-read.table("GSM1169466.out",sep="\t")


GSM2=c(1:53700)
GSM2

for (i in 1:20){
  k=1169465+i
  j=paste("GSM", k,".out", sep="")
 
  GSM1=read.table(j,sep="\t")
  GSM2=cbind(GSM2,GSM1[,3])
}

setwd("C:/Users/user/Documents/2021-1/bioinformatic/project/work")

for (i in 1:17){
  k=1169490+i
  j=paste("GSM", k,".out", sep="")
 
  GSM1=read.table(j,sep="\t")
  GSM2=cbind(GSM2,GSM1[,3])
}


GSM_f=cbind(GSM465[,3],GSM2)
dim(GSM_f)

GSM_f2=cbind(GSM_f[,1],GSM_f[,3:39])

write.csv(GSM_f2, file="project1_matrix_fin")

rownames(GSM_f2)=GSM465[,1]

vec=c()
for (m in 1:21){
  k=1169464+m
    
  j=paste("GSM", k, sep="")
  vec<-c(vec,j)
}

for (m in 1:17){
  k=1169490+m
    
  j=paste("GSM", k, sep="")
  vec<-c(vec,j)
}

colnames(GSM_f2)=vec
setwd("C:/Users/user/Documents/2021-1/bioinformatic/project")

coldata<-read.table("metadata.txt",sep="\t")
colnames(coldata)=c("sample","celltype")
##all
dds_f<-DESeqDataSetFromMatrix(countData=GSM_f2, colData=coldata, design=~celltype)
##DN1
setwd("C:/Users/user/Documents/2021-1/bioinformatic/project/plot/DN1_DN2")

GSM_DN1_DN2=cbind(GSM_f2[,26:27],GSM_f2[,35:36])
coldata1=rbind(coldata[26:27,],coldata[35:36,])
dds<-DESeqDataSetFromMatrix(countData=GSM_DN1_DN2, colData=coldata1, design=~celltype)
dds$celltype<-relevel(dds$celltype, ref="DN1")
dds<-DESeq(dds)


resultsNames(dds)
res<-results(dds, name="celltype_DN2_vs_DN1")
res<-lfcShrink(dds, coef="celltype_DN2_vs_DN1", type="apeglm")

descrip=c(rep(NA,nrow(normalized_count)))


normalized_count=counts(dds,normalized=TRUE)
normal=cbind(descrip, normalized_count)
df_n<-cbind(newColName=rownames(normal), normal)
rownames(df_n)<-1:nrow(df_n)
colnames(df_n)=c("NAME","DESCRIPTION",colnames(normalized_count))
write.table(df_n, file="normalized_count_DN1_DN2", quote=FALSE, sep="\t", row.names=FALSE)


plotMA(res, ylim=c(-2,2))
plotDispEsts(dds)
EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange',y='pvalue',title='DN1vs DN2',pCutoff=1e-5,pointSize=3.0,labSize=6.0)


vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$sample, vsd$celltype, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists, col=colors)

plotPCA(vsd, , intgroup="celltype")

res<-na.omit(res)
df = res[res$padj <0.05 & res$lfcSE >1,]

write.table(df, file="sign_gene_DN2vsDN1")

##DN2

setwd("C:/Users/user/Documents/2021-1/bioinformatic/project/plot/DN2_DN3")
GSM_DN2_DN3=cbind(GSM_f2[,27:28],GSM_f2[,36:37])
coldata2=rbind(coldata[27:28,],coldata[36:37,])
dds<-DESeqDataSetFromMatrix(countData=GSM_DN2_DN3, colData=coldata2, design=~celltype)
dds$celltype<-relevel(dds$celltype, ref="DN2")
dds<-DESeq(dds)



resultsNames(dds)
res<-results(dds, name="celltype_DN3_vs_DN2")
res<-lfcShrink(dds, coef="celltype_DN3_vs_DN2", type="apeglm")


descrip=c(rep(NA,nrow(normalized_count)))


normalized_count=counts(dds,normalized=TRUE)
normal=cbind(descrip, normalized_count)
df_n<-cbind(newColName=rownames(normal), normal)
rownames(df_n)<-1:nrow(df_n)
colnames(df_n)=c("NAME","DESCRIPTION",colnames(normalized_count))
write.table(df_n, file="normalized_count_DN2_DN3", quote=FALSE, sep="\t", row.names=FALSE)



plotMA(res, ylim=c(-2,2))
plotDispEsts(dds)
EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange',y='pvalue',title='DN2+ vs DN3',pCutoff=1e-5,pointSize=3.0,labSize=6.0)


vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$sample, vsd$celltype, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists, col=colors)

plotPCA(vsd, intgroup= "celltype")

res<-na.omit(res)
df = res[res$padj <0.05 & res$lfcSE >1,]

write.table(df, file="sign_gene_DN2vsDN3")
##DN3

setwd("C:/Users/user/Documents/2021-1/bioinformatic/project/plot/DN3_DN4")
GSM_DN3_DN4=cbind(GSM_f2[,28:29],GSM_f2[,37:38])
coldata3=rbind(coldata[28:29,],coldata[37:38,])
dds<-DESeqDataSetFromMatrix(countData=GSM_DN3_DN4, colData=coldata3, design=~celltype)
dds$celltype<-relevel(dds$celltype, ref="DN3")
dds<-DESeq(dds)



resultsNames(dds)
res<-results(dds, name="celltype_DN4_vs_DN3")
res<-lfcShrink(dds, coef="celltype_DN4_vs_DN3", type="apeglm")


normalized_count=counts(dds,normalized=TRUE)
normal=cbind(descrip, normalized_count)
df_n<-cbind(newColName=rownames(normal), normal)
rownames(df_n)<-1:nrow(df_n)
colnames(df_n)=c("NAME","DESCRIPTION",colnames(normalized_count))
write.table(df_n, file="normalized_count_DN3_DN4", quote=FALSE, sep="\t", row.names=FALSE)


plotMA(res, ylim=c(-2,2))
plotDispEsts(dds)
EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange',y='pvalue',title='DN3 vs DN4',pCutoff=1e-5,pointSize=3.0,labSize=6.0)


vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$sample, vsd$celltype, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists, col=colors)

plotPCA(vsd, intgroup="celltype")

res<-na.omit(res)
df = res[res$padj <0.05 & res$lfcSE >1,]

write.table(df, file="sign_gene_DN3vsDN4")
##DN4

setwd("C:/Users/user/Documents/2021-1/bioinformatic/project/plot/DN4_DP")
GSM_DN4_DP=cbind(GSM_f2[,22],GSM_f2[,25],GSM_f2[,29],GSM_f2[,31],GSM_f2[,34],GSM_f2[,38])
coldata4=rbind(coldata[22,],coldata[25,],coldata[29,],coldata[31,],coldata[34,],coldata[38,])
dds<-DESeqDataSetFromMatrix(countData=GSM_DN4_DP, colData=coldata4, design=~celltype)
dds$celltype<-relevel(dds$celltype, ref="DN4")
dds<-DESeq(dds)


resultsNames(dds)
res<-results(dds, name="celltype_DP_vs_DN4")
res<-lfcShrink(dds, coef="celltype_DP_vs_DN4", type="apeglm")



normalized_count=counts(dds,normalized=TRUE)
normal=cbind(descrip, normalized_count)
df_n<-cbind(newColName=rownames(normal), normal)
rownames(df_n)<-1:nrow(df_n)
colnames(df_n)=c("NAME","DESCRIPTION","GSM1169491", "GSM1169494", "GSM1169498","GSM1169500", 
"GSM1169503","GSM1169507")
write.table(df_n, file="normalized_count_DN4_DP", quote=FALSE, sep="\t", row.names=FALSE)



plotMA(res, ylim=c(-2,2))
plotDispEsts(dds)
EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange',y='pvalue',title='DP vs DN4',pCutoff=1e-5,pointSize=3.0,labSize=6.0)


vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$sample, vsd$celltype, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists, col=colors)

plotPCA(vsd, intgroup=c("sample", "celltype"))

res<-na.omit(res)
df = res[res$padj <0.05 & res$lfcSE >1,]

write.table(df, file="sign_gene_DN4vsDP")
##DP_CD8

setwd("C:/Users/user/Documents/2021-1/bioinformatic/project/plot/DP_CD8")

GSM_DP_CD8=cbind(GSM_f2[,22],GSM_f2[,24:25],GSM_f2[,31],GSM_f2[,33:34])
coldata6=rbind(coldata[22,],coldata[24:25,],coldata[31,],coldata[33:34,])
dds<-DESeqDataSetFromMatrix(countData=GSM_DP_CD8, colData=coldata6, design=~celltype)
dds$celltype<-relevel(dds$celltype, ref="DP")
dds<-DESeq(dds)


resultsNames(dds)

res<-results(dds, name="celltype_CD8._vs_DP")
res<-lfcShrink(dds, coef="celltype_CD8._vs_DP", type="apeglm")


normalized_count=counts(dds,normalized=TRUE)
normal=cbind(descrip, normalized_count)
df_n<-cbind(newColName=rownames(normal), normal)
rownames(df_n)<-1:nrow(df_n)
colnames(df_n)=c("NAME","DESCRIPTION","GSM1169491","GSM1169493","GSM1169494","GSM1169500","GSM1169502","GSM1169503")
write.table(df_n, file="normalized_count_DP_CD8", quote=FALSE, sep="\t", row.names=FALSE)



plotMA(res, ylim=c(-2,2))
plotDispEsts(dds)
EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange',y='pvalue',title='DP vs CD8',pCutoff=1e-5,pointSize=3.0,labSize=6.0)


vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$sample, vsd$celltype, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists, col=colors)

plotPCA(vsd, intgroup= "celltype")

res<-na.omit(res)
df = res[res$padj <0.05 & res$lfcSE >1,]

write.table(df, file="sign_gene_DPvsCD8+")
##DP_CD4
setwd("C:/Users/user/Documents/2021-1/bioinformatic/project/plot/DP_CD8")

GSM_DP_CD4=cbind(GSM_f2[,22:23],GSM_f2[,25],GSM_f2[,31:32],GSM_f2[,34])
coldata5=rbind(coldata[22:23,],coldata[25,],coldata[31:32,],coldata[34,])
dds<-DESeqDataSetFromMatrix(countData=GSM_DP_CD4, colData=coldata5, design=~celltype)
dds$celltype<-relevel(dds$celltype, ref="DP")
dds<-DESeq(dds)



resultsNames(dds)
res<-results(dds, name="celltype_CD4._vs_DP")
res<-lfcShrink(dds, coef="celltype_CD4._vs_DP", type="apeglm")

normalized_count=counts(dds,normalized=TRUE)
normal=cbind(descrip, normalized_count)
df_n<-cbind(newColName=rownames(normal), normal)
rownames(df_n)<-1:nrow(df_n)
colnames(df_n)=c("NAME","DESCRIPTION","GSM1169491","GSM1169492", "GSM1169494", "GSM1169500","GSM1169501","GSM1169503")
write.table(df_n, file="normalized_count_DP_CD4", quote=FALSE, sep="\t", row.names=FALSE)


plotMA(res, ylim=c(-2,2))
plotDispEsts(dds)
EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange',y='pvalue',title='DP vs CD4',pCutoff=1e-5,pointSize=3.0,labSize=6.0)


vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$sample, vsd$celltype, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists, col=colors)

plotPCA(vsd, intgroup="celltype")

res<-na.omit(res)
df = res[res$padj <0.05 & res$lfcSE >1,]

write.table(df, file="sign_gene_DPvsCD4+")
##th1
cts=cbind(GSM_f2[,1],GSM_f2[,4],GSM_f2[,5],GSM_f2[,6],GSM_f2[,7],GSM_f2[,14],GSM_f2[,18])

coldata_th1<-read.table("metadata1.txt", sep="\t")
colnames(coldata_th1)=c("sample","celltype")

colnames(cts)=c("GSM1169465","GSM1169468","GSM1169469","GSM1169470","GSM1169471","GSM1169478","GSM1169482")

dds<-DESeqDataSetFromMatrix(countData=cts, colData=coldata_th1, design=~celltype)
dds$celltype<-relevel(dds$celltype, ref="CD4+")
dds<-DESeq(dds)
resultsNames(dds)

res<-results(dds, name="celltype_Th1_vs_CD4.")
res<-lfcShrink(dds, coef="celltype_Th1_vs_CD4.", type="apeglm")

plotMA(res, ylim=c(-2,2))
plotDispEsts(dds)
EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange',y='pvalue',title='CD4+ vs Th1',pCutoff=1e-5,pointSize=3.0,labSize=6.0)


res<-na.omit(res)
res['pvalue']=res['pvalue'][,1]<0.05
res=res[res$pvalue==TRUE,]
res['padj']=res['padj'][,1]<0.05
res=res[res$padj==TRUE,]
rownames(res)
row=rownames(res)

write.table(rownames(res), file="sign_gene_th1vscd4")
save(row,file="significant_gene_th1_vs_cd4.rda")
