# Author: Christian Brion - 2022 - project IPF_DE
#
#-process ncRNA expression data: select only microRNA
#-perform differential gene expression analysis using limma (DEseq2 does not work on micro array data)
#-data visualization: volcano plot, MA plot, boxplot
#-second differential gene expression analysis using limma to check the effect of smoking (no effect)
#-clustering analysis on centered reduced DE genes using pheatmap
#-PCA analysis on all gene using prcomp
#

#initiation: install, load packages, define function, set work directory
library(limma)
library("pheatmap")

library(tidyverse)
library(data.table)
library(Matrix)
library(lattice)
library(fdrtool)
library(rpart)


colSdv<-function(df) {
  unlist(lapply(X = 1:ncol(df), function(x) {sd(df[,x])}))
}

rowMedian<-function(df) {
  unlist(lapply(X = 1:nrow(df), function(x) {median(df[x,])}))
}

rowSdv<-function(df) {
  unlist(lapply(X = 1:nrow(df), function(x) {sd(df[x,])}))
}

workdir <- "C:/Users/chris/Dropbox/Github/IPF_DE" #set work directory
setwd(workdir)

#====================================
# load datas
#====================================

sampleMD2 <- read.table(file = "sample_miMD2.txt",header = T,quote = "\"",sep = "\t",stringsAsFactors = T)
load("miRNA_data.RData")
colnames(data_mirna.matrix)<-substr(colnames(data_mirna.matrix),11,18)
data.matrix<-data_mirna.matrix
sum(colnames(data.matrix) %in% sampleMD2$Sample_title) #good
sum(colnames(data.matrix) == sampleMD2$Sample_title) #good
genepMD <- read.table(file = "miRNApMD.txt",header = T,quote = "\"",sep = "\t",stringsAsFactors = T,comment.char = "")

sum(rownames(data.matrix) %in% genepMD$ID)
genepMD <- genepMD[order(genepMD$ID),]
sum(rownames(data.matrix) == genepMD$ID)

#look at non-human ncRNA (before taking them out)
summary(genepMD$Species.Scientific.Name)
boxplot(rowMeans(data.matrix) ~ genepMD$Species.Scientific.Name)
boxplot(rowMeans(data.matrix) ~ genepMD$Species.Scientific.Name=="Homo sapiens")

#take out the non-human probs and non-microRNA probes
data.mat <- data.matrix[genepMD$Species.Scientific.Name=="Homo sapiens" &
                        !(genepMD$Sequence.Type %in% c("Oligonucleotide spike-in controls","5.8s rRNA","CDBox","HAcaBox","snoRNA","scaRna")),]
genepMD <- genepMD[genepMD$Species.Scientific.Name=="Homo sapiens" &
                     !(genepMD$Sequence.Type %in% c("Oligonucleotide spike-in controls","5.8s rRNA","CDBox","HAcaBox","snoRNA","scaRna")),]
summary(genepMD$Sequence.Type)

data.mat.miRNA<-data.mat
microRNApMD_human<-genepMD
save(data.mat.miRNA,microRNApMD_human,file = "microRNAexp_processeddata.RData")


#====================================
# limma DE
#====================================

#differential expression analysis
design_DF = model.matrix(~ sampleMD2$age + sampleMD2$gender + sampleMD2$final.diagnosis)
data.fit = lmFit(data.mat,design=design_DF)
data.fit.eb = eBayes(data.fit) #perform the multivariate differential expression analysis
data.fit.eb$coefficients[1:2,]

#Benjamini-Hochberg test correction (FDR)
FDR<-data.fit.eb$p.value
for (i in 2:ncol(FDR)) {
  FDR[,i]<-p.adjust(FDR[,i],method = "BH")
}
FDR<-FDR[,-1]

#prepare output test table
testoutput<-c("mean","age","Male","COP","DIP","IPF/UIP","NSIP","RB-ILD","UF")
testadjpv<-paste0("adjpv_",testoutput[-1])
testDF<-cbind(data.fit.eb$coefficients,FDR)
colnames(testDF)<-c(testoutput,testadjpv)
testDF[1:6,]

#all volcano plot
for (i in 1:8) {
  plot(testDF[,i+1],-log10(testDF[,i+9]),pch=16,cex=0.5,col=rgb(0,0,0,0.1),main=colnames(testDF)[i+1],xlab="log2FC",ylab="-log10adjpv")
  abline(h=-log10(0.01))
  legend("topleft",legend = sum(testDF[,i+1]<0 & testDF[,i+9] < 0.01),pch = NA)
  legend("topright",legend = sum(testDF[,i+1]>0 & testDF[,i+9] < 0.01),pch = NA)
}

svg(filename = "DEmicroRNA_IPF.svg",width = 6,height = 5)
i=5
plotNS<-testDF[abs(testDF[,i+1])<0.5 | testDF[,i+9] > 0.01,]
plotUP<-testDF[testDF[,i+1]>0.5 & testDF[,i+9] < 0.01,]
plotDOWN<-testDF[testDF[,i+1]< -0.5 & testDF[,i+9] < 0.01,]
plot(plotNS[,i+1],-log10(plotNS[,i+9]),pch=16,cex=0.7,col=rgb(0,0,0,0.2),main=colnames(testDF)[i+1],xlab="log2FC",ylab="-log10adjpv",xlim=c(-3,3),ylim=c(0,14))
points(plotUP[,i+1],-log10(plotUP[,i+9]),pch=16,cex=1,col="#ffa50099")
points(plotDOWN[,i+1],-log10(plotDOWN[,i+9]),pch=16,cex=1,col="#00bfff99")
abline(h=-log10(0.01),col="red")
abline(v=c(-0.5,0.5),lty=2)
legend("topleft",legend = nrow(plotDOWN),pch = NA)
legend("topright",legend = nrow(plotUP),pch = NA)
dev.off()

#MA plot
plot(testDF[,1],testDF[,6],pch=16,cex=0.5,col=rgb(0,0,0,0.1),main=colnames(testDF)[6],ylab="log2FC",xlab="average (log2)")
sum(testDF[,6]< -0.5 & testDF[,14] < 0.01)
sum(testDF[,6]> 0.5 & testDF[,14] < 0.01)

i=5
#data.mat[testDF[,i+1]>0 & testDF[,i+9] < 0.00001,][10,]
testDF[testDF[,i+1]< -0.5 & testDF[,i+9] < 0.01,][10,]

boxplot(data.mat[testDF[,i+1]< -0.5 & testDF[,i+9] < 0.01,][10,]~sampleMD2$final.diagnosis) #IPF effect for one miRNA
plot(data.mat[testDF[,i+1]< -0.5 & testDF[,i+9] < 0.01,][10,]~sampleMD2$smoke_impact_score) #smoking effect for the same miRNA

#===========================
#same analysis but subset with smoking score
# design_DF = model.matrix(~ sampleMD2$age + sampleMD2$gender + sampleMD2$final.diagnosis + sampleMD2$smoke_impact_score)
# data.fit = lmFit(data.mat[,!is.na(sampleMD2$smoke_impact_score)],design=design_DF)
# data.fit.eb = eBayes(data.fit)
# 
# data.fit.eb$coefficients[1:2,]
# 
# FDR<-data.fit.eb$p.value
# 
# for (i in 2:ncol(FDR)) {
#   FDR[,i]<-p.adjust(FDR[,i],method = "BH")
# }
# FDR<-FDR[,-1]
# 
# testoutput<-c("mean","age","Male","COP","DIP","IPF/UIP","NSIP","RB-ILD","UF","smoke_impact_score")
# testadjpv<-paste0("adjpv_",testoutput[-1])
# 
# testDF<-cbind(data.fit.eb$coefficients,FDR)
# colnames(testDF)<-c(testoutput,testadjpv)
# testDF[1:6,]
# 
# for (i in 1:9) {
#   plot(testDF[,i+1],-log10(testDF[,i+10]),pch=16,cex=0.5,col=rgb(0,0,0,0.1),main=colnames(testDF)[i+1],xlab="log2FC",ylab="-log10adjpv")
#   abline(h=-log10(0.01))
#   legend("topleft",legend = sum(testDF[,i+1]<0 & testDF[,i+10] < 0.01),pch = NA)
#   legend("topright",legend = sum(testDF[,i+1]>0 & testDF[,i+10] < 0.01),pch = NA)
# }
# 
# plot(testDF[,1],testDF[,6],pch=16,cex=0.5,col=rgb(0,0,0,0.1),main=colnames(testDF)[6],ylab="log2FC",xlab="average (log2)")

###not many genes significantly associated to smoking score


#====================================
# pheatmap heatmap
#====================================

#centered reduced data
data.mat.centered <- data.mat-rowMeans(data.mat)
data.mat.cr <- data.mat.centered/rowSdv(data.mat.centered)
data.mat.cr[1:10,1:10]

#QC mean vs sdv plot
plot(rowMeans(data.mat),rowSdv(data.mat.centered),pch=16,cex=0.5,col=rgb(0,0,0,0.2))

i=5
colnames(testDF)[i+1]
#only DE genes for IPF and only IPF & control samples
HMdata <- data.mat.cr[testDF[,i+9] < 0.01  & abs(testDF[,i+1]) > 0.5,sampleMD2$final.diagnosis %in% c("IPF/UIP","control")] # & abs(testDF[,i+1]) > 1
nrow(HMdata)
#clustering analysis
colset<-c(colorRampPalette(c("deepskyblue","deepskyblue"))(265),
          colorRampPalette(c("deepskyblue", "black", "yellow"))(300),
          colorRampPalette(c("yellow","yellow"))(150))
clustDATA<-pheatmap(HMdata, cluster_rows=T, show_rownames=FALSE, cluster_cols=T, show_colnames = FALSE, color = colset)
clustvis<-sampleMD2[sampleMD2$final.diagnosis %in% c("IPF/UIP","control"),][clustDATA$tree_col$order,]

#provide miRNA and samples with cluster ID
samplecut<-data.frame(clustDATA$tree_col$labels,cutree(clustDATA$tree_col,  k = 4))
colnames(samplecut)<-c("sampleID","clust")
summary(as.factor(samplecut$clust))
samplecut<-samplecut[clustDATA$tree_col$order,]
genecut<-data.frame(clustDATA$tree_row$labels,cutree(clustDATA$tree_row,  k = 2))
colnames(genecut)<-c("probID","clust")
summary(as.factor(genecut$clust))
genecut<-genecut[clustDATA$tree_row$order,]

# #save figure
# png("heatmap_miRNA.png",width = 800,height = 1000)
# pheatmap(HMdata, cluster_rows=T, show_rownames=FALSE, cluster_cols=T, show_colnames = FALSE, color = colset)
# dev.off()

# #clustering with the 20% most variable genes: not informative
# HMdata <- data.mat.cr[rowSdv(data.mat.centered)>quantile(rowSdv(data.mat.centered),probs = 0.80),]
# nrow(HMdata)
# colset<-c(colorRampPalette(c("deepskyblue","deepskyblue"))(200),
#           colorRampPalette(c("deepskyblue", "black", "yellow"))(300),
#           colorRampPalette(c("yellow","yellow"))(300))
# clustDATA<-pheatmap(HMdata, cluster_rows=T, show_rownames=FALSE, cluster_cols=T, show_colnames = FALSE, color = colset)
# clustvis<-sampleMD2[clustDATA$tree_col$order,]

#save data
write.table(testDF,"miRNATestvalue.txt",quote = F,sep = "\t",col.names = T,row.names = T)
write.table(genecut,"miRNAclust.txt",quote = F,sep = "\t",col.names = T,row.names = F)
write.table(samplecut,"samplemiRNAclust.txt",quote = F,sep = "\t",col.names = T,row.names = F)

#how the differential expressed miRNA correlate with age, smoke, clinical score
data.mat.onlyIPF<-data.mat[rownames(data.mat) %in% genecut$probID[genecut$clust %in% c(1)],sampleMD2$final.diagnosis=="IPF/UIP"]
see<-cor(t(data.mat.onlyIPF),sampleMD2[sampleMD2$final.diagnosis=="IPF/UIP",colnames(sampleMD2)%in%c("age","rin","st..george.s.total.score","fvc.pre.bronchodilator...predicted","dlco...predicted","smoke_impact_score")]
         , use = "complete.obs", method = "spearman")
boxplot(see,ylim=c(-1,1))
for (i in 1:6) {
  hist(see[,i],breaks=100,xlim=c(-1,1),main = colnames(see)[i])
}

#PCA visualization
data_PC = prcomp(t(data.mat),scale.=TRUE)
data_PC$percent <- data_PC$sdev/sum(data_PC$sdev)

colorPC<-rep("black",nrow(sampleMD2))
colorPC[sampleMD2$Sample_title %in% samplecut$sampleID[samplecut$clust%in%c(2,3)]]<-"orangered"
colorPC[sampleMD2$Sample_title %in% samplecut$sampleID[samplecut$clust==1]]<-"orange"
colorPC[sampleMD2$Sample_title %in% samplecut$sampleID[samplecut$clust==4]]<-"deepskyblue"

svg(filename = "microRNAEXP_PCA.svg",width = 6,height = 5)
plot(data_PC$x[,1:2],pch=16,col=colorPC,xlab=paste0("PC1 (",round(data_PC$percent[1]*100,1),"%)"),
     ylab=paste0("PC2(",round(data_PC$percent[2]*100,1),"%)"))
legend("topright",legend = c("cluster2/3: eavy impact","cluster1: ligh impact","cluster4: control "),pch=15,col=c("orangered","orange","deepskyblue"),bg = NA)
dev.off()