# Author: Christian Brion - 2022 - project IPF_DE
#
#-load and process expression data to match between gene and microRNA
#-spearman correlation analysis DE microRNA vs DE genes
#-identify 10 microRNA with high impact
#-correlation analysis microRNA high-impact vs nonDE genes: OFF target effect of micro RNA
#-use microRNA targets from mirDIP to identify ON/OFF target from interaction database
#

#initiation: install, load packages, define function, set work directory
library(corrplot)


colSdv<-function(df) {
  unlist(lapply(X = 1:ncol(df), function(x) {sd(df[,x])}))
}

rowMedian<-function(df) {
  unlist(lapply(X = 1:nrow(df), function(x) {median(df[x,])}))
}

rowSdv<-function(df) {
  unlist(lapply(X = 1:nrow(df), function(x) {sd(df[x,])}))
}

rowMin<-function(df) {
  unlist(lapply(X = 1:nrow(df), function(x) {min(df[x,])}))
}

colMin<-function(df) {
  unlist(lapply(X = 1:ncol(df), function(x) {min(df[,x])}))
}

rowMax<-function(df) {
  unlist(lapply(X = 1:nrow(df), function(x) {max(df[x,])}))
}

colMax<-function(df) {
  unlist(lapply(X = 1:ncol(df), function(x) {max(df[,x])}))
}

workdir <- "C:/Users/chris/Dropbox/Github/IPF_DE" #set work directory
setwd(workdir)

#===================
#load metadata

geneTestvalue <- read.table(file = "geneTestvalue.txt",header = T,quote = "\"",sep = "\t",stringsAsFactors = T)
sampleMD2 <- read.table(file = "sampleMD2.txt",header = T,quote = "\"",sep = "\t",stringsAsFactors = T)
load("geneexp_processeddata.RData")
sum(rownames(geneTestvalue)==genepMDunique$ID)

sample_miMD2 <- read.table(file = "sample_miMD2.txt",header = T,quote = "\"",sep = "\t",stringsAsFactors = T)
load("microRNAexp_processeddata.RData")

matching <- read.table(file = "matching.txt",header = T,quote = "\"",sep = "\t",stringsAsFactors = T)
geneclust <- read.table(file = "geneclust.txt",header = T,quote = "\"",sep = "\t",stringsAsFactors = T)
miRNAclus <- read.table(file = "miRNAclust.txt",header = T,quote = "\"",sep = "\t",stringsAsFactors = T)

#change to gene names
geneclust2<-merge(geneclust,genepMDunique[,colnames(genepMDunique)%in%c("ID","geneName")],by.x=1,by.y=1,all.x=T,all.y=F,sort = F)
microRNApMD_human$geneName<-gsub("_st","",microRNApMD_human$ID)
#microRNApMD_human$geneName<-gsub("-star","",microRNApMD_human$geneName)

miRNAclus2<-merge(miRNAclus,microRNApMD_human[,colnames(microRNApMD_human)%in%c("ID","geneName")],by.x=1,by.y=1,all.x=T,all.y=F,sort = F)
miRNAclus2$clust<-miRNAclus2$clust+100
summary(as.factor(geneclust2$clust))

#save cluster data in csv for easy copy paste on web base tool
allList<-rbind(geneclust2,miRNAclus2)
write.table(allList,file = "allList_clust.csv",sep=",",quote = F,row.names = F,col.names = T)

#===================
#process gene list and sample

#only DE genes and only matching "IPF/UIP" samples
data.mat.only<-data.mat[rownames(data.mat) %in% geneclust2$probID, colnames(data.mat) %in% sampleMD2$Sample_title[sampleMD2$Sample_title %in% matching$Sample_title & sampleMD2$final.diagnosis == "IPF/UIP"]]
matchgeneNAME<-merge(data.mat.only[,1],genepMDunique[,colnames(genepMDunique) %in% c("ID", "geneName")],by.x=0,by.y=1,all.x=T,all.y=F,sort=F)
sum(matchgeneNAME$Row.names==rownames(data.mat.only))
rownames(data.mat.only)<-matchgeneNAME$geneName
data.mat.miRNA.only<-data.mat.miRNA[rownames(data.mat.miRNA) %in% miRNAclus2$probID, colnames(data.mat.miRNA) %in% sample_miMD2$Sample_title[sample_miMD2$Sample_title %in% matching$Sample_title & sample_miMD2$final.diagnosis == "IPF/UIP"]]
rownames(data.mat.miRNA.only)<-gsub("_st","",rownames(data.mat.miRNA.only))
#rownames(data.mat.miRNA.only)<-gsub("-star","",rownames(data.mat.miRNA.only))

sum(colnames(data.mat.only)==colnames(data.mat.miRNA.only)) #control that the sample (columns) match perfectly between the two dataset

#===================
#correlation analysis DE microRNA vs DE genes (spearman)

M <- cor(t(data.mat.only),t(data.mat.miRNA.only),method = "spearman")

miRNAoi<-t(data.mat.miRNA.only)[,1]
x=1
spearmanpv<-function(x) { #function for calculating spearman pvalue (for lapply). need to find a better tool instead of doing it manually
  geneoi<-t(data.mat.only)[,x]
  test<-cor.test(geneoi,miRNAoi,method = "spearman")
  res<-test$p.value
  res
}

#calculating spearman pvalue manually (long process)
date()
pvalM<-M
for(i in 1:ncol(t(data.mat.miRNA.only))) {
  miRNAoi<-t(data.mat.miRNA.only)[,i]
  vect<-sapply(1:ncol(t(data.mat.only)),FUN = spearmanpv,simplify = T)
  pvalM[,i]<-vect
}
date() #3 minutes to run

#Benjamini-Hochberg test correction (FDR))
pvalMvect<-c(pvalM)
adjpvalM<-pvalM
adjpvalMvect<-p.adjust(pvalMvect,method = "BH")
sum(adjpvalMvect<0.01)
adjpvalM[]<-matrix(adjpvalMvect,nrow = nrow(adjpvalM),ncol = ncol(adjpvalM))
plot(pvalM[,1],adjpvalM[,1]) #good

sum(colMin(adjpvalM)<0.01)
sum(rowMin(adjpvalM)<0.01)

#compare number of positive/negative correlation for each microRNA
pos<-M>0
sign<-adjpvalM<0.01 
colmicroRNA<-rep("deepskyblue",ncol(M))
colmicroRNA[colnames(M) %in% miRNAclus2$geneName[miRNAclus2$clust==102]] <- "orange"
svg(filename = "ONtarget.svg",width = 6,height = 5)
plot(colSums(sign & pos),colSums(sign & !pos),pch=16,col=colmicroRNA,xlab = "# of ON target positive correlation",ylab = "# of ON target negative correlation")
abline(lm(colSums(sign & !pos)~colSums(sign & pos)))
abline(a=150,b=-1,col="red",lty=2)
miRNA_highimpact<-colnames(M)[colSums(sign)>150] #get the 10 DE miRNA with the most correlation
text(colSums(sign & pos)[colnames(M) %in% miRNA_highimpact],colSums(sign & !pos)[colnames(M) %in% miRNA_highimpact],
     labels = miRNA_highimpact,pos = 2,col=rgb(0,0,0,0.4))
dev.off()

M[rownames(M)=="MUC5B",names(colSums(sign)) %in% miRNA_highimpact]
M[rownames(M)=="CCDC19",names(colSums(sign)) %in% miRNA_highimpact]

#exemple of correlation
bestGene<-"CCDC19"
bestmiRNA<-"hsa-miR-34c-5p"
svg(filename = "correlationEG.svg",width = 6,height = 5)
plot(data.mat.miRNA.only[rownames(data.mat.miRNA.only)==bestmiRNA,],data.mat.only[rownames(data.mat.only)==bestGene,],
     xlab=bestmiRNA,ylab=bestGene)
dev.off()

#save 10 DE miRNA with the most correlation in csv for easy copy paste on web base tool
write.table(data.frame(miRNA_highimpact),"microRNA_highimpact.csv",sep = ',',row.names = F,col.names = F)

#select randomly 5 genes from each cluster to provide the corrplot
clustplot<-2
geneclust2<-geneclust2[order(geneclust2$clust),]
rand32<-c(sample((1:nrow(geneclust2))[geneclust2$clust==2],5,replace = F),
          sample((1:nrow(geneclust2))[geneclust2$clust%in%c(5,6,8)],5,replace = F),
          sample((1:nrow(geneclust2))[geneclust2$clust==4],5,replace = F),
          sample((1:nrow(geneclust2))[geneclust2$clust%in%c(1,7,3)],5,replace = F))
corrplot(M[as.character(geneclust2$geneName[rand32]),colnames(M) %in% miRNA_highimpact],
         p.mat = adjpvalM[as.character(geneclust2$geneName[rand32]),colnames(M) %in% miRNA_highimpact],
         sig.level = 0.10, addrect = 2)
# #to run only once as each run will provide a different plot
# svg(filename = "corrplot_upreg.svg",width = 6,height = 7)
# corrplot(M[as.character(geneclust2$geneName[rand32]),colnames(M) %in% miRNA_highimpact],
#          p.mat = adjpvalM[as.character(geneclust2$geneName[rand32]),colnames(M) %in% miRNA_highimpact],
#          sig.level = 0.10, addrect = 2)
# dev.off()


#===================
#correlation analysis microRNA high-impact vs nonDE genes: OFF target effect of micro RNA

#gene expression data only highly expressed gene unlikely to be linked to IPF (mean>5 and adjpv_IPF.UIP > 0.1)
geneOFF<-rownames(geneTestvalue)[geneTestvalue$mean>5 & geneTestvalue$adjpv_IPF.UIP > 0.1]
data.mat.other<-data.mat[rownames(data.mat) %in% geneOFF, colnames(data.mat) %in% sampleMD2$Sample_title[sampleMD2$Sample_title %in% matching$Sample_title & sampleMD2$final.diagnosis == "IPF/UIP"]]
matchgeneNAME<-merge(data.mat.other[,1],genepMDunique[,colnames(genepMDunique) %in% c("ID", "geneName")],by.x=0,by.y=1,all.x=T,all.y=F,sort=F)
sum(matchgeneNAME$Row.names==colnames(data.mat.other))
rownames(data.mat.other) <- matchgeneNAME$geneName

#miRNA expression data only for the 10 high impact miRNA
data.mat.miRNA.HI<-data.mat.miRNA.only[rownames(data.mat.miRNA.only) %in% miRNA_highimpact,]

#correlation analysis with pvalue and adj pvalue (same as above)
M_OFF <- cor(t(data.mat.other),t(data.mat.miRNA.HI),method = "spearman")
date()
pvalM_OFF<-M_OFF
spearmanpv_OFF<-function(x) {
  geneoi<-t(data.mat.other)[,x]
  test<-cor.test(geneoi,miRNAoi,method = "spearman")
  res<-test$p.value
  res
}
for(i in 1:ncol(t(data.mat.miRNA.HI))) {
  miRNAoi<-t(data.mat.miRNA.HI)[,i]
  vect<-sapply(1:ncol(t(data.mat.other)),FUN = spearmanpv_OFF,simplify = T)
  pvalM_OFF[,i]<-vect
  print(paste(i, date()))
}
date() #4 minutes 

pvalMvect_OFF<-c(pvalM_OFF)
adjpvalM_OFF<-pvalM_OFF
adjpvalMvect_OFF<-p.adjust(pvalMvect_OFF,method = "BH")
adjpvalM_OFF[]<-matrix(adjpvalMvect_OFF,nrow = nrow(adjpvalM_OFF),ncol = ncol(adjpvalM_OFF))
sum(adjpvalMvect_OFF<0.01)

sign_OFF<-adjpvalM_OFF<0.01
colSums(sign_OFF)
colSums(sign)[names(colSums(sign)) %in% miRNA_highimpact]

#compare number of ON/OFF target correlation for each microRNA
svg(filename = "ON-OFFtarget.svg",width = 6,height = 5)
plot(colSums(sign)[names(colSums(sign)) %in% miRNA_highimpact],colSums(sign_OFF),pch=16,col=colmicroRNA[names(colSums(sign)) %in% miRNA_highimpact],xlim = c(0,1000),ylim = c(0,1000),xlab = "# of ON target correlation",ylab = "# of OFF target correlation")
text(colSums(sign)[names(colSums(sign)) %in% miRNA_highimpact],colSums(sign_OFF),
     labels = names(colSums(sign_OFF)),pos = 2,col=rgb(0,0,0,0.4))
dev.off()

#===================
#process microRNA targets from mirDIP (http://ophid.utoronto.ca/mirDIP/index.jsp#r)

#import data
mirDIP_E_DOWN <- read.table(file = "mirDIP/mirDIP_E_DOWN.txt",header = T,quote = "\"",sep = "\t",stringsAsFactors = T,skip = 802)
mirDIP_E_UP <- read.table(file = "mirDIP/mirDIP_E_UP.txt",header = T,quote = "\"",sep = "\t",stringsAsFactors = T,skip = 1152)
mirDIP_microRNAhighimpact <- read.table(file = "mirDIP/mirDIP_microRNAhighimpact.txt",header = T,quote = "\"",sep = "\t",stringsAsFactors = T,skip = 20)

levels(mirDIP_microRNAhighimpact$MicroRNA)
miRNA_highimpact #microRNA names have been converted in mirDIP

#correct the microRNA names to unprocessed to match prob
mirDIP_E_DOWN$MicroRNA<-as.character(mirDIP_E_DOWN$MicroRNA)
mirDIP_E_DOWN$MicroRNA[mirDIP_E_DOWN$MicroRNA %in% c("hsa-miR-126-3p","hsa-miR-126-5p")] <- "hsa-miR-126"
mirDIP_E_DOWN$MicroRNA[mirDIP_E_DOWN$MicroRNA %in% c("hsa-miR-205-3p","hsa-miR-205-5p")] <- "hsa-miR-205"
mirDIP_E_DOWN$MicroRNA[mirDIP_E_DOWN$MicroRNA %in% c("hsa-miR-30a-3p","hsa-miR-30a-5p")] <- "hsa-miR-30a"
mirDIP_E_DOWN$MicroRNA[mirDIP_E_DOWN$MicroRNA %in% c("hsa-miR-30d-3p","hsa-miR-30d-5p")] <- "hsa-miR-30d"
mirDIP_E_DOWN$MicroRNA[mirDIP_E_DOWN$MicroRNA %in% c("hsa-miR-31-3p","hsa-miR-31-5p")] <- "hsa-miR-31"
mirDIP_E_DOWN$MicroRNA[mirDIP_E_DOWN$MicroRNA %in% c("hsa-miR-652-3p","hsa-miR-652-5p")] <- "hsa-miR-652"
mirDIP_E_DOWN$MicroRNA<-as.factor(mirDIP_E_DOWN$MicroRNA)
summary(mirDIP_E_DOWN$MicroRNA)

mirDIP_E_UP$MicroRNA<-as.character(mirDIP_E_UP$MicroRNA)
mirDIP_E_UP$MicroRNA[mirDIP_E_UP$MicroRNA %in% c("hsa-miR-126-3p","hsa-miR-126-5p")] <- "hsa-miR-126"
mirDIP_E_UP$MicroRNA[mirDIP_E_UP$MicroRNA %in% c("hsa-miR-205-3p","hsa-miR-205-5p")] <- "hsa-miR-205"
mirDIP_E_UP$MicroRNA[mirDIP_E_UP$MicroRNA %in% c("hsa-miR-30a-3p","hsa-miR-30a-5p")] <- "hsa-miR-30a"
mirDIP_E_UP$MicroRNA[mirDIP_E_UP$MicroRNA %in% c("hsa-miR-30d-3p","hsa-miR-30d-5p")] <- "hsa-miR-30d"
mirDIP_E_UP$MicroRNA[mirDIP_E_UP$MicroRNA %in% c("hsa-miR-31-3p","hsa-miR-31-5p")] <- "hsa-miR-31"
mirDIP_E_UP$MicroRNA[mirDIP_E_UP$MicroRNA %in% c("hsa-miR-652-3p","hsa-miR-652-5p")] <- "hsa-miR-652"
mirDIP_E_UP$MicroRNA<-as.factor(mirDIP_E_UP$MicroRNA)
summary(mirDIP_E_UP$MicroRNA)

mirDIP_microRNAhighimpact$MicroRNA<-as.character(mirDIP_microRNAhighimpact$MicroRNA)
mirDIP_microRNAhighimpact$MicroRNA[mirDIP_microRNAhighimpact$MicroRNA %in% c("hsa-miR-126-3p","hsa-miR-126-5p")] <- "hsa-miR-126"
mirDIP_microRNAhighimpact$MicroRNA[mirDIP_microRNAhighimpact$MicroRNA %in% c("hsa-miR-205-3p","hsa-miR-205-5p")] <- "hsa-miR-205"
mirDIP_microRNAhighimpact$MicroRNA[mirDIP_microRNAhighimpact$MicroRNA %in% c("hsa-miR-30a-3p","hsa-miR-30a-5p")] <- "hsa-miR-30a"
mirDIP_microRNAhighimpact$MicroRNA[mirDIP_microRNAhighimpact$MicroRNA %in% c("hsa-miR-30d-3p","hsa-miR-30d-5p")] <- "hsa-miR-30d"
mirDIP_microRNAhighimpact$MicroRNA[mirDIP_microRNAhighimpact$MicroRNA %in% c("hsa-miR-31-3p","hsa-miR-31-5p")] <- "hsa-miR-31"
mirDIP_microRNAhighimpact$MicroRNA[mirDIP_microRNAhighimpact$MicroRNA %in% c("hsa-miR-652-3p","hsa-miR-652-5p")] <- "hsa-miR-652"
mirDIP_microRNAhighimpact$MicroRNA<-as.factor(mirDIP_microRNAhighimpact$MicroRNA)
summary(mirDIP_microRNAhighimpact$MicroRNA)

#eliminate duplicates
mirDIP_microRNAhighimpact <- unique(mirDIP_microRNAhighimpact,by=c(1,4))

#filter for only miRNA high impact and eliminate duplicates
mirDIP_E_DOWN_highImpact <- mirDIP_E_DOWN[mirDIP_E_DOWN$MicroRNA%in%miRNA_highimpact,]
mirDIP_E_DOWN_highImpact <- unique(mirDIP_E_DOWN_highImpact,by=c(1,4))
mirDIP_E_UP_highImpact <- mirDIP_E_UP[mirDIP_E_UP$MicroRNA%in%miRNA_highimpact,]
mirDIP_E_UP_highImpact <- unique(mirDIP_E_UP_highImpact,by=c(1,4))
mirDIP_E_DOWN_highImpact$MicroRNA <- as.factor(as.character(mirDIP_E_DOWN_highImpact$MicroRNA))

mirDIP_E_ALL_highImpact <- rbind(mirDIP_E_UP_highImpact,mirDIP_E_DOWN_highImpact)
mirDIP_E_ALL_highImpact <- unique(mirDIP_E_ALL_highImpact,by=c(1,4))
mirDIP_E_ALL_highImpact$MicroRNA <- as.factor(as.character(mirDIP_E_ALL_highImpact$MicroRNA))

#percentage of ON target interaction compare to all
svg(filename = "mirDIP-ON-OFFtarget.svg",width = 11,height = 8)
b<-barplot(summary(mirDIP_E_ALL_highImpact$MicroRNA)/summary(mirDIP_microRNAhighimpact$MicroRNA)*100,ylim = c(0,15),col = "orange",ylab="mirDIP: percentage on target",names.arg = gsub("hsa-","",names(summary(mirDIP_E_ALL_highImpact$MicroRNA))))
text(b,summary(mirDIP_E_ALL_highImpact$MicroRNA)/summary(mirDIP_microRNAhighimpact$MicroRNA)*100,round(summary(mirDIP_E_ALL_highImpact$MicroRNA)/summary(mirDIP_microRNAhighimpact$MicroRNA)*100,2),pos = 3)
b<-barplot(summary(mirDIP_E_DOWN_highImpact$MicroRNA)/summary(mirDIP_microRNAhighimpact$MicroRNA)*100,ylim = c(0,15),col = "deepskyblue",add=T,xaxt='n',yaxt='n')
dev.off()
