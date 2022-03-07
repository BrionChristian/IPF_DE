# Author: Christian Brion - 2022 - project IPF_DE
#
#-import sample meta data from series matrix
#-parse sample information
#-compute smoking index (tot # of pack/#year w/o smoking)
#-import mRNA array prob meta data from GPL files using the txt file
#-parse prob meta data
#-import ncRNA prob meta data from GPL files using GEOquery
#-compare mRNA and ncRNA samples meta data
#-export datas
#

#initiation: install, load packages, define function, set work directory
#BiocManager::install("GEOquery",force = T)
#BiocManager::install("GEOmetadb",force = T)
library(GEOquery)
library(GEOmetadb)
library(Biobase)

library(tidyverse)
library(data.table)
library(Matrix)
library(lattice)
library(fdrtool)
library(rpart)

colSdv<-function(df) {
  unlist(lapply(X = 1:ncol(df), function(x) {sd(df[x,])}))
}

workdir <- "C:/Users/chris/Dropbox/Github/IPF_DE" #set work directory
setwd(workdir)

#====================================
# load and parse and explore meta datas
#====================================

#======================
#mRNA sample meta data

#import data
sampleMD <- as.data.frame(t(read.table(file = "GSE32537_series_matrix.txt",header = F,skip = 44,nrows = 87-44,quote = "\"",sep = "\t")),stringsAsFactors = F)
colnames(sampleMD) <- sampleMD[1,]
sampleMD <- sampleMD[-1,]
colnames(sampleMD) <- sub("!","",colnames(sampleMD))

#parse Sample_char columns
col_info <- grep("Sample_char",colnames(sampleMD))
sample_info_type <- unlist(strsplit(t(sampleMD[1,col_info]),": "))[(1:length(col_info))*2-1]
temp <- data.frame(lapply(col_info,function(x) {unlist(strsplit(t(sampleMD[,x]),": "))[(1:nrow(sampleMD))*2]}))
sampleMD[,col_info] <- temp
colnames(sampleMD)[col_info] <- sample_info_type
sampleMD[sampleMD=="--"] <- NA
sampleMD$`quit how many years ago`[sampleMD$`quit how many years ago`=="nonsmoker"] <- NA
sampleMD2 <- sampleMD[,c(1,2,col_info)]
for (i in 3:ncol(sampleMD2)) {
  if (sum(!is.na(as.numeric(sampleMD2[,i]))) > 0) {
    sampleMD2[,i] <- as.numeric(sampleMD2[,i])
  } else {
    sampleMD2[,i] <- as.factor(sampleMD2[,i])
  }
}

#calculate smoking score
sampleMD2$`smoking status`<- as.character(sampleMD2$`smoking status`)
sampleMD2$`smoking status`[is.na(sampleMD2$`smoking status`)]<-"NoData"
smoke_impact_score<-rep(0,nrow(sampleMD2))
smoke_impact_score[sampleMD2$`smoking status`=="NoData"]<-NA
smoke_impact_score[sampleMD2$`smoking status`=="current"]<-(sampleMD2$age[sampleMD2$`smoking status`=="current"]-18)*sampleMD2$`pack years`[sampleMD2$`smoking status`=="current"]
smoke_impact_score[sampleMD2$`smoking status`=="former"]<-(sampleMD2$age[sampleMD2$`smoking status`=="former"]-18-sampleMD2$`quit how many years ago`[sampleMD2$`smoking status`=="former"])*sampleMD2$`pack years`[sampleMD2$`smoking status`=="former"]/(sampleMD2$`quit how many years ago`[sampleMD2$`smoking status`=="former"]/1)
smoke_impact_score[smoke_impact_score<0] <- 0
smoke_impact_score<-log2(smoke_impact_score+1)
sampleMD2$smoke_impact_score <- smoke_impact_score

#quick visualization of the data
summary(sampleMD2$gender)
summary(sampleMD2$`final diagnosis`)
summary(sampleMD2$`tissue source`)
summary(sampleMD2$preservative)
summary(as.factor(paste(sampleMD2$gender,sampleMD2$`final diagnosis`)))
boxplot(sampleMD2$age ~ sampleMD2$gender)
boxplot(sampleMD2$age ~ sampleMD2$`final diagnosis`) #age and diagnosis are confounding
boxplot(sampleMD2$smoke_impact_score ~ sampleMD2$gender)
boxplot(sampleMD2$smoke_impact_score ~ sampleMD2$`final diagnosis`)
plot(sampleMD2$age,sampleMD2$smoke_impact_score)

svg(filename = "age_smoking_confounders.svg",width = 6,height = 8)
par(mfrow=c(2,1))
color<-c("deepskyblue",rep("orange",6))
boxplot(sampleMD2$age ~ sampleMD2$`final diagnosis`,xlab="final diagnosis",ylab="age",col=color)
boxplot(sampleMD2$smoke_impact_score ~ sampleMD2$`final diagnosis`,xlab="final diagnosis",ylab="smoking index",col=color)
dev.off()
par(mfrow=c(1,1))


#======================
#microRNA sample meta data

#import data
sample_miMD <- as.data.frame(t(read.table(file = "GSE32538_series_matrix.txt",header = F,skip = 44,nrows = 87-44,quote = "\"",sep = "\t")),stringsAsFactors = F)
colnames(sample_miMD) <- sample_miMD[1,]
sample_miMD <- sample_miMD[-1,]
colnames(sample_miMD) <- sub("!","",colnames(sample_miMD))

#parse Sample_char columns
col_info <- grep("Sample_char",colnames(sample_miMD))
sample_info_type <- unlist(strsplit(t(sample_miMD[1,col_info]),": "))[(1:length(col_info))*2-1]
temp <- data.frame(lapply(col_info,function(x) {unlist(strsplit(t(sample_miMD[,x]),": "))[(1:nrow(sample_miMD))*2]}))
sample_miMD[,col_info] <- temp
colnames(sample_miMD)[col_info] <- sample_info_type
sample_miMD[sample_miMD=="--"] <- NA
sample_miMD$`quit how many years ago`[sample_miMD$`quit how many years ago`=="nonsmoker"] <- NA
sample_miMD$`final diagnosis`[sample_miMD$`final diagnosis`=="Control"] <- "control"
sample_miMD2 <- sample_miMD[,c(1,2,col_info)]
for (i in 3:ncol(sample_miMD2)) {
  if (sum(!is.na(as.numeric(sample_miMD2[,i]))) > 0) {
    sample_miMD2[,i] <- as.numeric(sample_miMD2[,i])
  } else {
    sample_miMD2[,i] <- as.factor(sample_miMD2[,i])
  }
}

#calculate smoking score
sample_miMD2$`smoking status`<- as.character(sample_miMD2$`smoking status`)
sample_miMD2$`smoking status`[is.na(sample_miMD2$`smoking status`)]<-"NoData"
smoke_impact_score<-rep(0,nrow(sample_miMD2))
smoke_impact_score[sample_miMD2$`smoking status`=="NoData"]<-NA
smoke_impact_score[sample_miMD2$`smoking status`=="current"]<-(sample_miMD2$age[sample_miMD2$`smoking status`=="current"]-18)*sample_miMD2$`pack years`[sample_miMD2$`smoking status`=="current"]
smoke_impact_score[sample_miMD2$`smoking status`=="former"]<-(sample_miMD2$age[sample_miMD2$`smoking status`=="former"]-18-sample_miMD2$`quit how many years ago`[sample_miMD2$`smoking status`=="former"])*sample_miMD2$`pack years`[sample_miMD2$`smoking status`=="former"]/(sample_miMD2$`quit how many years ago`[sample_miMD2$`smoking status`=="former"]/1)
smoke_impact_score[smoke_impact_score<0] <- 0
smoke_impact_score<-log2(smoke_impact_score+1)
sample_miMD2$smoke_impact_score <- smoke_impact_score

#quick visualization of the data
summary(sample_miMD2$gender)
summary(sample_miMD2$`final diagnosis`)
summary(sample_miMD2$`tissue source`)
summary(sample_miMD2$preservative)
summary(as.factor(paste(sample_miMD2$gender,sample_miMD2$`final diagnosis`)))
boxplot(sample_miMD2$age ~ sample_miMD2$gender)
boxplot(sample_miMD2$age ~ sample_miMD2$`final diagnosis`) #age and diagnosis are confounding
boxplot(sample_miMD2$smoke_impact_score ~ sample_miMD2$gender)
boxplot(sample_miMD2$smoke_impact_score ~ sample_miMD2$`final diagnosis`)
plot(sample_miMD2$age,sample_miMD2$smoke_impact_score)

#compare mRNA and micro RNA (using Sample Title, as geo does not match)
sampleMD2$indID<-paste(sampleMD2$gender,sampleMD2$`final diagnosis`,sampleMD2$age,sampleMD2$`smoking status`,sampleMD2$smoke_impact_score,sampleMD2$`tissue source`,sampleMD2$rin,sep="-")
sample_miMD2$indID<-paste(sample_miMD2$gender,sample_miMD2$`final diagnosis`,sample_miMD2$age,sample_miMD2$`smoking status`,sample_miMD2$smoke_impact_score,sample_miMD2$`tissue source`,sample_miMD2$rin,sep="-")
sampleMD2$microRNA_data <- sampleMD2$Sample_title %in% sample_miMD2$Sample_title
sum(!sampleMD2$microRNA_data)
sample_miMD2$mRNA_data <- sample_miMD2$Sample_title %in% sampleMD2$Sample_title
sum(!sample_miMD2$mRNA_data)
matching <- merge(sampleMD2[,c(1:5,17)],sample_miMD2[,c(1:5,17)],by.x=1, all.x=T,by.y=1, all.y=T)
Not_matching <- matching[is.na(matching$Sample_geo_accession.x) | is.na(matching$Sample_geo_accession.y),]
matching <- merge(sampleMD2[,c(1:5,17)],sample_miMD2[,c(1:5,17)],by.x=1, all.x=F,by.y=1, all.y=F)
issue <- matching[matching$indID.x != matching$indID.y,][,c(6,11)] #the information about smoking are similar but not perfectly identical

#======================
#gene probe metadata

#import data
genepMD <- read.table(file = "GPL6244-17930.txt",header = T,skip = 12,quote = "\"",sep = "\t",comment.char = "",stringsAsFactors = F)
summary(as.factor(genepMD$category))

#parse gene name and function
x=1
parseGeneNames<-function(x) { #parsing function for lapply
  if(genepMD$gene_assignment[x]!="---") {
    geneID<-strsplit(genepMD$gene_assignment[x],split = " // ")[[1]][1]
    geneName<-strsplit(genepMD$gene_assignment[x],split = " // ")[[1]][2]
    funct<-strsplit(genepMD$gene_assignment[x],split = " // ")[[1]][3]
  } else if (genepMD$mrna_assignment[x]!="---") {
    geneID<-strsplit(genepMD$mrna_assignment[x],split = " // ")[[1]][1]
    geneName<-strsplit(genepMD$mrna_assignment[x],split = " // ")[[1]][2]
    funct<-strsplit(genepMD$mrna_assignment[x],split = " // ")[[1]][3]
  } else {
    geneID<-"none"
    geneName<-"none"
    funct<-"none"
  }
  res<-c(geneID,geneName,funct)
}
temp <- data.frame(t(data.frame(lapply(1:nrow(genepMD), parseGeneNames))))
colnames(temp) <- c("geneID","geneName","funct")
genepMD2 <- cbind(genepMD,temp)

#data exploration
sort(summary(as.factor(genepMD2$geneName)),decreasing = T)[1:20]
length(levels(as.factor(genepMD2$geneName)))
genepMDexp<-genepMD2[(genepMD2$geneName=="GOLGA6L10"),] #too look at why some gene have multiple entry

#======================
#mRNA probe metadata #long process directly using GEO
miRNApGEO <- getGEO(filename=  "GPL8786_family.soft.gz")
miRNApMD<- Table(miRNApGEO)


#======================
#save data
write.table(sampleMD2,"sampleMD2.txt",col.names = T,sep="\t",quote = F,row.names = F)
write.table(sample_miMD2,"sample_miMD2.txt",col.names = T,sep="\t",quote = F,row.names = F)
write.table(genepMD2,"genepMD.txt",col.names = T,sep="\t",quote = F,row.names = F)
write.table(miRNApMD,"miRNApMD.txt",col.names = T,sep="\t",quote = F,row.names = F)
write.table(matching,"matching.txt",col.names = T,sep="\t",quote = F,row.names = F)



