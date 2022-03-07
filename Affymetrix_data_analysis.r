# Author: Christian Brion - 2022 - project IPF_DE
#
#-import raw expression data from downloaded affymetrix .CEL files with affy/affyPLM 
#-quick data visualization for QC
#-normalization of the data using RMA function
#-quick post normalization data visualization for QC
#-export the data as a RData
#
#code inspired from https://wiki.bits.vib.be/index.php/Analyze_your_own_microarray_data_in_R/Bioconductor#Install_the_required_R_packages

#initiation: install, load packages, define function, set work directory
# BiocManager::install("affy",force = T) #for importing affymetrix Gene ST datas
# BiocManager::install("affyPLM ",force = T) #for importing affymetrix Gene ST datas
# BiocManager::install("limma",force = T)
# BiocManager::install("Biobase",force = T)
# BiocManager::install("Biostrings",force = T)
# BiocManager::install("genefilter",force = T)

library(affy)
library(affyPLM)
library(Biobase)
library(Biostrings)
library(genefilter)

library(Matrix)
library(lattice)
library(fdrtool)
library(rpart)

colSdv<-function(df) {
  unlist(lapply(X = 1:ncol(df), function(x) {sd(df[x,])}))
}

workdir <- "C:/Users/chris/Dropbox/Github/IPF_DE" #set work directory
setwd(workdir)

#=========first mRNA
# import CEL files containing raw probe-level data into an R AffyBatch object
data <- ReadAffy(celfile.path = "GSE32537_RAW/")

#rapid quality control
int <- intensity(data)
nrow(int)
ncol(int)
logint <- log2(int)
plot(colMeans(logint),colSdv(logint))

svg(filename = "array_mRNA_stdv_issue.svg",width = 6,height = 5)
barplot(colSdv(logint),ylab="standard dev",xlab="mRNA array samples")
dev.off()

data.rma <- rma(data) #normalizing data (Background correction, Log2 transformation, Quantile normalization, Probe normalization)
data.matrix = exprs(data.rma)

#rapid quality control
plot(colMeans(data.matrix),colSdv(data.matrix))

plot(colMeans(logint),colMeans(data.matrix))
plot(colSdv(logint),colSdv(data.matrix))

logint.PC = prcomp(t(logint),scale.=TRUE)
logint.PC$percent <- logint.PC$sdev/sum(logint.PC$sdev)
plot(logint.PC$x[,1:2])

data.PC = prcomp(t(data.matrix),scale.=TRUE)
data.PC$percent <- data.PC$sdev/sum(data.PC$sdev)
plot(data.PC$x[,1:2])

#MAplot(data,which=1) #too big?
MAplot(data.rma,which=1)

save(data.matrix,file = "gene_exp_data.RData")


#=========second microRNA
# import CEL files containing raw probe-level data into an R AffyBatch object
data_mirna <- ReadAffy(celfile.path = "GSE32538_RAW/")

#rapid quality control
int2 <- intensity(data_mirna)
nrow(int2)
ncol(int2)
logint2 <- log2(int2)
plot(colMeans(logint2),colSdv(logint2))

svg(filename = "array_microRNA_stdv_issue.svg",width = 6,height = 5)
barplot(colSdv(logint2),ylab="standard dev",xlab="microRNA array samples")
dev.off()

data_mirna.rma <- rma(data_mirna) #normalizing data (Background correction, Log2 transformation, Quantile normalization, Probe normalization)
data_mirna.matrix = exprs(data_mirna.rma)

#rapid quality control
plot(colMeans(data_mirna.matrix),colSdv(data_mirna.matrix))

plot(colMeans(logint2),colMeans(data_mirna.matrix))
plot(colSdv(logint2),colSdv(data_mirna.matrix))

logint2.PC = prcomp(t(logint2),scale.=TRUE)
logint2.PC$percent <- logint2.PC$sdev/sum(logint2.PC$sdev)
plot(logint2.PC$x[,1:2])

data_mirna.PC = prcomp(t(data_mirna.matrix),scale.=TRUE)
data_mirna.PC$percent <- data_mirna.PC$sdev/sum(data_mirna.PC$sdev)
plot(data_mirna.PC$x[,1:2])

MAplot(data_mirna,which=1)
MAplot(data_mirna.rma,which=1)

save(data_mirna.matrix,file = "miRNA_data.RData")
