# Author: Christian Brion - 2022 - project IPF_DE
#
#-load and prepare DE gene data for enrichment
#-enrichment analysis on all DE gene using clusterProfiler: GO(BP,MF,CC)
#-data visualization using enrichplot
#-kegg enrichment analysis on all DE gene and data visualization using pathview
#-GO)BP)enrichment analysis each four clusters
#
#inspired from https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

#initiation: install, load packages, define function, set work directory
# BiocManager::install("clusterProfiler")
# BiocManager::install("pathview")
# BiocManager::install("enrichplot",force = TRUE)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(pathview)


# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE,force = TRUE)
library(organism, character.only = TRUE)


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

#===================
#load metadata

sampleMD2 <- read.table(file = "sampleMD2.txt",header = T,quote = "\"",sep = "\t",stringsAsFactors = T)
load("geneexp_processeddata.RData")

#===================
#load clusterData and DE data

geneTestvalue <- read.table(file = "geneTestvalue.txt",header = T,quote = "\"",sep = "\t",stringsAsFactors = T)
geneclust <- read.table(file = "geneclust.txt",header = T,quote = "\"",sep = "\t",stringsAsFactors = T)
sampleGeneclust <- read.table(file = "sampleGeneclust.txt",header = T,quote = "\"",sep = "\t",stringsAsFactors = T)


#====================================
# general clusterProfiler enrichment
#====================================

keytypes(org.Hs.eg.db)

pdf("all_enrichment.pdf",width = 14,height = 14) #save all the plot in one PDF

# reading in data from deseq2

# we want the log2 fold change 
sum(rownames(geneTestvalue)==genepMDunique$ID) #good
original_gene_list <- geneTestvalue$IPF.UIP
names(original_gene_list) <- genepMDunique$geneName # name the vector

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# only keep DE genes
gene_list <- gene_list[names(gene_list) %in% genepMDunique$geneName[genepMDunique$ID %in% geneclust$probID]]

#go enrichment
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", #"SYMBOL"
             minGSSize = 3, 
             maxGSSize = 200, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH")

#enrichment visualization
require(DOSE)
dotplot(gse, showCategory=20, split=".sign",title = "all BP") + facet_grid(.~.sign)
gsepairsim<-pairwise_termsim(gse)
emapplot(gsepairsim, showCategory = 20)
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 6)
ridgeplot(gse) + labs(x = "enrichment distribution")

#go enrichment
gse <- gseGO(geneList=gene_list, 
             ont ="MF", 
             keyType = "SYMBOL", #"SYMBOL"
             minGSSize = 3, 
             maxGSSize = 200, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH")

#enrichment visualization
dotplot(gse, showCategory=20, split=".sign",title = "all MF") + facet_grid(.~.sign)
gsepairsim<-pairwise_termsim(gse)
emapplot(gsepairsim, showCategory = 20)
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 6)
ridgeplot(gse) + labs(x = "enrichment distribution")

#go enrichment
gse <- gseGO(geneList=gene_list, 
             ont ="CC", 
             keyType = "SYMBOL", #"SYMBOL"
             minGSSize = 3, 
             maxGSSize = 200, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH")

#enrichment visualization
dotplot(gse, showCategory=20, split=".sign",title = "all CC") + facet_grid(.~.sign)
gsepairsim<-pairwise_termsim(gse)
emapplot(gsepairsim, showCategory = 20)
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 6)
ridgeplot(gse) + labs(x = "enrichment distribution")

#=======================
#kegg pathway

kegg_organism = "hsa"
kegg_gene_list <- gene_list
temp <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism) #convert gene name into entrez, mandatory for kegg
kegg_gene_list<-kegg_gene_list[names(kegg_gene_list) %in% temp$SYMBOL] #lost 123 (1570-1447) genes
sum(names(kegg_gene_list)==temp$SYMBOL) #it match
kegg_gene_list2<-kegg_gene_list
names(kegg_gene_list2)<-temp$ENTREZID


#kegg enrichment
kk2 <- gseKEGG(geneList     = kegg_gene_list2,
               organism     = kegg_organism,
               minGSSize    = 3,
               maxGSSize    = 200,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               keyType       = "ncbi-geneid")

#enrichment visualization
dotplot(kk2, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
kk2pairsim<-pairwise_termsim(kk2)
emapplot(kk2pairsim)
cnetplot(kk2, categorySize="pvalue", foldChange=gene_list, showCategory = 4)
ridgeplot(kk2) + labs(x = "enrichment distribution")

#pathway visualization with pathview
setwd(paste0(workdir,"/pathview"))
#Neuroactive ligand-receptor interaction
dme <- pathview(gene.data=kegg_gene_list2, pathway.id="hsa04080", species = kegg_organism, 
                low = list(gene = "deepskyblue", cpd = "green"), mid = list(gene = "gray", cpd = "gray"), high = list(gene = "yellow", cpd = "red"))
#Neutrophil extracellular trap formation
dme <- pathview(gene.data=kegg_gene_list2, pathway.id="hsa04613", species = kegg_organism,
                low = list(gene = "deepskyblue", cpd = "green"), mid = list(gene = "gray", cpd = "gray"), high = list(gene = "yellow", cpd = "red"))
#NF-kappa B signaling pathway
dme <- pathview(gene.data=kegg_gene_list2, pathway.id="hsa04064", species = kegg_organism,
                low = list(gene = "deepskyblue", cpd = "green"), mid = list(gene = "gray", cpd = "gray"), high = list(gene = "yellow", cpd = "red"))
#VEGF signaling pathway
dme <- pathview(gene.data=kegg_gene_list2, pathway.id="hsa04370", species = kegg_organism,
                low = list(gene = "deepskyblue", cpd = "green"), mid = list(gene = "gray", cpd = "gray"), high = list(gene = "yellow", cpd = "red"))
#Hedgehog signaling pathway
dme <- pathview(gene.data=kegg_gene_list2, pathway.id="hsa04340", species = kegg_organism,
                low = list(gene = "deepskyblue", cpd = "green"), mid = list(gene = "gray", cpd = "gray"), high = list(gene = "yellow", cpd = "red"))
#TGF-beta signaling pathway
dme <- pathview(gene.data=kegg_gene_list2, pathway.id="hsa04350", species = kegg_organism,
                low = list(gene = "deepskyblue", cpd = "green"), mid = list(gene = "gray", cpd = "gray"), high = list(gene = "yellow", cpd = "red"))
setwd(workdir)

#=========
#some genes found on the web for control

GOI<-c("MUC5B", "SFTPA1", "SFTPA2B", "DKC1") # mucin , pulmonary surfactant proteins 1/2, dyskerin 

geneTestvalue$names<-genepMDunique$geneName
geneTestvalue[geneTestvalue$names %in% GOI,]

#====================================
# percluster clusterProfiler enrichment
#====================================

#separate DE genes per clusters
gene_list1 <- gene_list[names(gene_list) %in% genepMDunique$geneName[genepMDunique$ID %in% geneclust$probID[geneclust$clust %in% c(2)]]]
gene_list2 <- gene_list[names(gene_list) %in% genepMDunique$geneName[genepMDunique$ID %in% geneclust$probID[geneclust$clust %in% c(5,6,8)]]]
gene_list3 <- gene_list[names(gene_list) %in% genepMDunique$geneName[genepMDunique$ID %in% geneclust$probID[geneclust$clust %in% c(4)]]]
gene_list4 <- gene_list[names(gene_list) %in% genepMDunique$geneName[genepMDunique$ID %in% geneclust$probID[geneclust$clust %in% c(1,7,3)]]]

#go enrichment
gene_list_clust <- gene_list1
GOenrich <- enrichGO(gene=names(gene_list_clust),
                     OrgDb = organism, 
                     ont ="BP", 
                     keyType = "SYMBOL", #"SYMBOL"
                     minGSSize = 3, 
                     maxGSSize = 200, 
                     pvalueCutoff = 0.05, 
                     pAdjustMethod = "BH")

#enrichment visualization
dotplot(GOenrich, showCategory=20,title = "clustI BP")
cnetplot(GOenrich, categorySize="pvalue", foldChange=gene_list_clust, showCategory = 6)


#go enrichment
gene_list_clust <- gene_list2
GOenrich <- enrichGO(gene=names(gene_list_clust),
                     OrgDb = organism, 
                     ont ="BP", 
                     keyType = "SYMBOL", #"SYMBOL"
                     minGSSize = 3, 
                     maxGSSize = 200, 
                     pvalueCutoff = 0.05, 
                     pAdjustMethod = "BH")

#enrichment visualization
dotplot(GOenrich, showCategory=20,title = "clustII BP")
cnetplot(GOenrich, categorySize="pvalue", foldChange=gene_list_clust, showCategory = 6)

#go enrichment
gene_list_clust <- gene_list3
GOenrich <- enrichGO(gene=names(gene_list_clust),
                     OrgDb = organism, 
                     ont ="BP", 
                     keyType = "SYMBOL", #"SYMBOL"
                     minGSSize = 3, 
                     maxGSSize = 200, 
                     pvalueCutoff = 0.05, 
                     pAdjustMethod = "BH")

#enrichment visualization
dotplot(GOenrich, showCategory=20,title = "clustIII BP")
cnetplot(GOenrich, categorySize="pvalue", foldChange=gene_list_clust, showCategory = 6)

#go enrichment
gene_list_clust <- gene_list4
GOenrich <- enrichGO(gene=names(gene_list_clust),
                     OrgDb = organism, 
                     ont ="BP", 
                     keyType = "SYMBOL", #"SYMBOL"
                     minGSSize = 3, 
                     maxGSSize = 200, 
                     pvalueCutoff = 0.05, 
                     pAdjustMethod = "BH")

#enrichment visualization
dotplot(GOenrich, showCategory=20,title = "clustIV(down) BP")
cnetplot(GOenrich, categorySize="pvalue", foldChange=gene_list_clust, showCategory = 6)

Sys.sleep(60) #wait for the figure to load into the pdf (slow computer issue)
dev.off() #close pdf
