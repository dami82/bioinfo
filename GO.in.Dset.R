## R Script
## Retrieve gene expression levels (Log2FoldChange) of genes belonging 
## to a GO gene family of interest 
##
## Author: Damiano Fantini
## 29 oct 2015
## 

# load packages
source("http://bioconductor.org/biocLite.R")
library(GO.db)
library(org.Hs.eg.db)

# set WD for this project
setwd("../Documents/R/homemade_tools/leyo/")

# define keywords we want to search among FO terms
my_keyword<-c("protein acetylation","deacetylase", "deacetylation")

# We are using here a dataset downloaded from Oncomine
# The dataset comes as a couple of csv files, as the analyses for 
# upregulated genes and downregulated genes come separately
# We first have to prepare data by merging the two files in 
# one data frame. We have also to take care of the following problems:
# 1) probes -> symbols (one to many): remove all
# 2) duplicated symbols: check before starting the analysis and decide what to keep
# 3) Weird FoldChange column: transform in Log2FC (Oncomine provides values that are
#    >1 for upreg genes and < (-1) for downreg genes)
# 
# csv files are not provided, but you can easily download them from Oncomine

dataset_UP <- read.csv("Crabtree_Uterus_Gene_List_up.csv", header = FALSE, as.is=TRUE)
colnames(dataset_UP) <- dataset_UP[3,]
dataset_UP<- dataset_UP[-c(1,2,3),]

dataset_DOWN <- read.csv("Crabtree_Uterus_Gene_List_down.csv", header = FALSE, as.is=TRUE)
colnames(dataset_DOWN) <- dataset_DOWN[3,]
dataset_DOWN<- dataset_DOWN[-c(1,2,3),]

dataset_UP<- dataset_UP[as.numeric(dataset_UP$`Fold Change`)>0,]
dataset_DOWN<- dataset_DOWN[as.numeric(dataset_DOWN$`Fold Change`)<0,]
dataSet <- rbind(dataset_UP, dataset_DOWN)

dupli_genes <- dataSet[duplicated(dataSet$`Gene Symbol`),'Gene Symbol']
dataSet <- dataSet[!(dataSet$`Gene Symbol` %in% dupli_genes),]
dataSet <- dataSet[order(dataSet$`Gene Symbol`),]
log2FC<- sapply(as.numeric(dataSet$`Fold Change`), (function(fc){
  if ( fc >0) {
    log2(fc)
  } else {
    log2(1/(fc*(-1)))
  }
}))
dataSet <- cbind(dataSet,log2FC)
zScores <- scale(log2FC)
dataSet <- cbind(dataSet, zScores)

# DataSet is now ready to be analysed. The following line should return a volcano plot
plot(-log10(as.numeric(dataSet$`Q-value`))~((dataSet$log2FC)))

# Let's start working with the GO database
# We have to first extract data from GO.db
GOdb_keys <- keys(GO.db, keytype="GOID")
GOdb_list <- select(GO.db, keys=GOdb_keys, columns=c("GOID","TERM"), keytype = "GOID")

# Let's search for the terms of interest and retrieve only GOIDs of interest
# We use sapply as grep takes only one "patter" per time
my_goIDs <- sapply(my_keyword,(function(k_word){
  my_GOlist <- GOdb_list[grep(k_word, GOdb_list$TERM, value=FALSE),]
  my_GOlist$GOID 
}))

# We need a vector of GO IDs, but now my_goIDs is a list
my_goIDs <- unique(as.vector(unlist(my_goIDs)))

# Retireve human genes belonging to those terms
gene_list<-select(org.Hs.eg.db, keys=my_goIDs, columns=c("ENTREZID","SYMBOL","GOALL"), keytype = "GOALL")
# Let's remove NA
gene_list<-gene_list[!is.na(gene_list$ENTREZID),]
# And now let's save a txt file that includes all genes (SYMBOLS) we are looking for
write(as.character(unique(gene_list$SYMBOL)), "gene_shortlist.txt", sep="\n")

# Let's retrieve only DataSet genes included in our GO sene shortlist
GO_Set <- dataSet[dataSet$`Gene Symbol` %in% gene_list$SYMBOL,]
GO_Set <- GO_Set[order(as.numeric(GO_Set$log2FC)),]
rownames(GO_Set) <- GO_Set$`Gene Symbol`
GO_Set <- GO_Set[,4:ncol(GO_Set)]
# Now Go_Set contains all genes of interest (DE and non-DE genes)
# Let's save a csv copy
write.csv(GO_Set, "GOgenes_in_Dataset.csv", sep= ",")

# Let's remove what is non DE and write
GO_Set <- GO_Set[as.numeric(GO_Set$`Q-value`)<0.05 & abs(as.numeric(GO_Set$log2FC))>0.5,]
write.csv(GO_Set, "DE_GOgenes_in_Dataset.csv", sep= ",")
# Funny Volcano plot
# adjust offset based on the plot
setx_offset = 0.04
y_offset = 0.25
pdf("volcano.pdf", 6,8)
plot(-log10(as.numeric(dataSet$`Q-value`))~((dataSet$log2FC)), 
     pch=19, col="darkgrey", cex=0.5, main="Volcano Plot")
for (i in 1:nrow(GO_Set)) {
  my_x <- GO_Set[i,"log2FC"]
  my_y <- -log10(as.numeric(GO_Set[i,"Q-value"]))
  if (my_x > 0) { 
    my_col = "red" 
    my_adj = 0
    x_offset = setx_offset
  } else { 
    my_col = "darkgreen" 
    my_adj = 1
    x_offset = setx_offset * (-1)
  }  
  points(my_x,my_y,pch=19, col=my_col)
  text(my_x + x_offset, my_y + y_offset, rownames(GO_Set)[i],cex=0.86, col=my_col, adj= my_adj)
}
dev.off()
#Done! Success! :)

