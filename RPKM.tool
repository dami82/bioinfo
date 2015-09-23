#RPKM tool
#url::http://www.biotechworld.it/bioinf/2015/05/04/5/
#Damiano Fantini
#2015-09-21

#g.o.i
my_genes<- c("GAPDH","XRCC6","CUL4B", "LEF1", "TCF7L2")

#input data
setwd("R/homemade_tools/RPKM_tool/") # change to project directory... must contain raw reads csv file
list.files()
read_counts <- read.csv("counts.csv", header = T, as.is=T)

#load required libraries
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

#retrieve the total exon lenght of each gene of interest (feed symbol)
gene_lookup <- select(org.Hs.eg.db,keys=my_genes, keytype="SYMBOL",columns=c("ENTREZID","SYMBOL","GENENAME","CHR","CHRLOC"))
tx_lookup <- select(TxDb.Hsapiens.UCSC.hg19.knownGene, key=gene_lookup$ENTREZID, keytype="GENEID", columns= c("TXID","TXSTART","TXEND","EXONSTART","EXONEND","EXONID","TXSTRAND","EXONCHROM"))
gene_IDs <- unique(tx_lookup$GENEID)
results <- matrix(,ncol=length(my_genes), nrow=0)
colnames(results) <- my_genes

results <- rbind(results,
                 sapply (gene_IDs, (function(id){
                   ex_ids <- unique(tx_lookup[which(tx_lookup$GENEID == id),"EXONID"])
                   
                   sum( 
                     sapply(ex_ids, (function(ex){
                       exons <- tx_lookup[which(tx_lookup$EXONID==ex &
                                                  tx_lookup$GENEID==id),]
                       as.numeric(as.matrix(exons)[1,"EXONEND"])-as.numeric(as.matrix(exons)[1,"EXONSTART"])
                     })))
                 }))
)
rownames(results)[1] <- "all_exon_len"

#define treatment groups (CTRL vs. siRNA)
ctrls <- grep("CTRL",colnames(read_counts), ignore.case = T)
treat <- grep("CTRL",colnames(read_counts), ignore.case = T, invert = T)
treat <- treat[-which(colnames(read_counts)=="X")]

#Retrieve raw counts
sel_counts <- read_counts[which(read_counts$X %in% my_genes),]
results <- rbind(results,as.vector(
  sapply(1:ncol(results),(function(r){
    if (sum(sel_counts[which(sel_counts$X==colnames(results)[r]),ctrls]==0)>0 | sum(sel_counts[which(sel_counts$X==colnames(results)[r]),treat]==0)>0) {
      0
    } else {
      mean(as.numeric(as.character(sel_counts[which(sel_counts$X==colnames(results)[r]),c(ctrls, treat)])))
    }
  }))
))
rownames(results)[2] <- "counts"

#Attach total counts (across all samples)
rw_reads<-which(
  sapply(1:nrow(read_counts),(function(rrw){
    if(sum(read_counts[rrw,ctrls]>0) == length(ctrls) |
       sum(read_counts[rrw,treat]>0) == length(treat)) { T } else { F }
  }))
)
results <- rbind(results,
                 rep(sum(apply(read_counts[as.numeric(as.character(rw_reads)),c(ctrls, treat)],1,sum)),
                     ncol(results)))
rownames(results)[3] <- "tot_reads"

results <- rbind(results, sapply(1:ncol(results), (function(i){
  (10^9)*results["counts",i]/(results["tot_reads",i]*results["all_exon_len",i]) 
})))
rownames(results)[4] <- "RPKM"

barplot((results["RPKM",])+0.01,log="y")
write.csv(results, "rpkm.csv")
