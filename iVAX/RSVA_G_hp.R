##  Map RSVA -G class II epitope 9 mer to LK tree
library(treeio)
library(ggplot2)
library(ggtree)
library(phytools)
library(ape)
library(RColorBrewer)
## RSVA ML 
t1 <- read.tree(file = "RAxML/RSV_G/RSVA_G_date2/RAxML_bestTree.HRSVA_G.tree")
p1<-ggtree(t1)+
  geom_treescale(x=0.1, y=3,linesize= 1,,width=0.02)
p1

## add epitope
epitope<- read.table("iVAX/RSVA_G_analysis/HRSVA_G_hp_reorder.txt",  sep="\t", header=TRUE,check.names=FALSE, stringsAsFactor=F, fill=TRUE,row.names = 1)
##convert epitope as factor
for(i in 1:ncol(epitope)){
  epitope[,i] <- as.factor(epitope[,i])
}

#########################################################################
gheatmap(p1,epitope,offset = 0.018, width=2, font.size=2, colnames_position= "bottom", colnames_angle = 90, colnames_offset_y = 0, hjust = 0)+
  scale_fill_brewer(palette = "Set2") 
