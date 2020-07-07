## map epitope profile to a tree
library(treeio)
library(ggplot2)
library(ggtree)
library(phytools)
library(ape)
library(RColorBrewer)

## RSVA ML color with region and bootstrap
t1 <- read.tree(file = "RAxML/HRSVB_F/RAxML_bipartitions.HRSVB_F.tree")
p1<-ggtree(t1,branch.length="none") +
  geom_nodepoint(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) >= 70),color='black')
p1

##add metadata
metadata <- read.table("RAxML/HRSVB_F/metadata.csv", sep=",", header=TRUE,check.names=FALSE, stringsAsFactor=F, fill=TRUE)

p2<-p1 %<+% metadata+ aes(color=genotype_ref)

p2
epitope<- read.table("RAxML/HRSVB_F/HRSVB_F_classII_epimx_hpfull.txt",  sep="\t", header=TRUE,check.names=FALSE, stringsAsFactor=F, fill=TRUE,row.names = 1)
##convert epitope as factor
for(i in 1:ncol(epitope)){
  epitope[,i] <- as.factor(epitope[,i])
}


#########################################################################
gheatmap(p2,epitope,offset = 0.018, width=4, font.size=0, colnames_position= "top", colnames_angle = 90, colnames_offset_y = 0, hjust = 0)+
  scale_fill_viridis_d(option="D", name="epitope allel count")


