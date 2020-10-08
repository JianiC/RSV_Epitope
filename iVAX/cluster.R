library('tidyverse')
library('tidyr')
library('reshape2')
library('ggplot2')
library('varhandle')
library('RColorBrewer')
library('gridExtra')
library(dplyr)


## cluster analysis for RSVB_F 
## summarize the cluste File
sum1<-read.csv("ClusterReport/RSVB_F_cluster_test1.csv",header=FALSE)[-1:-10,]
sum2<-read.csv("ClusterReport/RSVB_F_cluster_test2.csv",header=FALSE)[-1:-10,]

## combine two files
sum<-rbind (sum1,sum2)

colnames(sum)<-c("File","seq","cluster.Address","cluster.AA","Hydrophobicity","Hits","Score1","Score2","cores")
## remove the empty seq rows
sum<- sum %>% filter(!cluster.AA == "")

## get the conserved cluster seq
n<-n_distinct(sum$seq) ## No. of sequences
## get the conserved cluster sequence
cluster_sum1<-sum %>%
  group_by(cluster.Address,cluster.AA) %>%
  dplyr::summarise(count = n())

cluster_sum1$coverage<-cluster_sum1$count/n
## filter the cluster sequence by coverage
conserved_cluster<-cluster_sum1%>%
  filter(coverage>=0.05)

write.csv(conserved_cluster,"ClusterReport/RSVB_F_conserved_cluster.csv",row.names = FALSE)