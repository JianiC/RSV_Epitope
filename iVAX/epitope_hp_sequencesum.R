## summarise epitope for heatmap generation
library(tidyverse)
## class I, try with a clean dataframe
df<-read.csv("ClusterReport/PRINT_CLUSTER_SUMMARY_HRSVA_F_CLEAN0627_PROTEIN_All_Proteins.csv",header=FALSE)[-1:-10,]
colnames(df)<-c("File","seq","cluster.Address","cluster.AA","Hydrophobicity","Hits","Score1","Score2","cores")
df<- df %>% filter(!cluster.AA == "")

## clacluate coverage
df<-df%>% distinct(seq,cluster.Address,cluster.AA)
df2<-df%>% 
  group_by(cluster.Address) %>%
  dplyr::summarise(count = n())

df3<-df%>%group_by(cluster.Address,cluster.AA) %>%
  dplyr::summarise(count2 = n())

n<-n_distinct(df$seq)
df2<- df2%>% mutate(coverage = count/as.numeric(n))
df3<- df3%>% mutate(coverage2 = count2/as.numeric(n))


df4<-left_join(df3,df2,by=c("cluster.Address"="cluster.Address"))


df4<-df4%>%filter(coverage>=0.01)

#df4<-df4 %>% group_by(address)%>%slice_max(coverage2, n = 4) ## only keep top 4 epitope
## create id within group
df4<-df4[order(-df4$coverage2),]
df4<-df4 %>% group_by(cluster.Address) %>% mutate(id = row_number(cluster.Address))

df4<-df4%>%filter(id<=4)
## reshape dataframe
df4<-df4%>%select(cluster.Address, cluster.AA,coverage,id)
df4<-dcast(df4,coverage+cluster.Address~id,value.var="cluster.AA")



######################################################
## prepare a function

epitope_sum<-function(datafile,subtype,protein){
  df<-read.csv(datafile)
  ## clacluate coverage
  df<-df%>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)
  df2<-df%>% 
    group_by(Frame.Start, Frame.Stop) %>%
    dplyr::summarise(count = n())
  
  df3<-df%>%group_by(Frame.Start, Frame.Stop,AA.Sequence) %>%
    dplyr::summarise(count2 = n())
  
  n<-n_distinct(df$Accession)
  df2<- df2%>% mutate(coverage = count/as.numeric(n))
  df3<- df3%>% mutate(coverage2 = count2/as.numeric(n))
  df4<-left_join(df3,df2,by=c("Frame.Start"="Frame.Start","Frame.Stop"="Frame.Stop"))
  df4<-df4%>%filter(coverage>=0.01)
  
  df4$address  = str_c(df4$Frame.Start,"-",df4$Frame.Stop)
  #df4<-df4 %>% group_by(address)%>%slice_max(coverage2, n = 4) ## only keep top 4 epitope
  ## create id within group
  df4<-df4[order(-df4$coverage2),]
  df4<-df4 %>% group_by(address) %>% mutate(id = row_number(address))
  df4<-df4%>%filter(id<=4)
  ## reshape dataframe
  df4<-df4%>%select(address, AA.Sequence,coverage,id)
  df4<-dcast(df4,coverage+address~id,value.var="AA.Sequence")
  df4$subtype=subtype
  df4$protein=protein
  return(df4)
}

RSVA_F_classI<-epitope_sum("iVAX/RSVA_F/RSVA_F_classI_analysis/RSVA_F_classI_epimax_filter.csv","RSV-A","F")
RSVB_F_classI<-epitope_sum("iVAX/RSVB_F/RSVB_F_classI_sum/HRSVB_F_classI_epimax_filter.csv","RSV-B","F")
RSVA_G_classI<-epitope_sum("iVAX/RSVA_G/RSVA_G_classI/RSVA_G_classI_epimax_filter.csv","RSV-A","G")
RSVB_G_classI<-epitope_sum("iVAX/RSVB_G/RSVB_G_classI_epimax_filter.csv","RSV-B","G")

classI_epitope<-rbind(RSVA_F_classI,RSVB_F_classI,RSVA_G_classI,RSVB_G_classI)
classI_epitope$class<-"MHC class I"

###############################################################################
## sum class II cluster epitope sequence 
cluster_sum<-function(clusterfile,subtype,protein){
  df<-read.csv(clusterfile,header=FALSE)[-1:-10,]
  colnames(df)<-c("File","seq","cluster.Address","cluster.AA","Hydrophobicity","Hits","Score1","Score2","cores")
  df<- df %>% filter(!cluster.AA == "")
  
  ## clacluate coverage
  df<-df%>% distinct(seq,cluster.Address,cluster.AA)
  df2<-df%>% 
    group_by(cluster.Address) %>%
    dplyr::summarise(count = n())
  
  df3<-df%>%group_by(cluster.Address,cluster.AA) %>%
    dplyr::summarise(count2 = n())
  
  n<-n_distinct(df$seq)
  df2<- df2%>% mutate(coverage = count/as.numeric(n))
  df3<- df3%>% mutate(coverage2 = count2/as.numeric(n))
  
  df4<-left_join(df3,df2,by=c("cluster.Address"="cluster.Address"))
  df4<-df4%>%filter(coverage>=0.01)
  
  #df4<-df4 %>% group_by(address)%>%slice_max(coverage2, n = 4) ## only keep top 4 epitope
  ## create id within group
  df4<-df4[order(-df4$coverage2),]
  df4<-df4 %>% group_by(cluster.Address) %>% mutate(id = row_number(cluster.Address))
  df4<-df4%>%filter(id<=4)
  ## reshape dataframe
  df4<-df4%>%select(cluster.Address, cluster.AA,coverage,id)
  df4<-dcast(df4,coverage+cluster.Address~id,value.var="cluster.AA")
  df4$subtype=subtype
  df4$protein=protein
  
  return(df4)
  
}

RSVA_F_classII<-cluster_sum("ClusterReport/PRINT_CLUSTER_SUMMARY_HRSVA_F_CLEAN0627_PROTEIN_All_Proteins.csv","RSV-A","F")
RSVB_F_classII<-cluster_sum("ClusterReport/RSVB_F_cluster_combine.csv","RSV-B","F")
RSVA_G_classII<-cluster_sum("ClusterReport/PRINT_CLUSTER_SUMMARY_HRSVA_G_CLEAN0819_PROTEIN_All_Proteins.csv","RSV-A","G")
RSVB_G_classII<-cluster_sum("ClusterReport/PRINT_CLUSTER_SUMMARY_HRSVB_G_CLEAN_0819_PROTEIN_All_Proteins.csv","RSV-B","G")

classII_epitope<-rbind(RSVA_F_classII,RSVB_F_classII,RSVA_G_classII,RSVB_G_classII)
classII_epitope$class<-"MHC class II"

write_csv(classI_epitope,"iVAX/RSV_classI_epitope_hp_sequence.csv")
write_csv(classII_epitope,"iVAX/RSV_classII_epitope_hp_sequence.csv")