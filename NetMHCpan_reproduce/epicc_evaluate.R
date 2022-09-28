library(reshape2)
library('tidyverse')
library('tidyr')
library('reshape2')
library('ggplot2')
library('varhandle')
library('RColorBrewer')
library('gridExtra')
library(dplyr)


## evaluate EPICC with IEDB
epicc_classI<-function(classI_vax,classI_strain,remove){
  df1_1<-read.table(classI_vax,sep=",", header=TRUE,check.names=FALSE, stringsAsFacto=F, fill=TRUE,row.names = 1)
  df1_2<-read.table(classI_strain,sep=",", header=TRUE,check.names=FALSE, stringsAsFacto=F, fill=TRUE,row.names = 1)

  data=df1_2+df1_2
  
  data <- select(data, -remove)
  data<-data[!(row.names(data) %in% remove), ]
  return(data)
}

RSVA_F_remove<-c("KT992094-VACCINE_D46_D53","RSVB_PDA_ANCESTRAL","AF035006-VACCINE_RA2CP","AF013255-VACCINE_CP52")
RSVA_F_classI_epicc<-epicc_classI("EpiCC/RSV_EpiCC_Data/RSVA_F_classI_vac_unique.csv",
                            "EpiCC/RSV_EpiCC_Data/RSVA_F_classI_strain_unique.csv",
                            RSVA_F_remove)
RSVA_F_classI_epicc$seq1<-row.names(RSVA_F_classI_epicc)

## filter epicc dataframe
epicc_test<-melt(RSVA_F_classI_epicc)
epicc_test<-epicc_test %>% separate(seq1, into = c("Accession1", "Subtype1", "year1", "country1", "date1"), sep = "-", extra = "merge")
epicc_test<-epicc_test %>% separate(variable, into = c("Accession2", "Subtype2", "year2", "country2", "date2"), sep = "-", extra = "merge")
epicc_test<-epicc_test %>%
  select(Accession1, Accession2, value)
iedb_accession<-scan("EpiCC/EpiCC_confrim/RSVA_F/iedb_accession.txt",what="",sep="\n")
epicc_test2<- epicc_test %>%
  subset(Accession1 %in% iedb_accession) %>%
  subset(Accession2 %in% iedb_accession)
  
## change back to distance matrix
epicc_test2<-dcast(epicc_test2, Accession1~Accession2, value.var="value")
epicc_test2<-epicc_test2 %>% column_to_rownames(., var = "Accession1")  ## set the first column as index

mds_epicc<-epicc_test2%>%
  #dist()%>%
  cmdscale(k=2) %>%
  as_tibble()

colnames(mds_epicc)<-c("Dim.1","Dim.2")
k<-kmeans(mds_epicc,centers=3,nstart = 123)
clust <- k$cluster %>% as.factor()
mds_epicc <- mds_epicc %>%
  mutate(groups = clust)

#write.csv(epicc_test,"EpiCC/EpiCC_confrim/RSVA_F/RSVA_F_epicc.csv")
mds_epicc$accession <-rownames(epicc_test2)

###############################################################################

## match with iedb dataframe
###############################################################################
iedb_test<-read.table("EpiCC/EpiCC_confrim/RSVA_F/RSVA_F_predict1_7_iedb.csv",sep=",", header=TRUE,check.names=FALSE, stringsAsFacto=F, fill=TRUE,row.names = 1)
iedb_test$seq1<-row.names(iedb_test)
iedb_test<-melt(iedb_test)

accession_dic<-read_csv("EpiCC/EpiCC_confrim/RSVA_F/RSVA_F_accession_dic.csv")

## try transform seq name by merge
iedb_test2<-merge(iedb_test, accession_dic, by.x="seq1", by.y="seq_num", all.x=TRUE)
iedb_test2<-merge(iedb_test2, accession_dic, by.x="variable", by.y="seq_num", all.x=TRUE)
iedb_test2<-iedb_test2 %>%
  select(accession.x, accession.y, value)

## change back to distance matrix
iedb_test2<-dcast(iedb_test2, accession.x~accession.y, value.var="value")
iedb_test2<-iedb_test2 %>% column_to_rownames(., var = "accession.x")  ## set the first 

mds_iedb<-iedb_test2%>%
  #dist()%>%
  cmdscale(k=2) %>%
  as_tibble()

k_iedb<-kmeans(mds_iedb,centers=3,nstart = 123)
clust_iedb <- k_iedb$cluster %>% as.factor()
mds_iedb <- mds_iedb %>%
  mutate(groups_iedb = clust_iedb)

colnames(mds_iedb)<-c("iedb_Dim.1","iedb_Dim.2","group_iedb")
mds_iedb$accession <-rownames(iedb_test2)

## merge two dataframe to check cluster grouping
iedb_group<- merge(mds_epicc,mds_iedb,by=c("accession"="accession"))

## visulization
ggplot(iedb_group, aes(x=iedb_Dim.1, y=iedb_Dim.2, color=groups)) + 
  geom_point(size=1.5, alpha=.75)+
  xlab("Dimension 1")+
  ylab("Dimension 2")+
  scale_color_brewer(palette = "Set2",name="T-cell epitope immune clusters deterrmined by EPICC")+
  theme_bw()+
  theme(legend.position = "bottom")+
  geom_encircle(aes(group=group_iedb),expand=0,linetype=2)+
  ggtitle("MHC class I epitope landscape estimated with IEDB")
 # scale_shape_manual(values = c(19,5))+


ggplot(iedb_group, aes(x=Dim.1, y=Dim.2, color=groups)) + 
  geom_point(size=1.5, alpha=.75)+
  xlab("Dimension 1")+
  ylab("Dimension 2")+
  scale_color_brewer(palette = "Set2",name="T-cell epitope immune clusters deterrmined by EPICC")+
  # scale_shape_manual(values = c(19,5))+
  theme_bw()+
  theme(legend.position = "bottom")+
  geom_encircle(aes(group=groups),expand=0,linetype=2)+
  ggtitle("MHC class I epitope landscape estimated with EPICC")
  
###########################
## for original datamatrix, plot heatmap
epicc_test2$seq1<-row.names(epicc_test2)
epicc_hp<-melt(epicc_test2)
library("viridis")
#install.packages("hrbrthemes")
library(hrbrthemes)
ggplot(epicc_hp, aes(seq1, variable, fill= value)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  theme(legend.position="bottom",
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  xlab(" ")+
  ylab("")+
  labs(fill= "EPICC epitope immune distance")+
  coord_fixed(ratio = 1)+
  ggtitle("Pairwise MHC class I epitope distance estimated with EPICC")


iedb_test2$seq1<-row.names(iedb_test2)
iedb_hp<-melt(iedb_test2)

ggplot(iedb_hp, aes(seq1, variable, fill= value)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  theme(legend.position="bottom",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  xlab("")+
  ylab("")+
  labs(fill= "IEDB epitope immune distance")+
  coord_fixed(ratio = 1)+
  ggtitle("Pairwise MHC class I epitope distance estimated with IEDB")
  #theme_ipsum()
############################################################################
## Eigen value decoposition
#install.packages("matlib")
library(matlib)
#ev_epicc<-eigen(epicc_test2)
#L_epicc <- ev_epicc$values
#V_epicc <- ev_epicc$vectors
library(RSpectra)
epicc_ev<-data.matrix(epicc_test2)
V_epicc<-eigs_sym(epicc_ev, k=210,opts = list(retvec = FALSE))

iedb_ev<-data.matrix(iedb_test2)
V_iedb<-eigs_sym(iedb_ev, k=210,opts = list(retvec = FALSE))

df <- do.call(rbind,mapply(cbind, V_epicc, V_iedb))
df <- as.data.frame(df, stringsAsFactors = FALSE)
#ev_iedb<-eigen(iedb_test2)
#L_iedb <- ev_iedb$values
ggplot(df, aes(x=V1, y=V2)) + 
  geom_point()+
  stat_smooth(method = lm)+
  theme_bw()+
  xlab("Eigenvalue of MHC class I epitope distance esitmate with EPICC")+
  ylab("Eigenvalue of MHC class I epitope distance esitmate with IEDB")

res <- cor.test(df$V1, df$V2, 
                method = "pearson")
res


