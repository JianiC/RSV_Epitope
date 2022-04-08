## RSVA G epicc visulize
library(reshape2)
library(ggplot2)
## read the raw dataframe
raw1<-read.table("EpiCC/RSV_EpiCC_Data/RSVA_G_classI_epicc_unique.csv", 
                 sep=",", header=TRUE,check.names=FALSE, stringsAsFacto=F, fill=TRUE,row.names = 1)
raw2<- read.table("EpiCC/RSV_EpiCC_Data/RSVA_G_classII_epicc_unique.csv", 
                  sep=",", header=TRUE,check.names=FALSE, stringsAsFacto=F, fill=TRUE,row.names = 1)

raw<-raw1+raw2



###################################################################  
## with cross T-cell immunity disimilarity build mds
## remove vaccine column and rows  
row.names.remove <- c("KT992094-VACCINE_D46_D53")

data <- select(raw, -row.names.remove) 
vac<-c("U63644-VACCINE_CPTS-248_404","AF013255-VACCINE_CP52","RSVA_PDA_ANCESTRAL","JX198138-A-1961-AUSTRALIA-7_16_1961-GA1")
data<-data[!(row.names(data) %in% row.names.remove), ]


mds<-data%>%
  #dist()%>%
  cmdscale(k=2) %>%
  as_tibble()
colnames(mds)<-c("Dim.1","Dim.2")

# plot MDS
ggscatter(mds, x = "Dim.1", y = "Dim.2", 
          #label = rownames(data),
          size = 3,
          repel = TRUE)
#################################################
## k-means clustering
## optimal number of k
## Elbow method
set.seed(123)
library(factoextra)
library(NbClust)
RSVA_G_wss<-fviz_nbclust(mds, kmeans, method = "wss")
#Average Silhouette Method
RSVA_G_sil<-fviz_nbclust(mds, kmeans, method = "silhouette")
#Gap Statistic Method
RSVA_G_gap<-fviz_nbclust(mds, kmeans, method = "gap_stat")

# K-means clustering
k<-kmeans(mds,centers=4,nstart = 123)
clust <- k$cluster %>%
  as.factor()
mds2 <- mds %>%
  mutate(groups = clust)
# Plot and color by groups
mds2$strain <-rownames(data)
taxa<-read.csv("EpiCC/RSVA_G/RSVA_G_ML/RSVA_G_epicc_taxa.csv")
mds_add<-left_join(mds2,taxa,by = c("strain" = "id"))

#install.packages("gghighlight")
library(gghighlight)
library(ggforce)
#install.packages("ggforce")


RSVA_G_mds<-ggplot(mds_add, aes(x=Dim.1, y=Dim.2,color=groups,shape=duplication,label=ifelse(strain %in% vac, label, ""))) + 
  geom_point(size=3)+
  xlab("Dimension 1")+
  ylab("Dimension 2")+
  scale_shape_manual(values = c(19,5))+
  scale_color_brewer(palette = "Set2",name="T-cell immunity clusters")+
  theme_light()+
  geom_encircle(aes(group=groups),expand=0,linetype=2)+
  geom_label_repel(fill = "white", xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),color="black",size=3)


## plot pairwise tree
## pair with a phylogenetic tree
## plot tree by side
library(treeio)
library(ggplot2)
library(ggtree)

library(ape)
library(RColorBrewer)
library("dplyr")
library("tibble")

t1 <- read.tree(file = "EpiCC/RSVA_G/RSVA_G_ML/RSVA_G_besttree_midpoint.nwk")
p1<-ggtree(t1) +
  geom_treescale(x=0, y=300,linesize= 1,,width=0.005)

library(tidyr)
mds2<-mds2 %>% 
  separate(strain, into = c("Accession", "Subtype", "year", "country", "date"), sep = "-", extra = "merge")
mds2$strain <-rownames(data)

taxa<-read.csv("EpiCC/RSVA_G/RSVA_G_ML/RSVA_G_epicc_taxa.csv")

taxa_group<-left_join(taxa,mds2,by = c("accession" = "Accession"))
taxa_group$groups<-as.factor(taxa_group$groups)
p2<-p1 %<+% taxa_group+ 
  geom_tippoint(aes(color=groups,shape=duplication), size=3, alpha=.75)+
  scale_shape_manual(values = c(19,5))+
  scale_color_brewer(palette = "Set2",name="T-cell immunity clusters")

p2
RSVA_G_ML<-p2
### clock signal of 
crossT_dist<-data%>%
  select(`JX198138-A-1961-AUSTRALIA-7_16_1961-GA1`)
crossT_dist$strain<-rownames(data)
crossT_dist<-crossT_dist %>% separate(strain, into = c("Accession", "Subtype", "year", "country", "date"), sep = "-", extra = "merge")
mds_add<-mds_add %>% 
  separate(strain, into = c("Accession", "Subtype", "year", "country", "date"), sep = "-", extra = "merge")
crossT_dist<-left_join(crossT_dist,mds_add,by=c("Accession"="Accession","Subtype"="Subtype",
                                             "year"="year","country"="country","date"="date"))

RSVA_G_immu<-ggplot(crossT_dist, aes(x=as.numeric(year), y=`JX198138-A-1961-AUSTRALIA-7_16_1961-GA1`,color=groups,shape=duplication)) + 
  geom_point(size=3)+
  scale_colour_brewer(palette = "Set2",name="T-cell immuno-clusters")+
  scale_shape_manual(values = c(19,5))+
  ylab("T cell immunity distance")+
  xlab("Isolated year")+
  theme_light()
######################################################################

## genetic hamming distance
ham_genetic<- read.table("EpiCC/RSVA_G/RSVA_G_epicc_subsample_genetic_hamming.csv", 
                         sep=",", header=TRUE,check.names=FALSE, stringsAsFacto=F, fill=TRUE,row.names = 1)

ham_dist<-ham_genetic%>%
  select(`JX198138-A-1961-Australia-7/16/1961-GA1`)


ham_dist$strain<-rownames(ham_dist)
ham_dist<-ham_dist %>% separate(strain, into = c("Accession", "Subtype", "year", "country", "date"), sep = "-", extra = "merge")
ham_dist<-left_join(ham_dist,mds_add,by=c("Accession"="Accession","Subtype"="Subtype",
                                       "year"="year"))
RSVA_G_ham<-ggplot(ham_dist, aes(x=as.numeric(year), y=`JX198138-A-1961-Australia-7/16/1961-GA1`,color=groups,shape=duplication)) + 
  geom_point(size=3)+
  scale_colour_brewer(palette = "Set2",name="T-cell immuno-clusters")+
  scale_shape_manual(values = c(19,5))+
  ylab("Genetic Hamming distance")+
  xlab("Isolated Year")+
  theme_light()

library(grid)
library(ggpubr)
ggarrange(RSVA_F_ML,RSVA_F_mds,RSVA_F_immu,RSVA_F_ham,
  RSVA_G_ML,RSVA_G_mds,RSVA_G_immu,RSVA_G_ham,
  labels = c("A", "B", "C","D","E","F","G","H"),
  ncol = 4, nrow = 2,
  common.legend = TRUE,
  widths = c(1.5,2,1.5,1.5))