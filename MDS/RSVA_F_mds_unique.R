## RSVA G epicc visulize
library(reshape2)
library(ggplot2)
## read the raw dataframe
raw1<-read.table("EpiCC/RSV_EpiCC_Data/RSVA_F_classI_epicc_unique.csv", 
                 sep=",", header=TRUE,check.names=FALSE, stringsAsFacto=F, fill=TRUE,row.names = 1)
raw2<- read.table("EpiCC/RSV_EpiCC_Data/RSVA_F_classII_epicc_unique.csv", 
                  sep=",", header=TRUE,check.names=FALSE, stringsAsFacto=F, fill=TRUE,row.names = 1)

raw<-raw1+raw2




###################################################################  
## with cross T-cell immunity disimilarity build mds
library("ggpubr")
library(dplyr)
## remove vaccine column and rows  
#data <- select(raw, -contains("VACCINE"))
data <- select(raw, -row.names.remove) 
vac<-c("U63644-VACCINE_CPTS-248_404","AF013255-VACCINE_CP52","RSVA_PDA_ANCESTRAL")
row.names.remove <- c("KT992094-VACCINE_D46_D53", "RSVB_PDA_ANCESTRAL","AF035006-VACCINE_RA2CP")
data<-data[!(row.names(data) %in% row.names.remove), ]


## create a function for GOF
r <- cmdscale(data, eig=TRUE)
plot(cumsum(r$eig) / sum(r$eig), 
     type="h", lwd=5, las=1, 
     xlab="Number of dimensions", 
     ylab=expression(R^2))

plot(r$eig, 
     type="h", lwd=5, las=1, 
     xlab="Number of dimensions", 
     ylab="Eigenvalues")


M1<-list()
M2<-list()
for(i in 1:3){
  #GOF<-list()
  GOF<-cmdscale(crossT_dis,i,eig=T)$GOF
  M1 <- append(M1,GOF[1])
  M2<-append(M2,GOF[2])
}

library (plyr)
m1 <- ldply (M1, data.frame)
m2<-ldply (M2, data.frame)

m<-cbind(m1,m2)
names(m)<-c("M1","M2")
m$dimensions <- seq.int(nrow(m))

GOF<-function(eurodist,subtype){
  M1<-list()
  M2<-list()
  for(i in 1:10){
    
    GOF<-cmdscale(eurodist,i,eig=T)$GOF
    M1 <- append(M1,GOF[1])
    M2<-append(M2,GOF[2])
  }
  m1 <- ldply (M1, data.frame)
  m2<-ldply (M2, data.frame)
  m<-cbind(m1,m2)
  names(m)<-c("M1","M2")
  m$dimensions <- seq.int(nrow(m))
  m$subtype<-subtype
  return(m)
  
}


RSVA_F_GOF<-GOF(data,"RSV-A F protein")

RSV_GOF<-rbind(RSVA_F_GOF,RSVA_G_GOF,RSVB_F_GOF,RSVB_G_GOF)


RSV_GOF2 <- melt(RSV_GOF, id.vars = c("subtype", "dimensions"),
                 variable.name = "method", 
                 value.name = "value")

ggplot(RSV_GOF2, aes(x=dimensions, y=value, color=method)) +
  geom_line()+
  geom_point()+
  scale_x_continuous(breaks=c(0,2,4,6,8,10))+
  scale_color_brewer(palette="Set2")+
  theme_light()+
  facet_wrap(.~subtype,scales="free",nrow=2)+
  ylab("Goodness-of-fit")+
  xlab("k dimensions")

###########################################################################


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
RSVA_F_wss<-fviz_nbclust(mds, kmeans, method = "wss")
#Average Silhouette Method
RSVA_F_sil<-fviz_nbclust(mds, kmeans, method = "silhouette")
#Gap Statistic Method
RSVA_F_nb<-fviz_nbclust(mds, kmeans, method = "gap_stat")
#install.packages("ggalt")
library(ggalt)
# K-means clustering
k<-kmeans(mds,centers=3,nstart = 123)
clust <- k$cluster %>%
  as.factor()
mds2 <- mds %>%
  mutate(groups = clust)
# Plot and color by groups
#install.packages("ggrepel")
#devtools::install_github("slowkow/ggrepel")
library(ggrepel)
mds2$strain <-rownames(data)
RSVA_F_mds<-ggplot(mds2, aes(x=Dim.1, y=Dim.2,color=groups,label=ifelse(strain %in% vac, strain, ""))) + 
  geom_point()+
  xlab("Dimension 1")+
  ylab("Dimension 2")+
  scale_color_brewer(palette = "Set2",name="T-cell immunity clusters")+
  theme_light()+
  geom_encircle(expand=0,linetype=2)+
  geom_label_repel(fill = "white", xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),color="black",size=2)
 
  #stat_ellipse(type = "norm", linetype = 2)

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

t1 <- read.tree(file = "EpiCC/RSVA_F/RSVA_F_ML/RSVA_F_besttree_midpoint.nwk")
p1<-ggtree(t1) +
  geom_treescale(x=0.03, y=3,linesize= 1,,width=0.005)

mds2$strain <-rownames(raw)
mds2<-mds2 %>% separate(strain, into = c("Accession", "Subtype", "year", "country", "date"), sep = "-", extra = "merge")

taxa<-read.csv("EpiCC/RSVA_F/RSVA_F_epicc_taxa.csv")

taxa_group<-left_join(taxa,mds2,by = c("accession" = "Accession"))
taxa_group$groups<-as.factor(taxa_group$groups)
p2<-p1 %<+% taxa_group+ 
  geom_tippoint(aes(color=groups), size=2, alpha=.75)+
  scale_color_brewer(palette = "Set2",name="T-cell immunity clusters")

RSVA_F_ML<-p2

### clock signal of 
crossT_dist<-data%>%
  select(`JX198138-A-1961-AUSTRALIA-7_16_1961-GA1`)
crossT_dist$strain<-rownames(crossT_dist)
crossT_dist<-crossT_dist %>% separate(strain, into = c("Accession", "Subtype", "year", "country", "date"), sep = "-", extra = "merge")
crossT_dist<-left_join(crossT_dist,mds2,by=c("Accession"="Accession","Subtype"="Subtype",
                                             "year"="year","country"="country","date"="date"))

RSVA_F_immu<-ggplot(crossT_dist, aes(x=as.numeric(year), y=`JX198138-A-1961-AUSTRALIA-7_16_1961-GA1`,color=groups)) + 
  geom_point()+
  scale_colour_brewer(palette = "Set2",name="T-cell immuno-clusters")+
  ylab("T cell immunity distance")+
  xlab("Isolated Year")+
  theme_light()
######################################################################

## genetic hamming distance
ham_genetic<- read.table("EpiCC/RSVA_F/RSVA_F_epicc_subsample_genetic_hamming.csv", 
                         sep=",", header=TRUE,check.names=FALSE, stringsAsFacto=F, fill=TRUE,row.names = 1)

ham_dist<-ham_genetic%>%
  select(`JX198138-A-1961-Australia-7/16/1961-GA1`)


ham_dist$strain<-rownames(ham_dist)
ham_dist<-ham_dist %>% separate(strain, into = c("Accession", "Subtype", "year", "country", "date"), sep = "-", extra = "merge")
ham_dist<-left_join(ham_dist,mds2,by=c("Accession"="Accession","Subtype"="Subtype",
                                       "year"="year"))
RSVA_F_ham<-ggplot(ham_dist, aes(x=as.numeric(year), y=`JX198138-A-1961-Australia-7/16/1961-GA1`,color=groups)) + 
  geom_point()+
  scale_colour_brewer(palette = "Set2",name="T-cell immuno-clusters")+
  ylab("Genetic Hamming distance")+
  xlab("Isolated Year")+
  theme_light()


##################################################################
## combine mds, ML tree, antigenetic and ham genetic change


