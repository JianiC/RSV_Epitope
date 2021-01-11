## RSVB F epicc visulize
library(reshape2)
## read the raw dataframe
raw1<-read.table("EpiCC/RSV_EpiCC_Data/RSVB_F_classI_epicc_share_raw.csv", 
                 sep=",", header=TRUE,check.names=FALSE, stringsAsFacto=F, fill=TRUE,row.names = 1)
raw2<- read.table("EpiCC/RSV_EpiCC_Data/RSVB_F_classII_epicc_share_raw.csv", 
                  sep=",", header=TRUE,check.names=FALSE, stringsAsFacto=F, fill=TRUE,row.names = 1)

raw<-raw1+raw2

## norm_function
norm<-function(df,n,m){
  x<-df[n,m]
  a<-df[n,n]
  y<-x/a
  return(y)
}

all_norm<-function(input){
  my_vector <- vector(mode="numeric")
  for(n in 1:nrow(input)){
    for(m in 1:ncol(input)){
      x<-norm(input,n,m)
      my_vector[ncol(input)*(n-1)+m]=x
    }
  }
  
  #print(is.numeric(ncol(df)))
  ## convert vector to df
  my_matrix <- matrix(my_vector, ncol=ncol(input), byrow=TRUE)
  output <- as.data.frame(my_matrix, stringsAsFactors=FALSE)
  strain<-colnames(input)
  colnames(output)<-strain
  row.names(output)<-strain
  return(output)
}


norm_crossT<-all_norm(raw)
#####################################################################################################

## select the vaccine column
library(dplyr)
library(tidyr)
vac<-norm_crossT%>%
  select(`U63644-VACCINE_CPTS-248_404`,`AF013255-VACCINE_CP52`,`KT992094-VACCINE_D46_D53`)

vac$taxa<-rownames(vac)

write.csv(vac, "EpiCC/RSVB_F/RSVB_F_vac.csv")
## after code the year and region

vac_meta<-read.csv("EpiCC/RSVB_F/RSVB_F_vac_code.csv")

vac_sum<-vac_meta%>% 
  group_by(year_code,region) %>%
  dplyr::summarise(
                   CP248=mean(`U63644.VACCINE_CPTS.248_404`),
                   CP52=mean(`AF013255.VACCINE_CP52`),)

vac_sum<-vac_sum %>% gather(vaccine,value,CP248:CP52)
vac_sum$meta<-paste(vac_sum$year_code,"\n",vac_sum$region)

RSVB_F<-vac_sum


ggplot(data = vac_sum,aes( y = value,x = interaction(region,year_code),
                           group=vaccine,colour = vaccine)) +   geom_point(size=2) + 
  geom_polygon(size = 1, alpha= 0.2)+coord_polar()+
  scale_x_discrete() +
  ylim(0,1.0)+
  theme_light()+
  xlab("WHO regoin ~ Isolation Year")+
  ylab("Normalized Cross T-cell Immunity")+
  scale_colour_brewer(palette = "Set2",name="Vaccine candidates")


###################################################################  
## with cross T-cell immunity disimilarity build mds
## remove vaccine column and rows  
norm_crossT_clean <- select(norm_crossT, -contains("VACCINE"))
norm_crossT_clean <- select(norm_crossT_clean, -contains("RSVA_PDA")) 

row.names.remove <- c("KT992094-VACCINE_D46_D53", "U63644-VACCINE_CPTS-248_404","AF013255-VACCINE_CP52","AF035006-VACCINE_RA2CP","RSVA_PDA_ANCESTRAL")
norm_crossT_clean<-norm_crossT_clean[!(row.names(norm_crossT_clean) %in% row.names.remove), ]




crossT_dis<-1-norm_crossT_clean 
mds<-crossT_dis%>%
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

fviz_nbclust(mds, kmeans, method = "wss")
#Average Silhouette Method
fviz_nbclust(mds, kmeans, method = "silhouette")
#Gap Statistic Method
fviz_nbclust(mds, kmeans, method = "gap_stat")

# K-means clustering
k<-kmeans(mds,centers=2,nstart = 123)
clust <- k$cluster %>%
  as.factor()
mds2 <- mds %>%
  mutate(groups = clust)
# Plot and color by groups
RSVB_F_mds<-ggplot(mds2, aes(x=Dim.1, y=Dim.2, color=groups)) + 
  geom_point()+
  scale_color_brewer(palette = "Set2",name="T-cell immunity clusters")+
  theme_light()

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

t1 <- read.tree(file = "EpiCC/RSVB_F/RSVB_F_subsample_epicc.nwk")
p1<-ggtree(t1) +
  geom_treescale(x=0.03, y=3,linesize= 1,,width=0.005)

mds2$strain <-rownames(crossT_dis)
mds2<-mds2 %>% separate(strain, into = c("Accession", "Subtype", "year", "country", "date"), sep = "-", extra = "merge")

taxa<-read.csv("EpiCC/RSVB_F/RSVB_F_epicc_subsample_taxa.csv")

taxa_group<-left_join(taxa,mds2,by = c("accession" = "Accession"))
taxa_group$groups<-as.factor(taxa_group$groups)
p2<-p1 %<+% taxa_group+ 
  geom_tippoint(aes(color=groups), size=2, alpha=.75)+
  scale_color_brewer(palette = "Set2",name="T-cell immunity clusters")

RSVB_F_ML<-p2

### clock signal of 
crossT_dist<-crossT_dis%>%
  select(`RSVB_PDA_ANCESTRAL`)
crossT_dist$strain<-rownames(crossT_dist)
crossT_dist<-crossT_dist %>% separate(strain, into = c("Accession", "Subtype", "year", "country", "date"), sep = "-", extra = "merge")
crossT_dist<-left_join(crossT_dist,mds2,by=c("Accession"="Accession","Subtype"="Subtype",
                                             "year"="year","country"="country","date"="date"))

RSVB_F_immu<-ggplot(crossT_dist, aes(x=as.numeric(year), y=`RSVB_PDA_ANCESTRAL`,color=groups)) + 
  geom_point()+
  scale_colour_brewer(palette = "Set2",name="T-cell immuno-clusters")+
  ylab("T cell immunity distance")+
  xlab("Isolated Year")+
  theme_light()
######################################################################

## genetic hamming distance
ham_genetic<- read.table("EpiCC/RSVB_F/RSVB_F_epicc_subsample_genetic_hamming.csv", 
                         sep=",", header=TRUE,check.names=FALSE, stringsAsFacto=F, fill=TRUE,row.names = 1)

ham_dist<-ham_genetic%>%
  select(`RSVB_pda_ancestral`)


ham_dist$strain<-rownames(ham_dist)
ham_dist<-ham_dist %>% separate(strain, into = c("Accession", "Subtype", "year", "country", "date"), sep = "-", extra = "merge")
ham_dist<-left_join(ham_dist,mds2,by=c("Accession"="Accession","Subtype"="Subtype",
                                       "year"="year"))
RSVB_F_ham<-ggplot(ham_dist, aes(x=as.numeric(year), y=`RSVB_pda_ancestral`,color=groups)) + 
  geom_point()+
  scale_colour_brewer(palette = "Set2",name="T-cell immuno-clusters")+
  ylab("Genetic Hamming distance")+
  xlab("Isolated Year")+
  theme_light()

##########################################
##################################################################

## sum RSVA and RSVB plot in the same page
ggarrange(RSVA_F_mds,RSVA_F_ML,RSVA_F_immu,RSVA_F_ham,
          RSVB_F_mds,RSVB_F_ML,RSVB_F_immu,RSVB_F_ham,
          labels = c("A", "B", "C","D","E","F","G","H"),
          ncol = 4, nrow = 2,
          common.legend = TRUE)
