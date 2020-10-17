# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name
## read the dataframe

data<-read.table("Test_data/smithetal_normshared_Oct15.csv",sep=",",row.names = 1,header=TRUE)

library(magrittr)
library(dplyr)
library(ggpubr)
library(purrr)

# compute MDS ## GOF determine diminsion?
mds<-data%>%
  dist() %>%
  cmdscale(k=3) %>%
  as_tibble()

colnames(mds)<-c("Dim.1","Dim.2","Dim.3")

# plot MDS
ggscatter(mds, x = "Dim.1", y = "Dim.2", 
          #label = rownames(data),
          size = 3,
          repel = TRUE)



############################################################
## find the ooptimal number of cluster
set.seed(123)

# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(mds, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

## other approaches
#install.packages("NbClust")
#install.packages("factoextra")
library(factoextra)
library(NbClust)
# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
fviz_nbclust(mds, kmeans, nstart = 25,  method = "gap_stat", nboot = 50,k.max = 15)+
  labs(subtitle = "Gap statistic method")

# Silhouette method
fviz_nbclust(mds, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# Elbow method
fviz_nbclust(mds, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")


NbClust(data = mds, diss = NULL, distance = "euclidean",min.nc = 2, max.nc = 15, method = kmeans)
###########################################################


# K-means clustering
clust <- kmeans(mds,10)$cluster %>%
  as.factor()
mds <- mds %>%
  mutate(groups = clust)
# Plot and color by groups
ggplot(mds, aes(x=Dim.1, y=Dim.2, color=groups)) + 
  geom_point()+
  scale_color_brewer(palette = "Set3",name="T-cell immuno-clusters") 



######################################
## plot tree by side
library(treeio)
library(ggplot2)
library(ggtree)
library(phytools)
library(ape)
library(RColorBrewer)
library("dplyr")
library("tibble")
mds$strain<-rownames(data)
t1 <- read.tree(file = "Test_data/Smith_tree/best_tree_smith_2004_HA1.nwk")
p1<-ggtree(t1) 
p1

#write.csv(mds,"Test_data/smith_HA1_1015.csv",row.names = FALSE)
cluster<-mds %>% select(strain,groups)
#cluster<-column_to_rownames(cluster, 'strain')
colnames(cluster)<-c("taxa","groups")
#cluster <- read.table("Test_data/mds_group_test.csv", sep=",", header=TRUE,check.names=FALSE, stringsAsFactor=F, fill=TRUE)
cluster$taxa<- gsub('_', '/', cluster$taxa)
cluster$groups<-as.factor(cluster$groups)
p2<-p1 %<+% cluster+ 
  geom_tippoint(aes(color=groups), size=2, alpha=.75)+
  scale_color_brewer(palette = "Set3",name="T-cell immuno-clusters") 
p2
##############################################
#calculate the distance between two points
# construct the mds table
library('tidyverse')
library('tidyr')
library(dplyr)

#calculate distance, eg(distance to BI/16190/68)
# x and y coord for reference point
X0=-0.18335288
Y0=-0.0782464947
Z0=	0.011120043
mds$dis<-sqrt((mds$Dim.1-X0)^2+(mds$Dim.2-Y0)^2+(mds$Dim.3-Z0)^2)
# seperate the season from strain name
mds<-mds%>%
  separate(strain,into=c("antigentic_group","extra","season"))

# a function to translate year 
foo <- function(x, year=1968){
  
  y <- ifelse(x >= year %% 100, 1900+x, 2000+x)
  return(y)
}

mds$test<-lapply(as.numeric(mds$season),foo)
library(ggplot2)
ggplot(mds, aes(x=as.numeric(test), y=dis,color=groups)) + 
  geom_point()+
  scale_colour_brewer(palette = "Set3",name="T-cell immuno-clusters")
####################################################################3
#use the relatness score
mds$relate<-data$BI_16190_68

ggplot(mds, aes(x=as.numeric(test), y=relate,color=groups)) + 
  geom_point()+
  scale_colour_brewer(palette = "Set3",name="T-cell immuno-clusters")
###################################################################
d <- dist(data) # euclidean distances between the rows
fit <- cmdscale (d,eig=TRUE, k=2) # k is the number of dim
fit # view results
head(fit)
# plot solution
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS", type="n")
text(x, y, labels = row.names(data), cex=.7)

## find the ooptimal number of cluster
et.seed(123)

# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(fit$points, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")



## cluster with k-means and plot
kmeans_clust <- kmeans(fit$points, 2)  # k-means wihth 32clusters.
plot(fit$points, type = "n", main="MDS with sammon() and clustered", xlab = "X-Dim", ylab="Y-Dim", pch=19)
text(fit$points, labels = rownames(data), col = kmeans_clust$cluster)
print(kmeans_clust)