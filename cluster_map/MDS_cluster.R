# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name
## read the dataframe

data<-read.table("Test_data/epicc_mds_test.csv",sep=",",row.names = 1,header=TRUE)
library(magrittr)
library(dplyr)
library(ggpubr)

# compute MDS
mds<-data%>%
  dist() %>%
  cmdscale() %>%
  as_tibble()

colnames(mds)<-c("Dim.1","Dim.2")

# plot MDS
ggscatter(mds, x = "Dim.1", y = "Dim.2", 
          #label = rownames(data),
          size = 3,
          repel = TRUE)



############################################################
## find the ooptimal number of cluster
et.seed(123)

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



###########################################################


# K-means clustering
clust <- kmeans(mds, 2)$cluster %>%
  as.factor()
mds <- mds %>%
  mutate(groups = clust)
# Plot and color by groups
ggscatter(mds, x = "Dim.1", y = "Dim.2", 
          #label = rownames(data),
          color = "groups",
          palette = "jco",
          size = 3, 
          ellipse = TRUE,
          ellipse.type = "convex",
          repel = TRUE)


##############################################
#calculate the distance between two points
# construct the mds table

mds$strain<-rownames(data)

#calculate distance, eg(north vaccine1 to north vaccine 2020)

dist()









d <- dist(data) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
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