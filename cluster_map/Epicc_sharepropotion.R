library('tidyverse')
library('tidyr')
library('reshape2')
library('ggplot2')
library('varhandle')
library('RColorBrewer')
library('gridExtra')
library(dplyr)

## read the EPICC shared raw
df<- read.table("Test_data/Smith_raw_shared.csv",  sep=",", header=TRUE,check.names=FALSE, stringsAsFactor=F, fill=TRUE,row.names = 1)


## test with the standrize function, test with column =2 , row =3
#x<-df[2,3]
#standarize with the Epitope sum score for each strain the column is a , row is b
#a<-df[2,2]
#b<-df[3,3]
# new score is y
#y<- x/((a+b)/2)

## wrap up into a function
epicc_stand<-function(df,n,m){
  x<-df[n,m]
  a<-df[n,n]
  b<-df[m,m]
  y<-x/((a+b)/2)
  return(y)
}
my_vector <- vector(mode="numeric")

#x<-epicc_stand(df,2,3)
#my_vector[1]<-x
# store the result to a new dataframe into a vector
## loop through the whole dataframe

for(n in 1:300){
  for(m in 1:300){
    x<-epicc_stand(df,n,m)
    my_vector[300*(n-1)+m]=x
  }
  
}

## covert vector to a dataframe
m <- matrix(my_vector, ncol=300, byrow=TRUE)
df2 <- as.data.frame(m, stringsAsFactors=FALSE)
strain<-colnames(df)
colnames(df2)<-strain
row.names(df2)<-strain

write.csv(df2, "smith_epicc_raw_shared_standarized1020.csv")
