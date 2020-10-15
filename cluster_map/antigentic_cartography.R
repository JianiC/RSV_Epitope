## installation
# devtools::install_github("YuLab-SMU/ggacmap")
#library(ggacmap)
#remotes::install_url("https://github.com/acorg/acmacs.r/releases/download/v4.6/acmacs.r_4.6_R_macOS-10.14.tgz", build = FALSE)
devtools::install_github("acorg/Racmacs")
library(Racmacs)
map_file <- system.file("extdata/h3map2004.ace", package = "Racmacs")
map <- read.acmap(map_file)
view(map)

path_to_titer_file <- system.file("extdata/h3map2004_hitable.csv", package = "Racmacs")

test<-read.table("iVAX/RSVB_F_analysis_old/test2.csv",sep=",",  header=TRUE,check.names=FALSE, stringsAsFactor=F, fill=TRUE,row.names = 1)

titer_table <- read.titerTable("iVAX/RSVB_F_analysis_old/test2.csv")

test<-read.table("iVAX/RSVB_F_analysis_old/test4.txt",sep="\t",  header=TRUE,check.names=FALSE, stringsAsFactor=F, fill=TRUE,row.names = 1)
test<-read.table("Test_data/smith-2004-data.txt",sep="\t",  header=TRUE,check.names=FALSE, stringsAsFactor=F, fill=TRUE,row.names = 1)

test<-read.table("Test_data/epicc_test.txt",sep="\t",  header=TRUE,check.names=FALSE, stringsAsFactor=F, fill=TRUE,row.names = 1)
map <- acmap(
  table = test
)

map <- optimizeMap(
  map                     = map,
  number_of_dimensions    = 2,
  number_of_optimizations = 10,
  minimum_column_basis    = "none"
)

view(map)


antigentic_distance<-mapDistances(map, optimization_number = NULL)
