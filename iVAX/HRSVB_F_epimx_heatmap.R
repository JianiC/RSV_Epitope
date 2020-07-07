library('tidyverse')
library('tidyr')
library('reshape2')
library('ggplot2')
library('varhandle')
library('RColorBrewer')
library('gridExtra')

#loop through directory and read files in r
filenames <- list.files(path = "./iVAX/HRSVB_F_classII_csv", pattern = "MergedFinalData_dir00*")
numfiles <- length(filenames) 

#### test ## change file name?
for (i in 1:numfiles){  
  name <- gsub(".csv","",filenames[i])
  name <- gsub("MergedFinalData","HRSVB_F_classII",name) 
  assign(name, read.csv(paste("./iVAX/HRSVB_F_classII_csv/", filenames[i], sep = ""), 
                        header = FALSE, stringsAsFactors = FALSE, sep = ","))
}

#to create name list and store dataframe name
dataname <- NULL
for (i in 1:numfiles){
  name[i] <- gsub(".csv","",filenames[i])
  name[i] <- gsub("MergedFinalData","HRSVB_F_classII",name[i]) 
  dataname <- rbind(dataname, name[i])
}

##combine mulitple dataframe
HRSVB_F_classII <-rbind(HRSVB_F_classII_dir001,HRSVB_F_classII_dir002,HRSVB_F_classII_dir003,HRSVB_F_classII_dir004,HRSVB_F_classII_dir005,HRSVB_F_classII_dir006,HRSVB_F_classII_dir007,HRSVB_F_classII_dir008)

********************************************************************
  #preclean-data function
  #to extract each sequence EMX (EpiMatrix) info
  #data reshaping

totalseq <- grep("File: HRSVB_F_CLEAN0627_PROTEIN", HRSVB_F_classII$V1, useBytes = TRUE)
endseq <- grep("Summarized Results", HRSVB_F_classII$V1, useBytes = TRUE)

HRSVB_F.epimx.data <- ""

for (i in 1:length(totalseq)) {
  epimx.data <- HRSVB_F_classII%>% slice((totalseq[i]):(endseq[i]-2))
  seqname <- gsub("File: HRSVB_F_CLEAN0627_PROTEIN - Sequence: ", "", epimx.data$V1[1])
  epimx.data <- epimx.data %>% slice(grep("Frame", V1, useBytes = T):nrow(epimx.data))
  colnames(epimx.data) <- c("Frame.Start", "AA.Sequence", "Frame.Stop", "DRB1*0101", "DRB1*0301", "DRB1*0401",  
                            "DRB1*0701", "DRB1*0801", "DRB1*0901", "DRB1*1101", "DRB1*1301", "DRB1*1501",  
                            "Hits", "NA")
  epimx.data <- epimx.data[-c(1:2),]
  epimx.data.reshape <- epimx.data %>% mutate(Seq.name = seqname) %>% 
    select(Seq.name, Frame.Start:Hits) %>%
    gather(HLA.Alleles, Z.score, "DRB1*0101":"DRB1*1501", factor_key = TRUE)
  HRSVB_F.epimx.data <- rbind(HRSVB_F.epimx.data, epimx.data.reshape)
}

#MORE data restructuring 
*************************************************************************************
HRSVB_F.epimx.data.2 <- HRSVB_F.epimx.data %>% filter(!Seq.name == "")
## remove duplicate from HRSVB_F.epimx.data.2

HRSVB_F.epimx.data.2 <-distinct(HRSVB_F.epimx.data.2)

##save the original epimx dataframe
write.csv(HRSVB_F.epimx.data.2,"iVAX/HRSVB_F_classII_epimx.csv")

#checkpoint
x <- data.frame(unique(HRSVB_F.epimx.data.2$Seq.name)) ## get all sequence name 
**************************************************************************************************
  #filter Z SCORE >= 1.64
HRSVB_F.epimx.data.3 <- HRSVB_F.epimx.data.2 %>% filter(Z.score >= 1.64)
write.csv(HRSVB_F.epimx.data.3,"iVAX/HRSVB_F_classII_epimx_filter.csv")
HRSVB_F.epimx.data.3$Z.score <- as.numeric(HRSVB_F.epimx.data.3$Z.score)


*************************************************************************************************
## prepare for heatmap
## get the count of record for each strain, each AA frame 
## try with grouped data

  
##  
HRSVB_F_heatmap1<-  HRSVB_F.epimx.data.3 %>% select(Seq.name, Frame.Start,HLA.Alleles)
HRSVB_F_heatmap2<-HRSVB_F_heatmap1 %>%
  group_by(Seq.name, Frame.Start) %>%
  summarize(count = n())

## reshape from value to column name

HRSVB_F_heatmap3<-dcast(HRSVB_F_heatmap2, Seq.name ~ Frame.Start) 
## sort by Frame start value
list<-str_sort(c(colnames(HRSVB_F_heatmap3)), numeric = TRUE)
list
HRSVB_F_heatmap4<-HRSVB_F_heatmap3[list]
write.csv(HRSVB_F_heatmap4,"iVAX/HRSVB_F_classII_epimx_hp.csv")

## create columns for empty start frame
col<-c(1:564)
unlist(col)
##covert to string
for (i in 1:564){
  col[i]<-toString(col[i])
}


for(i in 1:564){
  if (col[i] %in% colnames(HRSVB_F_heatmap4)) {
    x=col[i]
    message(sprintf(" %s is exist", x))
  } else {
    HRSVB_F_heatmap4[,col[i]]<-NA
  }
}

## replace 'NA' with 0
HRSVB_F_heatmap4[is.na(HRSVB_F_heatmap4)] <- 0
## sort again to generate a new datafrane

list2<-str_sort(c(colnames(HRSVB_F_heatmap3)), numeric = TRUE)
list



