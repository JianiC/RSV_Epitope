library('tidyverse')
library('tidyr')
library('reshape2')
library('ggplot2')
library('varhandle')
library('RColorBrewer')
library('gridExtra')
library(dplyr)
#############################################################################
#loop through directory and read files in r
filenames <- list.files(path = "./iVAX/HRSVB_F_classII_csv", pattern = "MergedFinalData_dir00*")
numfiles <- length(filenames) 

## change file name
for (i in 1:numfiles){  
  name <- gsub(".csv","",filenames[i])
  name <- gsub("MergedFinalData","HRSVB_F_classII",name) 
  assign(name, read.csv(paste("./iVAX/HRSVB_F_classII_csv/", filenames[i], sep = ""), 
                        header = FALSE, stringsAsFactors = FALSE, sep = ","))
}



#########################################################################
#preclean-data function
#to extract each sequence EMX (EpiMatrix) info
#data reshaping

for (i in 1:numfiles){
  HRSVB_F_classII_dir <- eval(parse(text=paste("HRSVB_F_classII_dir", i, sep="")))
  HRSVB_F.epimx.data <- rbind(HRSVB_F.epimx.data,next_df)
}













totalseq <- grep("File: HRSVB_F_CLEAN0627_PROTEIN", HRSVB_F_classII_dir003$V1, useBytes = TRUE)
endseq <- grep("Summarized Results", HRSVB_F_classII_dir003$V1, useBytes = TRUE)

HRSVB_F.epimx.data003 <- ""

for (i in 1:length(totalseq)) {
  epimx.data <- HRSVB_F_classII_dir003%>% slice((totalseq[i]):(endseq[i]-2))
  seqname <- gsub("File: HRSVB_F_CLEAN0627_PROTEIN - Sequence: ", "", epimx.data$V1[1])
  epimx.data <- epimx.data %>% slice(grep("Frame", V1, useBytes = T):nrow(epimx.data))
  colnames(epimx.data) <- c("Frame.Start", "AA.Sequence", "Frame.Stop", "DRB1*0101", "DRB1*0301", "DRB1*0401",  
                            "DRB1*0701", "DRB1*0801", "DRB1*0901", "DRB1*1101", "DRB1*1301", "DRB1*1501",  
                            "Hits", "NA")
  epimx.data <- epimx.data[-c(1:2),]
  epimx.data.reshape <- epimx.data %>% mutate(Seq.name = seqname) %>% 
    select(Seq.name, Frame.Start:Hits) %>%
    gather(HLA.Alleles, Z.score, "DRB1*0101":"DRB1*1501", factor_key = TRUE)
  HRSVB_F.epimx.data003 <- rbind(HRSVB_F.epimx.data003, epimx.data.reshape)
}


######################################################################################
## combine 8 dataframe 

HRSVB_F.epimx.data <- ""
for (i in 1:numfiles){
  next_df <- eval(parse(text=paste("HRSVB_F.epimx.data00", i, sep="")))
  HRSVB_F.epimx.data <- rbind(HRSVB_F.epimx.data,next_df)
}

######################################################################################
#MORE data restructuring 

HRSVB_F.epimx.data.2 <- HRSVB_F.epimx.data %>% filter(!Seq.name == "")
## remove duplicate from HRSVB_F.epimx.data.2

HRSVB_F.epimx.data.2 <-distinct(HRSVB_F.epimx.data.2)

##save the original epimx dataframe
write.csv(HRSVB_F.epimx.data.2,"iVAX/HRSVB_F_classII_epimx0812.csv")  
#checkpoint
x <- data.frame(unique(HRSVB_F.epimx.data.2$Seq.name)) ## get all sequence name   
HRSVB_F.epimx.data.2<-read.csv("iVAX/HRSVB_F_classII_epimx0812.csv")
#############################################################################3
## add missing sequence 

HRSVB_F_classII_missing <-read_csv("iVAX/HRSVB_F_classII_csv/HRSVB_F_classII_missing.csv")

totalseq <- grep("File: HRSVB_F_CLASSII_MISSING", HRSVB_F_classII_missing$`Class II EpiMatrix Report`, useBytes = TRUE)
endseq <- grep("Summarized Results", HRSVB_F_classII_missing$`Class II EpiMatrix Report`, useBytes = TRUE)

HRSVB_F.epimx.data_missing <- ""

for (i in 1:length(totalseq)) {
  epimx.data <- HRSVB_F_classII_missing%>% slice((totalseq[i]):(endseq[i]-2))
  seqname <- gsub("File: HRSVB_F_CLEAN0627_PROTEIN - Sequence: ", "", epimx.data$`Class II EpiMatrix Report`[1])
  epimx.data <- epimx.data %>% slice(grep("Frame", `Class II EpiMatrix Report`, useBytes = T):nrow(epimx.data))
  colnames(epimx.data) <- c("Frame.Start", "AA.Sequence", "Frame.Stop", "DRB1*0101", "DRB1*0301", "DRB1*0401",  
                            "DRB1*0701", "DRB1*0801", "DRB1*0901", "DRB1*1101", "DRB1*1301", "DRB1*1501",  
                            "Hits", "NA")
  epimx.data <- epimx.data[-c(1:2),]
  epimx.data.reshape <- epimx.data %>% mutate(Seq.name = seqname) %>% 
    select(Seq.name, Frame.Start:Hits) %>%
    gather(HLA.Alleles, Z.score, "DRB1*0101":"DRB1*1501", factor_key = TRUE)
  HRSVB_F.epimx.data_missing <- rbind(HRSVB_F.epimx.data_missing, epimx.data.reshape)
}
#MORE data restructuring 

HRSVB_F.epimx.data_missing <- HRSVB_F.epimx.data_missing %>% filter(!Seq.name == "")
## remove duplicate from HRSVB_F.epimx.data.2

HRSVB_F.epimx.data_missing <-distinct(HRSVB_F.epimx.data_missing)

HRSVB_F.epimx.data_missing<-HRSVB_F.epimx.data_missing%>%separate(Seq.name,into=c("File","Seq.name"),sep=" Sequence: ")
HRSVB_F.epimx.data_missing<-HRSVB_F.epimx.data_missing[c(2:8)]
HRSVB_F.epimx.data.2<-HRSVB_F.epimx.data.2[c(2:8)]


HRSVB_F.epimx.data.2_complete<-rbind(HRSVB_F.epimx.data.2,HRSVB_F.epimx.data_missing)


##save the original epimx dataframe
write.csv(HRSVB_F.epimx.data.2_complete,"iVAX/HRSVB_F_classII_epimx0815.csv")  
#checkpoint
x <- data.frame(unique(HRSVB_F.epimx.data.2_complete$Seq.name)) ## get all sequence name   

###########################################################################################
## filter with a clean sequence dataset
 
HRSVB_F.epimx.data.2_complete2<-HRSVB_F.epimx.data.2_complete %>% 
  separate(Seq.name, into = c("Accession", "Subtype", "Country", "date", "genotype"), sep = "-", extra = "merge")
seq <- scan("./RAxML/RAxML0708/RAxML_0708filter/HRSVB_F/HRSVB_F_taxa.txt", what="", sep="\n")

HRSVB_F.epimx.data.2_complete3 <- HRSVB_F.epimx.data.2_complete2 %>% 
  mutate(Type = ifelse(HRSVB_F.epimx.data.2_complete2$Accession %in% seq, "T", "F"))

HRSVB_F.epimx.data.2_complete4<-HRSVB_F.epimx.data.2_complete3%>% filter(HRSVB_F.epimx.data.2_complete3$Type =="T")
HRSVB_F.epimx.data.clean<-HRSVB_F.epimx.data.2_complete4[c(1:11)]
write_csv(HRSVB_F.epimx.data.clean,"iVAX/HRSVB_F_classII_epimx0815_clean.csv")
#####################################################################################  
#filter Z SCORE >= 1.64
HRSVB_F.epimx.data.3 <- HRSVB_F.epimx.data.clean %>% filter(Z.score >= 1.64)
#write.csv(HRSVB_F.epimx.data.3,"iVAX/HRSVB_F_classII_epimx0810_filter.csv")
#HRSVB_F.epimx.data.3$Z.score <- as.numeric(HRSVB_F.epimx.data.3$Z.score)

###################################################################################
#filter EPIBARS Hits >=4
HRSVB_F.epimx.data.3$Hits<- as.numeric(HRSVB_F.epimx.data.3$Hits)
HRSVB_F.epimx.data.4<- HRSVB_F.epimx.data.3 %>% filter(Hits >= 4)
HRSVB_F.epimx.data.4<- HRSVB_F.epimx.data.4[!(HRSVB_F.epimx.data.4$`AA.Sequence` == ""), ]
write.csv(HRSVB_F.epimx.data.4,"iVAX/HRSVB_F_classII_epimx0815_filter.csv")

#########################################################################################
## identify all of the potential epitope banner from the matix

## for each sequence, remove the duplicate from HLA allel 

HRSVB_F.epimx.data.5<-HRSVB_F.epimx.data.4 %>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)
## check
#x<-HRSVB_F.epimx.data.5%>% filter(Frame.Start==12)
#write_csv(x,"check2.csv")
### ignore the sequence information and get the counts of for each epitope banner

HRSVB_F.epimx.data.7<-HRSVB_F.epimx.data.5%>% 
  group_by(Frame.Start, AA.Sequence,Frame.Stop) %>%
  dplyr::summarise(count = n())

## calculate the coverage percentage
n<-n_distinct(HRSVB_F.epimx.data.clean$Accession)
HRSVB_F.epimx.data.7<- HRSVB_F.epimx.data.7%>% mutate(coverge = count/as.numeric(n))
## filter the interest epitope
epitope<-HRSVB_F.epimx.data.7%>%filter(coverge>=0.05)
epitope<-epitope[order(as.numeric(epitope$Frame.Start)),]
write_csv(epitope,"iVAX/HRSVB_epitope_raw0815.csv")

#########################################################################################
## get infor for sepecific epitope
## get the HLA allel for each epitope banner

epitope <- as.vector(epitope$AA.Sequence)
HRSVB_F.epimx.data.epitope<- HRSVB_F.epimx.data.4 %>% 
  mutate(Type = ifelse(HRSVB_F.epimx.data.4$AA.Sequence %in% epitope, "T", "F")) 
HRSVB_F.epimx.data.epitope<-HRSVB_F.epimx.data.epitope%>%
  filter(HRSVB_F.epimx.data.epitope$Type=="T")

HRSVB_F.epimx.data.epitope_HLA<-HRSVB_F.epimx.data.epitope%>% 
  group_by(Frame.Start, AA.Sequence,Frame.Stop,HLA.Alleles,Z.score) %>%
  dplyr::summarise(count = n())

HRSVB_F.epimx.data.epitope_HLA_list<-HRSVB_F.epimx.data.epitope_HLA %>% 
  group_by(AA.Sequence) %>% 
  mutate(HLAs = paste0(HLA.Alleles, collapse = ",")) %>%
  select("Frame.Start","AA.Sequence","HLAs")

HRSVB_F.epimx.data.epitope_HLA_list<-distinct(HRSVB_F.epimx.data.epitope_HLA_list)


### get the summarize Z score for each AA.Sequence
HRSVB_F.epimx.data.epitope_HLA_sum<-HRSVB_F.epimx.data.epitope_HLA%>%
  group_by(AA.Sequence)%>%
  summarise(sum_Zscore=sum(as.numeric(Z.score)),mean_Zscore=mean(as.numeric(Z.score)))


##merge all epitope infor dataframe 
epitope_join1<-merge(epitope, HRSVB_F.epimx.data.epitope_HLA_list,by.x="AA.Sequence",by.y="AA.Sequence")
epitope_join2<-merge(epitope_join1,HRSVB_F.epimx.data.epitope_HLA_sum,by.x="AA.Sequence",by.y="AA.Sequence")


write_csv(epitope_join2,"iVAX/HRSVB_epitope_HLA_sum")
##############################################################################
## prepare the dataframe that could map to each sequence
## start with the epitope filter data that ignore the HLA allels
HRSVB_F.epimx.data.epitope_seq<-HRSVB_F.epimx.data.epitope%>%
  select("Accession","AA.Sequence","Type")
HRSVB_F.epimx.data.epitope_seq<-distinct(HRSVB_F.epimx.data.epitope_seq)

## reshape, change the AA Sequence value to column name
HRSVB_F.epimx.data_seq<-dcast(HRSVB_F.epimx.data.epitope_seq, Accession ~ AA.Sequence) 
write_csv(HRSVB_F.epimx.data_seq,"iVAX/HRSVB_epitope_seq_hp")

#####################################################################################################
## prepare for the data that could be used for dollo model

## start with the epitope filter with Z score and hits
## remove the duplicates caused by HLA allels
HRSVB_F.epimx.data.dollo <-HRSVB_F.epimx.data.4%>%
  select("Accession","Frame.Start","AA.Sequence")
HRSVB_F.epimx.data.dollo<-distinct(HRSVB_F.epimx.data.dollo)
## sort 
HRSVB_F.epimx.data.dollo<-HRSVB_F.epimx.data.dollo[order(as.numeric(HRSVB_F.epimx.data.dollo$Frame.Start)),]
HRSVB_F.epimx.data.dollo$code<-'1'
HRSVB_F.epimx.data.dollo <-HRSVB_F.epimx.data.dollo%>%
  select("Accession","AA.Sequence","code")
HRSVB_F.epimx.data.dollo<-dcast(HRSVB_F.epimx.data.dollo, Accession ~ AA.Sequence) 
## replace 'NA' with 0
HRSVB_F.epimx.data.dollo[is.na(HRSVB_F.epimx.data.dollo)] <- 0

paste_noNA <- function(x,sep=", ") {
  gsub(", " ,sep, toString(x[!is.na(x) & x!="" & x!="NA"] ) ) }

sep=""

HRSVB_F.epimx.data.dollo$dollo <- apply( HRSVB_F.epimx.data.dollo[ , c(2:476) ] , 1 , paste_noNA , sep=sep)

HRSVB_F.epimx.data.dollo2<-HRSVB_F.epimx.data.dollo%>%
  select("Accession","dollo")
write_csv(HRSVB_F.epimx.data.dollo2,"iVAX/RSVB_F_sum/RSVB_F_classII_dollo.csv")




### For Babel add "ascertainment"
HRSVB_F.epimx.data.dollo2$Babel <- paste("0", HRSVB_F.epimx.data.dollo2$dollo, sep="")
HRSVB_F.epimx.data.Babel<-HRSVB_F.epimx.data.dollo2%>%
  select("Accession","Babel")
write_csv(HRSVB_F.epimx.data.Babel,"iVAX/RSVB_F_sum/RSVB_F_classII_Babel.csv")