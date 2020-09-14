## HRSVA-F
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
filenames <- list.files(path = "./iVAX/HRSVA_F_classII_csv", pattern = "MergedFinalData0*")
numfiles <- length(filenames) 

## change file name
for (i in 1:numfiles){  
  name <- gsub(".csv","",filenames[i])
  name <- gsub("MergedFinalData","HRSVA_F_classII",name) 
  assign(name, read.csv(paste("./iVAX/HRSVA_F_classII_csv/", filenames[i], sep = ""), 
                        header = FALSE, stringsAsFactors = FALSE, sep = ","))
}


HRSVA_F_miss<-read.csv("iVAX/CLASS2_PROTEIN_REPORT_HRSVA_F_IVAXMISSING0822_All_Proteins.csv",header = FALSE, stringsAsFactors = FALSE, sep = ",")
#########################################################################
#to extract each sequence EMX (EpiMatrix) info

EMX<-function(Filedir,epimx.data_output){
  
  totalseq <- grep("File: HRSVA_F_CLEAN0627_PROTEIN", Filedir$V1, useBytes = TRUE)
  endseq <- grep("Summarized Results", Filedir$V1, useBytes = TRUE)
  epimx.data_output <- ""
  
  for (i in 1:length(totalseq)){
    epimx.data <- Filedir%>% slice((totalseq[i]):(endseq[i]-2))
    seqname <- gsub("File: HRSVA_F_CLEAN0627_PROTEIN - Sequence: ", "", epimx.data$V1[1])
    epimx.data <- epimx.data %>% slice(grep("Frame", V1, useBytes = T):nrow(epimx.data))
    colnames(epimx.data) <- c("Frame.Start", "AA.Sequence", "Frame.Stop", "DRB1*0101", "DRB1*0301", "DRB1*0401",  
                              "DRB1*0701", "DRB1*0801", "DRB1*0901", "DRB1*1101", "DRB1*1301", "DRB1*1501",  
                              "Hits", "NA")
    epimx.data <- epimx.data[-c(1:2),]
    epimx.data.reshape <- epimx.data %>% mutate(Seq.name = seqname) %>%
      select(Seq.name, Frame.Start:Hits) %>%
      gather(HLA.Alleles, Z.score, "DRB1*0101":"DRB1*1501", factor_key = TRUE)
    epimx.data_output <- rbind(epimx.data_output, epimx.data.reshape)
    
  }
    
  return(epimx.data_output)

}


EMX2<-function(Filedir,epimx.data_output){
  
  totalseq <- grep("File: HRSVA_F_IVAXMISSING0822", Filedir$V1, useBytes = TRUE)
  endseq <- grep("Summarized Results", Filedir$V1, useBytes = TRUE)
  epimx.data_output <- ""
  
  for (i in 1:length(totalseq)){
    epimx.data <- Filedir%>% slice((totalseq[i]):(endseq[i]-2))
    seqname <- gsub("File: HRSVA_F_IVAXMISSING0822 - Sequence: ", "", epimx.data$V1[1])
    epimx.data <- epimx.data %>% slice(grep("Frame", V1, useBytes = T):nrow(epimx.data))
    colnames(epimx.data) <- c("Frame.Start", "AA.Sequence", "Frame.Stop", "DRB1*0101", "DRB1*0301", "DRB1*0401",  
                              "DRB1*0701", "DRB1*0801", "DRB1*0901", "DRB1*1101", "DRB1*1301", "DRB1*1501",  
                              "Hits", "NA")
    epimx.data <- epimx.data[-c(1:2),]
    epimx.data.reshape <- epimx.data %>% mutate(Seq.name = seqname) %>%
      select(Seq.name, Frame.Start:Hits) %>%
      gather(HLA.Alleles, Z.score, "DRB1*0101":"DRB1*1501", factor_key = TRUE)
    epimx.data_output <- rbind(epimx.data_output, epimx.data.reshape)
    
  }
  
  return(epimx.data_output)
  
}




HRSVA_F_epimx.data01<-EMX(HRSVA_F_classII01,output)
HRSVA_F_epimx.data02<-EMX(HRSVA_F_classII02,output)
HRSVA_F_epimx.data03<-EMX(HRSVA_F_classII03,output)
HRSVA_F_epimx.data04<-EMX(HRSVA_F_classII04,output)
HRSVA_F_epimx.data05<-EMX(HRSVA_F_classII05,output)
HRSVA_F_epimx.data06<-EMX(HRSVA_F_classII06,output)
HRSVA_F_epimx.data07<-EMX(HRSVA_F_classII07,output)
HRSVA_F_epimx.data08<-EMX(HRSVA_F_classII08,output)
HRSVA_F_epimx.data09<-EMX(HRSVA_F_classII09,output)
HRSVA_F_epimx.data10<-EMX(HRSVA_F_classII10,output)
HRSVA_F_epimx.data11<-EMX(HRSVA_F_classII11,output)
HRSVA_F_epimx.data12<-EMX(HRSVA_F_classII12,output)
HRSVA_F_epimx.data13<-EMX(HRSVA_F_classII13,output)
HRSVA_F_epimx.data14<-EMX2(HRSVA_F_miss,output)
######################################################################################
## combine 13 dataframe 

#HRSVA_F.epimx.data <- ""
for (i in 10:14){
  next_df <- eval(parse(text=paste("HRSVA_F_epimx.data", i, sep="")))
  HRSVA_F.epimx.data <- rbind(HRSVA_F.epimx.data,next_df)
}

######################################################################################
#MORE data restructuring 

HRSVA_F.epimx.data.2 <- HRSVA_F.epimx.data %>% filter(!Seq.name == "")
## remove duplicate from HRSVB_F.epimx.data.2

HRSVA_F.epimx.data.2 <-distinct(HRSVA_F.epimx.data.2)
x <- data.frame(unique(HRSVA_F.epimx.data.2$Seq.name)) ## get all sequence name  

write.csv(HRSVA_F.epimx.data.2,"iVAX/HRSVA_F_classII_epimx0831.csv")  
#HRSVA_F.epimx.data.2<-read.csv("iVAX/HRSVA_F_classII_epimx0812.csv")
#############################################################################

## filter with a clean sequence dataset

HRSVA_F.epimx.data.2<-HRSVA_F.epimx.data.2 %>% separate(Seq.name, into = c("Accession", "Subtype", "Country", "date", "genotype"), sep = "-", extra = "merge")
seq <- scan("./RAxML/RAxML0708/RAxML_0708filter/HRSVA_F/HRSVA_F_accession.txt", what="", sep="\n")

HRSVA_F.epimx.data.3 <- HRSVA_F.epimx.data.2 %>% 
  mutate(Type = ifelse(HRSVA_F.epimx.data.2$Accession %in% seq, "T", "F"))

HRSVA_F.epimx.data.3<-HRSVA_F.epimx.data.3%>% filter(HRSVA_F.epimx.data.3$Type =="T")
HRSVA_F.epimx.data.clean<-HRSVA_F.epimx.data.3[c(1:11)]
write_csv(HRSVA_F.epimx.data.clean,"iVAX/HRSVA_F_classII_epimx0831_clean.csv")
#####################################################################################  

#filter Z SCORE >= 1.64
HRSVA_F.epimx.data.3 <- HRSVA_F.epimx.data.clean %>% filter(Z.score >= 1.64)
#write.csv(HRSVB_F.epimx.data.3,"iVAX/HRSVB_F_classII_epimx0810_filter.csv")
#HRSVB_F.epimx.data.3$Z.score <- as.numeric(HRSVB_F.epimx.data.3$Z.score)

###################################################################################
#filter EPIBARS Hits >=4
HRSVA_F.epimx.data.3$Hits<- as.numeric(HRSVA_F.epimx.data.3$Hits)
HRSVA_F.epimx.data.4<- HRSVA_F.epimx.data.3 %>% filter(Hits >= 4)
HRSVA_F.epimx.data.4<- HRSVA_F.epimx.data.4[!(HRSVA_F.epimx.data.4$`AA.Sequence` == ""), ]
write.csv(HRSVA_F.epimx.data.4,"iVAX/HRSVA_F_classII_epimx0831_filter.csv")

##############################################################################################################
## Get the epitope profile for each eptiope frame start site
###1. coverage at each site
HRSVA_F.epimx.data.5<-HRSVA_F.epimx.data.4 %>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

HRSVA_F.epimx.data.5<-HRSVA_F.epimx.data.5%>% 
  group_by(Frame.Start, Frame.Stop) %>%
  dplyr::summarise(count = n())

n<-n_distinct(HRSVA_F.epimx.data.clean$Accession)
HRSVA_F.epimx.data.5<- HRSVA_F.epimx.data.5%>% mutate(coverge = count/as.numeric(n))

### 2. get the sepecific epitope profile, including coverage, HLA list and average z score
### coverage for each epitope
HRSVA_F.epimx.data.6<-HRSVA_F.epimx.data.4 %>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

HRSVA_F.epimx.data.6<-HRSVA_F.epimx.data.6%>% 
  group_by(Frame.Start, AA.Sequence,Frame.Stop) %>%
  dplyr::summarise(count = n())

## calculate the coverage percentage
HRSVA_F.epimx.data.6<- HRSVA_F.epimx.data.6%>% mutate(coverge = count/as.numeric(n))


### get the HLA list for each epitope
HRSVA_F.epimx.data.7<-HRSVA_F.epimx.data.4 %>% distinct(Frame.Start,AA.Sequence,Frame.Stop,HLA.Alleles)
HRSVA_F.epimx.data.epitope_HLA_list<-HRSVA_F.epimx.data.7 %>% 
  group_by(AA.Sequence) %>% 
  mutate(HLAs = paste0(HLA.Alleles, collapse = ",")) %>%
  select("Frame.Start","AA.Sequence","HLAs")

HRSVA_F.epimx.data.epitope_HLA_list<-distinct(HRSVA_F.epimx.data.epitope_HLA_list)

## get the average Z score for each epitope
HRSVA_F.epimx.data.8<-HRSVA_F.epimx.data.4 %>% distinct(Frame.Start,AA.Sequence,Frame.Stop,HLA.Alleles,Hits,Z.score)

HRSVA_F.epimx.data.8<-HRSVA_F.epimx.data.8%>%
  group_by(Frame.Start,AA.Sequence,Frame.Stop,Hits)%>%
  summarise(sum_Zscore=sum(as.numeric(Z.score)))

HRSVA_F.epimx.data.8<-HRSVA_F.epimx.data.8%>%mutate(mean.score=sum_Zscore/9)

##merge all of the epitope profile together

HRSVA_F_epitope_sum<-left_join(HRSVA_F.epimx.data.6,HRSVA_F.epimx.data.epitope_HLA_list,HRSVB_F.epimx.data.8,
                               by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))

HRSVA_F_epitope_sum<-left_join(HRSVA_F_epitope_sum,HRSVA_F.epimx.data.8,
                               by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))

## write these 445 epitope profile
write.csv(HRSVA_F_epitope_sum,"iVAX/HRSVA_F_classII_epitopesum_0831.csv",row.names=FALSE)

############################################################################################################################

## prepare the visulization for heatmap
## filter the epitpe site with > 5% coverge 

## filter the interest epitope
epitope_site<-HRSVA_F.epimx.data.5%>%filter(coverge>=0.05)
## from the epitope site to get the epitope profile
epitope_new <-left_join(epitope_site,HRSVA_F_epitope_sum,by=c("Frame.Start" = "Frame.Start"))

## sort the dataframe by Frame start and then number of counts for each epitope

epitope_new <-epitope_new[order( epitope_new$Frame.Start, -epitope_new$count.y),]
## write epitope for hp visulization 54 sites, 370 epitope
write.csv(epitope_new,"iVAXHRSVA_F_classII_epitopehp_0831.csv",row.names=FALSE)
## prepare heatmap for visulization

## transform the different epitope at the same location with different states
epitope_new$code<-epitope_new$Frame.Start

epitope_new$code<- make.unique(as.character(epitope_new$code))
epitope_new$code<-as.numeric(epitope_new$code)%%1 *10+1 


## if the presence of epitope < 1%
epitope_new<-epitope_new %>% mutate(code = ifelse(coverge.y<=0.01, "x",as.character(epitope_new$code)))
## transform to heatmap
## predominant site profile

site<-epitope_new%>%
  filter(code=="1")%>%
  select("Frame.Start","AA.Sequence")
site$epitope1<-paste(site$Frame.Start," ",site$AA.Sequence)
HRSVA_F.epimx.data.seq<-HRSVA_F.epimx.data.4 %>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

## add epitope infor
HRSVA_F.epimx.data.seq<-left_join(epitope_new,HRSVA_F.epimx.data.seq,
                                  by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))
## add epitope site infor
HRSVA_F.epimx.data.seq<-left_join(HRSVA_F.epimx.data.seq,site,
                                  by = c("Frame.Start" = "Frame.Start"))
HRSVA_F.epimx.data.hp<-HRSVA_F.epimx.data.seq%>%
  ungroup()%>%
  select("Accession","code","epitope1")

## reshape, change the AA Sequence value to column name
HRSVA_F.epimx.data.hp<-dcast(HRSVA_F.epimx.data.hp, Accession ~ epitope1,value.var = "code") 
write.csv(HRSVA_F.epimx.data.hp,"iVAX/HRSVA_F_classII_hp.csv",row.names=FALSE)




























#########################################################################################
## identify all of the potential epitope banner from the matix

## for each sequence, remove the duplicate from HLA allel 

HRSVA_F.epimx.data.5<-HRSVA_F.epimx.data.4 %>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

HRSVA_F.epimx.data.6<-HRSVA_F.epimx.data.5%>% 
  group_by(Frame.Start, AA.Sequence,Frame.Stop) %>%
  dplyr::summarise(count = n())

## calculate the coverage percentage
n<-n_distinct(HRSVA_F.epimx.data.clean$Accession)
HRSVA_F.epimx.data.6<- HRSVA_F.epimx.data.6%>% mutate(coverge = count/as.numeric(n))

## filter the interest epitope
epitope<-HRSVA_F.epimx.data.6%>%filter(coverge>=0.05)
epitope<-epitope[order(as.numeric(epitope$Frame.Start)),]
write_csv(epitope,"iVAX/HRSVA_epitope_raw0820.csv")

#########################################################################################
## get infor for sepecific epitope
## get the HLA allel for each epitope banner

epitope <- as.vector(epitope$AA.Sequence)
HRSVA_F.epimx.data.epitope<- HRSVA_F.epimx.data.4 %>% 
  mutate(Type = ifelse(HRSVA_F.epimx.data.4$AA.Sequence %in% epitope, "T", "F")) 


HRSVA_F.epimx.data.epitope<-HRSVA_F.epimx.data.epitope%>%
  filter(HRSVA_F.epimx.data.epitope$Type=="T")

HRSVA_F.epimx.data.epitope_HLA<-HRSVA_F.epimx.data.epitope%>% 
  group_by(Frame.Start, AA.Sequence,Frame.Stop,HLA.Alleles,Z.score) %>%
  dplyr::summarise(count = n())

HRSVA_F.epimx.data.epitope_HLA_list<-HRSVA_F.epimx.data.epitope_HLA %>% 
  group_by(AA.Sequence) %>% 
  mutate(HLAs = paste0(HLA.Alleles, collapse = ",")) %>%
  select("Frame.Start","AA.Sequence","HLAs")

HRSVA_F.epimx.data.epitope_HLA_list<-distinct(HRSVA_F.epimx.data.epitope_HLA_list)
##merge all epitope infor dataframe 
epitope_join1<-merge(epitope, HRSVA_F.epimx.data.epitope_HLA_list,by.x="AA.Sequence",by.y="AA.Sequence")
write_csv(epitope_join1,"iVAX/HRSVA_F_epitope_HLA_sum")

##############################################################################
## prepare the dataframe that could map to each sequence
## start with the epitope filter data that ignore the HLA allels
HRSVA_F.epimx.data.epitope_seq<-HRSVA_F.epimx.data.epitope%>%
  select("Accession","AA.Sequence","Type")
HRSVA_F.epimx.data.epitope_seq<-distinct(HRSVA_F.epimx.data.epitope_seq)

## reshape, change the AA Sequence value to column name
HRSVA_F.epimx.data_seq<-dcast(HRSVA_F.epimx.data.epitope_seq, Accession ~ AA.Sequence) 


write_csv(HRSVA_F.epimx.data_seq,"iVAX/HRSVA_F_epitope_seq_hp")

