## RSVB_G_classII epitope sum
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
filenames <- list.files(path = "./iVAX/RSV_G/RSVB_G_epimax_csv", pattern = "MergedFinalData_dir*")
numfiles <- length(filenames) 

## change file name
for (i in 1:numfiles){  
  name <- gsub(".csv","",filenames[i])
  name <- gsub("MergedFinalData","HRSVB_G_classII",name) 
  assign(name, read.csv(paste("./iVAX/RSV_G/RSVB_G_epimax_csv/", filenames[i], sep = ""), 
                        header = FALSE, stringsAsFactors = FALSE, sep = ","))
}
#####################################################################################################
#to extract each sequence EMX (EpiMatrix) info

EMX<-function(Filedir,epimx.data_output){
  
  totalseq <- grep("HRSVB_G_CLEAN_0819_PROTEIN ", Filedir$V1, useBytes = TRUE)
  endseq <- grep("Summarized Results", Filedir$V1, useBytes = TRUE)
  epimx.data_output <- ""
  
  for (i in 1:length(totalseq)){
    epimx.data <- Filedir%>% slice((totalseq[i]):(endseq[i]-2))
    seqname <- gsub("File: HRSVB_G_CLEAN_0819_PROTEIN - Sequence: ", "", epimx.data$V1[1])
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


HRSVB_G_epimx.data1<-EMX(HRSVB_G_classII_dir1,output)
HRSVB_G_epimx.data2<-EMX(HRSVB_G_classII_dir2,output)
HRSVB_G_epimx.data3<-EMX(HRSVB_G_classII_dir3,output)
HRSVB_G_epimx.data4<-EMX(HRSVB_G_classII_dir4,output)
HRSVB_G_epimx.data5<-EMX(HRSVB_G_classII_dir5,output)
HRSVB_G_epimx.data6<-EMX(HRSVB_G_classII_dir6,output)


######################################################################################
## combine 6 dataframe 

HRSVB_G.epimx.data <- ""
for (i in 1:6){
  next_df <- eval(parse(text=paste("HRSVB_G_epimx.data", i, sep="")))
  HRSVB_G.epimx.data <- rbind(HRSVB_G.epimx.data,next_df)
}

######################################################################################
#MORE data restructuring 

HRSVB_G.epimx.data.2 <- HRSVB_G.epimx.data %>% filter(!Seq.name == "")
## remove duplicate from HRSVB_F.epimx.data.2

HRSVB_G.epimx.data.2 <-distinct(HRSVB_G.epimx.data.2)
x <- data.frame(unique(HRSVB_G.epimx.data.2$Seq.name)) ## get all sequence name  

write.csv(HRSVB_G.epimx.data.2,"iVAX/HRSVB_G_classII_epimx1003.csv")  

#############################################################################


## seperate the sequence with duplication or without duplication

HRSVB_G.epimx.data.2<-HRSVB_G.epimx.data.2 %>% separate(Seq.name, into = c("Accession", "Subtype", "Country", "date", "genotype"), sep = "-", extra = "merge")

seq <- scan("clean_data/RSVB_G_deletion1004_accession.txt", what="", sep="\n")


HRSVB_G.epimx.data.3 <- HRSVB_G.epimx.data.2 %>% 
  mutate(Type = ifelse(HRSVB_G.epimx.data.2$Accession %in% seq, "T", "F"))

HRSVB_G.epimx.data.deletion<-HRSVB_G.epimx.data.3%>% filter(HRSVB_G.epimx.data.3$Type =="T")
HRSVB_G.epimx.data.nodel<-HRSVB_G.epimx.data.3%>% filter(HRSVB_G.epimx.data.3$Type =="F")


HRSVB_G.epimx.data.deletion<-HRSVB_G.epimx.data.deletion[c(1:11)]
HRSVB_G.epimx.data.nodel<-HRSVB_G.epimx.data.nodel[c(1:11)]

write_csv(HRSVB_G.epimx.data.deletion,"iVAX/HRSVB_G_epimx_deletion.csv")
write_csv(HRSVB_G.epimx.data.nodel,"iVAX/HRSVB_G_epimx_nodel.csv")

#####################################################################################  

#filter Z SCORE >= 1.64
HRSVB_G.epimx.data.deletion.2 <- HRSVB_G.epimx.data.deletion %>% filter(Z.score >= 1.64)
HRSVB_G.epimx.data.nodel.2 <- HRSVB_G.epimx.data.nodel %>% filter(Z.score >= 1.64)

#filter EPIBARS Hits >=4

HRSVB_G.epimx.data.deletion.2 <- HRSVB_G.epimx.data.deletion.2  %>% filter((as.numeric(Hits)) >= 4)
HRSVB_G.epimx.data.nodel.2 <- HRSVB_G.epimx.data.nodel.2  %>% filter((as.numeric(Hits)) >= 4)
write.csv(HRSVB_G.epimx.data.deletion.2 ,"iVAX/HRSVB_G_epimx_deletion_filter.csv")
write.csv(HRSVB_G.epimx.data.nodel.2,"iVAX/HRSVB_G_epimx_data_nodel.csv")

##############################################################################################################
##############################################################################################################
## Get the epitope profile for each eptiope frame start site
## with 6 nt deletion at the upstream

###1. coverage at each site
HRSVB_G.epimx.data.deletion.3<-HRSVB_G.epimx.data.deletion.2%>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

HRSVB_G.epimx.data.deletion.3<-HRSVB_G.epimx.data.deletion.3%>% 
  group_by(Frame.Start, Frame.Stop) %>%
  dplyr::summarise(count = n())

n<-n_distinct(HRSVB_G.epimx.data.deletion$Accession)
HRSVB_G.epimx.data.deletion.3<- HRSVB_G.epimx.data.deletion.3%>% mutate(coverge = count/as.numeric(n))

### 2. get the sepecific epitope profile, including coverage, HLA list and average z score
### coverage for each epitope
HRSVB_G.epimx.data.deletion.4<-HRSVB_G.epimx.data.deletion.2 %>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

HRSVB_G.epimx.data.deletion.4<-HRSVB_G.epimx.data.deletion.4%>% 
  group_by(Frame.Start, AA.Sequence,Frame.Stop) %>%
  dplyr::summarise(count = n())

## calculate the coverage percentage
HRSVB_G.epimx.data.deletion.4<- HRSVB_G.epimx.data.deletion.4%>% mutate(coverge = count/as.numeric(n))


### 3.get the HLA list for each epitope
HRSVB_G.epimx.data.deletion.5<-HRSVB_G.epimx.data.deletion.2 %>% distinct(Frame.Start,AA.Sequence,Frame.Stop,HLA.Alleles)
HRSVB_G.epimx.data.deletion.epitope_HLA_list<-HRSVB_G.epimx.data.deletion.5 %>% 
  group_by(AA.Sequence) %>% 
  mutate(HLAs = paste0(HLA.Alleles, collapse = ",")) %>%
  select("Frame.Start","AA.Sequence","HLAs")

HRSVB_G.epimx.data.deletion.epitope_HLA_list<-distinct(HRSVB_G.epimx.data.deletion.epitope_HLA_list)

## 4.get the average Z score for each epitope
HRSVB_G.epimx.data.deletion.6<-HRSVB_G.epimx.data.deletion.2 %>% distinct(Frame.Start,AA.Sequence,Frame.Stop,HLA.Alleles,Hits,Z.score)

HRSVB_G.epimx.data.deletion.6<-HRSVB_G.epimx.data.deletion.6%>%
  group_by(Frame.Start,AA.Sequence,Frame.Stop,Hits)%>%
  summarise(sum_Zscore=sum(as.numeric(Z.score)))

HRSVB_G.epimx.data.deletion.6<-HRSVB_G.epimx.data.deletion.6%>%mutate(mean.score=sum_Zscore/9)

##merge all of the epitope profile together

HRSVB_G.epimx.data.deletion_epitope_sum<-left_join(HRSVB_G.epimx.data.deletion.4,HRSVB_G.epimx.data.deletion.epitope_HLA_list,HRSVB_G.epimx.data.deletion.3,
                                  by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))

HRSVB_G.epimx.data.deletion_epitope_sum<-left_join(HRSVB_G.epimx.data.deletion_epitope_sum,HRSVB_G.epimx.data.deletion.6,
                                  by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))

## write these 445 epitope profile
write.csv(HRSVB_G.epimx.data.deletion_epitope_sum,"iVAX/HRSVB_G.epimx.data.deletion_epitope_sum1005.csv",row.names=FALSE)
##############################################################################################################################

## without deletion sum

###1. coverage at each site
HRSVB_G.epimx.data.nodel.3<-HRSVB_G.epimx.data.nodel.2%>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

HRSVB_G.epimx.data.nodel.3<-HRSVB_G.epimx.data.nodel.3%>% 
  group_by(Frame.Start, Frame.Stop) %>%
  dplyr::summarise(count = n())

n<-n_distinct(HRSVB_G.epimx.data.nodel$Accession)
HRSVB_G.epimx.data.nodel.3<- HRSVB_G.epimx.data.nodel.3%>% mutate(coverge = count/as.numeric(n))

### 2. get the sepecific epitope profile, including coverage, HLA list and average z score
### coverage for each epitope
HRSVB_G.epimx.data.nodel.4<-HRSVB_G.epimx.data.nodel.2 %>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

HRSVB_G.epimx.data.nodel.4<-HRSVB_G.epimx.data.nodel.4%>% 
  group_by(Frame.Start, AA.Sequence,Frame.Stop) %>%
  dplyr::summarise(count = n())

## calculate the coverage percentage
HRSVB_G.epimx.data.nodel.4<- HRSVB_G.epimx.data.nodel.4%>% mutate(coverge = count/as.numeric(n))


### 3.get the HLA list for each epitope
HRSVB_G.epimx.data.nodel.5<-HRSVB_G.epimx.data.nodel.2 %>% distinct(Frame.Start,AA.Sequence,Frame.Stop,HLA.Alleles)
HRSVB_G.epimx.data.nodel.epitope_HLA_list<-HRSVB_G.epimx.data.nodel.5 %>% 
  group_by(AA.Sequence) %>% 
  mutate(HLAs = paste0(HLA.Alleles, collapse = ",")) %>%
  select("Frame.Start","AA.Sequence","HLAs")

HRSVB_G.epimx.data.nodel.epitope_HLA_list<-distinct(HRSVB_G.epimx.data.nodel.epitope_HLA_list)

## 4.get the average Z score for each epitope
HRSVB_G.epimx.data.nodel.6<-HRSVB_G.epimx.data.nodel.2 %>% distinct(Frame.Start,AA.Sequence,Frame.Stop,HLA.Alleles,Hits,Z.score)

HRSVB_G.epimx.data.nodel.6<-HRSVB_G.epimx.data.nodel.6%>%
  group_by(Frame.Start,AA.Sequence,Frame.Stop,Hits)%>%
  summarise(sum_Zscore=sum(as.numeric(Z.score)))

HRSVB_G.epimx.data.nodel.6<-HRSVB_G.epimx.data.nodel.6%>%mutate(mean.score=sum_Zscore/9)

##merge all of the epitope profile together

HRSVB_G_nodel_epitope_sum<-left_join(HRSVB_G.epimx.data.nodel.4,HRSVB_G.epimx.data.nodel.epitope_HLA_list,HRSVB_G.epimx.data.nodel.3,
                                  by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))

HRSVB_G_nodel_epitope_sum<-left_join(HRSVB_G_nodel_epitope_sum,HRSVB_G.epimx.data.nodel.6,
                                  by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))

## write these 445 epitope profile
write.csv(HRSVB_G_nodel_epitope_sum,"iVAX/HRSVB_G_nodel_epitope_sum1005.csv",row.names=FALSE)
############################################################################################################################3
#################################################################################################################
## try to sum ON1 and GA together

#filter Z SCORE >= 1.64
HRSVB_G.epimx.data.3 <- HRSVB_G.epimx.data.2 %>% filter(Z.score >= 1.64)

#filter EPIBARS Hits >=4

HRSVB_G.epimx.data.4 <- HRSVB_G.epimx.data.3  %>% filter((as.numeric(Hits)) >= 4)

write.csv(HRSVB_G.epimx.data.4,"iVAX/HRSVB_G_epimx_filter.csv")
############################################################################################################################


###1. coverage at each site
HRSVB_G.epimx.data.5<-HRSVB_G.epimx.data.4%>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

HRSVB_G.epimx.data.5<-HRSVB_G.epimx.data.5%>% 
  group_by(Frame.Start, Frame.Stop) %>%
  dplyr::summarise(count = n())

n<-n_distinct(HRSVB_G.epimx.data.4$Accession)
HRSVB_G.epimx.data.5<- HRSVB_G.epimx.data.5%>% mutate(coverge = count/as.numeric(n))

### 2. get the sepecific epitope profile, including coverage, HLA list and average z score
### coverage for each epitope
HRSVB_G.epimx.data.6<-HRSVB_G.epimx.data.4 %>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

HRSVB_G.epimx.data.6<-HRSVB_G.epimx.data.6%>% 
  group_by(Frame.Start, AA.Sequence,Frame.Stop) %>%
  dplyr::summarise(count = n())

## calculate the coverage percentage
HRSVB_G.epimx.data.6<- HRSVB_G.epimx.data.6%>% mutate(coverge = count/as.numeric(n))


### 3.get the HLA list for each epitope
HRSVB_G.epimx.data.7<-HRSVB_G.epimx.data.4 %>% distinct(Frame.Start,AA.Sequence,Frame.Stop,HLA.Alleles)
HRSVB_G.epitope_HLA_list<-HRSVB_G.epimx.data.7 %>% 
  group_by(AA.Sequence) %>% 
  mutate(HLAs = paste0(HLA.Alleles, collapse = ",")) %>%
  select("Frame.Start","AA.Sequence","HLAs")

HRSVB_G.epitope_HLA_list<-distinct(HRSVB_G.epitope_HLA_list)

## 4.get the average Z score for each epitope
HRSVB_G.epimx.data.8<-HRSVB_G.epimx.data.4 %>% distinct(Frame.Start,AA.Sequence,Frame.Stop,HLA.Alleles,Hits,Z.score)

HRSVB_G.epimx.data.8<-HRSVB_G.epimx.data.8%>%
  group_by(Frame.Start,AA.Sequence,Frame.Stop,Hits)%>%
  summarise(sum_Zscore=sum(as.numeric(Z.score)))

HRSVB_G.epimx.data.8<-HRSVB_G.epimx.data.8%>%mutate(mean.score=sum_Zscore/9)

##merge all of the epitope profile together

HRSVB_G_epitope_sum<-left_join(HRSVB_G.epimx.data.6,HRSVB_G.epitope_HLA_list,HRSVB_G.epimx.data.5,
                               by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))

HRSVB_G_epitope_sum<-left_join(HRSVB_G_epitope_sum,HRSVB_G.epimx.data.8,
                               by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))

## write these 445 epitope profile
write.csv(HRSVB_G_epitope_sum,"iVAX/HRSVB_G_epitope_sum1005.csv",row.names=FALSE)
##############################################################################################################################3
## filter the epitpe site with > 5% coverge 

## filter the interest epitope
epitope_site<-HRSVB_G.epimx.data.5%>%filter(coverge>=0.05)
## from the epitope site to get the epitope profile
epitope_new <-left_join(epitope_site,HRSVB_G_epitope_sum,by=c("Frame.Start" = "Frame.Start"))

## sort the dataframe by Frame start and then number of counts for each epitope

epitope_new <-epitope_new[order(epitope_new$Frame.Start, -epitope_new$count.y),]
## write epitope for hp visulization 54 sites, 370 epitope
write.csv(epitope_new ,"iVAX/HRSVB_G_classII_epitopehp_1005.csv",row.names=FALSE)
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


HRSVB_G.epimx.data.seq<-HRSVB_G.epimx.data.2 %>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

## add epitope infor
HRSVB_G.epimx.data.seq<-left_join(epitope_new,HRSVB_G.epimx.data.seq,
                                  by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))
## add epitope site infor
HRSVB_G.epimx.data.seq<-left_join(HRSVB_G.epimx.data.seq,site,
                                  by = c("Frame.Start" = "Frame.Start"))
HRSVB_G.epimx.data.hp<-HRSVB_G.epimx.data.seq%>%
  ungroup()%>%
  select("Accession","code","epitope1")

## reshape, change the AA Sequence value to column name
HRSVB_G.epimx.data.hp<-dcast(HRSVB_G.epimx.data.hp, Accession ~ epitope1,value.var = "code") 
write.csv(HRSVB_G.epimx.data.hp,"iVAX/HRSVB_G_hp.csv",row.names=FALSE)


