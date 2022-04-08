## RSVB_F_classI epimax sum
## RSVB_F_classII epitope sum
library('tidyverse')
library('tidyr')
library('reshape2')
library('ggplot2')
library('varhandle')
library('RColorBrewer')
library('gridExtra')
library(dplyr)
#############################################################################

RSVA_F_classI<-read.csv("iVAX/CLASS_I_epitope/CLASS1_PROTEIN_REPORT_HRSVB_F_CLEAN0627_PROTEIN_All_Proteins_BY_9.csv", header = FALSE, stringsAsFactors = FALSE, sep = ",")
#####################################################################################################
#to extract each sequence EMX (EpiMatrix) info


  
totalseq <- grep("HRSVB_F_CLEAN0627_PROTEIN ", RSVA_F_classI$V1, useBytes = TRUE)
endseq <- grep("Summarized Results", RSVA_F_classI$V1, useBytes = TRUE)
epimx.data_output <- ""

epimx.data <- RSVA_F_classI%>% slice((totalseq[1]):(endseq[1]-2))
seqname <- gsub("File: HRSVB_F_CLEAN0627_PROTEIN - Sequence: ", "", epimx.data$V1[1])
epimx.data <- epimx.data %>% slice(grep("Frame", V1, useBytes = T):nrow(epimx.data))
colnames(epimx.data) <- c("Frame.Start", "AA.Sequence", "Frame.Stop","Hydro-", "A0101", "A0201", "A0301",  
                          "A2402", "B0702", "B4403", "Hits")
epimx.data <- epimx.data[-c(1:2),]
epimx.data.reshape <- epimx.data %>% mutate(Seq.name = seqname) %>% 
  dplyr::select(Seq.name, Frame.Start:Hits) %>%
  gather(HLA.Alleles, Z.score, "A0101":"B4403", factor_key = TRUE)

epimx.data_output <- rbind(epimx.data_output, epimx.data.reshape)






  
for (i in 1:length(totalseq)){
  epimx.data <- RSVA_F_classI%>% slice((totalseq[i]):(endseq[i]-2))
  seqname <- gsub("File: HRSVB_F_CLEAN0627_PROTEIN - Sequence: ", "", epimx.data$V1[1])
  epimx.data <- epimx.data %>% slice(grep("Frame", V1, useBytes = T):nrow(epimx.data))
  colnames(epimx.data) <- c("Frame.Start", "AA.Sequence", "Frame.Stop","Hydro-", "A0101", "A0201", "A0301",  
                            "A2402", "B0702", "B4403", "Hits")
  epimx.data <- epimx.data[-c(1:2),]
  epimx.data.reshape <- epimx.data %>% mutate(Seq.name = seqname) %>% 
    dplyr::select(Seq.name, Frame.Start:Hits) %>%
    gather(HLA.Alleles, Z.score, "A0101":"B4403", factor_key = TRUE)
  
  epimx.data_output <- rbind(epimx.data_output, epimx.data.reshape)
  
}
  
HRSVB_F_epimax<-epimx.data_output





######################################################################################

#MORE data restructuring 

HRSVB_F.epimx.data.2 <- HRSVB_F_epimax %>% filter(!Seq.name == "")
## remove duplicate from HRSVB_F.epimx.data.2

HRSVB_F.epimx.data.2 <-distinct(HRSVB_F.epimx.data.2)
x <- data.frame(unique(HRSVB_F.epimx.data.2$Seq.name)) ## get all sequence name  

write.csv(HRSVB_F.epimx.data.2,"iVAX/HRSVB_F_classI_epimx1020.csv")  

#####################################################################################  
HRSVB_F.epimx.data.2<-read.csv("iVAX/RSVB_F/RSVB_F_classI_sum/HRSVB_F_classI_epimx1020.csv")
#filter Z SCORE >= 1.64
HRSVB_F.epimx.data.3 <- HRSVB_F.epimx.data.2 %>% filter(Z.score >= 2.328)


#filter EPIBARS Hits >=4

HRSVB_F.epimx.data.3 <- HRSVB_F.epimx.data.3 %>% filter((as.numeric(Hits)) >= 1)

write.csv(HRSVB_F.epimx.data.3 ,"iVAX/HRSVB_F_classI_epimax_filter.csv")

#############################################################################


## seperate the sequence with duplication or without duplication

HRSVB_F.epimx.data.3<-HRSVB_F.epimx.data.3 %>% separate(Seq.name, into = c("Accession", "Subtype", "Country", "date", "genotype"), sep = "-", extra = "merge")

seq <- scan("RAxML/RSVB_F/RSVB_F_date2_accession.txt", what="", sep="\n")


HRSVB_F.epimx.data.3 <- HRSVB_F.epimx.data.3 %>% 
  mutate(Type = ifelse(HRSVB_F.epimx.data.3$Accession %in% seq, "T", "F"))
HRSVB_F.epimx.data.3<-HRSVB_F.epimx.data.3%>% filter(HRSVB_F.epimx.data.3$Type =="T")
##############################################################################################################
## Get the epitope profile for each eptiope frame start site
## with 6 nt deletion at the upstream

###1. coverage at each site
HRSVB_F.epimx.data.4<-HRSVB_F.epimx.data.3%>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

HRSVB_F.epimx.data.4<-HRSVB_F.epimx.data.4%>% 
  group_by(Frame.Start, Frame.Stop) %>%
  dplyr::summarise(count = n())

n<-n_distinct(HRSVB_F.epimx.data.3$Accession)
HRSVB_F.epimx.data.4<- HRSVB_F.epimx.data.4%>% mutate(coverge = count/as.numeric(n))

### 2. get the sepecific epitope profile, including coverage, HLA list and average z score
### coverage for each epitope
HRSVB_F.epimx.data.5<-HRSVB_F.epimx.data.3 %>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

HRSVB_F.epimx.data.5<-HRSVB_F.epimx.data.5%>% 
  group_by(Frame.Start, AA.Sequence,Frame.Stop) %>%
  dplyr::summarise(count = n())

## calculate the coverage percentage
HRSVB_F.epimx.data.5<- HRSVB_F.epimx.data.5%>% mutate(coverge = count/as.numeric(n))


### 3.get the HLA list for each epitope
HRSVB_F.epimx.data.6<-HRSVB_F.epimx.data.3 %>% distinct(Frame.Start,AA.Sequence,Frame.Stop,HLA.Alleles)
HRSVB_F.epimx_HLA_list<-HRSVB_F.epimx.data.6 %>% 
  group_by(AA.Sequence) %>% 
  mutate(HLAs = paste0(HLA.Alleles, collapse = ",")) %>%
  dplyr::select("Frame.Start","AA.Sequence","HLAs")

HRSVB_F.epimx_HLA_list<-distinct(HRSVB_F.epimx_HLA_list)

## 4.get the average Z score for each epitope
HRSVB_F.epimx.data.7<-HRSVB_F.epimx.data.3 %>% distinct(Frame.Start,AA.Sequence,Frame.Stop,HLA.Alleles,Hits,Z.score)

HRSVB_F.epimx.data.7<-HRSVB_F.epimx.data.7%>%
  group_by(Frame.Start,AA.Sequence,Frame.Stop,Hits)%>%
  summarise(sum_Zscore=sum(as.numeric(Z.score)))

HRSVB_F.epimx.data.7<-HRSVB_F.epimx.data.7%>%mutate(mean.score=sum_Zscore/9)

##merge all of the epitope profile together

HRSVB_F_epitope_sum<-left_join(HRSVB_F.epimx.data.5,HRSVB_F.epimx_HLA_list,HRSVB_F.epimx.data.6,by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))

HRSVB_F_epitope_sum<-left_join(HRSVB_F_epitope_sum,HRSVB_F.epimx.data.7,
                                                   by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))

## write these 445 epitope profile
write.csv(HRSVB_F_epitope_sum,"iVAX/RSVB_F/RSVB_F_classI_sum/HRSVB_F_classI_epitope_sum0810.csv",row.names=FALSE)
##############################################################################################################################3
## filter the epitpe site with > 5% coverge 

## filter the interest epitope
epitope_site<-HRSVB_F.epimx.data.4%>%filter(coverge>=0.05)
## from the epitope site to get the epitope profile
epitope_new <-left_join(epitope_site,HRSVB_F_epitope_sum,by=c("Frame.Start" = "Frame.Start"))

## sort the dataframe by Frame start and then number of counts for each epitope

epitope_new <-epitope_new[order(epitope_new$Frame.Start, -epitope_new$count.y),]
## write epitope for hp visulization 54 sites, 370 epitope
write.csv(epitope_new ,"iVAX/RSVB_F/RSVB_F_classI_sum/HRSVB_F_classI_epitope_hplist0813.csv",row.names=FALSE)
## prepare heatmap for visulization

## transform the different epitope at the same location with different states
epitope_new$code<-epitope_new$Frame.Start

epitope_new$code<- make.unique(as.character(epitope_new$code))
epitope_new$code<-as.numeric(epitope_new$code)%%1 *10+1 

epitope_new$code2<-as.character(epitope_new$code)
## if the presence of epitope < 1%
epitope_new<-epitope_new %>% mutate(code = ifelse(coverge.y<=0.01, "x",code2))
## transform to heatmap
## predominant site profile

site<-epitope_new%>%
  filter(code=="1")%>%
  dplyr::select("Frame.Start","AA.Sequence")
site$epitope1<-paste(site$Frame.Start," ",site$AA.Sequence)


HRSVB_F.epimx.data.seq<-HRSVB_F.epimx.data.3 %>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

## add epitope infor
HRSVB_F.epimx.data.seq<-left_join(epitope_new,HRSVB_F.epimx.data.seq,
                                  by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))
## add epitope site infor
HRSVB_F.epimx.data.seq<-left_join(HRSVB_F.epimx.data.seq,site,
                                  by = c("Frame.Start" = "Frame.Start"))
HRSVB_F.epimx.data.hp<-HRSVB_F.epimx.data.seq%>%
  ungroup()%>%
  dplyr::select("Accession","code","epitope1")

## reshape, change the AA Sequence value to column name
HRSVB_F.epimx.data.hp<-dcast(HRSVB_F.epimx.data.hp, Accession ~ epitope1,value.var = "code") 
write.csv(HRSVB_F.epimx.data.hp,"iVAX/RSVB_F/RSVB_F_classI_sum/HRSVB_F_classI_hp_0810.csv",row.names=FALSE)
#################################################################################################















##RSVB_F_classII

RSVB_F_classII<-read.csv("iVAX/RSVB_F_classsII_analysis/HRSVB_F_classII_epimx0815.csv")
#####################################################################################  

#filter Z SCORE >= 1.64
HRSVB_F.epimx.data.3 <- RSVB_F_classII %>% filter(Z.score >= 1.64)


#filter EPIBARS Hits >=4

HRSVB_F.epimx.data.3 <- HRSVB_F.epimx.data.3 %>% filter((as.numeric(Hits)) >= 3)

write.csv(HRSVB_F.epimx.data.3 ,"iVAX/HRSVB_F_classII_epimax_filter.csv")

#############################################################################


## seperate the sequence with duplication or without duplication

HRSVB_F.epimx.data.3<-HRSVB_F.epimx.data.3 %>% separate(Seq.name, into = c("Accession", "Subtype", "Country", "date", "genotype"), sep = "-", extra = "merge")

seq <- scan("RAxML/RAxML_date/RAxML_date2/RSVB_F_date2_accession.txt", what="", sep="\n")


HRSVB_F.epimx.data.3 <- HRSVB_F.epimx.data.3 %>% 
  mutate(Type = ifelse(HRSVB_F.epimx.data.3$Accession %in% seq, "T", "F"))
HRSVB_F.epimx.data.3<-HRSVB_F.epimx.data.3%>% filter(HRSVB_F.epimx.data.3$Type =="T")
##############################################################################################################
## Get the epitope profile for each eptiope frame start site
## with 6 nt deletion at the upstream

###1. coverage at each site
HRSVB_F.epimx.data.4<-HRSVB_F.epimx.data.3%>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

HRSVB_F.epimx.data.4<-HRSVB_F.epimx.data.4%>% 
  group_by(Frame.Start, Frame.Stop) %>%
  dplyr::summarise(count = n())

n<-n_distinct(HRSVB_F.epimx.data.3$Accession)
HRSVB_F.epimx.data.4<- HRSVB_F.epimx.data.4%>% mutate(coverge = count/as.numeric(n))

### 2. get the sepecific epitope profile, including coverage, HLA list and average z score
### coverage for each epitope
HRSVB_F.epimx.data.5<-HRSVB_F.epimx.data.3 %>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

HRSVB_F.epimx.data.5<-HRSVB_F.epimx.data.5%>% 
  group_by(Frame.Start, AA.Sequence,Frame.Stop) %>%
  dplyr::summarise(count = n())

## calculate the coverage percentage
HRSVB_F.epimx.data.5<- HRSVB_F.epimx.data.5%>% mutate(coverge = count/as.numeric(n))


### 3.get the HLA list for each epitope
HRSVB_F.epimx.data.6<-HRSVB_F.epimx.data.3 %>% distinct(Frame.Start,AA.Sequence,Frame.Stop,HLA.Alleles)
HRSVB_F.epimx_HLA_list<-HRSVB_F.epimx.data.6 %>% 
  group_by(AA.Sequence) %>% 
  mutate(HLAs = paste0(HLA.Alleles, collapse = ",")) %>%
  dplyr::select("Frame.Start","AA.Sequence","HLAs")

HRSVB_F.epimx_HLA_list<-distinct(HRSVB_F.epimx_HLA_list)

## 4.get the average Z score for each epitope
HRSVB_F.epimx.data.7<-HRSVB_F.epimx.data.3 %>% distinct(Frame.Start,AA.Sequence,Frame.Stop,HLA.Alleles,Hits,Z.score)

HRSVB_F.epimx.data.7<-HRSVB_F.epimx.data.7%>%
  group_by(Frame.Start,AA.Sequence,Frame.Stop,Hits)%>%
  summarise(sum_Zscore=sum(as.numeric(Z.score)))

HRSVB_F.epimx.data.7<-HRSVB_F.epimx.data.7%>%mutate(mean.score=sum_Zscore/9)

##merge all of the epitope profile together

HRSVB_F_epitope_sum<-left_join(HRSVB_F.epimx.data.5,HRSVB_F.epimx_HLA_list,HRSVB_F.epimx.data.6,by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))

HRSVB_F_epitope_sum<-left_join(HRSVB_F_epitope_sum,HRSVB_F.epimx.data.7,
                               by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))

## write these 445 epitope profile
write.csv(HRSVB_F_epitope_sum,"iVAX/HRSVB_F_classII_epitope_sum1021.csv",row.names=FALSE)
##############################################################################################################################3
## filter the epitpe site with > 5% coverge 

## filter the interest epitope
epitope_site<-HRSVB_F.epimx.data.4%>%filter(coverge>=0.05)
## from the epitope site to get the epitope profile
epitope_new <-left_join(epitope_site,HRSVB_F_epitope_sum,by=c("Frame.Start" = "Frame.Start"))

## sort the dataframe by Frame start and then number of counts for each epitope

epitope_new <-epitope_new[order(epitope_new$Frame.Start, -epitope_new$count.y),]
## write epitope for hp visulization 54 sites, 370 epitope
write.csv(epitope_new ,"iVAX/RSVB_F_classI_epitope_hp.csv",row.names=FALSE)
## prepare heatmap for visulization

## transform the different epitope at the same location with different states
epitope_new$code<-epitope_new$Frame.Start

epitope_new$code<- make.unique(as.character(epitope_new$code))
epitope_new$code<-as.numeric(epitope_new$code)%%1 *10+1 

epitope_new$code2<-as.character(epitope_new$code)
## if the presence of epitope < 1%
epitope_new<-epitope_new %>% mutate(code = ifelse(coverge.y<=0.01, "x",code2))
## transform to heatmap
## predominant site profile

site<-epitope_new%>%
  filter(code=="1")%>%
  dplyr::select("Frame.Start","AA.Sequence")
site$epitope1<-paste(site$Frame.Start," ",site$AA.Sequence)


HRSVB_F.epimx.data.seq<-HRSVB_F.epimx.data.3 %>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

## add epitope infor
HRSVB_F.epimx.data.seq<-left_join(epitope_new,HRSVB_F.epimx.data.seq,
                                  by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))
## add epitope site infor
HRSVB_F.epimx.data.seq<-left_join(HRSVB_F.epimx.data.seq,site,
                                  by = c("Frame.Start" = "Frame.Start"))
HRSVB_F.epimx.data.hp<-HRSVB_F.epimx.data.seq%>%
  ungroup()%>%
  dplyr::select("Accession","code","epitope1")

## reshape, change the AA Sequence value to column name
HRSVB_F.epimx.data.hp<-dcast(HRSVB_F.epimx.data.hp, Accession ~ epitope1,value.var = "code") 
write.csv(HRSVB_F.epimx.data.hp2,"iVAX/HRSVB_F_classII_hp_order.csv",row.names=FALSE)
order<-scan("iVAX/CLASS_II_epitope/RSVB_F_sum/RSVB_F_classII_columnorder.txt", what="", sep="\n")
HRSVB_F.epimx.data.hp2<-HRSVB_F.epimx.data.hp[order]
