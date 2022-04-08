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
filenames <- list.files(path = "./iVAX/RSV_G/RSVA_G_epimax_csv", pattern = "MergedFinalData_dir*")
numfiles <- length(filenames) 

## change file name
for (i in 1:numfiles){  
  name <- gsub(".csv","",filenames[i])
  name <- gsub("MergedFinalData","HRSVA_G_classII",name) 
  assign(name, read.csv(paste("./iVAX/RSV_G/RSVA_G_epimax_csv/", filenames[i], sep = ""), 
                        header = FALSE, stringsAsFactors = FALSE, sep = ","))
}
#####################################################################################################
#to extract each sequence EMX (EpiMatrix) info

EMX<-function(Filedir,epimx.data_output){
  
  totalseq <- grep("File: HRSVA_G_CLEAN0819_PROTEIN", Filedir$V1, useBytes = TRUE)
  endseq <- grep("Summarized Results", Filedir$V1, useBytes = TRUE)
  epimx.data_output <- ""
  
  for (i in 1:length(totalseq)){
    epimx.data <- Filedir%>% slice((totalseq[i]):(endseq[i]-2))
    seqname <- gsub("File: HRSVA_G_CLEAN0819_PROTEIN - Sequence: ", "", epimx.data$V1[1])
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


HRSVA_G_epimx.data1<-EMX(HRSVA_G_classII_dir1,output)
HRSVA_G_epimx.data2<-EMX(HRSVA_G_classII_dir2,output)
HRSVA_G_epimx.data3<-EMX(HRSVA_G_classII_dir3,output)
HRSVA_G_epimx.data4<-EMX(HRSVA_G_classII_dir4,output)
HRSVA_G_epimx.data5<-EMX(HRSVA_G_classII_dir5,output)
HRSVA_G_epimx.data6<-EMX(HRSVA_G_classII_dir6,output)
HRSVA_G_epimx.data7<-EMX(HRSVA_G_classII_dir7,output)
HRSVA_G_epimx.data8<-EMX(HRSVA_G_classII_dir8,output)
HRSVA_G_epimx.data9<-EMX(HRSVA_G_classII_dir9,output)

######################################################################################
## combine 9 dataframe 

HRSVA_G.epimx.data <- ""
for (i in 1:9){
  next_df <- eval(parse(text=paste("HRSVA_G_epimx.data", i, sep="")))
  HRSVA_G.epimx.data <- rbind(HRSVA_G.epimx.data,next_df)
}

######################################################################################
#MORE data restructuring 

HRSVA_G.epimx.data.2 <- HRSVA_G.epimx.data %>% filter(!Seq.name == "")
## remove duplicate from HRSVB_F.epimx.data.2

HRSVA_G.epimx.data.2 <-distinct(HRSVA_G.epimx.data.2)
x <- data.frame(unique(HRSVA_G.epimx.data.2$Seq.name)) ## get all sequence name  

write.csv(HRSVA_G.epimx.data.2,"iVAX/HRSVA_G_classII_epimx0927.csv")  

#############################################################################


## seperate the sequence with duplication or without duplication

HRSVA_G.epimx.data.2<-HRSVA_G.epimx.data.2 %>% separate(Seq.name, into = c("Accession", "Subtype", "Country", "date", "genotype"), sep = "-", extra = "merge")
seq <- scan("clean_data/RSVA_G_Duplication_accession.txt", what="", sep="\n")


HRSVA_G.epimx.data.3 <- HRSVA_G.epimx.data.2 %>% 
  mutate(Type = ifelse(HRSVA_G.epimx.data.2$Accession %in% seq, "T", "F"))

HRSVA_G.epimx.data.ON<-HRSVA_G.epimx.data.3%>% filter(HRSVA_G.epimx.data.3$Type =="T")
HRSVA_G.epimx.data.GA<-HRSVA_G.epimx.data.3%>% filter(HRSVA_G.epimx.data.3$Type =="F")


HRSVA_G.epimx.data.ON<-HRSVA_G.epimx.data.ON[c(1:11)]
HRSVA_G.epimx.data.GA<-HRSVA_G.epimx.data.GA[c(1:11)]

write_csv(HRSVA_G.epimx.data.ON,"iVAX/HRSVA_G_epimx_ON.csv")
write_csv(HRSVA_G.epimx.data.GA,"iVAX/HRSVA_G_epimx_GA.csv")
#####################################################################################  
#####################################################################################  

#filter Z SCORE >= 1.64
HRSVA_G.epimx.data.ON.2 <- HRSVA_G.epimx.data.ON %>% filter(Z.score >= 1.64)
HRSVA_G.epimx.data.GA.2 <- HRSVA_G.epimx.data.GA %>% filter(Z.score >= 1.64)

###################################################################################
#filter EPIBARS Hits >=4

HRSVA_G.epimx.data.ON.2 <- HRSVA_G.epimx.data.ON.2  %>% filter((as.numeric(Hits)) >= 4)
HRSVA_G.epimx.data.GA.2 <- HRSVA_G.epimx.data.GA.2  %>% filter((as.numeric(Hits)) >= 4)
write.csv(HRSVA_G.epimx.data.ON.2,"iVAX/HRSVA_G_epimx_ON_filter.csv")
write.csv(HRSVA_G.epimx.data.GA.2,"iVAX/HRSVA_G_epimx_GA_filter.csv")

##############################################################################################################
##############################################################################################################
## Get the epitope profile for each eptiope frame start site
## G duplication ON

###1. coverage at each site
HRSVA_G.epimx.data.ON.3<-HRSVA_G.epimx.data.ON.2%>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

HRSVA_G.epimx.data.ON.3<-HRSVA_G.epimx.data.ON.3%>% 
  group_by(Frame.Start, Frame.Stop) %>%
  dplyr::summarise(count = n())

n<-n_distinct(HRSVA_G.epimx.data.ON$Accession)
HRSVA_G.epimx.data.ON.3<- HRSVA_G.epimx.data.ON.3%>% mutate(coverge = count/as.numeric(n))

### 2. get the sepecific epitope profile, including coverage, HLA list and average z score
### coverage for each epitope
HRSVA_G.epimx.data.ON.4<-HRSVA_G.epimx.data.ON.2 %>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

HRSVA_G.epimx.data.ON.4<-HRSVA_G.epimx.data.ON.4%>% 
  group_by(Frame.Start, AA.Sequence,Frame.Stop) %>%
  dplyr::summarise(count = n())

## calculate the coverage percentage
HRSVA_G.epimx.data.ON.4<- HRSVA_G.epimx.data.ON.4%>% mutate(coverge = count/as.numeric(n))


### 3.get the HLA list for each epitope
HRSVA_G.epimx.data.ON.5<-HRSVA_G.epimx.data.ON.2 %>% distinct(Frame.Start,AA.Sequence,Frame.Stop,HLA.Alleles)
HRSVA_G.epimx.data.ON.epitope_HLA_list<-HRSVA_G.epimx.data.ON.5 %>% 
  group_by(AA.Sequence) %>% 
  mutate(HLAs = paste0(HLA.Alleles, collapse = ",")) %>%
  select("Frame.Start","AA.Sequence","HLAs")

HRSVA_G.epimx.data.ON.epitope_HLA_list<-distinct(HRSVA_G.epimx.data.ON.epitope_HLA_list)

## 4.get the average Z score for each epitope
HRSVA_G.epimx.data.ON.6<-HRSVA_G.epimx.data.ON.2 %>% distinct(Frame.Start,AA.Sequence,Frame.Stop,HLA.Alleles,Hits,Z.score)

HRSVA_G.epimx.data.ON.6<-HRSVA_G.epimx.data.ON.6%>%
  group_by(Frame.Start,AA.Sequence,Frame.Stop,Hits)%>%
  summarise(sum_Zscore=sum(as.numeric(Z.score)))

HRSVA_G.epimx.data.ON.6<-HRSVA_G.epimx.data.ON.6%>%mutate(mean.score=sum_Zscore/9)

##merge all of the epitope profile together

HRSVA_G_ON_epitope_sum<-left_join(HRSVA_G.epimx.data.ON.4,HRSVA_G.epimx.data.ON.epitope_HLA_list,HRSVA_G.epimx.data.ON.3,
                               by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))

HRSVA_G_ON_epitope_sum<-left_join(HRSVA_G_ON_epitope_sum,HRSVA_G.epimx.data.ON.6,
                               by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))

## write these 445 epitope profile
write.csv(HRSVA_G_ON_epitope_sum,"iVAX/HRSVA_G_ON_epitope_sum0927.csv",row.names=FALSE)

############################################################################################################################
## G  without duplication GA

###1. coverage at each site
HRSVA_G.epimx.data.GA.3<-HRSVA_G.epimx.data.GA.2%>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

HRSVA_G.epimx.data.GA.3<-HRSVA_G.epimx.data.GA.3%>% 
  group_by(Frame.Start, Frame.Stop) %>%
  dplyr::summarise(count = n())

n<-n_distinct(HRSVA_G.epimx.data.GA$Accession)
HRSVA_G.epimx.data.GA.3<- HRSVA_G.epimx.data.GA.3%>% mutate(coverge = count/as.numeric(n))

### 2. get the sepecific epitope profile, including coverage, HLA list and average z score
### coverage for each epitope
HRSVA_G.epimx.data.GA.4<-HRSVA_G.epimx.data.GA.2 %>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

HRSVA_G.epimx.data.GA.4<-HRSVA_G.epimx.data.GA.4%>% 
  group_by(Frame.Start, AA.Sequence,Frame.Stop) %>%
  dplyr::summarise(count = n())

## calculate the coverage percentage
HRSVA_G.epimx.data.GA.4<- HRSVA_G.epimx.data.GA.4%>% mutate(coverge = count/as.numeric(n))


### 3.get the HLA list for each epitope
HRSVA_G.epimx.data.GA.5<-HRSVA_G.epimx.data.GA.2 %>% distinct(Frame.Start,AA.Sequence,Frame.Stop,HLA.Alleles)
HRSVA_G.epimx.data.GA.epitope_HLA_list<-HRSVA_G.epimx.data.GA.5 %>% 
  group_by(AA.Sequence) %>% 
  mutate(HLAs = paste0(HLA.Alleles, collapse = ",")) %>%
  select("Frame.Start","AA.Sequence","HLAs")

HRSVA_G.epimx.data.GA.epitope_HLA_list<-distinct(HRSVA_G.epimx.data.GA.epitope_HLA_list)

## 4.get the average Z score for each epitope
HRSVA_G.epimx.data.GA.6<-HRSVA_G.epimx.data.GA.2 %>% distinct(Frame.Start,AA.Sequence,Frame.Stop,HLA.Alleles,Hits,Z.score)

HRSVA_G.epimx.data.GA.6<-HRSVA_G.epimx.data.GA.6%>%
  group_by(Frame.Start,AA.Sequence,Frame.Stop,Hits)%>%
  summarise(sum_Zscore=sum(as.numeric(Z.score)))

HRSVA_G.epimx.data.GA.6<-HRSVA_G.epimx.data.GA.6%>%mutate(mean.score=sum_Zscore/9)

##merge all of the epitope profile together

HRSVA_G_GA_epitope_sum<-left_join(HRSVA_G.epimx.data.GA.4,HRSVA_G.epimx.data.GA.epitope_HLA_list,HRSVA_G.epimx.data.GA.3,
                                  by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))

HRSVA_G_GA_epitope_sum<-left_join(HRSVA_G_GA_epitope_sum,HRSVA_G.epimx.data.GA.6,
                                  by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))

## write these 445 epitope profile
write.csv(HRSVA_G_GA_epitope_sum,"iVAX/HRSVA_G_GA_epitope_sum0927.csv",row.names=FALSE)
##############################################################################################################################3
################################################################################################################################
## prepare the visulization for heatmap for ON genotype
## filter the epitpe site with > 5% coverge 

## filter the interest epitope
ON_epitope_site<-HRSVA_G.epimx.data.ON.3%>%filter(coverge>=0.05)
## from the epitope site to get the epitope profile
ON_epitope_new <-left_join(ON_epitope_site,HRSVA_G_ON_epitope_sum,by=c("Frame.Start" = "Frame.Start"))

## sort the dataframe by Frame start and then number of counts for each epitope

ON_epitope_new <-ON_epitope_new[order(ON_epitope_new$Frame.Start, -ON_epitope_new$count.y),]
## write epitope for hp visulization 54 sites, 370 epitope
write.csv(ON_epitope_new ,"iVAX/HRSVA_G_classII_ON_epitopehp_0927.csv",row.names=FALSE)
## prepare heatmap for visulization

## transform the different epitope at the same location with different states
ON_epitope_new$code<-ON_epitope_new$Frame.Start

ON_epitope_new$code<- make.unique(as.character(ON_epitope_new$code))
ON_epitope_new$code<-as.numeric(ON_epitope_new$code)%%1 *10+1 


## if the presence of epitope < 1%
ON_epitope_new<-ON_epitope_new %>% mutate(code = ifelse(coverge.y<=0.01, "x",as.character(ON_epitope_new$code)))
## transform to heatmap
## predominant site profile

ON_site<-ON_epitope_new%>%
  filter(code=="1")%>%
  select("Frame.Start","AA.Sequence")
ON_site$epitope1<-paste(ON_site$Frame.Start," ",ON_site$AA.Sequence)


HRSVA_G.epimx.data.ON.seq<-HRSVA_G.epimx.data.ON.2 %>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

## add epitope infor
HRSVA_G.epimx.data.ON.seq<-left_join(ON_epitope_new,HRSVA_G.epimx.data.ON.seq,
                                  by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))
## add epitope site infor
HRSVA_G.epimx.data.ON.seq<-left_join(HRSVA_G.epimx.data.ON.seq,ON_site,
                                  by = c("Frame.Start" = "Frame.Start"))
HRSVA_G.epimx.data.ON.hp<-HRSVA_G.epimx.data.ON.seq%>%
  ungroup()%>%
  select("Accession","code","epitope1")

## reshape, change the AA Sequence value to column name
HRSVA_G.epimx.data.ON.hp<-dcast(HRSVA_G.epimx.data.ON.hp, Accession ~ epitope1,value.var = "code") 
write.csv(HRSVA_G.epimx.data.ON.hp,"iVAX/HRSVA_G_ON_hp.csv",row.names=FALSE)

#######################################################################################################
## prepare the visulizatiGA for heatmap for GA genotype
## filter the epitpe site with > 5% coverge 

## filter the interest epitope
GA_epitope_site<-HRSVA_G.epimx.data.GA.3%>%filter(coverge>=0.05)
## from the epitope site to get the epitope profile
GA_epitope_new <-left_join(GA_epitope_site,HRSVA_G_GA_epitope_sum,by=c("Frame.Start" = "Frame.Start"))

## sort the dataframe by Frame start and then number of counts for each epitope

GA_epitope_new <-GA_epitope_new[order(GA_epitope_new$Frame.Start, -GA_epitope_new$count.y),]
## write epitope for hp visulizatiGA 54 sites, 370 epitope
write.csv(GA_epitope_new ,"iVAX/HRSVA_G_classII_GA_epitopehp_0927.csv",row.names=FALSE)
## prepare heatmap for visulizatiGA

## transform the different epitope at the same locatiGA with different states
GA_epitope_new$code<-GA_epitope_new$Frame.Start

GA_epitope_new$code<- make.unique(as.character(GA_epitope_new$code))
GA_epitope_new$code<-as.numeric(GA_epitope_new$code)%%1 *10+1 


## if the presence of epitope < 1%
GA_epitope_new<-GA_epitope_new %>% mutate(code = ifelse(coverge.y<=0.01, "x",as.character(GA_epitope_new$code)))
## transform to heatmap
## predominant site profile

GA_site<-GA_epitope_new%>%
  filter(code=="1")%>%
  select("Frame.Start","AA.Sequence")
GA_site$epitope1<-paste(GA_site$Frame.Start," ",GA_site$AA.Sequence)


HRSVA_G.epimx.data.GA.seq<-HRSVA_G.epimx.data.GA.2 %>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

## add epitope infor
HRSVA_G.epimx.data.GA.seq<-left_join(GA_epitope_new,HRSVA_G.epimx.data.GA.seq,
                                     by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))
## add epitope site infor
HRSVA_G.epimx.data.GA.seq<-left_join(HRSVA_G.epimx.data.GA.seq,GA_site,
                                     by = c("Frame.Start" = "Frame.Start"))
HRSVA_G.epimx.data.GA.hp<-HRSVA_G.epimx.data.GA.seq%>%
  ungroup()%>%
  select("Accession","code","epitope1")

## reshape, change the AA Sequence value to column name
HRSVA_G.epimx.data.GA.hp<-dcast(HRSVA_G.epimx.data.GA.hp, Accession ~ epitope1,value.var = "code") 
write.csv(HRSVA_G.epimx.data.GA.hp,"iVAX/HRSVA_G_GA_hp.csv",row.names=FALSE)
#################################################################################################################
## try to sum ON1 and GA together

#filter Z SCORE >= 1.64
HRSVA_G.epimx.data.3 <- HRSVA_G.epimx.data.2 %>% filter(Z.score >= 1.64)

#filter EPIBARS Hits >=4

HRSVA_G.epimx.data.4 <- HRSVA_G.epimx.data.3  %>% filter((as.numeric(Hits)) >= 4)

write.csv(HRSVA_G.epimx.data.4,"iVAX/HRSVA_G_epimx_filter.csv")

############################################################################################################################


###1. coverage at each site
HRSVA_G.epimx.data.5<-HRSVA_G.epimx.data.4%>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

HRSVA_G.epimx.data.5<-HRSVA_G.epimx.data.5%>% 
  group_by(Frame.Start, Frame.Stop) %>%
  dplyr::summarise(count = n())

n<-n_distinct(HRSVA_G.epimx.data.4$Accession)
HRSVA_G.epimx.data.5<- HRSVA_G.epimx.data.5%>% mutate(coverge = count/as.numeric(n))

### 2. get the sepecific epitope profile, including coverage, HLA list and average z score
### coverage for each epitope
HRSVA_G.epimx.data.6<-HRSVA_G.epimx.data.4 %>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

HRSVA_G.epimx.data.6<-HRSVA_G.epimx.data.6%>% 
  group_by(Frame.Start, AA.Sequence,Frame.Stop) %>%
  dplyr::summarise(count = n())

## calculate the coverage percentage
HRSVA_G.epimx.data.6<- HRSVA_G.epimx.data.6%>% mutate(coverge = count/as.numeric(n))


### 3.get the HLA list for each epitope
HRSVA_G.epimx.data.7<-HRSVA_G.epimx.data.4 %>% distinct(Frame.Start,AA.Sequence,Frame.Stop,HLA.Alleles)
HRSVA_G.epitope_HLA_list<-HRSVA_G.epimx.data.7 %>% 
  group_by(AA.Sequence) %>% 
  mutate(HLAs = paste0(HLA.Alleles, collapse = ",")) %>%
  select("Frame.Start","AA.Sequence","HLAs")

HRSVA_G.epitope_HLA_list<-distinct(HRSVA_G.epitope_HLA_list)

## 4.get the average Z score for each epitope
HRSVA_G.epimx.data.8<-HRSVA_G.epimx.data.4 %>% distinct(Frame.Start,AA.Sequence,Frame.Stop,HLA.Alleles,Hits,Z.score)

HRSVA_G.epimx.data.8<-HRSVA_G.epimx.data.8%>%
  group_by(Frame.Start,AA.Sequence,Frame.Stop,Hits)%>%
  summarise(sum_Zscore=sum(as.numeric(Z.score)))

HRSVA_G.epimx.data.8<-HRSVA_G.epimx.data.8%>%mutate(mean.score=sum_Zscore/9)

##merge all of the epitope profile together

HRSVA_G_epitope_sum<-left_join(HRSVA_G.epimx.data.6,HRSVA_G.epitope_HLA_list,HRSVA_G.epimx.data.5,
                                  by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))

HRSVA_G_epitope_sum<-left_join(HRSVA_G_epitope_sum,HRSVA_G.epimx.data.8,
                                  by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))

## write these 445 epitope profile
write.csv(HRSVA_G_epitope_sum,"iVAX/HRSVA_G_epitope_sum0927.csv",row.names=FALSE)
##############################################################################################################################3

## filter the epitpe site with > 5% coverge 

## filter the interest epitope
epitope_site<-HRSVA_G.epimx.data.5%>%filter(coverge>=0.05)
## from the epitope site to get the epitope profile
epitope_new <-left_join(epitope_site,HRSVA_G_epitope_sum,by=c("Frame.Start" = "Frame.Start"))

## sort the dataframe by Frame start and then number of counts for each epitope

epitope_new <-epitope_new[order(epitope_new$Frame.Start, -epitope_new$count.y),]
## write epitope for hp visulization 54 sites, 370 epitope
write.csv(epitope_new ,"iVAX/HRSVA_G_classII_epitopehp_0927.csv",row.names=FALSE)
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


HRSVA_G.epimx.data.seq<-HRSVA_G.epimx.data.2 %>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

## add epitope infor
HRSVA_G.epimx.data.seq<-left_join(epitope_new,HRSVA_G.epimx.data.seq,
                                     by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))
## add epitope site infor
HRSVA_G.epimx.data.seq<-left_join(HRSVA_G.epimx.data.seq,site,
                                     by = c("Frame.Start" = "Frame.Start"))
HRSVA_G.epimx.data.hp<-HRSVA_G.epimx.data.seq%>%
  ungroup()%>%
  select("Accession","code","epitope1")

## reshape, change the AA Sequence value to column name
HRSVA_G.epimx.data.hp<-dcast(HRSVA_G.epimx.data.hp, Accession ~ epitope1,value.var = "code") 
write.csv(HRSVA_G.epimx.data.hp,"iVAX/HRSVA_G_hp.csv",row.names=FALSE)



##############################################################################
## RSVA G class I analysis


RSVA_G_classI<-read.csv("iVAX/CLASS_I_epitope/CLASS1_PROTEIN_REPORT_HRSVA_G_CLEAN0819_PROTEIN_All_Proteins_BY_9.csv"
                        ,header = FALSE, stringsAsFactors = FALSE, sep = ",")


totalseq <- grep("HRSVA_G_CLEAN0819_PROTEIN ", RSVA_G_classI$V1, useBytes = TRUE)
endseq <- grep("Summarized Results", RSVA_G_classI$V1, useBytes = TRUE)
epimx.data_output <- ""

epimx.data <- RSVA_G_classI%>% slice((totalseq[1]):(endseq[1]-2))
seqname <- gsub("File: HRSVA_G_CLEAN0819_PROTEIN - Sequence: ", "", epimx.data$V1[1])
epimx.data <- epimx.data %>% slice(grep("Frame", V1, useBytes = T):nrow(epimx.data))
colnames(epimx.data) <- c("Frame.Start", "AA.Sequence", "Frame.Stop","Hydro-", "A0101", "A0201", "A0301",  
                          "A2402", "B0702", "B4403", "Hits")
epimx.data <- epimx.data[-c(1:2),]
epimx.data.reshape <- epimx.data %>% mutate(Seq.name = seqname) %>% 
  dplyr::select(Seq.name, Frame.Start:Hits) %>%
  gather(HLA.Alleles, Z.score, "A0101":"B4403", factor_key = TRUE)

epimx.data_output <- rbind(epimx.data_output, epimx.data.reshape)



for (i in 1:length(totalseq)){
  epimx.data <- RSVA_G_classI%>% slice((totalseq[i]):(endseq[i]-2))
  seqname <- gsub("File: HRSVA_G_CLEAN0819_PROTEIN - Sequence: ", "", epimx.data$V1[1])
  epimx.data <- epimx.data %>% slice(grep("Frame", V1, useBytes = T):nrow(epimx.data))
  colnames(epimx.data) <- c("Frame.Start", "AA.Sequence", "Frame.Stop","Hydro-", "A0101", "A0201", "A0301",  
                            "A2402", "B0702", "B4403", "Hits")
  epimx.data <- epimx.data[-c(1:2),]
  epimx.data.reshape <- epimx.data %>% mutate(Seq.name = seqname) %>% 
    dplyr::select(Seq.name, Frame.Start:Hits) %>%
    gather(HLA.Alleles, Z.score, "A0101":"B4403", factor_key = TRUE)
  
  epimx.data_output <- rbind(epimx.data_output, epimx.data.reshape)
  
}
RSVA_G_epimax<-epimx.data_output
######################################################################################

#MORE data restructuring 

RSVA_G_epimx_data.2 <- RSVA_G_epimax %>% filter(!Seq.name == "")
## remove duplicate from HRSVB_F.epimx.data.2

RSVA_G.epimx.data.2 <-distinct(RSVA_G_epimx_data.2)
x <- data.frame(unique(RSVA_G.epimx.data.2$Seq.name)) ## get all sequence name  

write.csv(RSVA_G.epimx.data.2,"iVAX/RSVA_G_classI_epimx1120.csv")  

#####################################################################################  

RSVA_G.epimx.data.2 <- read.csv("iVAX/RSVA_G/RSVA_G_classI/RSVA_G_classI_epimx1120.csv")

#filter Z SCORE >= 1.64
RSVA_G.epimx.data.3 <- RSVA_G.epimx.data.2 %>% filter(Z.score >2.326)


#filter EPIBARS Hits >=4

RSVA_G.epimx.data.3 <- RSVA_G.epimx.data.3 %>% filter((as.numeric(Hits)) >= 1)

write.csv(RSVA_G.epimx.data.3 ,"iVAX/RSVA_G_classI_epimax_filter.csv")

#############################################################################


RSVA_G.epimx.data.3<-RSVA_G.epimx.data.3 %>% separate(Seq.name, into = c("Accession", "Subtype", "Country", "date", "genotype"), sep = "-", extra = "merge")
seq <- scan("clean_data/RSVA_G_accession.txt", what="", sep="\n")



RSVA_G.epimx.data.3 <- RSVA_G.epimx.data.3 %>% 
  mutate(Type = ifelse(RSVA_G.epimx.data.3$Accession %in% seq, "T", "F"))
RSVA_G.epimx.data.3<-RSVA_G.epimx.data.3%>% filter(RSVA_G.epimx.data.3$Type =="T")
##############################################################################################################
## Get the epitope profile for each eptiope frame start site
## with 6 nt deletion at the upstream

###1. coverage at each site
RSVA_G.epimx.data.4<-RSVA_G.epimx.data.3%>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

RSVA_G.epimx.data.4<-RSVA_G.epimx.data.4%>% 
  group_by(Frame.Start, Frame.Stop) %>%
  dplyr::summarise(count = n())

n<-n_distinct(RSVA_G.epimx.data.3$Accession)
RSVA_G.epimx.data.4<- RSVA_G.epimx.data.4%>% mutate(coverge = count/as.numeric(n))

### 2. get the sepecific epitope profile, including coverage, HLA list and average z score
### coverage for each epitope
RSVA_G.epimx.data.5<-RSVA_G.epimx.data.3 %>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

RSVA_G.epimx.data.5<-RSVA_G.epimx.data.5%>% 
  group_by(Frame.Start, AA.Sequence,Frame.Stop) %>%
  dplyr::summarise(count = n())

## calculate the coverage percentage
RSVA_G.epimx.data.5<- RSVA_G.epimx.data.5%>% mutate(coverge = count/as.numeric(n))


### 3.get the HLA list for each epitope
RSVA_G.epimx.data.6<-RSVA_G.epimx.data.3 %>% distinct(Frame.Start,AA.Sequence,Frame.Stop,HLA.Alleles)
RSVA_G.epimx_HLA_list<-RSVA_G.epimx.data.6 %>% 
  group_by(AA.Sequence) %>% 
  mutate(HLAs = paste0(HLA.Alleles, collapse = ",")) %>%
  dplyr::select("Frame.Start","AA.Sequence","HLAs")

RSVA_G.epimx_HLA_list<-distinct(RSVA_G.epimx_HLA_list)

## 4.get the average Z score for each epitope
RSVA_G.epimx.data.7<-RSVA_G.epimx.data.3 %>% distinct(Frame.Start,AA.Sequence,Frame.Stop,HLA.Alleles,Hits,Z.score)

RSVA_G.epimx.data.7<-RSVA_G.epimx.data.7%>%
  group_by(Frame.Start,AA.Sequence,Frame.Stop,Hits)%>%
  summarise(sum_Zscore=sum(as.numeric(Z.score)))

RSVA_G.epimx.data.7<-RSVA_G.epimx.data.7%>%mutate(mean.score=sum_Zscore/9)

##merge all of the epitope profile together

RSVA_G_epitope_sum<-left_join(RSVA_G.epimx.data.5,RSVA_G.epimx_HLA_list,RSVA_G.epimx.data.6,by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))

RSVA_G_epitope_sum<-left_join(RSVA_G_epitope_sum,RSVA_G.epimx.data.7,
                               by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))

## write these 445 epitope profile
write.csv(RSVA_G_epitope_sum,"iVAX/RSVA_G/RSVA_G_classI/RSVA_G_classI_epitope_sum0810.csv",row.names=FALSE)

##############################################################################################################################3
## filter the epitpe site with > 5% coverge 

## filter the interest epitope
epitope_site<-RSVA_G.epimx.data.4%>%filter(coverge>=0.05)
## from the epitope site to get the epitope profile
epitope_new <-left_join(epitope_site,RSVA_G_epitope_sum,by=c("Frame.Start" = "Frame.Start"))

## sort the dataframe by Frame start and then number of counts for each epitope

epitope_new <-epitope_new[order(epitope_new$Frame.Start, -epitope_new$count.y),]
## write epitope for hp visulization 54 sites, 370 epitope
write.csv(epitope_new ,"iVAX/RSVA_G/RSVA_G_classI/RSVA_G_classI_epitope_hplist_0813.csv",row.names=FALSE)
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


RSVA_G.epimx.data.seq<-RSVA_G.epimx.data.3 %>% distinct(Accession, Frame.Start,AA.Sequence,Frame.Stop)

## add epitope infor
RSVA_G.epimx.data.seq<-left_join(epitope_new,RSVA_G.epimx.data.seq,
                                  by = c("Frame.Start" = "Frame.Start", "AA.Sequence" = "AA.Sequence"))
## add epitope site infor
RSVA_G.epimx.data.seq<-left_join(RSVA_G.epimx.data.seq,site,
                                  by = c("Frame.Start" = "Frame.Start"))

RSVA_G.epimx.data.hp<-RSVA_G.epimx.data.seq%>%
  ungroup()%>%
  dplyr::select("Accession","code","epitope1")

## reshape, change the AA Sequence value to column name
RSVA_G.epimx.data.hp<-dcast(RSVA_G.epimx.data.hp, Accession ~ epitope1,value.var = "code") 
write.csv(RSVA_G.epimx.data.hp,"iVAX/RSVA_G/RSVA_G_classI/RSVA_G_classI_hp_0810.csv",row.names=FALSE)
#################################################################################################







