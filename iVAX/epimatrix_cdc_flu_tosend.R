if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('tidyr')) install.packages('tidyr'); library('tidyr')
if (!require('reshape2')) install.packages('reshape2'); library('reshape2')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('varhandle')) install.packages('varhandle'); library('varhandle')
if (!require('RColorBrewer')) install.packages('RColorBrewer'); library('RColorBrewer')
if (!require('gridExtra')) install.packages('gridExtra'); library('gridExtra')

#loop through directory and read files in r
#"^" same starting pattern
filenames <- list.files(path = "../Result/", pattern = "H3VAX.*csv$")
numfiles <- length(filenames) 

for (i in 1:numfiles){  
    name <- gsub(".csv","",filenames[i]) #optional: change dataf rame name
    name <- gsub("CLASS2_PROTEIN_REPORT_H3VAXSTRAINS_NH_SH_164SEQ_APR22_V3_","emx_class2_H3Vac",name) 
    name <- gsub("All_Proteins","",name) 
    name <- gsub("_",".",name)  
    assign(name, read.csv(paste("../Result/", filenames[i], sep = ""), 
                          header = FALSE, stringsAsFactors = FALSE, sep = ","))
}

#to create name list and store
dataname <- NULL
for (i in 1:numfiles){
    name[i] <- gsub(".csv","",filenames[i])
    name[i] <- gsub("CLASS2_PROTEIN_REPORT_H3VAXSTRAINS_NH_SH_164SEQ_APR22_V3_","emx_class2_H3Vac",name[i])
    name[i] <- gsub("All_Proteins","",name[i])
    name[i] <- gsub("_",".",name[i])
    dataname <- rbind(dataname, name[i])
}


#preclean-data function
#to extract each sequence EMX (EpiMatrix) info
#data reshaping
totalseq <- grep("File: H3VAXSTRAINS_NH_SH_164SEQ_APR22_V3", emx.class2.H3Vac$V1, useBytes = TRUE)
endseq <- grep("Summarized Results", emx.class2.H3Vac$V1, useBytes = TRUE)

cdcH3flu.epimx.data <- ""

for (i in 1:length(totalseq)) {
    epimx.data <- emx.class2.H3Vac %>% slice((totalseq[i]):(endseq[i]-2))
    seqname <- gsub("File: H3VAXSTRAINS_NH_SH_164SEQ_APR22_V3 - Sequence: ", "", epimx.data$V1[1])
    epimx.data <- epimx.data %>% slice(grep("Frame", V1, useBytes = T):nrow(epimx.data))
    colnames(epimx.data) <- c("Frame.Start", "AA.Sequence", "Frame.Stop", "DRB1*0101", "DRB1*0301", "DRB1*0401",  
                              "DRB1*0701", "DRB1*0801", "DRB1*0901", "DRB1*1101", "DRB1*1301", "DRB1*1501",  
                              "Hits", "NA")
    epimx.data <- epimx.data[-c(1:2),]
    epimx.data.reshape <- epimx.data %>% mutate(Seq.name = seqname) %>% 
        select(Seq.name, Frame.Start:Hits) %>%
        gather(HLA.Alleles, Z.score, "DRB1*0101":"DRB1*1501", factor_key = TRUE)
    cdcH3flu.epimx.data <- rbind(cdcH3flu.epimx.data, epimx.data.reshape)
}

#MORE data restructuring

cdcH3flu.epimx.data.2 <- cdcH3flu.epimx.data %>% filter(!Seq.name == "")

#checkpoint
x <- data.frame(unique(cdcH3flu.epimx.data.2$Seq.name))
x.2 <- x %>% separate(unique.cdcH3flu.epimx.data.2.Seq.name., into = c("Uniq.ID", "Protein", "Strain", "VacRegion", "YearRange", "Subtype"), sep = "-") %>%
    separate(Strain, into = c("Flutype", "Area", "Strain.ID", "Year"), sep = "_")


#INCLUDE ALL epitopes

cdcH3flu.epimx.data.3 <- cdcH3flu.epimx.data.2 %>%
    separate(Seq.name, into = c("Uniq.ID", "Protein", "Strain", "VacRegion", "YearRange", "Subtype"), sep = "-") %>% mutate(Strain.Name = Strain) %>%
    separate(Strain, into = c("Flutype", "Area", "Strain.ID", "Year"), sep = "_") %>%
    select(Uniq.ID, Strain.Name, VacRegion, YearRange, Area, Year, Frame.Start:Z.score)


#to get total # of epitope for each seq
seq.by.epitope <- aggregate(AA.Sequence ~ Uniq.ID, data = cdcH3flu.epimx.data.3, FUN = length)
seq.by.epitope <- seq.by.epitope %>% mutate(total.epitopes.per.alleles = AA.Sequence/9) %>% arrange(desc(total.epitopes.per.alleles))

#to count epitope by seq
epitope.by.seq <- aggregate(Uniq.ID ~ AA.Sequence, data = cdcH3flu.epimx.data.3, FUN = length)
epitope.by.seq <- epitope.by.seq %>% mutate(count.per.alleles = Uniq.ID/9) %>% arrange(desc(count.per.alleles))

#TABULATE EPITOPE ACCORDING TO VACCINE REGION AND YEAR
cdcH3flu.epimx.data.3a <- cdcH3flu.epimx.data.3 %>% spread(HLA.Alleles, Z.score) %>% mutate(ROY = YearRange) %>%
  separate(ROY, into = c("NH_range", "SH_range"), sep = "__") %>% separate(NH_range, into = c("NH_start", "NH_end"), sep = "_") %>%
  separate(SH_range, into = c("SH_start", "SH_end"), sep = "_") %>% mutate(VacLabel = paste0(VacRegion,"/",NH_start,"/",SH_start, "/", Strain.Name)) %>%
  arrange(VacLabel) %>% mutate(AA.seq.label = paste0(Frame.Start,"/",AA.Sequence))

a <- as.data.frame(table(cdcH3flu.epimx.data.3a$AA.seq.label, cdcH3flu.epimx.data.3a$VacLabel))
a <- a %>% arrange(Var1) %>% spread(Var2, Freq)

#to calculate sum of each row
a$row_sum <- apply(a[,c(2:36)], 1, sum)
a <- a %>% arrange(desc(row_sum))

#export output file for visualization
write.csv(a, "H3_vac_all_epitopes.csv")

#--------------------------------------------------------------------------------------------------------------------------------

#filter Z SCORE >= 1.64
cdcH3flu.epimx.data.4 <- cdcH3flu.epimx.data.3 %>% filter(Z.score >= 1.64)
cdcH3flu.epimx.data.4$Z.score <- as.numeric(cdcH3flu.epimx.data.4$Z.score)

#to get total # of epitope for each seq
seq.by.epitope.filter <- aggregate(AA.Sequence ~ Uniq.ID, data = cdcH3flu.epimx.data.4, FUN = length)
seq.by.epitope.filter <- seq.by.epitope.filter %>% arrange(desc(AA.Sequence))

#to count epitope by seq
epitope.by.seq.filter <- aggregate(Uniq.ID ~ AA.Sequence + VacRegion, data = cdcH3flu.epimx.data.4, FUN = length)
epitope.by.seq.filter <- epitope.by.seq.filter %>% arrange(desc(Uniq.ID)) 


#filter EPIBARS
cdcH3flu.epimx.data.5 <- ""
# for (i in 1:length(epitope.by.seq.filter$AA.Sequence)) {
#     temp <- cdcH3flu.epimx.data.4 %>% filter(AA.Sequence == epitope.by.seq.filter$AA.Sequence[i])
#     cdcH3flu.epimx.data.5 <- rbind(cdcH3flu.epimx.data.5, temp)
# }

cdcH3flu.epimx.data.4$Hits <- as.numeric(cdcH3flu.epimx.data.4$Hits)

cdcH3flu.epimx.data.5 <- cdcH3flu.epimx.data.4 %>% filter(Hits >= 4)

# #to get total # of epitope for each seq
# seq.by.epitope.filter.epibar <- aggregate(AA.Sequence ~ Uniq.ID, data = cdcH3flu.epimx.data.5, FUN = length)
# seq.by.epitope.filter.epibar <- seq.by.epitope.filter.epibar %>% arrange(desc(AA.Sequence))

#to count epitope by seq
epitope.by.seq.filter.epibar <- aggregate(Uniq.ID ~ AA.Sequence + VacRegion, data = cdcH3flu.epimx.data.5, FUN = length)
epitope.by.seq.filter.epibar <- epitope.by.seq.filter.epibar %>% arrange(desc(Uniq.ID)) 

#TABULATE EPITOPE ACCORDING TO VACCINE REGION AND YEAR
cdcH3flu.epimx.data.5a <- cdcH3flu.epimx.data.5 %>% spread(HLA.Alleles, Z.score) %>% mutate(ROY = YearRange) %>%
  separate(ROY, into = c("NH_range", "SH_range"), sep = "__") %>% separate(NH_range, into = c("NH_start", "NH_end"), sep = "_") %>%
  separate(SH_range, into = c("SH_start", "SH_end"), sep = "_") %>% mutate(VacLabel = paste0(VacRegion,"/",NH_start,"/",SH_start, "/", Strain.Name)) %>%
  arrange(VacLabel) %>% mutate(AA.seq.label = paste0(Frame.Start,"/",AA.Sequence))

b <- as.data.frame(table(cdcH3flu.epimx.data.5a$AA.seq.label, cdcH3flu.epimx.data.5a$VacLabel))
b <- b %>% arrange(Var1) %>% spread(Var2, Freq)

#to calculate sum of each row
b$row_sum <- apply(b[,c(2:36)], 1, sum)
b <- b %>% arrange(desc(row_sum))

#export output file 
write.csv(b, "H3_vac_epibars_epitopes_8May2020.csv")


