## calculate the immunogic potential at different site
## load single strain from RSVA_F_class I as a test
df<-read.csv("iVAX/RSVB_F/RSVB_F_classI_sum/HRSVB_F_classI_epimx1020.csv")
df1<-df%>%filter(Accession=="MG642043")
df1<-df[grepl("^MG642043", df$Seq.name), ]
df1$site=df1$Frame.Start
## assgin the Frame start as epitope postion and assign for each site
total_df<-df1
for (i in 1:8){
  next_df<-df1
  next_df$site=as.numeric(df1$site)+i
  total_df<-rbind(total_df,next_df)
  
}

## sum the Z score for each site
epitope_score<-total_df%>%group_by(site)%>%
  summarise(epitope_score=sum(Z.score))
epitope_score<-epitope_score%>%select(epitope_score)
#epitope_score<-epitope_score%>% filter((site>=26 && site <=99) || (site >=147 && site <=506))
## write out a dataframe
write.csv(epitope_score,"iVAX/RSVB_F_classI_epitope_site.txt",row.names=FALSE,col.names = FALSE)
### wrap into function
epitope_site<-function(epimatrix,seq){
  df<-read.csv(epimatrix)
  df1<-df%>%filter(Accession==seq)
  df1$site=df1$Frame.Start
  ## assgin the Frame start as epitope postion and assign for each site
  total_df<-df1
  for (i in 1:8){
    next_df<-df1
    next_df$site=as.numeric(df1$site)+i
    total_df<-rbind(total_df,next_df)
  }
  epitope_score<-total_df%>%group_by(site)%>%
    summarise(epitope_score=sum(Z.score))
  epitope_score<-epitope_score%>%select(epitope_score)
  return(epitope_score)
}

RSVA_F_classI_score<-epitope_site("iVAX/RSVA_F/RSVA_F_classI_analysis/RSVA_F_classI_epimax_clean.csv","KU316120")
RSVA_F_classII_score<-epitope_site('iVAX/RSVA_F/RSVA_F_classII_analysis/HRSVA_F_classII_epimx0831_clean.csv',"KU316120")
write.csv(RSVA_F_classII_score,"iVAX/RSVA_F_classII_epitope_site.txt",row.names=FALSE,col.names = FALSE)
write.csv(RSVA_F_classI_score,"iVAX/RSVA_F_classI_epitope_site.txt",row.names=FALSE,col.names = FALSE)
RSVB_F_classII_score<-epitope_site("iVAX/HRSVB_F_classII_csv/RSVB_F_classsII_analysis/HRSVB_F_classII_epimx0815.csv","MG642043")
RSVB_F_classI_score<-epitope_site("iVAX/RSVA_F/RSVA_F_classI_analysis/RSVA_F_classI_epimax_clean.csv","MG642043")