## RSV F protein mds

## parse the EPICC unique output to generat final data, remove is a list with sequence name need to be removed
epicc_unique<-function(classI_vax,classI_strain,classII_vax,classII_strain,remove){
  df1_1<-read.table(classI_vax,sep=",", header=TRUE,check.names=FALSE, stringsAsFacto=F, fill=TRUE,row.names = 1)
  df1_2<-read.table(classI_strain,sep=",", header=TRUE,check.names=FALSE, stringsAsFacto=F, fill=TRUE,row.names = 1)
  df2_1<-read.table(classII_vax,sep=",", header=TRUE,check.names=FALSE, stringsAsFacto=F, fill=TRUE,row.names = 1)
  df2_2<-read.table(classII_strain,sep=",", header=TRUE,check.names=FALSE, stringsAsFacto=F, fill=TRUE,row.names = 1)
  data=df1_2+df1_2+df2_1+df2_2
  
  data <- select(data, -remove)
  data<-data[!(row.names(data) %in% remove), ]
  return(data)
}

## peform mds and k-means return the table fill with strain name and, coordinates and group
## further calculate T-cell immun distance and genetic distance
mds<-function(data,k=n,x,tmrca,ham,tmrca2,label_file){
  mds<-data%>%
    #dist()%>%
    cmdscale(k=2) %>%
    as_tibble()
  
  colnames(mds)<-c("Dim.1","Dim.2")
  ## k-means
  k<-kmeans(mds,centers=x,nstart = 123)
  clust <- k$cluster %>%
    as.factor()
  mds <- mds %>%
    mutate(groups = clust)
  mds$strain <-rownames(data)
  
  ## calculte T-cell epitope immune distance
  crossT_dist<-data%>%select(one_of(tmrca))
  crossT_dist$strain<-rownames(crossT_dist)
  crossT_dist<-crossT_dist %>% separate(strain, into = c("Accession", "Subtype", "year", "country", "date"), sep = "-", extra = "merge")
  crossT_dist$strain2<-rownames(crossT_dist)
  mds<-left_join(mds,crossT_dist,by=c("strain"="strain2"))
  
  ## add genetic hamming distance
  ham_dist<- read.table(ham, sep=",", header=TRUE,check.names=FALSE, stringsAsFacto=F, fill=TRUE,row.names = 1)
  
  ham_dist<-ham_dist%>%select(one_of(tmrca2))
  #names(ham_dist)[names(ham_dist) == "tmrca"] <- "ham_dist"
  ham_dist$strain<-rownames(ham_dist)
  ham_dist<-ham_dist %>% separate(strain, into = c("Accession", "Subtype", "year", "country", "date"), sep = "-", extra = "merge")
  mds<-left_join(mds,ham_dist,by=c("Accession"="Accession","Subtype"="Subtype","year"="year"))
  
  ## add label column from external dataframe
  label<-read.csv(label_file)
  mds<-left_join(mds,label,by=c("strain"="id"))
  
  
  return(mds)
  
}

## function to find the optimize number of cluster

plot_mds<-function(mds,label_list,color){
  p<-ggplot(mds, aes(x=Dim.1, y=Dim.2, color=groups,label=ifelse(strain %in% label_list, label, ""))) + 
    geom_point(size=1.5, alpha=.75)+
    xlab("Dimension 1")+
    ylab("Dimension 2")+
    scale_color_brewer(palette = color,name="T-cell epitope immune clusters")+
    theme_light()+
    geom_encircle(expand=0,linetype=2)+
    geom_label_repel(fill = "white", xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),color="black",size=2)

  
  return(p)
  
}


plot_mds_G<-function(mds,label_list,color){
  p<-ggplot(mds, aes(x=Dim.1, y=Dim.2, color=groups,shape=duplication,label=ifelse(strain %in% label_list, label, ""))) + 
    geom_point(size=1.5, alpha=.75)+
    xlab("Dimension 1")+
    ylab("Dimension 2")+
    scale_color_brewer(palette = color,name="T-cell epitope immune clusters")+
    scale_shape_manual(values = c(19,5))+
    theme_light()+
    geom_encircle(aes(group=groups),expand=0,linetype=2)+
    geom_label_repel(fill = "white", xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),color="black",size=2)
  
  
  return(p)
  
}
## function to plot the tree


plot_tree<-function(tree,mds,taxa_file,color){
  t<-read.tree(file = tree)
  p1<-ggtree(t) +
    geom_treescale(x=0, y=300,linesize= 1,,width=0.005)
  
  taxa<-read.csv(taxa_file)
  taxa_group<-left_join(taxa,mds,by = c("accession" = "Accession"))
  taxa_group$groups<-as.factor(taxa_group$groups)
  taxa_group<-taxa_group%>% select(taxa,groups)
  p2<-p1 %<+% taxa_group+ 
    geom_tippoint(aes(color=groups), size=1.5, alpha=.75)+
    scale_color_brewer(palette = color,name="T-cell epitope immune clusters")
  
  return(p2)
}


plot_tree_G<-function(tree,mds,taxa_file,color){
  t<-read.tree(file = tree)
  p1<-ggtree(t) +
    geom_treescale(x=0, y=300,linesize= 1,,width=0.005)
  
  taxa<-read.csv(taxa_file)
  taxa_group<-left_join(taxa,mds,by = c("accession" = "Accession"))
  taxa_group$groups<-as.factor(taxa_group$groups)
  
  taxa_group<-taxa_group%>% select(taxa,groups,duplication.x)
  p2<-p1 %<+% taxa_group+ 
    geom_tippoint(aes(color=groups,shape=duplication.x), size=1.5, alpha=.75)+
    scale_shape_manual(values = c(19,5))+
    scale_color_brewer(palette = color,name="T-cell epitope immune clusters")
  
  return(p2)
}

plot_immunedist<-function(mds,tmrca,color){

  p<-ggplot(mds, aes(x=as.numeric(year), y=get(tmrca),color=groups)) + 
    geom_point(size=1.5, alpha=.75)+
    scale_colour_brewer(palette = color,name="T-cell epitope immune clusters")+
    ylab("T-epitope immune distance")+
    xlab("Isolated Year")+
    theme_light()
  
}

plot_immunedist_G<-function(mds,tmrca,color){
  
  p<-ggplot(mds, aes(x=as.numeric(year), y=get(tmrca),color=groups,shape=duplication)) + 
    geom_point(size=1.5, alpha=.75)+
    scale_colour_brewer(palette = color,name="T-cell epitope immune clusters")+
    scale_shape_manual(values = c(19,5))+
    ylab("T-epitope immune distance")+
    xlab("Isolated Year")+
    theme_light()
  
}


plot_hamdist<-function(mds,tmrca2,color){
  ggplot(mds, aes(x=as.numeric(year), y=get(tmrca2),color=groups)) + 
    geom_point(size=1.5, alpha=.75)+
    scale_colour_brewer(palette = color,name="T-cell immuno-clusters")+
    ylab("Genetic Hamming distance")+
    xlab("Isolated Year")+
    theme_light()
  
}

plot_hamdist_G<-function(mds,tmrca2,color){
  ggplot(mds, aes(x=as.numeric(year), y=get(tmrca2),color=groups,shape=duplication)) + 
    geom_point(size=1.5, alpha=.75)+
    scale_colour_brewer(palette = color,name="T-cell immuno-clusters")+
    scale_shape_manual(values = c(19,5))+
    ylab("Genetic Hamming distance")+
    xlab("Isolated Year")+
    theme_light()
  
}
#######################################################################################################################
## test with RSVA_F
RSVA_F_remove<-c("KT992094-VACCINE_D46_D53","RSVB_PDA_ANCESTRAL","AF035006-VACCINE_RA2CP","AF013255-VACCINE_CP52")
RSVA_F_unique<-epicc_unique("EpiCC/RSV_EpiCC_Data/RSVA_F_classI_vac_unique.csv",
                            "EpiCC/RSV_EpiCC_Data/RSVA_F_classI_strain_unique.csv",
                            "EpiCC/RSV_EpiCC_Data/RSVA_F_classII_vac_unique.csv",
                            "EpiCC/RSV_EpiCC_Data/RSVA_F_classII_strain_unique.csv",RSVA_F_remove)

RSVA_F_mds<-mds(RSVA_F_unique,k=2,3,"RSVA_PDA_ANCESTRAL","EpiCC/RSVA_F/RSVA_F_epicc_subsample_genetic_hamming.csv","JX198138-A-1961-Australia-7/16/1961-GA1","EpiCC/RSVA_F/RSVA_F_epicc_label.csv")
RSVA_F_label<-c("U63644-VACCINE_CPTS-248_404","AF013255-VACCINE_CP52","RSVA_PDA_ANCESTRAL")
p_RSVA_F_mds<-plot_mds(RSVA_F_mds,RSVA_F_label,"RdPu")
p_RSVA_F_tree<-plot_tree("EpiCC/RSVA_F/RSVA_F_ML/RSVA_F_besttree_midpoint.nwk",RSVA_F_mds,"EpiCC/RSVA_F/RSVA_F_epicc_taxa.csv","RdPu")
p_RSVA_F_immudist<-plot_immunedist(RSVA_F_mds,"RSVA_PDA_ANCESTRAL", "RdPu")
p_RSVA_F_hamdist<-plot_hamdist(RSVA_F_mds,"JX198138-A-1961-Australia-7/16/1961-GA1","RdPu")


## test with RSVB_F
RSVB_F_remove<-c("KT992094-VACCINE_D46_D53","RSVA_PDA_ANCESTRAL","AF035006-VACCINE_RA2CP","U63644-VACCINE_CPTS-248_404")
RSVB_F_unique<-epicc_unique("EpiCC/RSV_EpiCC_Data/RSVB_F_classI_vac_unique.csv",
                            "EpiCC/RSV_EpiCC_Data/RSVB_F_classI_strain_unique.csv",
                            "EpiCC/RSV_EpiCC_Data/RSVB_F_classII_vac_unique.csv",
                            "EpiCC/RSV_EpiCC_Data/RSVB_F_classII_strain_unique.csv",RSVB_F_remove)

RSVB_F_mds<-mds(RSVB_F_unique,k=2,2,"RSVB_PDA_ANCESTRAL","EpiCC/RSVB_F/RSVB_F_epicc_subsample_genetic_hamming.csv","RSVB_pda_ancestral","EpiCC/RSVB_F/RSVB_F_epicc_label.csv")
RSVB_F_label<-c("AF013255-VACCINE_CP52","RSVB_PDA_ANCESTRAL")
p_RSVB_F_mds<-plot_mds(RSVB_F_mds,RSVB_F_label,"Blues")
p_RSVB_F_tree<-plot_tree("EpiCC/RSVB_F/RSVB_F_subsample_epicc.nwk",RSVB_F_mds,"EpiCC/RSVB_F/RSVB_F_epicc_subsample_taxa.csv","Blues")
p_RSVB_F_immudist<-plot_immunedist(RSVB_F_mds,"RSVB_PDA_ANCESTRAL","Blues")
p_RSVB_F_hamdist<-plot_hamdist(RSVB_F_mds,"RSVB_pda_ancestral","Blues")

######################3
## test with G protein
RSVA_G_remove<-c("KT992094-VACCINE_D46_D53","AF013255-VACCINE_CP52")
RSVA_G_unique<-epicc_unique("EpiCC/RSV_EpiCC_Data/RSVA_G_classI_vac_unique.csv",
                            "EpiCC/RSV_EpiCC_Data/RSVA_G_classI_strain_unique.csv",
                            "EpiCC/RSV_EpiCC_Data/RSVA_G_classII_vac_unique.csv",
                            "EpiCC/RSV_EpiCC_Data/RSVA_G_classII_strain_unique.csv",RSVA_G_remove)
RSVA_G_mds<-mds(RSVA_G_unique,k=2,3,"JX198138-A-1961-AUSTRALIA-7_16_1961-GA1","EpiCC/RSVA_G/RSVA_G_epicc_subsample_genetic_hamming.csv","JX198138-A-1961-Australia-7/16/1961-GA1","EpiCC/RSVA_G/RSVA_G_epicc_label.csv")
RSVA_G_label<-c("U63644-VACCINE_CPTS-248_404","JX198138-A-1961-AUSTRALIA-7_16_1961-GA1")
p_RSVA_G_mds<-plot_mds_G(RSVA_G_mds,RSVA_G_label,"OrRd")
p_RSVA_G_tree<-plot_tree_G("EpiCC/RSVA_G/RSVA_G_ML/RSVA_G_besttree_midpoint.nwk",RSVA_G_mds,"EpiCC/RSVA_G/RSVA_G_ML/RSVA_G_epicc_taxa.csv","OrRd")
p_RSVA_G_immudist<-plot_immunedist_G(RSVA_G_mds,"JX198138-A-1961-AUSTRALIA-7_16_1961-GA1","OrRd")
p_RSVA_G_hamdist<-plot_hamdist_G(RSVA_G_mds,"JX198138-A-1961-Australia-7/16/1961-GA1","OrRd")

##RSVB_G
RSVB_G_remove<-c("KT992094-VACCINE_D46_D53","U63644-VACCINE_CPTS-248_404")
RSVB_G_unique<-epicc_unique("EpiCC/RSV_EpiCC_Data/RSVB_G_classI_vac_unique.csv",
                            "EpiCC/RSV_EpiCC_Data/RSVB_G_classI_strain_unique.csv",
                            "EpiCC/RSV_EpiCC_Data/RSVB_G_classII_vac_unique.csv",
                            "EpiCC/RSV_EpiCC_Data/RSVB_G_classII_strain_unique.csv",RSVB_G_remove)
RSVB_G_mds<-mds(RSVB_G_unique,k=2,2,"KU316116-B-1977-USA-7_16_77-GB3","EpiCC/RSVB_G/RSVB_G_epicc_subsample_genetic_hamming.csv","KU316116-B-1977-USA-7/16/77-GB3","EpiCC/RSVB_G/RSVB_G_epicc_label.csv")
RSVB_G_label<-c("AF013255-VACCINE_CP52","KU316116-B-1977-USA-7_16_77-GB3")
p_RSVB_G_mds<-plot_mds_G(RSVB_G_mds,RSVB_G_label,"Greens")
p_RSVB_G_tree<-plot_tree_G("EpiCC/RSVB_G/RSVB_G_subsample_epicc.nwk",RSVB_G_mds,"EpiCC/RSVB_G/RSVB_G_epicc_subsample_taxa.csv","Greens")
p_RSVB_G_immudist<-plot_immunedist_G(RSVB_G_mds,"KU316116-B-1977-USA-7_16_77-GB3","Greens")
p_RSVB_G_hamdist<-plot_hamdist_G(RSVB_G_mds,"KU316116-B-1977-USA-7/16/77-GB3","Greens")
#####################################################

ggarrange(p_RSVA_F_mds,p_RSVA_F_tree,p_RSVA_F_immudist,p_RSVA_F_hamdist,
          p_RSVB_F_mds,p_RSVB_F_tree,p_RSVB_F_immudist,p_RSVB_F_hamdist,
          p_RSVA_G_mds,p_RSVA_G_tree,p_RSVA_G_immudist,p_RSVA_G_hamdist,
          p_RSVB_G_mds,p_RSVB_G_tree,p_RSVB_G_immudist,p_RSVB_G_hamdist,
          #labels = c("A", "B", "C","D","E","F","G","H"),
          ncol = 4, nrow = 4,
          widths = c(2,1.5,1.5,1.5),
          common.legend = TRUE)

