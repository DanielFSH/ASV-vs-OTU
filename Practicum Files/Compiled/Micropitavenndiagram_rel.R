library(VennDiagram)
library(tidyverse)

asvs<- read.csv('relative-asv.csv')
otus<-read.csv('relative-otu.csv')
micropita20<-read.csv('Micropita20_rel.csv')
micropita50<-read.csv('Micropita50_rel.csv')
micropita100<-read.csv('Micropita100_rel.csv')
Microtable<-as.data.frame(asvs$gid_16s)


#20
for (i in 1:length(asvs$allergy)) {
  a=asvs$gid_16s%in%micropita20$ï..REP.ASV
  b=otus$gid_16s%in%micropita20$REP.OTU
  c=asvs$gid_16s%in%micropita20$EXT.ASV
  d=otus$gid_16s%in%micropita20$EXT.OTU
  e=asvs$gid_16s%in%micropita20$DIV.ASV
  f=otus$gid_16s%in%micropita20$DIV.OTU
  g=asvs$gid_16s%in%micropita20$DIST.ASV
  h=otus$gid_16s%in%micropita20$DIST.OTU
  j=asvs$gid_16s%in%micropita20$DISCR.ASV
  k=otus$gid_16s%in%micropita20$DISCR.OTU
  l=asvs$gid_16s%in%micropita20$PCA.ASV
  m=otus$gid_16s%in%micropita20$PCA.OTU
  n=asvs$gid_16s%in%micropita20$SPCA.ASV
  o=otus$gid_16s%in%micropita20$SPCA.OTU
  Microtable$Repasv20[i]<-ifelse(a[i]==TRUE,1,0)
  Microtable$Repotu20[i]<-ifelse(b[i]==TRUE,1,0)
  Microtable$Extasv20[i]<-ifelse(c[i]==TRUE,1,0)
  Microtable$Extotu20[i]<-ifelse(d[i]==TRUE,1,0)
  Microtable$Divasv20[i]<-ifelse(e[i]==TRUE,1,0)
  Microtable$Divotu20[i]<-ifelse(f[i]==TRUE,1,0)
  Microtable$Distasv20[i]<-ifelse(g[i]==TRUE,1,0)
  Microtable$Distotu20[i]<-ifelse(h[i]==TRUE,1,0)
  Microtable$Discrasv20[i]<-ifelse(j[i]==TRUE,1,0)
  Microtable$Discrotu20[i]<-ifelse(k[i]==TRUE,1,0)
  Microtable$Pcaasv20[i]<-ifelse(l[i]==TRUE,1,0)
  Microtable$Pcaotu20[i]<-ifelse(m[i]==TRUE,1,0)
  Microtable$Spcaasv20[i]<-ifelse(n[i]==TRUE,1,0)
  Microtable$Spcaotu20[i]<-ifelse(o[i]==TRUE,1,0)
}
#50
for (i in 1:length(asvs$allergy)) {
  a=asvs$gid_16s%in%micropita50$ï..REP.ASV
  b=otus$gid_16s%in%micropita50$REP.OTU
  c=asvs$gid_16s%in%micropita50$EXT.ASV
  d=otus$gid_16s%in%micropita50$EXT.OTU
  e=asvs$gid_16s%in%micropita50$DIV.ASV
  f=otus$gid_16s%in%micropita50$DIV.OTU
  g=asvs$gid_16s%in%micropita50$DIST.ASV
  h=otus$gid_16s%in%micropita50$DIST.OTU
  j=asvs$gid_16s%in%micropita50$DISCR.ASV
  k=otus$gid_16s%in%micropita50$DISCR.OTU
  l=asvs$gid_16s%in%micropita50$PCA.ASV
  m=otus$gid_16s%in%micropita50$PCA.OTU
  n=asvs$gid_16s%in%micropita50$SPCA.ASV
  o=otus$gid_16s%in%micropita50$SPCA.OTU
  Microtable$Repasv50[i]<-ifelse(a[i]==TRUE,1,0)
  Microtable$Repotu50[i]<-ifelse(b[i]==TRUE,1,0)
  Microtable$Extasv50[i]<-ifelse(c[i]==TRUE,1,0)
  Microtable$Extotu50[i]<-ifelse(d[i]==TRUE,1,0)
  Microtable$Divasv50[i]<-ifelse(e[i]==TRUE,1,0)
  Microtable$Divotu50[i]<-ifelse(f[i]==TRUE,1,0)
  Microtable$Distasv50[i]<-ifelse(g[i]==TRUE,1,0)
  Microtable$Distotu50[i]<-ifelse(h[i]==TRUE,1,0)
  Microtable$Discrasv50[i]<-ifelse(j[i]==TRUE,1,0)
  Microtable$Discrotu50[i]<-ifelse(k[i]==TRUE,1,0)
  Microtable$Pcaasv50[i]<-ifelse(l[i]==TRUE,1,0)
  Microtable$Pcaotu50[i]<-ifelse(m[i]==TRUE,1,0)
  Microtable$Spcaasv50[i]<-ifelse(n[i]==TRUE,1,0)
  Microtable$Spcaotu50[i]<-ifelse(o[i]==TRUE,1,0)
}
#100
for (i in 1:length(asvs$allergy)) {
  a=asvs$gid_16s%in%micropita100$ï..REP.ASV
  b=otus$gid_16s%in%micropita100$REP.OTU
  c=asvs$gid_16s%in%micropita100$EXT.ASV
  d=otus$gid_16s%in%micropita100$EXT.OTU
  e=asvs$gid_16s%in%micropita100$DIV.ASV
  f=otus$gid_16s%in%micropita100$DIV.OTU
  g=asvs$gid_16s%in%micropita100$DIST.ASV
  h=otus$gid_16s%in%micropita100$DIST.OTU
  j=asvs$gid_16s%in%micropita100$DISCR.ASV
  k=otus$gid_16s%in%micropita100$DISCR.OTU
  l=asvs$gid_16s%in%micropita100$PCA.ASV
  m=otus$gid_16s%in%micropita100$PCA.OTU
  n=asvs$gid_16s%in%micropita100$SPCA.ASV
  o=otus$gid_16s%in%micropita100$SPCA.OTU
  Microtable$Repasv100[i]<-ifelse(a[i]==TRUE,1,0)
  Microtable$Repotu100[i]<-ifelse(b[i]==TRUE,1,0)
  Microtable$Extasv100[i]<-ifelse(c[i]==TRUE,1,0)
  Microtable$Extotu100[i]<-ifelse(d[i]==TRUE,1,0)
  Microtable$Divasv100[i]<-ifelse(e[i]==TRUE,1,0)
  Microtable$Divotu100[i]<-ifelse(f[i]==TRUE,1,0)
  Microtable$Distasv100[i]<-ifelse(g[i]==TRUE,1,0)
  Microtable$Distotu100[i]<-ifelse(h[i]==TRUE,1,0)
  Microtable$Discrasv100[i]<-ifelse(j[i]==TRUE,1,0)
  Microtable$Discrotu100[i]<-ifelse(k[i]==TRUE,1,0)
  Microtable$Pcaasv100[i]<-ifelse(l[i]==TRUE,1,0)
  Microtable$Pcaotu100[i]<-ifelse(m[i]==TRUE,1,0)
  Microtable$Spcaasv100[i]<-ifelse(n[i]==TRUE,1,0)
  Microtable$Spcaotu100[i]<-ifelse(o[i]==TRUE,1,0)
}


#TEMPORARY, the "list" is manually changed from N=20,50,100
#make the Venn diagrams for 20

venn.diagram(
  x = list(
    Microtable %>% filter(Repasv20==1) %>% select(`asvs$gid_16s`) %>% unlist(), 
    Microtable %>% filter(Repotu20==1) %>% select(`asvs$gid_16s`) %>% unlist()
  ),
  category.names = c("Representative ASVs" , "Representative OTUs"),
  filename = 'Rep20venn_diagramm.png',
  output = TRUE,
  # Circles
  lwd = 2,
  col='#21908dff')
venn.diagram(
  x = list(
    Microtable %>% filter(Extasv20==1) %>% select(`asvs$gid_16s`) %>% unlist(), 
    Microtable %>% filter(Extotu20==1) %>% select(`asvs$gid_16s`) %>% unlist()
  ),
  category.names = c("Extreme ASVs" , "Extreme OTUs"),
  filename = 'Ext20venn_diagramm.png',
  output = TRUE)
venn.diagram(
  x = list(
    Microtable %>% filter(Divasv20==1) %>% select(`asvs$gid_16s`) %>% unlist(), 
    Microtable %>% filter(Divotu20==1) %>% select(`asvs$gid_16s`) %>% unlist()
  ),
  category.names = c("Diverse ASVs" , "Diverse OTUs"),
  filename = 'Div20venn_diagramm.png',
  output = TRUE)
venn.diagram(
  x = list(
    Microtable %>% filter(Distasv20==1) %>% select(`asvs$gid_16s`) %>% unlist(), 
    Microtable %>% filter(Distotu20==1) %>% select(`asvs$gid_16s`) %>% unlist()
  ),
  category.names = c("Distinct ASVs" , "Distinct OTUs"),
  filename = 'Disct20venn_diagramm.png',
  output = TRUE)
venn.diagram(
  x = list(
    Microtable %>% filter(Discrasv20==1) %>% select(`asvs$gid_16s`) %>% unlist(), 
    Microtable %>% filter(Discrotu20==1) %>% select(`asvs$gid_16s`) %>% unlist()
  ),
  category.names = c("Discriminant ASVs" , "Discriminant OTUs"),
  filename = 'Discr20venn_diagramm.png',
  output = TRUE)
venn.diagram(
  x = list(
    Microtable %>% filter(Pcaasv20==1) %>% select(`asvs$gid_16s`) %>% unlist(), 
    Microtable %>% filter(Pcaotu20==1) %>% select(`asvs$gid_16s`) %>% unlist()
  ),
  category.names = c("PCA ASVs" , "PCA OTUs"),
  filename = 'Pca20venn_diagramm.png',
  output = TRUE)
venn.diagram(
  x = list(
    Microtable %>% filter(Spcaasv20==1) %>% select(`asvs$gid_16s`) %>% unlist(), 
    Microtable %>% filter(Spcaotu20==1) %>% select(`asvs$gid_16s`) %>% unlist()
  ),
  category.names = c("SPCA ASVs" , "SPCA OTUs"),
  filename = 'Spca20venn_diagramm.png',
  output = TRUE)

#make the venn diagrams for 50
venn.diagram(
  x = list(
    Microtable %>% filter(Repasv50==1) %>% select(`asvs$gid_16s`) %>% unlist(), 
    Microtable %>% filter(Repotu50==1) %>% select(`asvs$gid_16s`) %>% unlist()
  ),
  category.names = c("Representative ASVs" , "Representative OTUs"),
  filename = 'Rep50venn_diagramm.png',
  output = TRUE)
venn.diagram(
  x = list(
    Microtable %>% filter(Extasv50==1) %>% select(`asvs$gid_16s`) %>% unlist(), 
    Microtable %>% filter(Extotu50==1) %>% select(`asvs$gid_16s`) %>% unlist()
  ),
  category.names = c("Extreme ASVs" , "Extreme OTUs"),
  filename = 'Ext50venn_diagramm.png',
  output = TRUE)
venn.diagram(
  x = list(
    Microtable %>% filter(Divasv50==1) %>% select(`asvs$gid_16s`) %>% unlist(), 
    Microtable %>% filter(Divotu50==1) %>% select(`asvs$gid_16s`) %>% unlist()
  ),
  category.names = c("Diverse ASVs" , "Diverse OTUs"),
  filename = 'Div50venn_diagramm.png',
  output = TRUE)
venn.diagram(
  x = list(
    Microtable %>% filter(Distasv50==1) %>% select(`asvs$gid_16s`) %>% unlist(), 
    Microtable %>% filter(Distotu50==1) %>% select(`asvs$gid_16s`) %>% unlist()
  ),
  category.names = c("Distinct ASVs" , "Distinct OTUs"),
  filename = 'Disct50venn_diagramm.png',
  output = TRUE)
venn.diagram(
  x = list(
    Microtable %>% filter(Discrasv50==1) %>% select(`asvs$gid_16s`) %>% unlist(), 
    Microtable %>% filter(Discrotu50==1) %>% select(`asvs$gid_16s`) %>% unlist()
  ),
  category.names = c("Discriminant ASVs" , "Discriminant OTUs"),
  filename = 'Discr50venn_diagramm.png',
  output = TRUE)
venn.diagram(
  x = list(
    Microtable %>% filter(Pcaasv50==1) %>% select(`asvs$gid_16s`) %>% unlist(), 
    Microtable %>% filter(Pcaotu50==1) %>% select(`asvs$gid_16s`) %>% unlist()
  ),
  category.names = c("PCA ASVs" , "PCA OTUs"),
  filename = 'Pca50venn_diagramm.png',
  output = TRUE)
venn.diagram(
  x = list(
    Microtable %>% filter(Spcaasv50==1) %>% select(`asvs$gid_16s`) %>% unlist(), 
    Microtable %>% filter(Spcaotu50==1) %>% select(`asvs$gid_16s`) %>% unlist()
  ),
  category.names = c("SPCA ASVs" , "SPCA OTUs"),
  filename = 'Spca50venn_diagramm.png',
  output = TRUE)

#make the Venn diagrams for 100
venn.diagram(
  x = list(
    Microtable %>% filter(Repasv100==1) %>% select(`asvs$gid_16s`) %>% unlist(), 
    Microtable %>% filter(Repotu100==1) %>% select(`asvs$gid_16s`) %>% unlist()
  ),
  category.names = c("Representative ASVs" , "Representative OTUs"),
  filename = 'Rep100venn_diagramm.png',
  output = TRUE)
venn.diagram(
  x = list(
    Microtable %>% filter(Extasv100==1) %>% select(`asvs$gid_16s`) %>% unlist(), 
    Microtable %>% filter(Extotu100==1) %>% select(`asvs$gid_16s`) %>% unlist()
  ),
  category.names = c("Extreme ASVs" , "Extreme OTUs"),
  filename = 'Ext100venn_diagramm.png',
  output = TRUE)
venn.diagram(
  x = list(
    Microtable %>% filter(Divasv100==1) %>% select(`asvs$gid_16s`) %>% unlist(), 
    Microtable %>% filter(Divotu100==1) %>% select(`asvs$gid_16s`) %>% unlist()
  ),
  category.names = c("Diverse ASVs" , "Diverse OTUs"),
  filename = 'Div100venn_diagramm.png',
  output = TRUE)
venn.diagram(
  x = list(
    Microtable %>% filter(Distasv100==1) %>% select(`asvs$gid_16s`) %>% unlist(), 
    Microtable %>% filter(Distotu100==1) %>% select(`asvs$gid_16s`) %>% unlist()
  ),
  category.names = c("Distinct ASVs" , "Distinct OTUs"),
  filename = 'Disct100venn_diagramm.png',
  output = TRUE)
venn.diagram(
  x = list(
    Microtable %>% filter(Discrasv100==1) %>% select(`asvs$gid_16s`) %>% unlist(), 
    Microtable %>% filter(Discrotu100==1) %>% select(`asvs$gid_16s`) %>% unlist()
  ),
  category.names = c("Discriminant ASVs" , "Discriminant OTUs"),
  filename = 'Discr100venn_diagramm.png',
  output = TRUE)
venn.diagram(
  x = list(
    Microtable %>% filter(Pcaasv100==1) %>% select(`asvs$gid_16s`) %>% unlist(), 
    Microtable %>% filter(Pcaotu100==1) %>% select(`asvs$gid_16s`) %>% unlist()
  ),
  category.names = c("PCA ASVs" , "PCA OTUs"),
  filename = 'Pca100venn_diagramm.png',
  output = TRUE)
venn.diagram(
  x = list(
    Microtable %>% filter(Spcaasv100==1) %>% select(`asvs$gid_16s`) %>% unlist(), 
    Microtable %>% filter(Spcaotu100==1) %>% select(`asvs$gid_16s`) %>% unlist()
  ),
  category.names = c("SPCA ASVs" , "SPCA OTUs"),
  filename = 'Spca100venn_diagramm.png',
  output = TRUE)

#The following code generates the tsv files used for LEfSe analysis
Microtable=Microtable %>% 
  rename(gid_16s=`asvs$gid_16s`)
#A function, that iterates through Microtable, and makes a text file with the needed
#format for running LeFSE
strat_names<-colnames(Microtable)
rows.asv<-read.csv("to_subsample_asv.csv")
rows.otu<-read.csv("to_subsample_otu.csv")

For.lefse<-function(dataframe, column=2, asv=TRUE){
  strategy<-Microtable[c(1,column)]
  subsample<-merge(strategy,dataframe, by= 'gid_16s')
  #need to transpose the dataframe for LefSE
  ids <- subsample$gid_16s

  # transpose all but the first column (name)
  subsample_t <- as.data.frame(t(subsample[,-1]))
  if(asv==TRUE){
    subsample_t <- cbind(c(strat_names[column],rows.asv$ï..gid_16s), subsample_t) 
  }
  else{
    subsample_t <- cbind(c(strat_names[column],rows.otu$ï..gid_16s), subsample_t)
  }
  colnames(subsample_t) <- c("gid_16s",ids)
  
  #exporting the csv file, I need to add the | manually after to each
  write.table(subsample_t, 
              file = paste0(strat_names[column], ".txt"),
              row.names=FALSE, sep="\t", quote = FALSE)
}
for (i in (1:21)*2) {
  For.lefse(dataframe = asvs, column = i, asv = TRUE)
}
for (i in (1:21)*2+1) {
  For.lefse(dataframe = otus, column = i, asv = FALSE)
}


#For metagenomic data
#I use the metadata to extract the 16s and wgs sample ids,
#as well as the allergy variable,
#I make use of the "microtable" made for the venn_diagrams
kareliametadata<-read.csv('kareliametadata.csv')
sample_ids<-subset(kareliametadata, select = c(gid_16s.y,gid_wgs.y, allergy))
#dropping the y from the columns
sample_ids=sample_ids %>% 
  rename(
    gid_wgs=gid_wgs.y,
    gid_16s=gid_16s.y
  )
#below are the wgs_sampleid of the samples I could run metaphlan on
#generated by doing a merge of all 
#available metaphlan outputs and isolating the wgs column
wgs_samples<- read.csv('wgs_samples.csv')
wgs_samples=wgs_samples %>% 
  rename(
    gid_wgs=ï..wgs_sampleid
  )
#I then merge the two data frames
metaphlan_samples<-merge(x=wgs_samples,y=sample_ids,by="gid_wgs",all.x=TRUE)
#I then merge with Microtable to make the dataframe used for subsampling


metaphlan_subsamples<-merge(x=metaphlan_samples,
                            y=Microtable,by="gid_16s",all.x=TRUE)



subsetter<-function(dataframe=metaphlan_subsamples, column=4, allergy=FALSE){
  new_df=subset(dataframe, dataframe[column]==1, select = gid_wgs)
  names(new_df)<-NULL
  data_colname<-colnames(dataframe)
  if(allergy==TRUE){
    new_df=subset(dataframe, dataframe[column]==1, select = c(gid_wgs, allergy))
    write.csv(new_df, 
                file = paste0(data_colname[column], "_allergy.csv"),
                row.names=FALSE, sep="\t", quote = FALSE)
  }
  else{
  write.table(new_df, 
              file = paste0(data_colname[column], ".txt"),
              row.names=FALSE, sep="\t", quote = FALSE)
  }
}

for (i in 4:45) {
  subsetter(column = i)
}


#need file with allergies to insert the corresponding class for LEfSe
for (i in 4:45) {
  subsetter(column = i, allergy = TRUE)
}

colnames(Microtable)


