#Applying inclusion exclusion criteria
#load in the libraries used
largemetadata<- read.csv('Metadata.csv')
library(dplyr)
library(tidyverse)
# Adding completeness column based on whether there is 16s or wgs data:
largemetadata<-largemetadata %>% 
  mutate(Complete=case_when(gid_wgs == ''&gid_16s=='' ~ "Both Missing", 
                            gid_wgs == ''& gid_16s!='' ~ "Wgs Missing",
                            gid_wgs != ''&gid_16s=='' ~ "16S Missing",
                            gid_wgs!='' & gid_16s!=''~"Complete")
)
#finding number of distinct karelia
#only karelia cohort
karelia<-largemetadata[largemetadata$cohort=="karelia",]
#only complete samples
kareliacomplete<-karelia[karelia$Complete=="Complete",]
#only 1 sample per subject
karelia_distinct<-kareliacomplete%>%distinct(subjectID,.keep_all = TRUE)


#Import DIABIMMUNE_Karelia_metadata available from DIABIMMUNE database
metadata <- metadata %>% 
  mutate(allergy = if_else(allergy_milk == TRUE|allergy_egg==TRUE|
                             allergy_peanut==TRUE|allergy_dustmite==TRUE|
                             allergy_cat==TRUE|allergy_dog==TRUE|
                             allergy_birch==TRUE|allergy_timothy, 1, 0))
#making all missing allergy instances 0
metadata$allergy[is.na(metadata$allergy)] <- min(metadata$allergy, na.rm = T)
table(metadata$allergy)


#reanming the sampleID column to do a join
metadata<-metadata %>% 
  rename(
    sampleID = SampleID
  )
sample_metadata = merge(x=metadata,y=karelia_distinct,by="sampleID", all.y = TRUE)
write.csv(sample_metadata,"kareliametadata.csv", row.names = FALSE)
#The gid_16s column was used to manually create the manifest file in excel, 
#so that only these samples were processed using QIIME 2