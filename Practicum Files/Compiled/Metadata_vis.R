largemetadata<- read.csv('Metadata.csv')
library(dplyr)
library(tidyverse)
# Adding column based on other column:
largemetadata<-largemetadata %>% mutate(Complete=
                     case_when(gid_wgs == ''&gid_16s=='' ~ "Both Missing", 
                               gid_wgs == ''& gid_16s!='' ~ "Wgs Missing",
                               gid_wgs != ''&gid_16s=='' ~ "16S Missing",
                               gid_wgs!='' & gid_16s!=''~"Complete")
)
largemetadata
ages<-seq(26, 1166, by = 60)
ages
length(ages)
metatest<-largemetadata
post<-(head(metatest,1))
agegroup<-metatest[metatest$age_at_collection<=ages[2],]
agegroup


for (i in 1:20){
  agegroup<-metatest[metatest$age_at_collection<=ages[i],]
  metatest<-metatest[metatest$age_at_collection>ages[i],]
  agegroup<-agegroup %>% distinct(subjectID, .keep_all = TRUE)
  post <- rbind(post, agegroup)
  }
post
#need to drop the duplicated first row and last row
post<-post %>% distinct()
ages
agelabels=c("26-86","87-146","147-206","207-266","267-326","327-386",
            "387-446","447-506","507-566","567-626","627-686","687-746",
            "747-806","807-866","867-926","927-986","987-1046","1047-1106"
            ,"1107-1166")
library("ggplot2")
post$group <- cut(post$age_at_collection, breaks = ages,labels =agelabels ,right = TRUE)
ggplot(post, aes(x = group, fill = Complete)) + 
  geom_bar()+ggtitle("All cohorts")



metatest<-largemetadata
abx<-post[post$cohort == 'abx', ]
karelia<-post[post$cohort =='karelia',]
t1d<-post[post$cohort =='t1d',]
ggplot(abx, aes(x = group, fill = Complete)) + 
  geom_bar()+ggtitle("Antibiotics Cohort")
ggplot(karelia, aes(x = group, fill = Complete)) + 
  geom_bar()+ggtitle("Three Country Cohort")
ggplot(t1d, aes(x = group, fill = Complete)) + 
  geom_bar()+ggtitle("Type 1 Diabetes Cohort")
metatest
complete_only<-metatest[metatest$Complete=="Complete",]
complete_only<-complete_only%>%distinct(subjectID,.keep_all = TRUE)
complete_only

complete_only$group <- cut(complete_only$age_at_collection, breaks = ages,labels =agelabels ,right = TRUE)
ggplot(complete_only, aes(x = group, fill = Complete)) + 
  geom_bar()+ggtitle("All cohorts")

abx1<-complete_only[complete_only$cohort == 'abx', ]
karelia1<-complete_only[complete_only$cohort =='karelia',]
t1d1<-complete_only[complete_only$cohort =='t1d',]
ggplot(abx1, aes(x = group, fill = Complete)) + 
  geom_bar()+ggtitle("Antibiotics Cohort")
ggplot(karelia1, aes(x = group, fill = Complete)) + 
  geom_bar()+ggtitle("Three Country Cohort")
ggplot(t1d1, aes(x = group, fill = Complete)) + 
  geom_bar()+ggtitle("Type 1 Diabetes Cohort")

write.csv(karelia1,"karelia.csv", row.names = FALSE)


karelia1




##First histograms *
#onlt samples with both mgs, and 16s from same day
metadata2<- read.csv('Metadata(complete only).csv')
abx2<-metadata2[metadata2$cohort == 'abx', ]
karelia2<-metadata2[metadata2$cohort =='karelia',]
t1d2<-metadata2[metadata2$cohort =='t1d',]
hist(abx2$age_at_collection)
hist(karelia2$age_at_collection)
hist(t1d2$age_at_collection)

karelia2

#metadata from kaerlia only
#downloaded metadata from diabimmune project, the karelia cohort only
metadata2
metadata$SampleID
metadata<-metadata %>% 
  rename(
    sampleID = SampleID
  )
#making allergy column
metadata$allergy
metadata <- metadata %>% 
  mutate(allergy = if_else(allergy_milk == TRUE|allergy_egg==TRUE|
                           allergy_peanut==TRUE|allergy_dustmite==TRUE|
                             allergy_cat==TRUE|allergy_dog==TRUE|
                             allergy_birch==TRUE|allergy_timothy, 1, 0))
summary(metadata$allergy)
metadata$allergy[is.na(metadata$allergy)] <- min(metadata$allergy, na.rm = T)
table(metadata$allergy)

df = merge(x=metadata,y=karelia1,by="sampleID", all.y = TRUE)
ggplot(df, aes(x =  group, fill = Complete)) + 
  geom_bar()+ggtitle("Three Country Cohort")
ggplot(df, aes(x =  delivery, fill = Complete)) + 
  geom_bar()+ggtitle("Three Country Cohort")
ggplot(df, aes(x =  gender, fill = Complete)) + 
  geom_bar()+ggtitle("Three Country Cohort")
ggplot(df, aes(x =  Exclusive_breast_feeding, fill = Complete)) + 
  geom_bar()+ggtitle("Three Country Cohort")
ggplot(df, aes(x =  allergy_milk, fill = Complete)) + 
  geom_bar()+ggtitle("Three Country Cohort")

summary(df)
write.csv(df,"kareliametadata.csv", row.names = FALSE)
 
#finding number of distinct karelia

karelia3<-largemetadata[largemetadata$cohort=="karelia",]
kareliacomplete<-karelia3[karelia3$Complete=="Complete",]
karelia_distinct<-kareliacomplete%>%distinct(subjectID,.keep_all = TRUE)

