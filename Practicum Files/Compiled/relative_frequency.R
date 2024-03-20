library(tidyverse)
#Code to take pseudo-asv and pseudo-otu table into relative frequency tables
pseudo.asv<-read.csv("pseudo-asv.csv")
relative.asv<-(pseudo.asv[1,-1])
for (i in 1:199) {
  new.row<-pseudo.asv[i,-1]/sum(pseudo.asv[i,-1])
  relative.asv<-rbind(relative.asv, new.row)
}
relative.asv<-relative.asv[-1,]
relative.asv<-cbind(pseudo.asv$ï..OTU.ID,relative.asv)

pseudo.otu<-read.csv("pseudo-otu.csv")
relative.otu<-(pseudo.otu[1,-1])
for (i in 1:199) {
  new.row<-pseudo.otu[i,-1]/sum(pseudo.otu[i,-1])
  relative.otu<-rbind(relative.otu, new.row)
}
relative.otu<-relative.otu[-1,]
relative.otu<-cbind(pseudo.otu$ï..OTU.ID,relative.otu)

#use allergy.csv, generated using kareliametadata (file is in R-code directory)
#code that was used to generate allergy variable is in Metadata_vis.R
#code used to generate allergy.csv is below
kareliametadata<-read.csv('kareliametadata.csv')
kareliametadata=kareliametadata %>% 
  rename(
    gid_16s=gid_16s.x
  )
allergy <- kareliametadata[c("gid_16s","allergy")]
#write.csv(allergy,"allergy.csv", row.names = FALSE)


#imported allergy.csv
relative.asv=relative.asv %>% 
  rename(
    gid_16s=`pseudo.asv$ï..OTU.ID`
  )
relative.otu=relative.otu %>% 
  rename(
    gid_16s=`pseudo.otu$ï..OTU.ID`
  )

relative.asv<-merge(allergy, relative.asv)
relative.otu<-merge(allergy,relative.otu)

write.csv(relative.asv,"relative-asv.csv", row.names = FALSE)
write.csv(relative.otu,"relative-otu.csv", row.names = FALSE)




