library(tidyverse)
#First open the asv files
asvs<- read.csv('relative_asvT.csv')
otus<-read.csv('relative-otu.csv')

#need to import the LEfSe LDA files
#Disc20asv<-read.table('LEfSe/Discasv20_LDA.txt')

#Need to do a join somehow
asvs[1:5,1]
#Need to make the features consistent between tables
#First remove periods
#Use \\| to eascape the | character, turn them to periods
#Replace .__ with ___
asvs$ï..gid_16s<-gsub("\\.", "_", as.character(asvs$ï..gid_16s))
asvs$ï..gid_16s<-gsub("\\|", ".", as.character(asvs$ï..gid_16s))
asvs$ï..gid_16s<-gsub("\\.__", "___", as.character(asvs$ï..gid_16s))
asvs=asvs %>% 
  rename(
    feature=ï..gid_16s
  )
Discasv20_LDA=Discasv20_LDA %>% 
  rename(
    feature=V1
  )
sum(asvs$feature %in% Discasv20_LDA$feature)
trial=merge(x = Discasv20_LDA, y = asvs, by = "feature", all.x = TRUE)
#Note that the LDA file has more features, this is because of ASVS which show up 
#at only at a genus level in the ASV generate extra features at each level in the LDA file



#asv<-read.csv('to_subsample_asv.csv')
#taxa_asv<-asv[2:327,1]
#taxa_asv


#gsub replaces the first option with the second option
#format is D_0__kingdom.D_1__phylum.D_2__class.D_3__order.D_4__family.D_5__genus
#when its unknown
trial$feature
taxa_asv<-as.character(trial$feature[2:492])
#the following replaces all .D_O with .D_O__Other
#taxa1_asv<-gsub(".D_O",".D_O__Other",as.character(asv$ï..gid_16s[2:327]))
#Then it collapses the repeated .D_O__Other to a single one
#i=1
#while(i<4){
  #i=i+1
  #taxa1_asv<-gsub(".D_O__Other.D_O__Other",".D_O__Other", taxa1_asv)
#}
taxa1_asv<-strsplit(taxa_asv,".D_[12345]_")
taxa1_asv
taxa_asv
unique(sapply(taxa1_asv,length))
taxa2_asv<-sapply(taxa1_asv,function(s){
  len<-length(s)
  if(len==1){
    s1<-s
  }else if(len==2){
    s1<-paste0(s[1],".1_",s[2])
  }else if(len==3){
    s1<-paste0(s[1],".1_",s[2],".2_",s[3])
  }else if(len==4){
    s1<-paste0(s[1],".1_",s[2],".2_",s[3],".3_",s[4])
  }else if(len==5){
    s1<-paste0(s[1],".1_",s[2],".2_",s[3],".3_",s[4],".4_",s[5])
  }else  
    s1<-paste0(s[1],".1_",s[2],".2_",s[3],".3_",s[4],".4_",s[5],".5_",s[6])
  return(s1)})
taxa2_asv

#possible problems ".Ambiguous_taxa"
#and "uncultured.something"
taxa1_asv
taxa3_asv<-sapply(taxa1_asv,function(s){
  len<-length(s)
  s[s=="_"]<-"_Other"
  s<-gsub("Ambiguous_taxa","Unclasified",gsub("_","",gsub("[.]","-",s)))
  return(paste(s,collapse = "."))})
taxa3_asv<-c("allergy",taxa3_asv)
taxa3_asv
#past<-taxa3_asv
past
trial$otu1<-taxa2
#trial$otu2<- gsub("_Other","_Unclasified",taxa2)
trial$asv3<- taxa3
b<-as.data.frame(cbind(trial[1:5],taxa3_asv))
#b looks as desired
b
#which(taxa$otu2 %in% alltaxa)
#hertaxaB<-subset(taxa, H2r>0.5 & abs(Kurtosis)<2)
#dim(hertaxaB); head(hertaxaB)

#ummary(taxa$H2r)
#taxa$H2r_nrm<-taxa$H2r/max(taxa$H2r)
#summary(taxa$H2r_nrm)

write.table(taxa3_asv[2:492],"guide2.txt",quote = F,row.names = F,col.names = F,sep = "\t")


#for below, sink starts the text file?,
#Cat outputs each line as a line in the text file
#\t is tab 


#The next line adds every row of the taxa table to the file,
#and adds the text "ring_alpha 1 ?" after it, where ? is hte H2R_nrm value
write.table(cbind(taxa$otu3,"ring_alpha","1",taxa$H2r_nrm),"circular_plots/annot_3.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)

#It then does the same thing but with "ring_color 1 #AAAA00"
write.table(cbind(taxa$otu3,"ring_color","1","#AAAA00"),"circular_plots/annot_3.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)

sink(file = "circular_plots/annot_3.txt",append = T)
cat("Significant OTUs\tclade_marker_color\t#FF0000
Significant OTUs\tclade_marker_size\t20
Significant OTUs\tclass_label\tOTUs with significant H2r
")
sink()
write.table(cbind(taxa$otu3[taxa$H2rpval<0.05/253],"class","Significant OTUs"),"circular_plots/annot_3.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)

##Anotate classes
class<-unique(sapply(strsplit(taxa$otu3[taxa$H2rpval<0.05/253],"[.]"),function(x)return(x[5])))
class<-class[!is.na(class)]
class


c<-which(b$V3>=1 & !is.na(b$V3))
d<-which(b$V3<=0 & !is.na(b$V3))




### ASV results  - circular plot annotation
sink(file = "annot_asv20.txt")
cat("title\tASVs
title_font_size\t13
start_rotation\t270
total_plotted_degrees\t360	
*\tbranch_thickness\t0.916126068736
*\tannotation_background_separation\t0.001
*\tbranch_color_from_ancestor\t0
*\tbranch_bracket_depth\t0.75
*\tring_edge_width\t1\t0.0
annotation_font_size\t7
annotation_legend_font_size\t10
ring_internal_separator_thickness\t1\t0.5
ring_width\t1\t0.5
ring_height\t1\t0.75
")
sink()
sink(file = "annot_asv20.txt",append = T)
cat("Discriminative\tclade_marker_color\t#556B2F
Discriminative\tclade_marker_size\t20
Discriminative\tclass_label\tDiscriminative
")
sink()
#the over expressed
write.table(cbind(b$taxa3_asv[c],"ring_color","2","#556B2F"),"annot_asv.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)
write.table(cbind(b$taxa3_asv[c],"ring_shape","2","^"),"annot_asv.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)
#the under expressed
write.table(cbind(b$taxa3_asv[d],"ring_color","3","#DD0000"),"annot_asv.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)
write.table(cbind(b$taxa3_asv[d],"ring_shape","3","v"),"annot_asv.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)









#write.table(cbind(class,"clade_marker_color","#DD0000"),"circular_plots/annot_3.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)
write.table(cbind(class,"annotation",class),"circular_plots/annot_3.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)
write.table(cbind(class,"annotation_background_color","#DD0000"),"circular_plots/annot_3.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)
##FFD700


#Im thinking, run file through LefSe, download the LDA effect sizes table,
#if the row has a 1, it is over expressed in the subsample,
#if its 0 it is underexpressed in the subsample
#else its not significantly different
#so specify this and add it to the annotation file somehow?
#like if its 0 do a down arrow for this clade with this color
#if its 1 do an up arrow for this clade with this color
#do a diff ring level per each subsample method?
#do this per each type of subsample????
#Do a diff annotation file for each size 20,50,100!



