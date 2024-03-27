setwd("/Volumes/garcia/GEM")

source("Sent_13_10_2014_Osvaldo_GWAS/alltaxa.R")


taxa<-read.table('combined_output.txt', header=T, sep = "\t", quote = '')
head(taxa); dim(taxa)
## Put the taxa in the right format
tail(as.character(taxa$otu))

taxa$otu

#gsub finds a match for a pattern, replaces it with a different value
#Here, gsub replaces "_Other" with "_t__Other", for the OTUs with unspecified genus or species
taxa1<-gsub("_Other","_t__Other",as.character(taxa$otu))
taxa1
#The structure of the taxa is k_..._p__..._c__..._f__..._g__...
#We then split the string by clade group, grab all the ... from above
taxa1<-strsplit(taxa1,"_[pcofgt]_")
taxa1
#The following checks the "lenght" of each row
unique(sapply(taxa1,length))
#since the max is 6, no rows specify more than a genus (no species)
#Then we reconnect each split string, but the new format is:
#k__kingdom.p_phylum.c__class.o__order.f__family.g__genus
taxa2<-sapply(taxa1,function(s){
  len<-length(s)
  if(len==2){
    s1<-paste0(s[1],".p_",s[2])
  }else if(len==3){
    s1<-paste0(s[1],".p_",s[2],".c_",s[3])
  }else if(len==4){
    s1<-paste0(s[1],".p_",s[2],".c_",s[3],".o_",s[4])
  }else if(len==5){
    s1<-paste0(s[1],".p_",s[2],".c_",s[3],".o_",s[4],".f_",s[5])
  }else  
    s1<-paste0(s[1],".p_",s[2],".c_",s[3],".o_",s[4],".f_",s[5],".g_",s[6])
  return(s1)})
taxa2
#The next line checks if there are dupes in the otu list, all are unique
sum(duplicated(taxa$otu))

taxa$otu[224]
#Using the split strings in taxa1, it first replaces any strings that
#are solely '_' with "Other"
#the following replaces periods with -, then removes _, and then replaces "Other
#with "Unclassified"
#Lastly, it pastes the substrings in each row with a . separating each component
taxa3<-sapply(taxa1,function(s){
  len<-length(s)
  s[s=="_"]<-"_Other"
  s<-gsub("Other","Unclasified",gsub("_","",gsub("[.]","-",s)))
  return(paste(s,collapse = "."))})

taxa3

taxa$otu1<-taxa2
taxa$otu2<- gsub("_Other","_Unclasified",taxa2)
taxa$otu3<- taxa3


which(taxa$otu2 %in% alltaxa)
hertaxaB<-subset(taxa, H2r>0.5 & abs(Kurtosis)<2)
dim(hertaxaB); head(hertaxaB)

summary(taxa$H2r)
taxa$H2r_nrm<-taxa$H2r/max(taxa$H2r)
summary(taxa$H2r_nrm)
taxa3
#The final one to use was taxa3!
write.table(taxa3,"guide.txt",quote = F,row.names = F,col.names = F,sep = "\t")

### Heritable taxa results  - circular plot annotation
sink(file = "circular_plots/annot_3.txt")
cat("title\tHeritable Taxa
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
write.table(cbind(taxa$otu3,"ring_alpha","1",taxa$H2r_nrm),"circular_plots/annot_3.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)
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
#write.table(cbind(class,"clade_marker_color","#DD0000"),"circular_plots/annot_3.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)
write.table(cbind(class,"annotation",class),"circular_plots/annot_3.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)
write.table(cbind(class,"annotation_background_color","#DD0000"),"circular_plots/annot_3.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)
##FFD700


### GWAS' results  - circular plot annotation
names(taxa)
log_res<-read.csv("circular_plots/imputed_exwc_1098_Taxa_GWAS_sum_log_nonlow.csv")
tp_res<-read.csv("circular_plots/imputed_exwc_1098_Taxa_GWAS_sum_tplog_nonlow.csv")
counts<-read.table("Sent_13_10_2014_Osvaldo_GWAS/MB_exm_1098_ALL_taxa_166.txt",header=F)
counts<-counts[,3:ncol(counts)]
names(log_res)

names(counts)<-alltaxa
zp<-apply(counts,2,function(x)sum(x==0)/length(x))
zp_df<-data.frame(otu=names(zp),prop=zp)
head(zp_df)

gwas_info<-merge(merge(log_res,tp_res,by=c("taxa_id","taxa","taxa_short"),suffixes = c("_log","_tp"),sort=F),zp_df,by.x="taxa",by.y="otu",sort=F)

taxa_gwas<-merge(taxa,gwas_info,by.x="otu2",by.y="taxa",sort=F)
head(taxa_gwas); dim(taxa_gwas)

##Taxa indexes (on discovery) that need to be replicated
ids_disc1<-read.table("/Users/garcia/Documents/scinet/FDR/genomescan_illumina_201501/OTUs_replication_heritable.txt",header=T)
ids_for_rep<-unique(c(ids_disc1$Osvaldo_taxa_number))
ids_for_rep<-ids_for_rep[order(ids_for_rep)]
length(ids_for_rep)
alltaxa[ids_for_rep]

#check if all otus are in the list
taxa_gwas$otu2[taxa_gwas$taxa_id%in%ids_for_rep] %in% alltaxa[ids_for_rep]

###  Guide file
write.table(taxa_gwas$otu3,"circular_plots/guide_gwas.txt",quote = F,row.names = F,col.names = F,sep = "\t")

### Anotation file
sink(file = "circular_plots/annot_3_gwas.txt")
cat("title\tGWAS results
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
write.table(cbind(taxa_gwas$otu3,"ring_alpha","1",taxa_gwas$prop),"circular_plots/annot_3_gwas.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)
write.table(cbind(taxa_gwas$otu3,"ring_color","1","#8B1A1A"),"circular_plots/annot_3_gwas.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)
sink(file = "circular_plots/annot_3_gwas.txt",append = T)
cat("Significant OTUs\tclade_marker_color\t#66CD00
Significant OTUs\tclade_marker_size\t30
Significant OTUs\tclass_label\tGWAS significant OTUs
")
sink()
write.table(cbind(taxa_gwas$otu3[taxa_gwas$taxa_id%in%ids_for_rep],"class","Significant OTUs"),"circular_plots/annot_3_gwas.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)
##rings for GWAS models
sink(file = "circular_plots/annot_3_gwas.txt",append = T)
cat("ring_internal_separator_thickness\t2\t1
ring_internal_separator_thickness\t3\t1
")
sink()
### log models
id0<-which(taxa_gwas$lambda_log>=0.95 & taxa_gwas$lambda_log<=1.06 & !is.na(taxa_gwas$lambda_log))
write.table(cbind(taxa_gwas$otu3[id0],"ring_color","2","#556B2F"),"circular_plots/annot_3_gwas.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)
id1<-which((taxa_gwas$lambda_log<0.95 | taxa_gwas$lambda_log>1.06) & !is.na(taxa_gwas$lambda_log))
write.table(cbind(taxa_gwas$otu3[id1],"ring_color","2","#CD6600"),"circular_plots/annot_3_gwas.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)
id2<-which(is.na(taxa_gwas$lambda_log))
if(length(id2)!=0)write.table(cbind(taxa_gwas$otu3[id2],"ring_color","2","#FFFFFF"),"circular_plots/annot_3_gwas.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)

### two-part models
id0<-which(taxa_gwas$lambda_tp>=0.95 & taxa_gwas$lambda_tp<=1.06 & !is.na(taxa_gwas$lambda_tp))
write.table(cbind(taxa_gwas$otu3[id0],"ring_color","3","#556B2F"),"circular_plots/annot_3_gwas.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)
id1<-which((taxa_gwas$lambda_tp<0.95 | taxa_gwas$lambda_tp>1.06) & !is.na(taxa_gwas$lambda_tp))
write.table(cbind(taxa_gwas$otu3[id1],"ring_color","3","#CD6600"),"circular_plots/annot_3_gwas.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)
id2<-which(is.na(taxa_gwas$lambda_tp))
if(length(id2)!=0)write.table(cbind(taxa_gwas$otu3[id2],"ring_color","3","#FFFFFF"),"circular_plots/annot_3_gwas.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)


sink(file = "circular_plots/annot_3_gwas.txt",append = T)
cat("ring_internal_separator_thickness\t4\t1
")
sink()
### maximum p-value
  ##ids to be considered in the plots (lambda in range)
  id0a<-taxa_gwas$lambda_log>=0.95 & taxa_gwas$lambda_log<=1.06 & !is.na(taxa_gwas$lambda_log)
  id0b<-taxa_gwas$lambda_tp>=0.95 & taxa_gwas$lambda_tp<=1.06 & !is.na(taxa_gwas$lambda_tp)
taxa_gwas1<-taxa_gwas
taxa_gwas1$max_log[!id0a]<-NA
taxa_gwas1$max_tp[!id0b]<-NA
taxa_gwas1[,c("max_log","max_tp")]
max_pvals<-cbind(otu3=taxa_gwas1$otu3,
                 do.call(rbind,
                         apply(taxa_gwas1[,c("max_log","max_tp")],1,
                               function(x){
                                 if(all(is.na(x))) return(data.frame(mod="NA",mx=NA))
                                 if(is.na(x[1])) return(data.frame(mod="tp",mx=x[2]))
                                 if(is.na(x[2])) return(data.frame(mod="log",mx=x[1]))
                                 mx<-max(x,na.rm=T);if(mx==x[1])return(data.frame(mod="log",mx)) else return(data.frame(mod="tp",mx))})))

#max_pvals<-cbind(otu3=taxa_gwas$otu3[taxa_gwas$taxa_id%in%ids_for_rep],do.call(rbind,apply(taxa_gwas[taxa_gwas$taxa_id%in%ids_for_rep,c("max_log","max_tp")],1,function(x){mx<-max(x,na.rm=T);if(mx==x[1])return(data.frame(mod="log",mx)) else return(data.frame(mod="tp",mx))})))
max_pvals$otu3<-as.character(max_pvals$otu3)
max_pvals$mx1<-max_pvals$mx/sd(max_pvals$mx,na.rm=T)

write.table(cbind(max_pvals$otu3[max_pvals$mod=="log"],"ring_height","4",max_pvals$mx[max_pvals$mod=="log"]),"circular_plots/annot_3_gwas.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)
write.table(cbind(max_pvals$otu3[max_pvals$mod=="log"],"ring_color","4","#008B8B"),"circular_plots/annot_3_gwas.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)

write.table(cbind(max_pvals$otu3[max_pvals$mod=="tp"],"ring_height","4",max_pvals$mx[max_pvals$mod=="tp"]),"circular_plots/annot_3_gwas.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)
write.table(cbind(max_pvals$otu3[max_pvals$mod=="tp"],"ring_color","4","#00008B"),"circular_plots/annot_3_gwas.txt",quote = F,row.names = F,col.names = F,sep = "\t",append = T)





