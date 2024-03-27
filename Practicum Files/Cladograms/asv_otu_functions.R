library(tidyverse)
library(RColorBrewer)
#Functions to read tsv files and generate the annotation files

#NOTE:Changed the path to accomodate the folder being moved!

#The first function should read in the right tsv file,
#using fill so that it will grab all columns
#Type will be "asv" or "otu"
#Subsamplingmethod: Diverse, Distinct, Discriminant, Extreme, Representative, PCA, SPCA
#Size: 20, 50, 100
subsampler<-function(samplingmethod,type,size){
  #First open the asv files
  featuretable<- read.csv(paste0('relative_',type,'T.csv'))
  #need to import the LEfSe LDA files
  samplingtable<-read.table(paste0('LEfSe/',samplingmethod,type,size,'.res'),fill=TRUE,
                   col.names = c('feature', 'V2', samplingmethod,'V4','V5'))
  #Need to add a condition if the first row with values in cols V3,V4 is not 1 or 0
  for(i in 1:nrow(samplingtable[samplingmethod])){
    if (samplingtable[i,3]!=1 & samplingtable[i,3]!=0){
      samplingtable[i,3]<-'-'
    }
  }
  #Need to do a join somehow
  #Need to make the features consistent between tables
  #First remove periods
  #Use \\| to eascape the | character, turn them to periods
  #Replace .__ with ___
  featuretable$ï..gid_16s<-gsub("\\.", "_", as.character(featuretable$ï..gid_16s))
  featuretable$ï..gid_16s<-gsub("\\|", ".", as.character(featuretable$ï..gid_16s))
  featuretable$ï..gid_16s<-gsub("\\.__", "___", as.character(featuretable$ï..gid_16s))
  featuretable=featuretable %>% 
    rename(
      feature=ï..gid_16s
    )
  #The next step is to join it to the original asv/otu table
  #doing a left join with the asv table
  mergedtable=merge(x = featuretable, y = samplingtable, by = "feature", all.x = TRUE)
  methodtable=cbind(mergedtable[1],mergedtable[201:204])
  #Need to remove "allergy" as a clade
  return(methodtable[2:286,])
}
#the following has to become a function, that only runs once!
make_guide<-function(subsample_table,type){
#im going to use one of the above to remake the guide!
  features1<-strsplit(as.character(subsample_table$feature),".D_[12345]_")
  features1
#possible problems ".Ambiguous_taxa"
#and "uncultured.something"
  guide<-sapply(features1, function(s){
    len<-length(s)
    s[s=="_"]<-"_Other"
    s<-gsub("UnclassifiedUnclassified","Unclassified",
            gsub("Unclassifiedbacterium","Unclassified",
              gsub("Unclassifiedorganism","Unclassified",      
                 gsub("uncultured","Unclassified",
                      gsub("Ambiguoustaxa","Unclassified",
                           gsub("_","",gsub("[.]","-",s)))))))
    return(paste(s,collapse = "."))})
  guide<-as.data.frame(guide)
  #maybe make into a data frame instead of adding it back to an existing data frame?
  #let it be its own thing?
  #feature_table$guide<-features2
  #Write the guide to the current directory
  write.table(guide,paste0("guide",type,".txt"),quote = F,row.names = F,col.names = F,sep = "\t")
  #Return the guide to use for annotation
  return(guide)
  ###this will remake the guide for ASVs and OTUs at each sample size 20,50,100
  #maybe add another var to prevent it from being remade?
}
#The next step is to 'sink' or append the right changes to the annotation file
annot_append<-function(guide,subsample_table, level,color, color2,samplingmethod, type, size){
  #the following has to be part of the appending function! otherwise it'll get annoying
  sink(file = paste0("annot_",type,size,".txt"),append = T)
  cat(paste0("ring_internal_separator_thickness\t",level,"\t0.5
ring_separator_color\t",level,"\t#888888
ring_label\t", level,"\t",samplingmethod,"
ring_label_color\t",level,"\t",color,"
" ))
  sink()
  #identify the 1s and 0s
  over<-which(subsample_table[3]==1 & !is.na(subsample_table[3]))
  under<-which(subsample_table[3]==0 & !is.na(subsample_table[3]))
  #The overexpressed
  write.table(cbind(guide[over,],"ring_color", level,color),paste0("annot_",type,size,".txt"), quote = F,row.names = F,col.names = F,sep = "\t",append = T)
  #the following changes the transparency
  write.table(cbind(guide[under,],"ring_alpha",level,"0.75"),paste0("annot_",type,size,".txt"),quote = F,row.names = F,col.names = F,sep = "\t",append = T)
  #write.table(cbind(guide[over,],"ring_shape",level,"^"),paste0("annot_",type,size,".txt"),quote = F,row.names = F,col.names = F,sep = "\t",append = T)
  #the under expressed
  #write.table(cbind(guide[under,],"ring_color",level,color2),paste0("annot_",type,size,".txt"),quote = F,row.names = F,col.names = F,sep = "\t",append = T)
  #temp change with new colors
  write.table(cbind(guide[under,],"ring_color",level,color),paste0("annot_",type,size,".txt"),quote = F,row.names = F,col.names = F,sep = "\t",append = T)
  write.table(cbind(guide[under,],"ring_alpha",level,"0.40"),paste0("annot_",type,size,".txt"),quote = F,row.names = F,col.names = F,sep = "\t",append = T)
  #write.table(cbind(guide[under,],"ring_shape",level,"v"),paste0("annot_",type,size,".txt"),quote = F,row.names = F,col.names = F,sep = "\t",append = T)
  #Need an expression for if they're both empty
  if(identical(over, integer(0))&identical(over,integer(0))){
    over<-c(1)
    write.table(cbind(guide[over,],"ring_color",level,"#ffffff"),paste0("annot_",type,size,".txt"),quote = F,row.names = F,col.names = F,sep = "\t",append = T)
  }
}
#mistakes<-subsampler("Rep","otu","50")
#i need to remove allergy from the list
#the following is a function to change the short names to long names

#guide<-make_guide(mistakes,"otu")
#annot_append(guide,mistakes,6,fig_colors[6],fig_colors[6+7],"Representative","otu","50")

#over<-which(mistakes[3]==1 & !is.na(mistakes[3]))
#under<-which(mistakes[3]==0 & !is.na(mistakes[3]))
#if(identical(over, integer(0))&identical(over,integer(0))){
 # over<-mistakes[1,1]
  #write.table(cbind(guide[over,],"ring_color",level,"#ffffff"),paste0("annot_",type,size,".txt"),quote = F,row.names = F,col.names = F,sep = "\t",append = T)
#}
#Below is text im playing around with

#to write the text in the ring, need the folowing:
#ring_internal_separator_thickness	1	0.5
#ring_separator_color	1	#888888
#ring_label	1	Stool (prevalence)
#ring_label_color	1	#0000FF

#annotation_font_size\t7
#annotation_legend_font_size\t5

#ring_width\t1\t0.8
#ring_height\t1\t0.5


#samplingmethod,"\tclade_marker_color\t",color,"
#",samplingmethod,"\tclade_marker_size\t20
#",samplingmethod,"\tclass_label\t",samplingmethod,

#branch_bracket_width	0.0 (makes everyting square)
#branch_bracket_depth	0.75 (makes clades longer)
#clade_separation	.75 (it spreads out the clades, good)
#annotation_background_offset 2 (This made space between the rings and the clades I think)

#need to use grep to get sampling method, color, and level from just
#asv or otu and the sample size
#Need to make a list of chosen colors for the figure:
#library(RColorBrewer)
#par(mar=c(3,4,2,2))
#display.brewer.all()
#colorRampPalette(brewer.pal(9,"Blues"))(100)

#removed the following
# *\tannotation_background_separation\t0.001

fig_colors<-brewer.pal(8,"Dark2")
#fig_colors<-c("#cc00ff","#0044ff","#05f2ea","#00ff1c","#ffbb00","#ff0000","#854604",
#               "#a88ca7","#9199c4","#82aeb0","#82b087","#ffe9ab","#ff7d7d","#7d6552")
#There is an error in how the files are being read in, the PCA has the first instance of it

make_annotation<-function(type, size){
  #first get a list sampling methods
  allfiles<-list.files(path = "C:/Daniel/UofT/Practicum Files/Article/Cladograms/LEfSe",
                       pattern = NULL, all.files = TRUE,)
  annotation_type<-grep(paste0("*",type,size,".res"),allfiles, value = TRUE)
  method<-gsub(paste0(type,size,".res"),"",annotation_type)
  #should output this for method: "Discr" "Dist"  "Div"   "Ext"   "Pca"   "Rep"   "Spca" 
  #Need to change these into their long form
  #return(method)
  for (i in 1:7){
    #first it will subsample
    subsample_table<-subsampler(samplingmethod=method[i], type = type, size = size)
    #then it will make the guide once, and write the annotation file (in case some annotations are dependent on the guide)
    if(i==1){
      guide<-make_guide(subsample_table,type)
      ### ASV results  - circular plot annotation
      sink(file = paste0("annot_",type,size,".txt"))
      cat(paste0("title\t",toupper(type)," n=",size,"
title_font_size\t13
start_rotation\t90
total_plotted_degrees\t330
clade_marker_size\t5
clade_marker_edge_width\t0.3
clade_marker_color\t#176603
clade_separation\t0.5
annotation_background_offset\t0.1
annotation_legend_font_size\t5	
annotation_font_size\t8
*\tbranch_thickness\t0.5
*\tbranch_color_from_ancestor\t0
*\tbranch_bracket_depth\t0.75
*\tring_edge_width\t1\t0.0
ring_label_font_size	1	10
ring_label_font_size	2	10
ring_label_font_size	3	10
ring_label_font_size	4	10
ring_label_font_size	5	10
ring_label_font_size	6	10
ring_label_font_size	7	10
internal_label	2	Phy.
internal_label	3	Classes
internal_label	4	Orders
internal_label	5	Families
internal_label	6	Genera
internal_label_font_size\t8
internal_labels_rotation	90
"))
      sink()
      barheights=guide
    }
    #then it will make the annotations to the data frame
    annot_append(guide,subsample_table,level=i,fig_colors[i],fig_colors[i+7],method[i],type,size)
    barheights=cbind(barheights,subsample_table[3])
  }
  #guides <- list(guide,guide2)
  #temporary return both things to test somet stuff out
  #return(guides)
  ##if error delete from here
  ####
  ####
  
  #Need to add the values for each row
  #Need to make all cols into numeric
  for (j in 2:8) {
    barheights[j]<-as.numeric(unlist(barheights[j]))
  }
  #Make all the 0s into 1s?
  barheights[barheights==0]<-1
  barheights['Heights']<-rowSums(barheights[2:8], na.rm = TRUE)
  
  #Need to list all the rows with a sum value>0,
  #identify the ones with bars!
  bars<-which(barheights['Heights']>0)
  sink(file = paste0("annot_",type,size,".txt"),append = T)
  cat(paste0("ring_internal_separator_thickness\t8\t0.5
ring_separator_color\t8\t#000000
ring_label\t8\tDiference expressed
ring_label_color\t8\t",fig_colors[8],"
ring_color\t8\t",fig_colors[8],"
" ))
  sink()
  #the following will append it to the right file!
  write.table(cbind(barheights[bars,1],"ring_height",8,barheights[bars,9]),paste0("annot_",type,size,".txt"),quote = F,row.names = F,col.names = F,sep = "\t",append = T)
  #next need to annotate certain ones
  #Okay, now we isolate the clades with bars greater than 1
  annot<-which(barheights['Heights']>1)
  #Need to split the string so we can isolate the last clade of each asv/otu
  clades<-strsplit(barheights[annot,1],split='[.]')
  #If the clade is a genus, want it on the legend,
  #if its a lower clade, want it on the ring
  annotations<-c()
  for(i in 1:length(clades)){
    if(length(clades[[i]])>5){
      annotations<-c(annotations,paste0("*:",tail(clades[[i]],1)))
    }else{
      annotations<-c(annotations,tail(clades[[i]],1))
    }
  }
  annotations<-as.data.frame(annotations)
  #Then I need to write the annotation
  #each needs 3 lines!
  #the annotation text
  #Anaerococcus	annotation	Anaerococcus
  #The annotation background color
  #Anaerococcus	annotation_background_color	k
  #the annotation text size
  #Anaerococcus	annotation_font_size	14
  write.table(cbind(barheights[annot,1],"annotation",annotations),paste0("annot_",type,size,".txt"),quote = F,row.names = F,col.names = F,sep = "\t",append = T)
  write.table(cbind(barheights[annot,1],"annotation_background_color","k"),paste0("annot_",type,size,".txt"),quote = F,row.names = F,col.names = F,sep = "\t",append = T)
  #write.table(cbind(barheights[annot,1],"annotation_font_size",5),paste0("annot_",type,size,".txt"),quote = F,row.names = F,col.names = F,sep = "\t",append = T)
  return(barheights)
}

#Okay all the annotations look good?
#still need to make the dark shades a  bit transparent
#lefs<-make_annotation("otu","50")
types<-c("asv","otu")
sizes<-c("20","50","100")
#use a for loop to make all the annotation files
for (i in 1:2){
  for (j in 1:3){
    make_annotation(types[i],sizes[j])
  }
}
#great now to use Graphlan to check them



#Need to add the values for each row
#Need to make all cols into numeric
for (j in 2:8) {
  asv20methods[j]<-as.numeric(unlist(asv20methods[j]))
}
#Make all the 0s into 1s?
asv20methods[asv20methods==0]<-1
asv20methods['Sums']<-rowSums(asv20methods[2:8], na.rm = TRUE)

#Next need to "sink" something,
#Need to list all the rows with a value>0,
#What does the text need to say?
#The structure should be, the genus, ring_8?
#t645951833	ring_height	8	0.44330825
#identify the ones with bars!
bars<-which(asv20methods['Sums']>0)
cbind(asv20methods[bars,1],8,asv20methods[bars,9])

#the following will add bars to the right OTUs
write.table(cbind(asv20methods[bars,1],8,asv20methods[bars,9]),paste0("annot_",type,size,".txt"),quote = F,row.names = F,col.names = F,sep = "\t",append = T)

#And then give text for the height of the bar!

#time to make the large function
#It should use subsampler to grab the right file


#It will only look for the right "type" at the right sample "size"



#classtype will be 'relative_asvT.csv' or relative-otuT.csv', the transposed relative asv or otu files



#The next step is to 'sink' or append the right changes to the annotation file
#old version
second<-function(guide,temptable, level,color, annot_type,samplingmethod){
  #the following has to be part of the appending function! otherwise it'll get annoying
  sink(file = annot_type,append = T)
  cat(paste0(samplingmethod,"\tclade_marker_color\t",color,"
",samplingmethod,"\tclade_marker_size\t20
",samplingmethod,"\tclass_label\t",samplingmethod,"
"))
  sink()
  #identify the 1s and 0s
  c<-which(temptable[3]==1 & !is.na(temptable[3]))
  d<-which(temptable[3]==0 & !is.na(temptable[3]))
  #The overexpressed
  write.table(cbind(guide[c,6],"ring_color", level,color),annot_type,quote = F,row.names = F,col.names = F,sep = "\t",append = T)
  write.table(cbind(guide[c,6],"ring_shape",level,"^"),annot_type,quote = F,row.names = F,col.names = F,sep = "\t",append = T)
  #the under expressed
  write.table(cbind(guide[d,6],"ring_color",level,color),annot_type,quote = F,row.names = F,col.names = F,sep = "\t",append = T)
  write.table(cbind(guide[d,6],"ring_shape",level,"v"),annot_type,quote = F,row.names = F,col.names = F,sep = "\t",append = T)
}



?grep()
getwd()
a<-list.files(path = "C:/Daniel/UofT/Practicum Files/Article/Cladograms/LEfSe",
           pattern = NULL, all.files = TRUE,)
b<-grep("*asv*",a, value = TRUE)
c<-gsub("asv20_LDA.txt","",b)
c

?grep()



### ASV results  - circular plot annotation
sink(file = "annot_asv20.txt")
cat("title\tASVs
title_font_size\t13
start_rotation\t98
total_plotted_degrees\t355
clade_marker_size\t5
clade_marker_edge_width\t0.3
clade_marker_color\t#176603
*\tbranch_thickness\t0.5
*\tannotation_background_separation\t0.001
*\tbranch_color_from_ancestor\t0
*\tbranch_bracket_depth\t0.75
*\tring_edge_width\t1\t0.0
annotation_font_size\t7
annotation_legend_font_size\t10
ring_internal_separator_thickness\t1\t0.5
ring_width\t1\t0.8
ring_height\t1\t0.5
")
sink()
annot_append(lol,hi,3,"#556B2F","Discriminant","asv","20")
annot_append(lol,hi,4,"#34EBE1","Distinct","asv","20")




#regular expressions
#maybe color clades that are subsampled in all subsampling methods
#grep
?grep()



#goals for the day:
#remake the guide with the 327 asvs, make a function that will do all these steps automatically (done)
#finish function that will sink the right lines into the tsv file(done)
#make a function that can read the right file in, and then do all steps to make the annotation file!(done)


#for the extra annotation, where we change the color of a clade if its highlighted by many, 
#need a function that takes the V3 col from the tsvs, and checks the row to count the number of 1s or 0s,
#if its over a certain number, annotate that row??

round(a[6])


#First idea, use the guide column 1, and append to it subsample_table[3] each time, then return it?


#HMP example is great, use the rings!!! label them by subsample method, square of dark if overexpressed,
#square with light coloring if underexpressed!
#frees up the label/legend for actual clades and stuff such as the diff phyla or order!!!



#things to do:
#Histograms (the bars are present, need to change the color?)
#change pallette (maybe done?)
#Annotate Families (which ones!!!)
#Idea:Use the guide, split on the .,
#grab unique combination of first and second components? 
#OR annotate those where sum col is >1??
#lets try this!
bars<-which(barheights['Heights']>0)

write.table(cbind(barheights[bars,1],"ring_height",8,barheights[bars,9]),paste0("annot_",type,size,".txt"),quote = F,row.names = F,col.names = F,sep = "\t",append = T)


#Okay, now we isolate the clades with bars greater than 1
bars<-which(a['Heights']>1)
d<-a[bars,]
c<-c()
b<-strsplit(a[bars,1],split='[.]')
for(i in 1:length(b)){
  c<-c(c,tail(b[[i]],1))
  }
c<-paste0("*:",c)
as.data.frame(c)
c[!is.na(c)]
as.data.frame(c)

#^that became this:
bars<-which(barheights['Heights']>1)
annotations<-c()
clades<-strsplit(barheights[bars,1],split='[.]')
for(i in 1:length(clades)){
  annotations<-c(annotations,tail(clades[[i]],1))
}
annotations<-as.data.frame(annotations)
#But it wasnt working for now
#Okay it works but looks awful
#Next idea, if its a genus, put a letter and make a legend,
#If its any higher clade, annotate it

j=1
c<-c()
l<-c()
for(i in 1:length(b)){
  if(length(b[[i]])>5 && j<=26){
    c<-c(c,letters[j])
    l<-c(l,tail(b[[i]],1))
    j<-j+1
  }else if(length(b[[i]])>5 && j>26){
    c<-c(c,LETTERS[j-26])
    l<-c(l,tail(b[[i]],1))
    j<-j+1
  }else{
    c<-c(c,tail(b[[i]],1))
  }
}
c
as.data.frame(l)
#which became
j=1
legends<-c()
for(i in 1:length(clades)){
  if(length(clades[[i]])>5 && j<=26){
    annotations<-c(annotations,letters[j])
    legends<-c(legends,tail(b[[i]],1))
    j<-j+1
  }else if(length(clades[[i]])>5 && j>26){
    annotations<-c(annotations,LETTERS[j-26])
    legends<-c(legends,tail(b[[i]],1))
    j<-j+1
  }else{
    annotations<-c(annotations,tail(b[[i]],1))
  }
}
annotations<-as.data.frame(annotations)
legends<-as.data.frame(legends)
#Trying something else for now

#Then I need to write the legend 

#Then I need to write the annotation
#each needs 3 lines!
#the annotation text
#Anaerococcus	annotation	Anaerococcus
#The annotation background color
#Anaerococcus	annotation_background_color	k
#the annotation text size
#Anaerococcus	annotation_font_size	14
write.table(cbind(c,"annotation",c),paste0("annot_",type,size,".txt"),quote = F,row.names = F,col.names = F,sep = "\t",append = T)
write.table(cbind(c,"annotation_background_color","k"),paste0("annot_",type,size,".txt"),quote = F,row.names = F,col.names = F,sep = "\t",append = T)
write.table(cbind(c,"annotation_font_size",14),paste0("annot_",type,size,".txt"),quote = F,row.names = F,col.names = F,sep = "\t",append = T)

#the above should work? replace the c with something

#LEfSe(download it, do in compute node?)
#change the green, 
#try playing with transparency 
#change the title
#metaphlan 2.8 on the alliance
#module load cc n

#Okay so the annotations are done,
#Tmrw, play with the colors, the size, and the transparancy
#get LEfSe to work!

#if its a family 



#The LEfSe table comes out at .res, does it work?

lefs<-read.table(paste0('LEfSe/Repotu50.res'),fill=TRUE,
                          col.names = c('feature', 'V2', 'samplingmethod','V4','V5'))
#it will work no problem it seems! :D

al<-list.files(path = "C:/Daniel/UofT/Practicum Files/Article/Cladograms/LEfSe",
                       pattern = NULL, all.files = TRUE,)
an<-grep(paste0("*asv20.res*"),al, value = TRUE)
gsub(paste0("asv","20",".res"),"",an)
subsampler("Rep","otu","50")



#some taxa are "unknown" or "ambiguous" how to deal with those?
allfiles<-list.files(path = "C:/Daniel/UofT/Practicum Files/Article/Cladograms/LEfSe",
                     pattern = NULL, all.files = TRUE,)
annotation_type<-grep(paste0("*asv20.res"),allfiles, value = TRUE)
method<-gsub(paste0("asv20.res"),"",annotation_type)
size<-c("20","50","100")
type<-c("asv","otu")
newnames<-c("Discriminant","Distinct","Extreme","Diverse","PCA","Representative","SPCA")
#Changing the file names
for (i in 1:7) {
  for (j in 1:3) {
    for (k in 1:2) {
      file.rename(paste0(method[i],type[k],size[j],".res"),
                  paste0(newnames[i],type[k],size[j],".res")) 
    }
  }
}

#July 10th Fixing some of the odd text in the guides
# unculture, uncilturedbacterium ,unculturedorganism, Ambiguoustaxa, Unknownfamily

#Okay so the make guide wansnt replacing Ambiguous_taxa, due to the underscore being removed
#The following does replace it, but in some cases Unclassified is there twice
#Need other stuff to replace with unclassified
a<-read.delim("LEfSe/Discriminantasv20.res", header = FALSE)
taxas<-strsplit(as.character(a$V1),".D_[12345]_")
guide<-sapply(taxas, function(s){
  len<-length(s)
  s[s=="_"]<-"_Other"
  s<-gsub("Unclassifiedbacterium","Unclassified",
     gsub("Unclassifiedorganism","Unclassified",      
     gsub("uncultured","Unclassified",
     gsub("Ambiguoustaxa","Unclassified",
     gsub("_","",gsub("[.]","-",s))))))
  return(paste(s,collapse = "."))})
guide<-as.data.frame(guide)

#The above works looks like, will replace make_guide, but old version is below

make_guide<-function(subsample_table,type){
  #im going to use one of the above to remake the guide!
  features1<-strsplit(as.character(subsample_table$feature),".D_[12345]_")
  features1
  #possible problems ".Ambiguous_taxa"
  #and "uncultured.something"
  guide<-sapply(features1,function(s){
    len<-length(s)
    s[s=="_"]<-"_Other"
    s<-gsub("Ambiguous_taxa","Unclasified",gsub("_","",gsub("[.]","-",s)))
    return(paste(s,collapse = "."))})
  guide<-as.data.frame(guide)
  #maybe make into a data frame instead of adding it back to an existing data frame?
  #let it be its own thing?
  #feature_table$guide<-features2
  #Write the guide to the current directory
  write.table(guide,paste0("guide",type,".txt"),quote = F,row.names = F,col.names = F,sep = "\t")
  #Return the guide to use for annotation
  return(guide)
  ###this will remake the guide for ASVs and OTUs at each sample size 20,50,100
  #maybe add another var to prevent it from being remade?
}


#Next need to adjust, make sure it works with Metaphlan results

