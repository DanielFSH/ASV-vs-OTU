#Venn Diagrams figure! or Upsetplots!
library(UpSetR)
#following the example given on website
Subsamples <- read.csv("Subsampleselection.csv")
#samplingmethods<-names(Subsamples[2:43])
#upset(Subsamples, nsets = 14)
#Okay to get meaningful info, need only some columns.
#maybe the same method?
#Okay, so this loop makes a table with only Subsamples[rep]
for (j in seq(2,14, by=2)) {
methodtype<-c()
l<-j
for(i in 1:3){
  methodtype<-c(methodtype,l,l+1)
  l<-l+14
}
toupset<-Subsamples[methodtype]
print(upset(toupset, nsets = 6))
}
rm(i,j,l)

#Next build a larger loop which will print each type of upset plot
Subsamples[2:3]

##use grep
names(Subsamples)<- gsub("Rep","Representative",names(Subsamples))
names(Subsamples)
#compare across methods in the same samples size
#maybe for supplemental

#representative by batches
#I need a function that replaces the shorthand names with the long names
#can make it into one loop? yes im gonna neeed a loop with the names
shorttolong<-function(shortform){
  names(shortform)<- gsub("Rep","Representative",names(shortform))
  names(shortform)<- gsub("Ext","Extrme",names(shortform))
  names(shortform)<- gsub("Div","Diverse",names(shortform))
  names(shortform)<- gsub("Disc","Discriminant",names(shortform))
  names(shortform)<- gsub("Dist","Distinct",names(shortform))
  return(shortform)
}


for (j in seq(2,14, by=2)) {
  methodtype<-c()
  l<-j
  for(i in 1:3){
    methodtype<-c(methodtype,l,l+1)
    l<-l+14
  }
  toupset<-Subsamples[methodtype]
  print(upset(toupset, nsets = 6))
}


#Okay need to make two functions, one to rename cols from the
#shortnames to the long names, and one to make the upset plots

shorttolong<-function(shortform,short,long){
  names(shortform)<- gsub(short,long,names(shortform))
  return(shortform)
}
Upsetplotter<-function(shortform,short,long){
  newnames<-shorttolong(shortform,short,long)
  grouped<-grep(long,names(newnames))
  print(upset(newnames[grouped], nsets = 6))
}
slnames<-c("Rep","Ext","Div","Dist","Disc","Pca","Spca",
         "Representative","Extrme","Diverse","Distinct","Discriminant",
         "PCA","SPCA")
for (i in 1:7) {
  Upsetplotter(Subsamples,slnames[i],slnames[i+7])
  
}
shorttolong(Subsamples,"Pca","PCA")

#make the text larger
#try to make labels perpendicular