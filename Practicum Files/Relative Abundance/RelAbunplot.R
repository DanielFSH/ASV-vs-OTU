library(qiime2R)
library(ggplot2)
library(dplyr)
#this worked
asv_phylo<-qza_to_phyloseq(features="asv-karelia-comp-table.qza",
                   taxonomy="asv-karelia-taxonomy.qza", 
                   tree="asv_rooted_tree.qza",
                   metadata = "karelia-metadata.tsv")
#The following gets the relative abundance counts
asv_rel_abund <- transform_sample_counts(asv_phylo, function(x) {x / sum(x)})
#Trying to get the average relative abundance


#GPf = tax_glom(GlobalPatterns, "Family") %>% transform_sample_counts(function(x) {x * 100/sum(x)})
# The taxa orders are the same in tax_table() and otu_table()
# Use 'rowMeans(otu_table(GPf))' to calculate per-row average in the OTU table
asv_avgrel = data.frame(Order = tax_table(asv_rel_abund)[,"Order"], Mean = rowMeans(otu_table(asv_rel_abund)), row.names = NULL)
asv_avgrel = asv_avgrel[order(-asv_avgrel$Mean),]
head(asv_avgrel)
#ggplot(df, aes(fill=Phylum, y="Mean", x="")) + 
 # geom_bar(position="stack", stat="identity")
#The following adds a "Type" column to differentiate ASVs from OTUs
asv_avgrel$Type<-"ASV"

#Repeat for OTU
otu_phylo<-qza_to_phyloseq(features="otu-karelia-comp-table.qza",
                           taxonomy="otu-karelia-taxonomy.qza", 
                           tree="otu_rooted_tree.qza",
                           metadata = "karelia-metadata.tsv")
#The following gets the relative abundance counts
otu_rel_abund <- transform_sample_counts(otu_phylo, function(x) {x / sum(x)})
#Trying to get the average relative abundance


#GPf = tax_glom(GlobalPatterns, "Family") %>% transform_sample_counts(function(x) {x * 100/sum(x)})
# The taxa orders are the same in tax_table() and otu_table()
# Use 'rowMeans(otu_table(GPf))' to calculate per-row average in the OTU table
otu_avgrel = data.frame(Order = tax_table(otu_rel_abund)[,"Order"], Mean = rowMeans(otu_table(otu_rel_abund)), row.names = NULL)
otu_avgrel = otu_avgrel[order(-otu_avgrel$Mean),]

head(otu_avgrel)
#names(otu_avgrel)<-c("Class","Mean")

#add the otu column
otu_avgrel$Type<-"OTU"

#I then rbidn the two datasets whihc have the same columns
df<-rbind(asv_avgrel,otu_avgrel)
#Plot the figure
ggplot(df, aes(fill=Order, y=Mean, x=Type)) + 
  geom_bar(position="stack", stat="identity")+
  labs(y = "Relative Abundance", title = "Order Relative Abundance")

#Testing successful, can now make figure 
ggplot(df, aes(x = Type, y = "Abundance", fill = Phylum)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(y = "Relative Abundance", title = "Phylum Relative Abundance")





psmelt(df) %>%
  ggplot(aes(x = Type, y = Mean, fill = Phylum)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(y = "Relative Abundance", title = "Phylum Relative Abundance")

nonMatch_Uniquedf1 <- asv_avgrel %>% 
  filter(!asv_avgrel$Order %in% otu_avgrel$Order) #For df1 values not in df2
nonMatch_Uniquedf2 <- otu_avgrel %>% 
  filter(!otu_avgrel$Order %in% asv_avgrel$Order) #For df2 values not in df1

ggplot(nonMatch_Uniquedf1, aes(x="", y=Mean, fill=Order)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  labs(y = "", title = "Distinct clades present in ASV")

nonMatch_Uniquedf1
nonMatch_Uniquedf2
#adding a row to nonMath_uniquedf1 to display in the figure

heh<-data.frame(NA,0,"OTU")
names(heh)<-c("Order","Mean","Type")
heh2 <- rbind(nonMatch_Uniquedf1, heh)

#theres no OTUs not in ASVs, plot a bar plot of nonMatch
ggplot(df, aes(fill=Order, y=Mean, x=Type)) + 
  geom_bar(position="stack", stat="identity")+
  labs(y = "Relative Abundance", title = "Order Relative Abundance")
#previous version
ggplot(nonMatch_Uniquedf1, aes(x="", y=Mean, fill=Order)) +
  geom_bar(position="stack", stat="identity")+
  labs(y = "Relative Abundance", title = "Relative Adbundance of Distinct Order Clades")

ggplot(heh2, aes(x=Type, y=Mean, fill=Order)) +
  geom_bar(position="stack", stat="identity")+
  labs(y = "Relative Abundance", title = "Relative Adbundance of Distinct Order Clades")
nonMatch_Uniquedf1
#the above makes the desire plot


#to get correct average relative abundances need to sum them together
asv_sum<-aggregate(asv_avgrel$Mean, by=list(Order=asv_avgrel$Order), FUN=sum)
#Renaming columns
colnames(asv_sum)<-c("Order", "Mean")
#now need to add the NA relative abundances
asv_sum<-rbind(asv_sum,c("Unclassified",sum(asv_avgrel[is.na(asv_avgrel),2])))
nonMatch_Uniquedf3 <- asv_sum %>% 
  filter(!asv_sum$Order %in% otu_avgrel$Order) #For df1 values not in df2

ggplot(nonMatch_Uniquedf3, aes(x="", y=as.numeric(Mean), fill=Order)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  geom_label(aes(label = round(as.numeric(Mean)*100, 3)),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE, size=3)+
  theme_void()+
  labs(y = "Average Relative Abundance",
       title = "Distinct clades present in ASV")+
  scale_fill_viridis_d()






#okay need to add percents, need to re add x or y axis label
  
geom_label(aes(label = round(as.numeric(Mean)*100, 3)),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE)
  
  geom_text(aes(x = c(1, 1, 1, 1, 1, 1, 1, 1,1,1,1), # Solution for part 2,
                y = c(1,2,3,4,5,6,7,8,9,10,1), 
                label=round(as.numeric(Mean)*100, 3)))+




unique(asv_avgrel$Order)
unique(otu_avgrel$Order)
#the code above works to get the average of each 

#properly installed qiime2R
taxa_names("")