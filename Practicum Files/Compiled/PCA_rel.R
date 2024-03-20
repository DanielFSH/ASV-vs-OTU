set.seed(123)
BiocManager::valid()
library(mixOmics)
library(tidyverse)
library(umap)
#Remove the allergy and otu columns before running the pca
asvs<- read.csv('relative-asv.csv')
otus<-read.csv('relative-otu.csv')

#removing the allergy and otuid columns
asv2 <-subset( asvs, select = -c(gid_16s, allergy) )
otu2<-subset( otus, select = -c(gid_16s, allergy) )

asvs.pca <- pca(asv2)     # 1 Run the method
plotIndiv(asvs.pca, title="PCA for ASVs")    # 2 Plot the samples
groups1<-asvs$gid_16s
plotIndiv(asvs.pca, group = groups1, 
          legend = FALSE)
#otus pca
otus.pca <- pca(otu2)     # 1 Run the method
plotIndiv(otus.pca, title = "PCA for OTUs")    # 2 Plot the samples
groups2<-otus$gid_16s
plotIndiv(otus.pca, group = groups2,  legend = FALSE)
#interesting results

#SPCA
asvs.spca<-spca(asv2, ncomp = 2)
plotIndiv(asvs.spca, title = "SPCA for ASVS")    # 2 Plot the samples
otus.spca<-spca(otu2, ncomp = 2)
plotIndiv(otus.spca, title="PCA for OTUs")    # 2 Plot the samples


#re-do with pseudo count
pseudootus<-read.csv('pseudo-otu.csv')
pseudootu2<-subset( pseudootus, select = -c(ï..OTU.ID) )
potus.spca<-spca(pseudootu2, ncomp = 2, logratio = 'CLR')
plotIndiv(potus.spca, title="PCA for OTUs")    # 2 Plot the samples

pseudoasvs<-read.csv('pseudo-asv.csv')
pseudoasv2<-subset( pseudoasvs, select = -c(ï..OTU.ID) )
pasvs.spca<-spca(pseudoasv2, ncomp = 2, logratio = 'CLR')
plotIndiv(pasvs.spca, title="SPCA for ASVs")    # 2 Plot the samples


#pseudootu3<-apply(pseudootu2,2,function(x){x/sum(x)})
#potus2.spca<-spca(pseudootu2, ncomp = 2, logratio = "CLR")
#plotIndiv(potus2.spca, title="SPCA for OTUs")    # 2 Plot the samples



#How to organzie the pca results
#N=20,50 or 100
#N=20
#n=floor(N/3)
#m=(N-2*n)/2
#ASVs 
PCA.subsampling<-function(PCA.components, N=20){
  n=floor(N/3)
  m=(N-2*n)/2
  a=as.data.frame(PCA.components$variates$X)
  OrderedPCA1=order(a$PC1)
  OrderedPCA2=order(a$PC2)
  PCA=c(OrderedPCA1[1:n],OrderedPCA1[(199-n+1):199])
  nonduplicate=OrderedPCA2[OrderedPCA2%in%PCA==FALSE]
  L=length(nonduplicate)
  PCA=c(PCA, nonduplicate[1:m], nonduplicate[(L-m+1):L])
  PCA.samples<-c()
  for (i in PCA) {
    PCA.samples<-c(PCA.samples, asvs$gid_16s[i])
  } 
  PCA.samples<-as.data.frame(PCA.samples)
  return(PCA.samples)
}
PCA.subsample.20<-PCA.subsampling(PCA.components = asvs.pca)
PCA.subsample.20<-cbind(PCA.subsample.20,PCA.subsampling(PCA.components = otus.pca))
PCA.subsample.20<-cbind(PCA.subsample.20,PCA.subsampling(PCA.components = asvs.spca))
PCA.subsample.20<-cbind(PCA.subsample.20,PCA.subsampling(PCA.components = otus.spca))

PCA.subsample.50<-PCA.subsampling(PCA.components = asvs.pca, N=50)
PCA.subsample.50<-cbind(PCA.subsample.50,PCA.subsampling(PCA.components = otus.pca, N=50))
PCA.subsample.50<-cbind(PCA.subsample.50,PCA.subsampling(PCA.components = asvs.spca, N=50))
PCA.subsample.50<-cbind(PCA.subsample.50,PCA.subsampling(PCA.components = otus.spca, N=50))

PCA.subsample.100<-PCA.subsampling(PCA.components = asvs.pca, N=100)
PCA.subsample.100<-cbind(PCA.subsample.100,PCA.subsampling(PCA.components = otus.pca, N=100))
PCA.subsample.100<-cbind(PCA.subsample.100,PCA.subsampling(PCA.components = asvs.spca, N=100))
PCA.subsample.100<-cbind(PCA.subsample.100,PCA.subsampling(PCA.components = otus.spca, N=100))
 
colnames(PCA.subsample.20) <- c("PCA-ASV", "PCA-OTU", "SPCA-ASV","SPCA-OTU")
colnames(PCA.subsample.50) <- c("PCA-ASV", "PCA-OTU", "SPCA-ASV","SPCA-OTU")
colnames(PCA.subsample.100) <- c("PCA-ASV", "PCA-OTU", "SPCA-ASV","SPCA-OTU")

write.csv(PCA.subsample.20, "PCA.20.csv")
write.csv(PCA.subsample.50, "PCA.50.csv")
write.csv(PCA.subsample.100, "PCA.100.csv")




#PLSDA
#which outcome? Impute the missing allergy variables
X <- asv2
Y <- asvs$allergy
table(is.na(Y))
asv.plsda<-plsda(X,Y)
summary(asv.plsda)


#UMAP
asvs.umap = umap(asv2)
otus.umap = umap(otu2)
umap_asv <-asvs.umap$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2")
#I can add metadata here
umap_asv2 <- cbind(umap_asv, asvs$ï..OTUID)
umap_asv2 <- cbind(umap_asv2, asvs$allergy)

umap_asv2 %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color= asvs$allergy))+
  geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot")
umap_asv2 %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2,
  )) +
  geom_point(size=3, alpha=0.5)+
  facet_wrap(~asvs$allergy)+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle="UMAP plot")+
  theme(legend.position="bottom")


#OTU UMAP
umap_otu <-otus.umap$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2")
#I can add metadata here
umap_otu2 <- cbind(umap_otu, otus$ï..OTUID)
umap_otu2 <- cbind(umap_otu2, otus$allergy)

umap_otu2 %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color= otus$allergy))+
  geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot")


umap_otu2 %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2,
  )) +
  geom_point(size=3, alpha=0.5)+
  facet_wrap(~otus$allergy)+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle="UMAP plot")+
  theme(legend.position="bottom")

#for Micropita
c=1:5
1:5%in%c(1,5)
#clustering algorithmS?
