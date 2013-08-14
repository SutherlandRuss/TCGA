# script to identify hypermutated samples

library("scales")
plot(mutationsPerIndiv[,2],mutationsSilentPerIndiv)
abline(0,1, col="red")

cor(mutationsPerIndiv[,2], mutationsSilentPerIndiv,method= "pearson")
abline(0,cor(mutationsPerIndiv[,2], mutationsSilentPerIndiv,method= "pearson"),col="blue")

#barplot of SNV and indel frequency ordered by SNV frequency

snvIndelBarPlotPoints<- cbind(mutationsPerIndiv[,2],mutationsSilentPerIndiv)[order(mutationsPerIndiv[,2], decreasing =TRUE),]
barplot(mutationsPerIndiv[order(mutationsPerIndiv[,2], decreasing=TRUE),2], border="gray")

silentVsNonSilentRatio<- sapply(seq(1:length(intersectSamples)), function(x) (snvIndelBarPlotPoints[x,2]/sum(snvIndelBarPlotPoints[x,]))*100)
#1  indelSNVratio<- sapply(seq(1:length(uniqueSamples)), function(x) ((snvIndelBarPlotPoints[x,1]+1)/(snvIndelBarPlotPoints[x,2]+1)))
#2 indelSNVratio<- apply(snvIndelBarPlotPoints+1,1, function(x) x[1]/x[2])
indelSNVratio<- (snvIndelBarPlotPoints[,2])/(snvIndelBarPlotPoints[,1])


#plotbycancerType
snvIndelBarPlotPoints<-snvIndelBarPlotPoints[which(cancerType=="rectum"),]
indelSNVratio<-indelSNVratio[which(cancerType=="rectum")]
silentVsNonSilentRatio<-silentVsNonSilentRatio[which(cancerType=="rectum")]



#plot
par(mar=c(3,5,5,5))
barplot(t(snvIndelBarPlotPoints), col=c("blue", "red"),border=c("blue","red"), beside=FALSE,axisnames=FALSE, main= "non-silent mutation frequency and Indel:SNV ratio")# red = silent mutations and blue=nonsilent mutations
par(new=T)
plot(indelSNVratio, axes=F, pch=20, cex=0.7,xlab="",ylab="") # a plot of the percentage of indels amongst all non-silent mutations
axis(4, pretty(c(0, max(indelSNVratio))), pos= length(indelSNVratio)+2)
mtext("mutation frequency per sample", side=2, line=0, adj=0.5, padj=-3.5)
mtext("ratio of indels to SNV mutations",side=4,line=0, adj=0.5, padj=2.5)
mtext("samples ordered by mutation frequency",side=1, adj=0.5,padj=2.0)
legend(0.8*length(indelSNVratio),0.5, c("SNV","indel"), lty=c(1,1),lwd=c(2.5,2.5),col=c("blue","red"))
length(which(silentVsNonSilentRatio>20))
#hypermutatedSampleNames<-rownames(silentAndNonSilentMutations[which(silentVsNonSilentRatio>20),])# I need this because the names are ordered differently in 
# colnames(mutM) and rownames(silentAndNonSilentMutations)
#mutMHMRemoved<- mutM[,-match(hypermutatedSampleNames,colnames(mutM))]# use this to generate new plots
hyperIndex<-match(hypermutatedSampleNames,colnames(mutM))

#sampleColors[match(hypermutatedSampleNames, colnames(mutM))]<- "green" #



#plotting the mutation frequency of each sample


silentMutations<-scores[which(scores$Variant_Classification=="Silent"),]

silentAndNonSilentMutations <-cbind(table(silentMutations$sampleID), table(mutations$sampleID))
silentAndNonSilentMutations<-silentAndNonSilentMutations[order(silentAndNonSilentMutations[,2], decreasing=TRUE),]
hypermutatedSampleNames<-rownames(silentAndNonSilentMutations[which(silentVsNonSilentRatio>20),])# I need this because the names are ordered differently in 
hyperIndex<-match(hypermutatedSampleNames,colnames(mutM))


plot(silentAndNonSilentMutations[,2], log="y", ylim=c(1, max(silentAndNonSilentMutations[,2])), col="red", pch=20,cex=0.5, xaxt="n",
     main="total number of mutations per sample \n ordered from largest to smallest number of non-silent mutations")
points(silentAndNonSilentMutations[,1], col="blue", pch=20, cex=0.5)
axis(1,at=1:length(dimnames(silentAndNonSilentMutations)[[1]]), labels=(dimnames(silentAndNonSilentMutations)[[1]][1:length(dimnames(silentAndNonSilentMutations)[[1]])]), las=2, cex.axis=0.3)
legend(240,3000,c("non-silent mutations","silent mutations"),cex=1, pch=20, col=c("red","blue"))

silentMutationsForPlot<-sort(table(silentMutations$sampleID), decreasing=TRUE)
mutationsForPlot<- sort(table(mutations$sampleID), decreasing=TRUE)#identifying the hypermutated samples.
plot(mutationsForPlot, log = "y",ylim=c(1, max(mutationsForPlot[1])), col= "red", pch=20, xaxt="n",axis(1,at = mutationsForPlot,labels= dimnames(mutationsForPlot)[[1]]))
points(silentMutationsForPlot, col="blue", pch =20)




## mds plot of sample similarities BEFORE processing using the networ_informed_clustering_function script
# mutationMatrixLogicalOrdered is used as the input to the count match algorithm

t<- countMatch1(mutM)
u<- compDiss(t,1,mutM)
v<-cmdscale(u,k=2, eig=TRUE)
par(mar=c(5,5,5,25))
par(xpd=TRUE)
plot(v$points[,1],v$points[,2], pch = 19, col=alpha(smpclrs[,1],0.5), cex=1.0, cex.lab= 1.5, cex.main = 1.5, cex.axis=1.5, xlab="principal co-ordinate 1", ylab="principal co-ordinate 2", main ="mds pco plot of colorectal cancer samples \n before network processing of mutation matrix")
points(v$points[hyperIndex,1], v$points[hyperIndex,2],pch=1, cex = 1.5)
#legend(0.5, 0.11, c("Illumina", "Solid","hypermutated"), cex=1., pch=c(19,19,1),col=c("red","blue","black"))
legend(0.75, max(v$points[,2]+0.009), unique(smpclrs[,2]), cex=1.2,pt.cex=1.0, pch=c(rep(19, length(unique(smpclrs[,2])))),col=seq_along(1:length(unique(smpclrs[,1]))))

#start my loop to create all possible colorings of the MDS plot.look tongiht at how to structure a tapply function
# iterate over all possible colourings

##function of the code above##

MDSplot<- function(mdsPoints,smplclrs,hyperMutIndex){
  par(mar=c(5,5,5,25))
  par(xpd=TRUE)
  plot(mdsPoints$points[,1],mdsPoints$points[,2], pch = 19, col=alpha(smplclrs[,1],0.5), cex=1.2, cex.lab= 1.5, cex.main = 1.5, cex.axis=1.5, xlab="principal co-ordinate 1", ylab="principal co-ordinate 2", main ="mds pco plot of colorectal cancer samples \n before network processing of mutation matrix")
  points(mdsPoints$points[hyperIndex,1], mdsPoints$points[hyperIndex,2],pch=1, cex = 1.9)
  #legend(0.5, 0.11, c("Illumina", "Solid","hypermutated"), cex=1., pch=c(19,19,1),col=c("red","blue","black"))
  legend(0.72, max(v$points[,2]+0.009), c(unique(smplclrs[,2]),"hypermutated"), cex=1.5,pt.cex=1.2, pch=c(rep(19, length(unique(smplclrs[,2]))),1),col=c(seq_along(1:length(unique(smplclrs[,1]))),1))
}

for (i in seq(35,length(names(stratification)))){
  pdf(paste(i,"_",names(stratification[i]),".pdf", sep=""), width=15, height=10)
  MDSplot(v,stratification[[i]],hyperIndex)
  dev.off()
}
##some of the stratification variables have NA values. Investigate this.

plot(v$points[,1],v$points[,2], type = "n", xlab="principal co-ordinate 1", ylab="principal co-ordinate 2", main ="mds pco plot of coloectal cancer samples \n before network processing of mutation matrix")
text(v$points[,1],v$points[,2], labels = as.character(1:nrow(v$points)), col = sampleColors, cex = 0.5, xlab="principal co-ordinate 1", ylab="principal co-ordinate 2", main ="mds pco plot of coloectal cancer samples \n before network processing of mutation matrix")
points(v$points[hyperIndex,1], v$points[hyperIndex,2],pch=1, cex = 2)

plot(e[,1],e[,2], type = "n", xlab="principal co-ordinate 1", ylab="principal co-ordinate 2", main ="mds pco plot of coloectal cancer samples \n before network processing of mutation matrix")
text(e[,1],e[,2], labels = as.character(1:nrow(e)), col = sampleColors, cex = 0.5, xlab="principal co-ordinate 1", ylab="principal co-ordinate 2", main ="mds pco plot of coloectal cancer samples \n before network processing of mutation matrix")
points(e[hyperIndex,1], e[hyperIndex,2],pch=1, cex = 2)


## mds plot of sample similarities after processing using the networ_informed_clustering_function script
# c= the compDiss derived disimilarity matrix of samples
e<-cmdscale(c,k=2)
par(mar=c(5,5,5,5))
plot(e[,1],e[,2], col=sampleColors, pch =19,cex=0.8, xlab="principal co-ordinate 1", ylab="principal co-ordinate 2", main ="mds pco plot of colorectal cancer samples \n after network processing of mutation matrix")
points(e[hyperIndex,1], e[hyperIndex,2],pch=1, cex = 1.5)
legend(0.4, 0.15, c("Illumina", "Solid","hypermutated"), cex=1, pch=c(19,19,1),col=c("red","blue","black"))

##AFTER REMOVING HYPERMUTATED SAMPLES

## mds plot of sample similarities BEFORE processing using the networ_informed_clustering_function script hypermutated samples removed
# c= the compDiss derived disimilarity matrix of samples
j<- countMatch2(mutMHMRemoved)
k<- compDiss(j,1,mutMHMRemoved)
l<-cmdscale(k,2)
plot(l[,1],l[,2], col=sampleColors, cex=0.5, xlab="principal co-ordinate 1", ylab="principal co-ordinate 2", main ="mds pco plot of coloectal cancer samples \n before network processing of mutation matrix\n samples over 20% indels removed")
legend(0.5, 0.08, c("colon", "rectum"), cex=1, pch=1,col=c("red","blue"))

## mds plot of sample similarities AFTER processing using the networ_informed_clustering_function script hypermutated samples removed
# c= the compDiss derived disimilarity matrix of samples
m<- netInf(mutMHMRemoved,geneAdj)
n<- countMatch1(m)
p<- compDiss(n,1,mutMHMRemoved)
q<-cmdscale(p,2)
plot(q[,1],q[,2], col=sampleColors, cex=0.5, xlab="principal co-ordinate 1", ylab="principal co-ordinate 2", main ="mds pco plot of coloectal cancer samples \n after network processing of mutation matrix\n samples over 20% indels removed")
legend(0.5, 0.08, c("colon", "rectum"), cex=1, pch=1,col=c("red","blue"))


#do we now want to use partition around medoids basd on the number of clusters I think we have from the MDS plot?



#table giving the crosstab of two metadata variables
 crossT<-table(cbind(cancerType[,2], cancerType[,1], seqTech[,1])[,2],cbind(cancerType[,2], cancerType[,1], seqTech[,1])[,3])
 
#identify the genes driving the difference between sequencing technologies using fishers exact test.
# Generate a table for each gene (row of mutMatLogicalOrdered) and the sequencing machine label from seqTech.

table(mutMatLogicalOrdered[1,],seqTech[,1])

crossTB<- lapply(seq(1:length(mutMatLogicalOrdered[,1])), function(x) table(mutMatLogicalOrdered[x,], seqTech[,1]))
names(crossTB)<- rownames(mutMatLogicalOrdered)

FtestsResult<- matrix(NA, nrow=length(mutMatLogicalOrdered[,1]), ncol=1, dimnames=list(rownames(mutMatLogicalOrdered),"p_value"))
Ftests<- lapply(seq(1:length(crossTB)), function(x) fisher.test(matrix(crossTB[[x]],2,2))[[1]])
Ftests<-unlist(Ftests)
names(Ftests)<- rownames(mutMatLogicalOrdered)
#Ftests<-Ftests[order(Ftests,decreasing=FALSE)]
Ftests[which(Ftests<0.05/(length(Ftests)))]# this lets me know if there is a bias between sequencing technologies in the number of individuals carrying at least one mutation in a gene.
crossTB[which(Ftests<0.05/(length(Ftests)))]


#test to see if there are more total mutatuions identified by Illumina than Solid
fisher.test(table(mutations[,1], mutations$Sequencer))

genesBySequencer<-table(mutations$Hugo_Symbol,mutations$Sequencer)

#testforsequencer bias
#not sure which test to use as they give different p-values
binom.test(genesBySequencer[1,], p= 1/(sum(genesBySequencer[1,])/2))
binom.test(genesBySequencer[2,], n=sum(genesBySequencer[2,]))


#non-metricMDS doesn't work properly yet. sample165 and 187 are identical...
isoMDS(u, y = cmdscale(u, 2), k = 2, maxit = 50, trace = TRUE,
       +        tol = 1e-3, p = 2)

