# The code needed to create a gene adjacency list object, representing the neighborhood structure of gene interations within a network, and a logical
# mutation matrix of m genes and n individuals, with TRUE>=1 non-silent mutation in gene m and individual m and FALSE indicating otherwise.
#
#The geneAdjIndex object, the mutMatLogicalOrdered object, and the cancerType object are used as arguments in the network informed clustering function.
#

library(igraph)

# The data file in vcf-like format.
#scores.path <- "/Users/Russ/Dropbox/PhD/tumour_classifier_data/colorectal_somatic_mutations/combined"
#scores.path <- "C:/Users/rds/Dropbox/PhD/tumour_classifier_data/colorectal_somatic_mutations/combined"
scores.path <- "C:/Users/rsutherland/Dropbox/PhD/tumour_classifier_data/colorectal_somatic_mutations/combined"

scores.file <- "colorectalcancer.maf"
#scores.file <- basename("C:/Users/rsutherland/Dropbox/PhD/tumour_classifier_data/colorectal_somatic_mutations/combined/colorectalcancer.maf")

dataFile<-file.path(scores.path,scores.file)

#The function to read the data in to the scores table
createScoreTable<- function(dataFile){
  scores<- read.table(dataFile , header=TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
  scores<- unique(scores)
  #Splitting the TCGA barcode in to the minimum  number of fields to uniquely identify individuals
  samples<- strsplit(scores$Tumor_Sample_Barcode, split = "-", fixed = TRUE)
  # Unique sample IDs
  #sampleID<- matrix(sapply(1:length(samples),function(x) paste0(samples[[x]][1], samples[[x]][2], samples[[x]][3])))
  sampleID<- unlist(lapply(samples,function(x) paste(x[1:3],collapse="")), use.names=FALSE)
  
  #adding the sampleID to the geneScore dataframe
  scores <-cbind(scores, sampleID)
  return(scores)
}

scores<-createScoreTable(dataFile)

# in order to analyse each sequencer's sample separately
#muts<-split(scores, scores$Sequencer)
#scores<- muts[[1]]
#scores<- muts[2]

#the list of unique samples
uniqueSamples<- unique(scores$sampleID)
#attach unique sampleIDs to each mutation

# gives me the matrix of "functional mutations"
# requires a raw "scores" table as input, like you would get from the createScoreTable function
getFuncMutations<-function(scoresT){
  mutations<-scoresT[which(scores$Variant_Classification!="Silent"),]
  colnames(mutations)<- colnames(scoresT)
  return (mutations)
}


mutations<-scores[which(scores$Variant_Classification!="Silent"),]
colnames(mutations)<- colnames(scores)



# nonsilent mutations per individual
mutationsTypePerIndiv<-table(mutations$Variant_Classification, mutations$sampleID)
mutationsSilentPerIndiv<-sapply(seq(1:length(uniqueSamples)), function(x) sum(mutationsTypePerIndiv[1:4,x]))
names(mutationsSilentPerIndiv)<-colnames(mutationsTypePerIndiv)

mutationsPerIndiv <-aggregate(mutations[,1], by=list(mutations$sampleID), length)
mutationsPerGene <-aggregate(mutations[,1], by=list(mutations$Hugo_Symbol), length)


# The matrix of "non-silent mutations"
mutations<-scores[which(scores$Variant_Classification!="Silent"),]
colnames(mutations)<- colnames(scores)

mutationTable<- table(mutations$Hugo_Symbol, mutations$sampleID)
ranksGenesTest<- rownames(mutationTable)

#####################################################################################################################################

mutationMatrix<-matrix(mutationTable,length(mutationTable[,1]),length(mutationTable[1,]))# This is my mutation matrix to be used to calculate sample distances

rownames(mutationMatrix)<- rownames(mutationTable)
colnames(mutationMatrix)<- colnames(mutationTable)
mutationMatrixBinary<- mutationMatrix

mutationMatrixBinary[which(mutationMatrix>1)]<-1# reassigns genes carrying more than one mutation to 1.
# all genes carrying more than one mutation in each individual ahve been set to one. This matrix
# now indicates if an individual carries at least one mutation per gene.

# logical mutation matrix all 1s converted to TRUE and all 0s converted to FALSE.
mutationMatrixLogical <- matrix(as.logical(mutationMatrixBinary),length(ranksGenesTest),length(uniqueSamples))
rownames(mutationMatrixLogical)<- rownames(mutationTable)
colnames(mutationMatrixLogical)<- colnames(mutationTable)



###########################################################################################################################
# load network
#########################################################################################################################

#network.path  <- "C:/Users/rds/Dropbox/PhD/PINA/"
network.path  <- "C:/Users/rsutherland/Dropbox/PhD/PINA/"
#network.path  <- "/Users/Russ/Dropbox/PhD/PINA/"

network.name  <- "pina101212_min2_noUBC"
network.file  <- paste0(network.name,".simple")

#Network file
network<- scan(file = file.path(network.path, network.file), skip = 1L, what = list(nameA = character(), nameB = character(),nPubs = numeric()), sep = "\t")
network2<-cbind(network[[1]],network[[2]])

networkGeneList<-(sort.int(unique(c(network$nameA,network$nameB))))# list of genes in network
indexRankedGenesTest<-seq(1:length(ranksGenesTest))

# Define the network
g<- graph.edgelist(as.matrix(network2), directed=FALSE)


#position of each gene in ranksGenesTest in the networkGeneList
indexMutatedGeneInNetwork<- sapply(seq(1:length(ranksGenesTest)), function(x) match(ranksGenesTest[x], networkGeneList ))
indexMappedGenes<- matrix(which(indexMutatedGeneInNetwork!= "NA"))

#A list of the genes which map to the network
MappedGenes<- networkGeneList[indexMutatedGeneInNetwork[which(indexMutatedGeneInNetwork!= "NA")]]


geneNeighbourhoods<- sapply(1:length(MappedGenes), 
                            function(x) V(g)[unlist(neighborhood(g, 1,MappedGenes[x]))])

geneNeighbourhoodsNames<- sapply(1:length(MappedGenes), function(x) V(g)$name[unlist(geneNeighbourhoods[x])])
# The names of the genes in each gene's neighbourhood in the network. Ordered as in MappedGenes and geneNeighbourhoods.
# I can now match the names back to the column names in genesVsSamplesTest and sum those columns to collapse the information.

#########################################################################################################################################################
#
###Adding genes from the network which do not occur in mutationMatrixLogical to mutationMatrixLogical2.

#The genes carrying functional mutations that map to the network
networkAndMutated<-intersect(networkGeneList,ranksGenesTest)
#The genes in the network that do not carry mutations and are not present in ranksTestGenes.
networkGAdd<-setdiff(networkGeneList,intersect(networkGeneList,ranksGenesTest))

# a FALSE matrix for the networkGAdd genes to be appended to mutationMatrixLogical.
networkGAddF<-matrix(rep(FALSE,(length(uniqueSamples)*length(networkGAdd))),length(networkGAdd),length(uniqueSamples))
rownames(networkGAddF)<-networkGAdd

# The logical matrix for all genes carrying at least 1 mutation in one sample and genes within the network that do not
# carry any mutations. This is important to include for the network based transfromation of the data that follows.
mutationMatrixLogical2<-rbind(mutationMatrixLogical,networkGAddF)


################################################################################################################################
##
#The two variables used as arguments in the network informed modification of the logical mutation matrix.
#mutMatLogicalOrdered and geneAdjIndex

# The logical mutation matrix to be input to the network informed score modification function.
# TRUE==At least one non-silent mutation within gene i in individual j. FALSE== No non silent mutation within gene i and individual j.
mutMatLogicalOrdered<- mutationMatrixLogical2[order(rownames(mutationMatrixLogical2), decreasing=FALSE),]

# The gene adjacency list for each gene mapping to the network.
# Elements in each indicate an index position in the rownames(mutMatLogicalordered) representing an adjacent gene in the network.
geneAdjIndex<-lapply(seq_len(length(geneNeighbourhoodsNames)), function(y) unlist(lapply(seq_len(length(geneNeighbourhoodsNames[[y]])), function(x) which(geneNeighbourhoodsNames[[y]][x]==rownames(mutMatLogicalOrdered)))))

################################################################################################################################



##Variables to indicate the cancerType of each sample and the colour associated with each cancer type
seqTech<-unique(cbind(mutations$Sequencer,as.character(mutations$sampleID)))
seqTech<- sub(" ","", seqTech)
seqTech<- seqTech[order(seqTech[,2], decreasing = FALSE),]
cancerType<-unique(cbind(mutations$Cancer_type,as.character(mutations$sampleID)))
cancerType<- cancerType[order(cancerType[,2], decreasing = FALSE),]


colourSamples<-function(metadata){
  
  sampleColors <- vector("character", length(metadata[,1]))
  
  sampleColors[which(metadata[,1]=="IlluminaHiSeq")]<-"red"
  sampleColors[which(metadata[,1]=="SOLID")]<-"blue"
  return(sampleColors)
}
sampleColors<-colourSamples(seqTech)

# a function to generate a color able for the samples, you just need to supply the variable name you wish to use as a categorical variable.
# This will return an numeric vector with samples labelled according to their metadata group using integers. This integer is used by the col attribute in the plot function to generate point colors.
sampstr<-function(metadata, variableName){
  met<-colorectal[which(rownames(colorectal)==variableName),]# the metadata variable
  groups<-unique(colorectal[which(rownames(colorectal)==variableName),])
  colrs<-sapply(tst, function(x) which(groups==x))# the colour indeces to be used to color plots points
  return(cbind(colrs,met))
}

smpclrs<-sampstr(metadata,"histological_type")


tst<-colorectal[which(rownames(colorectal)=="histological_type"),]
ret<-sapply(tst, function(x) which(groups==x))

#after generating the tables from the create metadata structure program
metadata<-colorectal[,match(colnames(mutMatLogicalOrdered), colnames(colorectal))]

#checks that the metadata table contains only the records for the ids I have scores for. I have 9 NAs
setdiff(colnames(metadata),colnames(mutMatLogicalOrdered))

#makes sure to match the metadata samples to those int mutMatLogicalOrdered
matchingSamplesIndex<-which(colnames(mutMatLogicalOrdered)%in%colnames(metadata))
mutMatLogicalOrdered<- mutMatLogicalOrdered[,matchingSamplesIndex]
metadata<- metadata[,matchingSamplesIndex]
# the metadata and mutMatLogicalOrdered matrices now have the same samples and I can run logistic regression.
#then calculate the dissimilarity between samples based on the mutMatLogicalOrdered

# 
#
### I need to filter each of the matrices I use in other scripts so that they only use data for the samples thatare present in both the metadata and the sequence data files




library(cluster)
medoidsModel<-pam(c, k=3)


#for the logistic regression


help(glm)