
library(igraph)
library(plyr)
library(scales)
###########################################################################################################################
# load network
#########################################################################################################################

#network.path  <- "C:/Users/rds/Dropbox/PhD/PINA/"
#network.path  <- "C:/Users/rsutherland/Dropbox/PhD/PINA/"
network.path  <- "/Users/Russ/Dropbox/PhD/PINA/"

network.name  <- "pina101212_min2_noUBC"
network.file  <- paste0(network.name,".simple")
##########################################################################################################################
##load seq data
##########################################################################################################################
# The data file in vcf-like format.
scores.path <- "/Users/Russ/Dropbox/PhD/tumour_classifier_data/colorectal_somatic_mutations/combined"
#scores.path <- "C:/Users/rds/Dropbox/PhD/tumour_classifier_data/colorectal_somatic_mutations/combined"
#scores.path <- "C:/Users/rsutherland/Dropbox/PhD/tumour_classifier_data/colorectal_somatic_mutations/combined"

scores.file <- "colorectalcancer.maf"
#scores.file <- basename("C:/Users/rsutherland/Dropbox/PhD/tumour_classifier_data/colorectal_somatic_mutations/combined/colorectalcancer.maf")

###########################################################################################################################
#load metadata
###########################################################################################################################

#The basename for the clinical files
#colonClinical.path <- "C:/Users/rds/Dropbox/PhD/tumour_classifier_data/colorectal_somatic_mutations/coloncancer/Clinical_17_12_12/Biotab/"
#colonClinical.path <- "C:/Users/rsutherland/Dropbox/PhD/tumour_classifier_data/colorectal_somatic_mutations/coloncancer/Clinical_17_12_12/Biotab/"
colonClinical.path <- "/Users/Russ/Dropbox/PhD/tumour_classifier_data/colorectal_somatic_mutations/coloncancer/Clinical_17_12_12/Biotab/"

#The basename for the rectal clinical files
#rectumClinical.path <- "C:/Users/rds/Dropbox/PhD/tumour_classifier_data/rectum_adenocarcinoma/Clinical_18_02_2013/Biotab/"
#rectumClinical.path <- "C:/Users/rsutherland/Dropbox/PhD/tumour_classifier_data/rectum_adenocarcinoma/Clinical_18_02_2013/Biotab/"
rectumClinical.path <- "/Users/Russ/Dropbox/PhD/tumour_classifier_data/rectum_adenocarcinoma/Clinical_18_02_2013/Biotab/"

############################################################################################################################
#loading metadata
############################################################################################################################

clinical.files <- list.files(colonClinical.path, full.names = TRUE)
names(clinical.files) <- basename(clinical.files)

colonClinicalTables <- lapply(clinical.files, read.table, header=TRUE,as.is =TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings=c("NA","[Not Available]"))



clinicalDataTable<- function(rectumClinical.path){
  #the names of the rectumClinical files
  rectumClinical.files<- list(list.files(rectumClinical.path))
  rectumClinicalTables<- lapply(seq(1:length(rectumClinical.files[[1]])), function(x) cbind(read.table(file.path(rectumClinical.path,rectumClinical.files[[1]][x]) , header=TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings=c("NA", "[Not Available]"))))
  names(rectumClinicalTables)<- unlist(rectumClinical.files)
  
  rectumClinicalTablesSamples<- lapply(rectumClinicalTables, function(x) strsplit(x[,1], split="-", fixed =TRUE))
  #variables for the table
  TablesForTable<- which(sapply(rectumClinicalTablesSamples, function(x)length(x[[1]]))==3)
  #modified table list based on the data in tables with duplicate SampleIDs. ommitted tables are not needed, drug table may be included at a later date.
  TablesForTable<- TablesForTable[-c(3,6,8)]# This needs to be changed manually
  
  #Identify the tables which use individual IDs as opposed to tissue sample ids which can come from the same patient.
  rectumClinicalTablesSamples<-lapply(seq(1:length(TablesForTable)),function(x) rectumClinicalTablesSamples[[TablesForTable[x]]])
  #get the appropriate Clinical Tables
  #no factor types here yet
  rectumClinicalTables<-lapply(seq(1:length(TablesForTable)),function(x)rectumClinicalTables[[TablesForTable[x]]])
  #get the sample names for each row of each metadata table
  rectumClinicalTablesSampleIDs1<- lapply(rectumClinicalTablesSamples, function(x) unlist(lapply(x, function(y) paste(y[1:3],collapse="")), use.names=FALSE))
  # assign sample ids to rectumClinicalTables. Now each table should have the same sampleIDS and I can look for uniques etc...
  #no factor tyes here yet
  rectumClinicalTables<-lapply(seq(1:length(rectumClinicalTables)), function(x) data.frame(rectumClinicalTables[[x]],rectumClinicalTablesSampleIDs1[[x]], stringsAsFactors=FALSE))
  
  #rename the newly added sample IDs column
  colnames(rectumClinicalTables[[1]])[length(rectumClinicalTables[[1]][1,])]<-"SampleIds"
  #assign new sampleid column names to all tables of rectumClinical Samples
  for( i in seq_along(rectumClinicalTables)){
    colnames(rectumClinicalTables[[i]])[length(rectumClinicalTables[[i]][1,])]<- "SampleIds"
  }
  #get the list of unique samples across all files
  uniqueSamples<-as.list(as.character(unique(unlist(lapply(rectumClinicalTables, function(x) x$SampleIds), use.names=FALSE))))
  
  
  # find the number of variables across all of my tables
  vnumber<- function(tabList){
    vnum<-0
    for(i in seq_along(tabList)){
      vnum<-vnum+((dim(tabList[[i]])[2])-2)# the -2 is because there are two id columns in each table that I don't want
    }
    return(vnum)
  }
  
  vn<-vnumber(rectumClinicalTables)
  # find the variable names
  #vars<-lapply(rectumClinicalTables, function(x) colnames(x)[c(-1,-(length(x[1,])))])
  #varNames<-unlist(vars)
  
  # create a table of all data from all tables based on unique sampleIDs
  
  #The join_all function acts like a mysql join, to combine all of the listed tables based on matching "sampleIds" variables.
  clinicalData<- join_all(rectumClinicalTables, by="SampleIds", type="full", match="first")
  #disease type matrix
  tumourType<- matrix("rectum", ncol=1, nrow=length(clinicalData[,1]))
  # Add disease identifier to the clinicalData
  #no factors at this point
  clinicalData<-data.frame(clinicalData,tumourType, stringsAsFactors=FALSE)
  #change the SampleIds names to sampleID
  colnames(clinicalData)[1]<-"sampleID"
  return(clinicalData)
}


# function to combine metadata tables and ensures that variables remain in the correct order so that columns match netween the two joined tables
combineMetadata<- function(table1,table2){
  table1<-table1[,order(match(colnames(table1), colnames(table2)))]
  combinedT<-as.data.frame(rbind(table1,table2), stringsAsFactors=FALSE)
  rownames(combinedT)<-combinedT[,1]
  return(combinedT)
}


#load colon and rectum data
colon<-clinicalDataTable(colonClinical.path)
rectum<-clinicalDataTable(rectumClinical.path)
# order the colon columns the same way as they rectum columns
#colon<-colon[,order(match(colnames(colon),colnames(rectum)))]
#colorectal<- t(rbind(colon,rectum))# my metadata table to combine with the genescores data
#colnames(colorectal)<-colorectal[which(rownames(colorectal)=="SampleIds"),]

#combine the colon and rectum metadata tables
colorectal<- combineMetadata(colon,rectum)

#colorectal<-as.data.frame(colorectal) I don't need this to be a datframe 

#####################################################################################################################################################################
##loading seqdata
#####################################################################################################################################################################

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
  scores$sampleID<- as.character(scores$sampleID)
  
  return(scores)
}

scores<-createScoreTable(dataFile)
scoresSamples<- unique(scores$sampleID)
intersectSamples<-(intersect(scoresSamples,rownames(colorectal)))

extractIntersectSamples<- function(dataNames,intersectList,data){
  if(class(dataNames)!= "character"||class(intersectList)!="character"){
    stop("Both SeqInfoNames and MetadataNames must be of character type")}else{
      SampleMatchI<- which(dataNames%in%intersectList)
      NewSeqInfo <- data[SampleMatchI,]
      #NewSeqInfo<- NewSeqInfo[order(NewSeqInfo$sampleID),]# it is better if I order the sample names after uniqueing
      return(NewSeqInfo)
    }
}

scores<-extractIntersectSamples(scores$sampleID,intersectSamples,scores)# this produces the same result as in the original code, but much faster
scores<- scores[order(scores$sampleID),]

#the metadata must be a matrix
metadata2<-extractIntersectSamples(rownames(colorectal), intersectSamples, colorectal)# this produced the same result as in the original code but much faster
metadata2<- metadata2[order(rownames(metadata2)),]# the rownames(sampleIDS) are in the same order as in seqTech, cancerType and the sampleID variable of scores. Everything is now in the same order
#metadata_tab2<-cbind(metadata_tab,extractSamples(seqTech[,2],colnames(metadata), seqTech[,1]))
## I need to check that this new variable is the same as the metadata_tab table, to make sure the above fuinction is working. Once I know it is working I can delete the other code and supplement it with the function.
#identical(metadata_tab[,!c(105,106)],metadata_tab2[,!c(105,106,107)])# the two dataframes are identical and the matrices match back to the metadata matrix. I should use metadata_tab2 as the basis for my analysis.





getFuncMutations<-function(scoresT){
  mutations<-scoresT[which(scores$Variant_Classification!="Silent"),]
  colnames(mutations)<- colnames(scoresT)
  return (mutations)
}

# The matrix of "non-silent mutations"
mutations<-scores[which(scores$Variant_Classification!="Silent"),]
colnames(mutations)<- colnames(scores)





# nonsilent mutations per individual
mutationsTypePerIndiv<-table(mutations$Variant_Classification, mutations$sampleID)
mutationsSilentPerIndiv<-sapply(seq(1:length(intersectSamples)), function(x) sum(mutationsTypePerIndiv[1:4,x]))
names(mutationsSilentPerIndiv)<-colnames(mutationsTypePerIndiv)

#mutationsPerIndiv <-aggregate(mutations[,1], by=list(mutations$sampleID), length)
#mutationsPerGene <-aggregate(mutations[,1], by=list(mutations$Hugo_Symbol), length)
#the functions below are faster



#######################################################################################################################
## Plot mean number of mutations oer gene per individual
#######################################################################################################################
mutationsPerIndiv <-tapply(mutations$Hugo_Symbol, mutations$sampleID,length)
mutationsPerGene <-tapply(mutations$sampleID, mutations$Hugo_Symbol,length)

test<- tapply(mutations$sampleID, mutations[,c("Hugo_Symbol","sampleID")], length)

#gives me the mutations per individual per gene. I can now split this across a factor
mutsPerIndivPerGene<-as.matrix(table(mutations$sampleID,mutations$Hugo_Symbol))
muts
#meanMutsPerIndivPerGene<- sapply(seq(1:dim(mutsPerIndivPerGene)[2]), function(x) mean(mutsPerIndivPerGene[,x]))

#converts a table Object to a matrix object
tableToMatrix<- function(tableData){
  data<-as.matrix(tableData)
  class(data)<- "matrix"
  return(data)
}


matrixMutsPerIndivPerGene<- tableToMatrix(mutsPerIndivPerGene1)

# Returns the mean number of mutations carried in each Individual for each gene.
getPerGenePerSamp<- function(PerIndpGene){
  #perIndivPerGene<-as.matrix(table(mutDataFrame$sampleID,mutDataFrame$Hugo_Symbol))
  meanPerIndivPerGene<- sapply(seq(1:dim(PerIndpGene)[2]), function(x) mean(PerIndpGene[,x]))
  return(meanPerIndivPerGene)
}

#test5<- getPerGenePerSamp(mutsPerIndivPerGene)#identical to meanMutsPerIndivPerGene

stratMeanMutPIPG<-by(matrixMutsPerIndivPerGene, stratification$seqTech[,2], getPerGenePerSamp)# Using this I should simply be able to plug it in to getPerGenePerSamp
#variables to order the mutations per individual per gene graph across chromosomes and bp
maxMutPos<-tapply(mutations$Start_Position, mutations$Hugo_Symbol, max)
geneChrNum<-mutations$Chr[match(ranksGenesTest,mutations$Hugo_Symbol)]



plot(stratMeanMutPIPG[[1]], col=alpha(1,0.5), log="y")# The beginings of my plot according to sequencer Type
points(stratMeanMutPIPG[[2]],col=alpha(2,0.5))

##################################################################################################################################################################################
##################################################################################################################################################################################


##extra metadata_variables
##Variables to indicate the cancerType of each sample and the colour associated with each cancer type
seqTech<-unique(cbind(mutations$Sequencer,as.character(mutations$sampleID)))
seqTech<- sub(" ","", seqTech)
seqTech<- seqTech[order(seqTech[,2], decreasing = FALSE),]
cancerType<-unique(cbind(mutations$Cancer_type,as.character(mutations$sampleID)))
cancerType<- cancerType[order(cancerType[,2], decreasing = FALSE),]

metadata2<-cbind.data.frame(metadata2, cancerType=cancerType[,1], seqTech=seqTech[,1], stringsAsFactors=FALSE)# adds the cancerType and seqTech metadata to the metadata object

# The matrix of "non-silent mutations"
#mutations<-scores[which(scores$Variant_Classification!="Silent"),]
#colnames(mutations)<- colnames(scores)

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
mutationMatrixLogical <- matrix(as.logical(mutationMatrixBinary),length(ranksGenesTest),length(intersectSamples))
rownames(mutationMatrixLogical)<- rownames(mutationTable)
colnames(mutationMatrixLogical)<- colnames(mutationTable)


#Begin again from here tomorrow. runcreate metadata and then this script
######################################################################################################################################
#I need this to reassign all of the NA values to "[Not Available]" so that the which function in sampstr will work
metadata3<-metadata2# in doing this, the metadata2 is no longer a dataframe with the correct data types.
metadata3[is.na(metadata3)]<-"[Not Available]"

sampstr<-function(metadata, variableName){
  met<-metadata[,which(colnames(metadata)==variableName)]# the metadata variable
  groups<-unique(metadata[,which(colnames(metadata)==variableName)])
  colrs<-sapply(met, function(x) which(groups==x))# the colour indeces to be used to color plots points
  return(cbind.data.frame(colrs,met,stringsAsFactors=FALSE))
}

t1<-unique(metadata3[,which(colnames(metadata3)=="histological_type")])

t2<- unique(metadata3[,which(colnames(metadata3)=="histological_type")])

smpclrs<-sampstr(metadata3,"histological_type")# This now has the correct number of samples


stratification<-lapply(colnames(metadata3), function(x) sampstr(metadata3,x))# integers representing metadata classes for coloring of samples in the MDS plots
names(stratification)<- colnames(metadata3)
#################################################################################################################################################################################
#network data
#################################################################################################################################################################################

#Network file
network<- scan(file = file.path(network.path, network.file), skip = 1L, what = list(nameA = character(), nameB = character(),nPubs = numeric()), sep = "\t")
network2<-cbind(network[[1]],network[[2]])

networkGeneList<-(sort.int(unique(c(network$nameA,network$nameB))))# list of genes in network
indexRankedGenesTest<-seq(1:length(ranksGenesTest))# do I need this?

# Define the network
g<- graph.edgelist(as.matrix(network2), directed=FALSE)


#position of each gene in ranksGenesTest in the networkGeneList
#indexMutatedGeneInNetwork<- sapply(seq(1:length(ranksGenesTest)), function(x) match(ranksGenesTest[x], networkGeneList ))
#indexMappedGenes<- matrix(which(indexMutatedGeneInNetwork!= "NA"))
#not needed for MappedGenes


#A list of the genes which map to the network
MappedGenes<- intersect(networkGeneList, ranksGenesTest)


geneNeighbourhoods<- sapply(1:length(MappedGenes), 
                            function(x) V(g)[unlist(neighborhood(g, 1,MappedGenes[x]))])

geneNeighbourhoodsNames<- sapply(1:length(MappedGenes), function(x) V(g)$name[unlist(geneNeighbourhoods[x])])



#########################################################################################################################################################
#
###Adding genes from the network which do not occur in mutationMatrixLogical to mutationMatrixLogical2.

#The genes carrying functional mutations that map to the network
networkAndMutated<-intersect(networkGeneList,ranksGenesTest)
#The genes in the network that do not carry mutations and are not present in ranksTestGenes.
networkGAdd<-setdiff(networkGeneList,MappedGenes)

# a FALSE matrix for the networkGAdd genes to be appended to mutationMatrixLogical.
networkGAddF<-matrix(rep(FALSE,(length(intersectSamples)*length(networkGAdd))),length(networkGAdd),length(intersectSamples))
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
## run the network informed clustering function here to assign mutMatLogicalOrdered and geneAdjIndex to mutM and geneAdj.
################################################################################################################################
################################################################################################################################
