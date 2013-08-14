#####
## A program to create a clinical file for the cancer of choice from the Cancer Genome Atlas.
library("plyr")
##load the data

#The basename for the clinical files
#colonClinical.path <- "C:/Users/rds/Dropbox/PhD/tumour_classifier_data/colorectal_somatic_mutations/coloncancer/Clinical_17_12_12/Biotab/"
colonClinical.path <- "C:/Users/rsutherland/Dropbox/PhD/tumour_classifier_data/colorectal_somatic_mutations/coloncancer/Clinical_17_12_12/Biotab/"
#colonClinical.path <- "/Users/Russ/Dropbox/PhD/tumour_classifier_data/colorectal_somatic_mutations/coloncancer/Clinical_17_12_12/Biotab/"

#the names of the clinical files
#clinical.files<-list(list.files("C:/Users/rds/Dropbox/PhD/tumour_classifier_data/colorectal_somatic_mutations/coloncancer/Clinical_17_12_12/Biotab"))
clinical.files <- list.files(colonClinical.path, full.names = TRUE)
names(clinical.files) <- basename(clinical.files)
#clinical.files<-list(list.files("/Users/Russ/Dropbox/PhD/tumour_classifier_data/colorectal_somatic_mutations/coloncancer/Clinical_17_12_12/Biotab"))

# Loads all of the available clinical datafiles for the cancer of choice in to a list of lists
colonClinicalTables <- lapply(clinical.files, read.table, header=TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)


#The basename for the rectal clinical files
#rectumClinical.path <- "C:/Users/rds/Dropbox/PhD/tumour_classifier_data/rectum_adenocarcinoma/Clinical_18_02_2013/Biotab/"
rectumClinical.path <- "C:/Users/rsutherland/Dropbox/PhD/tumour_classifier_data/rectum_adenocarcinoma/Clinical_18_02_2013/Biotab/"
#rectumClinical.path <- "/Users/Russ/Dropbox/PhD/tumour_classifier_data/rectum_adenocarcinoma/Clinical_18_02_2013/Biotab/"




clinicalDataTable<- function(rectumClinical.path){
  #the names of the rectumClinical files
  rectumClinical.files<- list(list.files(rectumClinical.path))
  rectumClinicalTables<- lapply(seq(1:length(rectumClinical.files[[1]])), function(x) cbind(read.table(file.path(rectumClinical.path,rectumClinical.files[[1]][x]) , header=TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)))
  names(rectumClinicalTables)<- unlist(rectumClinical.files)
  
  rectumClinicalTablesSamples<- lapply(rectumClinicalTables, function(x) strsplit(x[,1], split="-", fixed =TRUE))
  #variables for the table
  TablesForTable<- which(sapply(rectumClinicalTablesSamples, function(x)length(x[[1]]))==3)
  #modified table list based on the data in tables with duplicate SampleIDs. ommitted tables are not needed, drug table may be included at a later date.
  TablesForTable<- TablesForTable[-c(3,6,8)]# This needs to be changed manually
  
  #Identify the tables which use individual IDs as opposed to tissue sample ids which can come from the same patient.
  rectumClinicalTablesSamples<-lapply(seq(1:length(TablesForTable)),function(x) rectumClinicalTablesSamples[[TablesForTable[x]]])
  #get the appropriate Clinical Tables
  rectumClinicalTables<-lapply(seq(1:length(TablesForTable)),function(x)rectumClinicalTables[[TablesForTable[x]]])
  #get the sample names for each row of each metadata table
  rectumClinicalTablesSampleIDs1<- lapply(rectumClinicalTablesSamples, function(x) unlist(lapply(x, function(y) paste(y[1:3],collapse="")), use.names=FALSE))
  # assign sample ids to rectumClinicalTables. Now each table should have the same sampleIDS and I can look for uniques etc...
  rectumClinicalTables<-lapply(seq(1:length(rectumClinicalTables)), function(x) cbind(rectumClinicalTables[[x]],rectumClinicalTablesSampleIDs1[[x]]))
  
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
  clinicalData<-cbind(clinicalData,tumourType)
  #change the SampleIds names to sampleID
  colnames(clinicalData)[1]<-"sampleID"
  return(clinicalData)
}


# function to combine metadata tables and ensures that variables remain in the correct order so that columns match netween the two joined tables
combineMetadata<- function(table1,table2){
  table1<-table1[,order(match(colnames(table1), colnames(table2)))]
  combinedT<-t(rbind(table1,table2))
  colnames(combinedT)<-combinedT[1,]
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







setwd("C:/users/rds/Desktop")

write.table(clinicalData,file="rectumClinicalData.txt", append=FALSE, sep="\t", eol="\n", na="NA", dec=".", row.names=TRUE, col.names=TRUE)

# do this after I have merged all of the tables together within each disease. This way I can add a disease identifier column.I don't need to do this
test1<- rbind(colonClinicalTables[[1]][1], rectumClinicalTables[[1]][1])
test2<- merge(colonClinicalTables[[1]], colonClinicalTables[[2]], by=intersect(names(colonClinicalTables[[1]][1]),names(colonClinicalTables[[1]][2])))




sapply(rectumClinicalTables, function(x) duplicated(x[,length(x[1,])]))
sapply(rectumClinicalTables, function(x) length(x[,1]))
