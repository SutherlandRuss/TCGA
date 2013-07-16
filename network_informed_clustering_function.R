#' @param mutM    A logical matrix genes are rows and samples are columns. 
#'                TRUE if the gene is mutated in the sample and FALSE otherwise.
#' @param geneAdj A list of genes that are adjacent to.....
#' @param geneID  An lexicographically ordered list of genes.
#' 
#' @return mutNet The logical matrix MutM updated using network information.
#'
#' Requires packages "igraph" "pam" and "fpc".
#'
#' The length(geneAdj)==length(geneID)==nrows(mutM)
#' The rownames(mutM)==geneID==names(geneAdj)
#'
#'  
#'  
setwd("C:/Users/rsutherland/Dropbox/PhD/tumour_classifier_data/colorectal_somatic_mutations/combined/1014_06_2013")
#setwd("C:/Users/rds/Dropbox/PhD/tumour_classifier_data/colorectal_somatic_mutations/combined/2428_06_2013")



mutM<-mutMatLogicalOrdered
geneAdj<-geneAdjIndex
isColon<- cancerType[,1]
names(isColon) <- cancerType[,2]

if(FALSE) {
  geneID = c("A", "B", "C", "D", "E", "F")
  mutM = cbind(S1 = c(T, F, F, F, T, F), 
               S2 = c(F, T, F, F, T, F),
               S3 = c(T, F, F, F, F, T),
               S4 = c(F, F, F, T, F, F))
  rownames(mutM) = geneID
  geneAdj = list(A = c(2L, 5L), 
                 B = c(1L, 3L, 4L), 
                 C = c(2L, 4L), 
                 D = c(2L, 3L),
                 E = c(1L), 
                 F = integer())
}

netInf <- function(mutM, geneAdj){
  ## iterate over each sample/column and add the network information
  for(s in seq_len(ncol(mutM))) {
    ## select the functional/true genes
    selG <- which(mutM[, s])
    if(length(selG) == 0L)
      next
    ## select all adjacent genes
    adjG <- unique.default(unlist(geneAdj[selG], use.names = FALSE))
    if(length(adjG) > 0L)
      mutM[adjG, s] <- TRUE  
  }
  
  return(mutM)
}

##############################################################################
#using the above functions and plotting a comparison between before and after network inference

#result<-func1(mutM,geneAdj)
#resultdist<-dist(t(result),method= "binary")
#edist <- dist(t(result))

#countM0<-rowSums(mutM)
#countM1<-rowSums(result)
#smoothScatter(countM0,countM1)

#countS0<-colSums(mutM)
#countS1<-colSums(result)
#smoothScatter(countS0,countS1)
######################################################################


#pos = (which(mutM) - 1L) %% nrow(mutM) + 1L
#col = (which(mutM) - 1L) %/% nrow(mutM) + 1L

##split(pos, col)
#l = lapply(split(geneAdj[pos], col), unlist, use.names = FALSE)

#newcol = rep.int(seq_along(l), sapply(l, length))
#newtrue = unlist(l, use.names = FALSE)

#newidx = (newcol - 1) * nrow(mutM) + newtrue
#a = mutM
#a[newidx] <- TRUE
#a
#m
#identical(a,m)

#####################################################################################



#counts the number of TRUE matches between observations.
countMatch1<- function(mutM){
  nc<-ncol(mutM)
  countM<- matrix(NA,nc,nc)
  for(s in seq_len(nc-1)){
    for(t in (s+1):nc){
      countM[t, s]<-length(which(mutM[which(mutM[,s]),t]))
      #countM[s,t]<-sum(mutM[,s]*(mutM[,s]==mutM[,t]))
    }
  }
  return(countM)
}

#function to compute the disimilarity between
#a binary matrix of n observations and m variables
countMatch2<- function(mutM){
  nc<-ncol(mutM)
  countM<- matrix(NA,nc,nc)
  for(s in seq_len(nc-1)){
    a <- mutM[,s,drop=TRUE]
    countM[(s+1):nc, s] <- colSums((mutM[,(s+1):nc, drop = FALSE]==a) * a)
  }
  return(countM)
}

# computes the jaccard similarity between all samples. Either this should be run or 
#countMatch AND compDiss.

jaccardSim<- function(mutM){
  nc<-ncol(mutM)
  jaccard<- matrix(NA,nc,nc)
  for(s in seq_len(nc-1)){
    for(t in (s+1):nc){
      jaccard[t, s]<-clujaccard(mutM[,s], mutM[,t])
    }
  } 
  jaccard<-jaccard[!is.na(jaccard)]
  return(jaccard)
}



#compute the dissimilarity matrix
compDiss<-function(countM, type, mutM){
  if(type==1){
    distM<- as.dist(1-countM/max(countM, na.rm=TRUE)) 
  }else if(type==2){
    distM<- as.dist(max(countM, na.rm =TRUE)-countM)
  }#add in jaccards similarity here
  attr(distM, "Labels") = colnames(mutM)
  return(distM)
}


#function for applying various cluster analyses including hierarchical
#and partition around medoid methods
grouping<- function(distM,isColon,cmethod,cutN,title){
  if(cmethod=="pam"){
    clustr<- pam(distM,k=cutN,diss=TRUE)
    clustCut<-clustr$clustering
    plot(clustr, main =title)
  }else{
    clustr<-hclust(distM,method=cmethod)
    clustCut<-cutree(clustr,k=cutN)
    plot(clustr,main=title,labels = c("o", "--------")[(isColon[clustr$labels] == "colon") + 1])
  }
  names(clustCut)<-colnames(mutM)
  groupT<-table(as.integer(clustCut), isColon[names(clustCut)])
  enrichment<-fisher.test(groupT)
  return(list(groupT,enrichment))
}


##example
a<-netInf(mutM,geneAdj)
b<-countMatch1(a)
c<-compDiss(b,1,mutM)
d<-grouping(c,isColon,"ward",5,"dissimilarity=(number of samples - matchcount) \n cluster method= ward  k=2")
##for jaccards as d
d<-grouping(b[!is.na(b)],isColon,"pam",2,"jaccards and two groups")

