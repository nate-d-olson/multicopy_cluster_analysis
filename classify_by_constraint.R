#title: "Exploring Cluster Genome Counts"
#author: "Ethan Ertel"
#date: "December 6, 2015"

#If the lpSolve library is not installed, use install.packages("lpSolve")
library(lpSolve)
library(readr)

#parse inputs
args=commandArgs(trailingOnly = T)
trainFile = args[1]
testFile = args[2]
outFile = args[3]
taxa99train <- read_csv(trainFile)
taxa99test <- read_csv(testFile)

# These are indices of the unique entries under each category
species <- unique(taxa99train$species)
clusters <- unique(taxa99train$cluster)
genus <- unique(taxa99train$genus)

#The matrix classMat contains the number of 16S copies per taxon present in each identified cluster
classMat <- matrix(0,length(species),length(clusters))

#add counts for each cluster and species combination to the matrix
for(i in c(1:dim(taxa99train)[1])){
  classMat[which(species==taxa99train$species[i]),which(clusters==taxa99train$cluster[i])] <- 1 + classMat[which(species==taxa99train$species[i]),which(clusters==taxa99train$cluster[i])]
}

#testClusters is the right hand side of the linear equation; it is a summary of the number of counts mapped to a cluster id in the taxa99test data set
testClusters = rep(0,length(clusters))
for(i in c(1:length(taxa99test$cluster))){
  idx = which(clusters == taxa99test$cluster[i])
  testClusters[idx] = testClusters[idx] + taxa99test$count[i]
}

#lpObject is the linear programming result of the problem [classMat][x]=[testClusters]
lpObject<- lp(direction = "min",objective.in = rep(1,length(species)),const.mat =t(classMat),const.dir = ">=",const.rhs = testClusters)
outMat <- data.frame(genus,rep(0,length(genus)))

#classifications are pooled to genus level and then output
for(i in 1:length(species)){ #taxa99test$genus[
  pos <- which(taxa99test$species==species[i])[1]
  gns <- taxa99test$genus[pos]
  outMat[which(genus==gns),2] <- outMat[which(genus==gns),2] + lpObject$solution[which(species==taxa99test$species[pos])]
}
write_csv(x = outMat,outFile)