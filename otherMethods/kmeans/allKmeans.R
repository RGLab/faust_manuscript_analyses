library(ggridges)
library(viridis)
library(ggplot2)
library(cowplot)
library(tidyr)
library(readxl)
library(dplyr)
library(stringr)
library(flowCore)
library(flowWorkspace)
library(Rphenograph)

gs <- load_gs(file.path(normalizePath("../.."),"dataSets","publication_gating_sets","CITN-09","whole_blood_samples","citn09_fresh_whole_blood_gs"))
set.seed(123)
pmData <- flowCore::parameters(gh_pop_get_data(gs[[1]],"root"))@data
firstSample <- TRUE
for (sampleName in sampleNames(gs)){
    print(sampleName)
    subjectID <- str_extract(sampleName,"\\d\\d\\d-\\d\\d-\\d\\d\\d")
    rawData <- flowCore::exprs(gh_pop_get_data(gs[[sampleName]],"L"))
    colnames(rawData) <- as.character(pmData[match(colnames(rawData),pmData$name),"desc"])
    colnames(rawData)[which(is.na(colnames(rawData)))] <- pmData[which(is.na(colnames(rawData))),"name"]
    sampleID <- rep(sampleName,nrow(rawData))
    if (firstSample) {
        all_train_data <- rawData[,-which(colnames(rawData)=="Time")]
        sampleLookup <- sampleID
        firstSample <- FALSE
    }
    else {
        all_train_data <- rbind(all_train_data,rawData[,-which(colnames(rawData)=="Time")])
        sampleLookup <- append(sampleLookup,sampleID)
    }
}
saveRDS(sampleLookup,"./sampleLookup.rds")
activeChannels <- c("CD278 ICOS","CD3","CD127","CD197 CCR7","CD279 PD1","CD8","CD4","CD28","CD25","HLA DR","CD45RA")
train_data <- all_train_data[,activeChannels]


print("kmeans 100 start ")
kmeans100Result <- kmeans(x=train_data,centers=100,iter.max=10000,algorithm="Lloyd")
kmeans100ClusterLabels <- as.numeric(kmeans100Result$cluster)
saveRDS(kmeans100ClusterLabels,"./kmeans_out_k_100_all.rds")
print("kmeans 100 end ")

print("kmeans 400 start ")
kmeans400Result <- kmeans(x=train_data,centers=400,iter.max=10000,algorithm="Lloyd")
kmeans400ClusterLabels <- as.numeric(kmeans400Result$cluster)
saveRDS(kmeans400ClusterLabels,"./kmeans_out_k_400_all.rds")
print("kmeans 400 end ")
