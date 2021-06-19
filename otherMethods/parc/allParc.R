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
library(reticulate)
parc <- import("parc",convert=FALSE)

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

print("Starting parc")
parcResult <- parc$PARC(r_to_py(train_data))
parcResult$run_PARC()
print("parc complete")
parcClusterLabels <- unlist(py_to_r(parcResult$labels))
saveRDS(parcClusterLabels,"./parc_out_all.rds")
