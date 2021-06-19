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
library(FlowSOM)

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

ff <- flowFrame(train_data)
fSOM <- FlowSOM::ReadInput(ff, transform = FALSE, scale = FALSE)
fSOM <- FlowSOM::BuildSOM(
                     fSOM,
                     colsToUse = seq(length(colnames(train_data))),
                     xdim = 10,
                     ydim = 10
                 )
fSOM <- FlowSOM::BuildMST(fSOM)
FlowSOM_out <- fSOM$map$mapping[, 1]
saveRDS(FlowSOM_out,"./FlowSOM_out_grid_10_all.rds")

ff20 <- flowFrame(train_data)
fSOM20 <- FlowSOM::ReadInput(ff20, transform = FALSE, scale = FALSE)
fSOM20 <- FlowSOM::BuildSOM(
                     fSOM20,
                     colsToUse = seq(length(colnames(train_data))),
                     xdim = 20,
                     ydim = 20
                 )
fSOM20 <- FlowSOM::BuildMST(fSOM20)
FlowSOM20_out <- fSOM20$map$mapping[, 1]
saveRDS(FlowSOM20_out,"./FlowSOM_out_grid_20_all.rds")
