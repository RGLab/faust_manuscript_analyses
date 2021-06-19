library(faust)
library(flowWorkspace)
library(xtable)

gs <- load_gs(file.path(normalizePath("../../.."),
                        "dataSets",
                        "publication_gating_sets",
                        "CITN-09",
                        "whole_blood_samples",
                        "citn09_fresh_whole_blood_gs"))
xtable(gh_pop_get_data(gs[[1]],"L")@parameters@data)
activeChannelsIn <- c("CD278 ICOS","CD3","CD127","CD197 CCR7","CD279 PD1",
                      "CD8","CD4","CD28","CD25","HLA DR","CD45RA")

channelBoundsIn <- structure(c(-500, 3000, -500, 3000, -500, 3000, 500, 2250, -250, 
3500, 500, 3000, 500, 3000, -500, 3000, 500, 2250, -500, 3000, 
1000, 3000), .Dim = c(2L, 11L), .Dimnames = list(c("Low", "High"
), c("CD278 ICOS", "CD3", "CD127", "CD197 CCR7", "CD279 PD1", 
"CD8", "CD4", "CD28", "CD25", "HLA DR", "CD45RA")))
xtable(t(channelBoundsIn))
supervisedList <- list(`CD279 PD1` = list(actionType = "Preference", action = 2))

faust(
    gatingSet=gs,
    activeChannels=activeChannelsIn,
    channelBounds=channelBoundsIn,
    startingCellPop="L",
    experimentalUnit="name",
    imputationHierarchy="VISIT",
    debugFlag = TRUE,
    seedValue = 1234,
    nameOccuranceNum = 7,
    annotationsApproved = TRUE,
    supervisedList = supervisedList,
    drawAnnotationHistograms = FALSE,
    threadNum=11,
    plottingDevice = "pdf"
)


