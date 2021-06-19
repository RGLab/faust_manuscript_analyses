library(flowWorkspace)
library(faust)
library(xtable)

gs <- load_gs(file.path(normalizePath("../.."),
                        "dataSets",
                        "publication_gating_sets",
                        "CITN-07",
                        "citn07_longitudinal_gs"))

startNode <- "45+"
activeChannelsIn <- c("CD123", "CD4", "CD14", "CD11C", "CD56", "CD8", "CD16", "CD3", "CD122", "CD19", "HLA DR")
channelBoundsIn <- structure(
    c(1500, 3000, 100, 2500, 1, 3000, 1, 3000, 1, 3000, 1, 3000,
      1000, 3000, 1, 3000, 1000, 3000, 1501, 3000, 1, 3500),
    .Dim = c(2L, 11L),
    .Dimnames = list(
        c("Low", "High"),
        c("CD123", "CD4", "CD14", "CD11C", "CD56", "CD8", "CD16", "CD3", "CD122", "CD19", "HLA DR")
    )
)
xtable(t(channelBoundsIn))

faust(
    gatingSet = gs,
    experimentalUnit="name",
    imputationHierarchy="genVisit",
    activeChannels = activeChannelsIn,
    channelBounds = channelBoundsIn,
    startingCellPop = startNode,
    debugFlag = TRUE,
    threadNum = 11,
    seedValue = 12345,
    supervisedList = list("CD4" = list(actionType = "Preference", action = 2)),
    annotationsApproved = TRUE,
    nameOccuranceNum = 50,
    drawAnnotationHistograms = FALSE,
    plottingDevice = "pdf"
)
