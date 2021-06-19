library(faust)
library(flowWorkspace)
library(xtable)

gs <- load_gs(file.path(normalizePath("../.."),
                        "dataSets",
                        "publication_gating_sets",
                        "melanoma_facs",
                        "melanoma_facs_gs"))
#
#Define faust parameters and run faust
#
activeChannelsIn <- c("CD3","CD4","CD11b","CD33","HLA-DR",
                      "CD56","CD45RO","CD11c","CD16","CD14",
                      "CD19")
channelBoundsIn <- structure(c(-Inf, Inf, -20, Inf, -Inf, Inf, -Inf, Inf, -Inf, 
Inf, -Inf, Inf, -Inf, Inf, -Inf, Inf, 85, Inf, -Inf, Inf, -Inf, 
Inf), .Dim = c(2L, 11L), .Dimnames = list(c("Low", "High"), c("CD3", 
"CD4", "CD11b", "CD33", "HLA-DR", "CD56", "CD45RO", "CD11c", 
"CD16", "CD14", "CD19")))
xtable(t(channelBoundsIn))


faust(
    gatingSet=gs,
    activeChannels=activeChannelsIn,
    channelBounds=channelBoundsIn,
    startingCellPop="life",
    experimentalUnit="name",
    depthScoreThreshold = 0.01,
    selectionQuantile = 0.0,
    debugFlag = TRUE,
    threadNum = 11,
    seedValue = 1234,
    annotationsApproved = TRUE,
    drawAnnotationHistograms=FALSE,
    plottingDevice = "pdf"
)

