library(flowCore)
library(flowWorkspace)
library(CytoML)
library(faust)

gs <- load_gs(file.path(normalizePath("../../.."),
                        "dataSets",
                        "publication_gating_sets",
                        "CITN-09",
                        "cryopreserved_pbmcs",
                        "citn09_cryo_pbmcs_gs"))
activeChannelsIn <- c("CD3", "CD127", "CD197 CCR7", "CD279 PD1", "CD8", "CD4", "CD28", 
"CD25", "HLA DR", "CD45RA")

cbIN <- structure(c(1200, 2690.337, 750, 1750, 1500, 2250, 1400, 2700, 
1800, 2800, -698.515, 3007.58, 500, 1750, -500, 2278.365, -798.624, 
3504.31, 1500, 2750), .Dim = c(2L, 10L), .Dimnames = list(c("Low", 
"High"), c("CD3", "CD127", "CD197 CCR7", "CD279 PD1", "CD8", 
"CD4", "CD28", "CD25", "HLA DR", "CD45RA")))

faust(
    gatingSet=gs,
    activeChannels=activeChannelsIn,
    channelBounds=cbIN,
    startingCellPop="L",
    projectPath=normalizePath("."),
    experimentalUnit="name",
    depthScoreThreshold = 0.01,
    selectionQuantile = 1.0,
    debugFlag = TRUE,
    seedValue = 12345,
    annotationsApproved = FALSE,
    drawAnnotationHistograms = FALSE,
    threadNum=16,
    plottingDevice = "pdf"
)
