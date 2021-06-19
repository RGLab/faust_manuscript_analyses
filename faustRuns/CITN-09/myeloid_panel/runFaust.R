library(faust)
library(flowWorkspace)
library(xtable)

gs <- load_gs(file.path(normalizePath("../../.."),
                        "dataSets",
                        "publication_gating_sets",
                        "CITN-09",
                        "myeloid_panel",
                        "citn09_myeloid_gs"))

cbList <- list(
    imputation01 = structure(c(1000, Inf, 1000, Inf, 1000, Inf, 
                               1000, Inf, 2000, Inf, -1000, 2100,
                               -Inf, Inf, 1000, Inf, 1000, 
                               Inf, 1750, Inf, -1000, 3750),
                             .Dim = c(2L, 11L),
                             .Dimnames = list(
                                 c("Low", "High"),
                                 c("CD11B", "CD20", "CD14", "CD11C", "CD56", 
                                   "CD33", "CD16", "CD3", "CD15", "CD19", "HLA DR"))),
    imputation02 = structure(c(-1500, Inf, -1000, Inf, 1000, Inf,
                               -5000, Inf, -2000, Inf, 2200, Inf, 
                               1000, Inf, -3000, Inf, 1000, Inf, -1000,
                               Inf, 1000, 3750),
                             .Dim = c(2L, 11L),
                             .Dimnames = list(
                                 c("Low", "High"),
                                 c("CD11B", "CD20", "CD14", "CD11C",
                                   "CD56", "CD33", "CD16", "CD3", "CD15", "CD19", "HLA DR"))),
    imputation03 = structure(c(-1500, Inf, 0, Inf, 500, Inf, -5000, Inf,
                               -2000, Inf, 0, Inf, 1000, Inf, -5000, Inf, 0, Inf, 
                               -1000, Inf, -500, 3750),
                             .Dim = c(2L, 11L),
                             .Dimnames = list(
                                 c("Low", "High"),
                                 c("CD11B", "CD20", "CD14", "CD11C", "CD56", 
                                   "CD33", "CD16", "CD3", "CD15", "CD19", "HLA DR"))),
    imputation04 = structure(c(0, Inf, -2500, Inf, 1000, Inf, 1000, Inf, 2300, Inf, -20000, Inf, 
                               1000, Inf, 2000, Inf, 0, Inf, -1000, Inf, -2500, 3750),
                             .Dim = c(2L, 11L),
                             .Dimnames = list(
                                 c("Low", "High"),
                                 c("CD11B", "CD20", "CD14", "CD11C", "CD56",
                                   "CD33", "CD16", "CD3", "CD15", "CD19", "HLA DR"))))
                                     
reportMat <- Reduce(rbind,cbList)
rownames(reportMat) <- paste0(c("IH01_","IH01_","IH02_","IH02_","IH03_","IH03_","IH04_","IH04_"),rownames(reportMat))
xtable(t(reportMat),digits=0)

activeChannelsIn <- c("CD11B","CD20","CD14","CD11C","CD56","CD33","CD16","CD3","CD15","CD19","HLA DR")
supervisedList <- list(
    CD33 = list(actionType = "Preference", action = 2),
    `HLA DR` = list(actionType = "Preference", action = 2),
    CD15 = list(actionType = "Preference", action = 1)
)

faust(
    gatingSet=gs,
    activeChannels=activeChannelsIn,
    channelBounds=cbList,
    startingCellPop="45+",
    experimentalUnit="name",
    imputationHierarchy="imputationID",
    debugFlag = TRUE,
    seedValue = 12345,
    threadNum=11,
    annotationsApproved=TRUE,
    supervisedList=supervisedList,
    drawAnnotationHistograms=FALSE,
    nameOccuranceNum=9,
    plottingDevice = "pdf"
)
