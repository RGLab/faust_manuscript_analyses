library(faust)
library(flowWorkspace)
library(openCyto)
library(flowCore)
#
#read in data and create gating set
#
proj_path <- file.path(normalizePath("../.."),
                       "dataSets",
                       "publication_gating_sets",
                       "melanoma_cytof")
allF <- list.files(file.path(proj_path,"FlowRepository_FR-FCM-ZYKP_files"))
ffNames <- allF[grepl(".fcs",allF)]
ffRef <- ffIn <- read.FCS(filename=file.path(proj_path,"FlowRepository_FR-FCM-ZYKP_files",ffNames[1]))
pmData <- ffRef@parameters@data
activeChannels <- pmData[-which(is.na(pmData[,"desc"])),"desc"]
ffList <- list()
for (ffn in ffNames) {
    ffIn <- read.FCS(filename=file.path(proj_path,"FlowRepository_FR-FCM-ZYKP_files",ffn))
    rData <- flowCore::exprs(ffIn)
    colnames(rData)[!is.na(pmData[,"desc"])] <- pmData[!is.na(pmData[,"desc"]),"desc"]
    for (channel in activeChannels) {
        rData[,channel] <- asinh(rData[,channel,drop=TRUE]/5)
    }
    ffTr <- flowCore::flowFrame(rData)
    ffList <- append(ffList,ffTr)
    names(ffList)[length(ffList)] <- ffn
}
fs <- as(ffList, "flowSet")
gsIn <- GatingSet(fs)
#
#subset to unstimulated samples
#
gs <- gsIn[grepl("US",sampleNames(gsIn))]
#
#Replicate gating strategy reported in figure 1. First, gate cells.
#
gs_add_gating_method(
    gs=gs,
    alias = "cells",
    pop = "+",
    dims = "140Ce_Beads",
    parent = "root",
    gating_method = "boundary",
    gating_args = "min=-Inf,max = 6"
)
#
#Next, add a singlet gate to each sample, translating the template as needed.
#
singletGateTemplate <- new("ellipsoidGate",
                           mean = c(`191Ir_DNA1` = 4.0,
                                    `193Ir_DNA2` = 4.6),
                           cov = structure(c(0.025, 0.015,
                                             0.015,0.02),
                                           .Dim = c(2L, 2L),
                                           .Dimnames = list(c("191Ir_DNA1", "193Ir_DNA2"),
                                                            c("191Ir_DNA1", "193Ir_DNA2"))),
                           distance = 7,
                           parameters = new("parameters", 
    .Data = list(new("unitytransform", .Data = function () 
    NULL, parameters = "191Ir_DNA1", transformationId = "defaultUnityTransform"), 
        new("unitytransform", .Data = function () 
        NULL, parameters = "193Ir_DNA2", transformationId = "defaultUnityTransform"))), 
    filterId = "singlets")
singletGatingList <- list()
for (sampleName in sampleNames(gs)) {
    singletGatingList <- append(singletGatingList,list(singletGateTemplate))
    names(singletGatingList)[length(singletGatingList)] <- sampleName
}
singletGatingList[["000-US_01_normalized.fcs"]]@mean <- c(4.1,4.7)
singletGatingList[["039-US_01_normalized.fcs"]]@mean <- c(4.6,5.2)
singletGatingList[["097-US_01_normalized.fcs"]]@mean <- c(4.4,5.0)
singletGatingList[["103-US_01_normalized.fcs"]]@mean <- c(4.1,4.7)
singletGatingList[["146-US_cct_normalized.fcs"]]@mean <- c(3.85,4.45)
singletGatingList[["162-US_01_normalized.fcs"]]@mean <- c(4.7,5.3)
singletGatingList[["165-US_01_normalized.fcs"]]@mean <- c(4.4,5.0)
singletGatingList[["191-US_01_normalized.fcs"]]@mean <- c(4.7,5.3)
singletGatingList[["203-US_01_normalized.fcs"]]@mean <- c(4.2,4.8)
singletGatingList[["240-US_cct_normalized.fcs"]]@mean <- c(4.2,4.8)
singletGatingList[["246-US_01_normalized.fcs"]]@mean <- c(4.4,5.0)
singletGatingList[["257-US_cct_normalized.fcs"]]@mean <- c(4.3,4.9)
singletGatingList[["266_US_01_normalized.fcs"]]@mean <- c(4.15,4.75)
singletGatingList[["287-US_01_normalized.fcs"]]@mean <- c(4.7,5.3)
singletGatingList[["343-US_01_normalized.fcs"]]@mean <- c(4.4,5.0)
singletGatingList[["346-US_01_normalized.fcs"]]@mean <- c(3.9,4.5)
singletGatingList[["427-US_01_normalized.fcs"]]@mean <- c(4.3,4.9)
singletGatingList[["475-US_cct_normalized.fcs"]]@mean <- c(4.15,4.75)
singletGatingList[["489-US_01_normalized.fcs"]]@mean <- c(4.5,5.1)
singletGatingList[["598-US_01_normalized.fcs"]]@mean <- c(4.1,4.7)
singletGatingList[["632-US_01_normalized.fcs"]]@mean <- c(4.1,4.7)
singletGatingList[["669-US_01_normalized.fcs"]]@mean <- c(4.05,4.65)
singletGatingList[["685-US_01_normalized.fcs"]]@mean <- c(4.2,4.8)
singletGatingList[["702-US_01_normalized.fcs"]]@mean <- c(4.2,4.8)
singletGatingList[["711_US_01_normalized.fcs"]]@mean <- c(4.4,5.0)
singletGatingList[["712-US_01_normalized.fcs"]]@mean <- c(4.25,4.85)
singletGatingList[["748-US_cct_normalized.fcs"]]@mean <- c(4.4,5.0)
singletGatingList[["809-US_01_normalized.fcs"]]@mean <- c(4.8,5.4)
singletGatingList[["843-US_01_normalized.fcs"]]@mean <- c(4.4,5.0)
singletGatingList[["883_US_01-concat_normalized.fcs"]]@mean <- c(4.4,5.0)
singletGatingList[["929-US_01-concat_normalized.fcs"]]@mean <- c(4.4,5.0)
singletGatingList[["999-US_01_normalized.fcs"]]@mean <- c(4.05,4.65)
singletGatingList[["C334-US_01_normalized.fcs"]]@mean <- c(4.225,4.825)
singletGatingList[["C335-US_01_normalized.fcs"]]@mean <- c(4.225,4.825)
singletGatingList[["C351-US_01_normalized.fcs"]]@mean <- c(3.85,4.45)
singletGatingList[["C360-US_01_normalized.fcs"]]@mean <- c(4.15,4.75)
singletGatingList[["C374-US_cct_normalized.fcs"]]@mean <- c(3.5,4.1)
singletGatingList[["C392-US_cct_normalized.fcs"]]@mean <- c(4.25,4.85)
singletGatingList[["C394-US_01_normalized.fcs"]]@mean <- c(4.2,4.8)
singletGatingList[["C403-US_01_normalized.fcs"]]@mean <- c(4.4,5.0)
singletGatingList[["C412-US_01_normalized.fcs"]]@mean <- c(4.25,4.85)
singletGatingList[["C56-US_01_normalized.fcs"]]@mean <- c(4.2,4.8)
singletGatingList[["P239-US_cct_normalized.fcs"]]@mean <- c(4.15,4.75)
singletGatingList[["P245-US_01_normalized.fcs"]]@mean <- c(4.2,4.8)
singletGatingList[["P262-US_cct_normalized.fcs"]]@mean <- c(4.4,5.0)
for (sampleName in sampleNames(gs)) {
    gs_pop_add(gs[sampleName],singletGatingList[[sampleName]], parent = "cells", name = "singlets")
    recompute(gs[sampleName])
}
#
#intact gate
#
gs_add_gating_method(
    gs,
    alias = "intact",
    pop = "+",
    dims = "Event_length",
    parent = "singlets",
    gating_method = "boundary",
    gating_args = "min=-Inf,max = 20"
)
#
#live gate
#
gs_add_gating_method(
    gs,
    alias = "live",
    pop = "-",
    dims = "115In_Dead",
    parent = "intact",
    gating_method = "gate_quantile",
    gating_args = "probs = 0.925"
)
#
#specify the starting node for the faust anlaysis here.
#
startingNode <- "live"

#
#specify markers and bounds
#
activeChannels <- c("145Nd_CD4","146Nd_CD8","141Pr_CD25","153Eu_CD45RA","165Ho_CD127",
                    "169Tm_CCR7","154Sm_CD3","157Gd_HLA-DR","170Er_PD-1","155Gd_CD28")

channelBoundsIn <- structure(c(0.75, 3.75, 1, 4, 0.25, 1.5, 1, 5, 0, 5.318, 0, 2.75, 
0, 6.229, 2, 6.25, 0.25, 2, 0.25, 4.535), .Dim = c(2L, 10L), .Dimnames = list(
    c("Low", "High"), c("145Nd_CD4", "146Nd_CD8", "141Pr_CD25", 
    "153Eu_CD45RA", "165Ho_CD127", "169Tm_CCR7", "154Sm_CD3", 
    "157Gd_HLA-DR", "170Er_PD-1", "155Gd_CD28")))

#
#run faust.
#
faust(
    gatingSet=gs,
    activeChannels=activeChannels,
    channelBounds=channelBoundsIn,
    startingCellPop=startingNode,
    depthScoreThreshold=0.01,
    selectionQuantile=1.0,
    threadNum=16,
    debugFlag=TRUE,
    seedValue=12345,
    annotationsApproved=FALSE,
    plottingDevice = "pdf"
)



