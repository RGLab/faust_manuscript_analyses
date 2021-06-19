library(ggplot2)
library(cowplot)
library(flowWorkspace)
library(tidyr)
library(readxl)
library(dplyr)
library(stringr)
library(viridis)
set.seed(71823)
#
#load subject metadata
#
dataPath <- file.path(normalizePath("../.."),"faustRuns","CITN-09","whole_blood_analysis")
gs <- load_gs(file.path(normalizePath("../.."),
                        "dataSets",
                        "publication_gating_sets",
                        "CITN-09",
                        "whole_blood_samples",
                        "citn09_fresh_whole_blood_gs"))
meta_data_all <- as.data.frame(pData(gs))
meta_data_all$subjectID <- meta_data_all$publicationID
meta_data <- meta_data_all[which(meta_data_all$VISIT=="C01"),]
bioData <- meta_data[,c("subjectID","tumorViralStatus","OverallResponse")]
#
#load faust count matrix. merge metadata and derive variables.
#
tcellCountDF <- as.data.frame(readRDS(file.path(dataPath,"faustData","faustCountMatrix.rds")))
allTcellPanelPops <- colnames(tcellCountDF)
tcellPops <- allTcellPanelPops[grepl("CD3\\+",allTcellPanelPops)]
tcellCountDF$parentCount <- apply(tcellCountDF[,tcellPops],1,sum)
tcellCountDF$sampleName <- rownames(tcellCountDF)
tcellCountDF$subjectID <- as.character(sapply(rownames(tcellCountDF),function(x){str_extract(x,"Subject-\\d\\d")}))
tcellCountDF$visitSTR <- as.character(sapply(rownames(tcellCountDF),function(x){str_extract(x,"C\\d\\d|EOT")}))
tcellAnalysisDF <- left_join(tcellCountDF,bioData,by="subjectID")
#
#time dynamics
#
brightPops <- tcellPops[grepl("PD1Bright",tcellPops)]
brightCD8 <- brightPops[grepl("CD8\\+",brightPops)]
brightCD8 <- brightCD8[grepl("CD4-",brightCD8)]
allPlotDF <- tcellAnalysisDF[,c(brightCD8,"parentCount","visitSTR","OverallResponse","subjectID")]
newORStr <- tcellAnalysisDF[,"OverallResponse"]
newORStr[which(newORStr == "PD")] <- "Non-Responder"
newORStr[which(newORStr == "SD")] <- "Non-Responder"
newORStr[which(newORStr == "PR")] <- "Responder"
newORStr[which(newORStr == "CR")] <- "Responder"
allPlotDF$resType <- newORStr
allPlotDF$brightCounts <- apply(allPlotDF[,brightCD8],1,sum)
plotDF <- data.frame(
    pct=(100*(allPlotDF[,"brightCounts"]/allPlotDF[,"parentCount"])),
    visit=as.factor(allPlotDF[,"visitSTR"]),
    Response=allPlotDF[,"resType"],
    subjectID=allPlotDF[,"subjectID"],
    stringsAsFactor=FALSE
)
p <- ggplot(plotDF,aes(x=visit,y=pct))+
    geom_line(aes(group=subjectID,color=Response),size=1.0)+
    scale_color_manual(values=c("#000000","#000000"))+
    theme_classic(base_size = 14)+
    scale_linetype_manual(values=c("solid", "dashed"))+
    theme(
        legend.position="none",
        legend.text=element_text(size=10),
        plot.title=element_text(size=14)
    )+
    scale_y_continuous(labels=function(x) paste0(x,"%"))+
    ggtitle("")+
    xlab("Immunotherapy Cycle")+
    ylab("Aggregate CD8+ PD-1 Bright CD3+ CD4-\nCells % of CD3+ by sample")+
    guides(linetype=guide_legend(keywidth = 5, keyheight = 1))
saveRDS(p,file.path(normalizePath("."),"artifacts","longitudinal_figure_panel_a.rds"))
