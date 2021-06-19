library(ggplot2)
library(ggpubr)
library(rstatix)
library(cowplot)
library(lme4)
library(tidyr)
library(readxl)
library(dplyr)
library(stringr)
library(viridis)
library(flowWorkspace)
library(faust)
library(xtable)
library(grid)
library(gridExtra)
#library(devtools)
#devtools::install_github("clauswilke/gridtext")
library(gridtext)
library(rlang)
set.seed(71823)
source("./helperFunctions.R")

#
#global plot parameters 
#
global_plot_font_size <- 6
global_plot_label_size <- 10
global_plot_point_size <- 2
citn09DataPath <- file.path(normalizePath("../.."),"faustRuns","CITN-09","whole_blood_analysis")
#
#logging function
#
paperLog <- function(logExprs,logDesc) {
    print("")
    print("****************************************************************")
    print("Paper log")
    print(logDesc)
    print(logExprs)
    print("****************************************************************")
    print("")
}
#
#load previously derived metadata from the gating set and subset 
#to baseline for modeling
#
gs <- load_gs(file.path(normalizePath("../.."),
                        "dataSets",
                        "publication_gating_sets",
                        "CITN-09",
                        "whole_blood_samples",
                        "citn09_fresh_whole_blood_gs"))
meta_data_all <- as.data.frame(pData(gs))
meta_data_all$subjectID <- meta_data_all$publicationID
meta_data <- meta_data_all[which(meta_data_all$VISIT=="C01"),]
#suppress warnings -- coercing from character "NA" to numeric NA 
suppressWarnings(meta_data[,"PD1_Total"] <- as.numeric(meta_data[,"PD1_Total"]))
suppressWarnings(meta_data[,"productiveClonality"] <- as.numeric(meta_data[,"productiveClonality"]))
#
#load faust count matrix. derive metadata.
#
allTcellCountDF <- as.data.frame(readRDS(file.path(citn09DataPath,"faustData","faustCountMatrix.rds")))
paperLog(quantile(1-(allTcellCountDF[,"0_0_0_0_0"]/apply(allTcellCountDF,1,sum))),"Annotation quantiles")
allTcellPanelPops <- colnames(allTcellCountDF)
paperLog((length(allTcellPanelPops)-1),"Total number of cell populations")
tcellPops <- names(which(sapply(allTcellPanelPops,function(x){grepl("CD3\\+",x)})))
paperLog(length(tcellPops),"Number of tcell populations")
allTcellCountDF$sampleName <- rownames(allTcellCountDF)
allTcellCountDF$subjectID <- str_extract(rownames(allTcellCountDF),"Subject-\\d\\d")
allTcellCountDF$visitSTR <- str_extract(rownames(allTcellCountDF),"C\\d\\d|EOT")
#
#get exact total CD3 counts
#
sel_meta_data <- data.frame(
    subjectID=as.factor(allTcellCountDF[,"subjectID"]),
    sampleName=allTcellCountDF[,"sampleName"]
)
cd3CountDF <- getCountsForTargetMarkers(
    projectPath=citn09DataPath,
    referencePhenotype=tcellPops[1],
    targetMarkers=c("CD3","CD4","CD8"),
    metaDataDF=sel_meta_data,
    markersInParentPhenotype=c("CD3")
)
cd3ParentDF <- cd3CountDF[,c("sampleName","parentCount")]
allTcellCountDF <- left_join(allTcellCountDF,cd3ParentDF,by=c("sampleName"))
#
#subset to baseline for modeling and derive CR/PR merged responder type.
#
tcellCountDF <- allTcellCountDF[which(allTcellCountDF[,c("visitSTR")]=="C01"),]
tcellAnalysisDF <- inner_join(tcellCountDF,meta_data,by="subjectID")
paperLog(length(table(tcellAnalysisDF[,"subjectID"])),"number of patients across trial.")
paperLog(table(tcellAnalysisDF[,"subjectID"]),"number of samples per patient.")
paperLog(dim(tcellAnalysisDF),"Dimension pre-second-subset")
tcellSubDF <- tcellAnalysisDF[which(tcellAnalysisDF[,c("visitSTR")]=="C01"),]
paperLog(dim(tcellAnalysisDF),"Dimension post-second-subset (expect same)")
tcellResType <- rep(0,nrow(tcellSubDF))
tcellResType[which(tcellSubDF$OverallResponse %in% c("CR","PR"))] <- 1
tcellSubDF$modelRT <- tcellResType
paperLog(table(tcellSubDF$modelRT),"response rate.")
#
#specify the GLMM
#
tcellSafeMod <- function(dataSet) {
    out <- tryCatch(
    {
        m <- glmer(cbind(childCount,(parentCount-childCount)) ~ resType + (1|subjectID),
                   data=dataSet,
                   family="binomial",
                   control=glmerControl(
                       optimizer="bobyqa", optCtrl = list(maxfun = 1e9),
                       check.conv.singular = .makeCC(action = "warning",  tol = 1e-4))
                   )
        rv <- c(coefficients(summary(m))[2,4])
        return(rv)
    },
    error=function(cond){
        message("Error!")
        message(cond)
        return(NA)
    },
    warning=function(cond){
        message("Warning!")
        message(cond)
        return(NA)
    },
    finally={
        message("done!")
    }
    )
    return(out)
}
#
#fit model to all tcell populations
#
mePvalVec <- c()
for (cellPop in tcellPops) {
    mdf <- data.frame(resType=as.factor(tcellSubDF[,"modelRT"]),
                      subjectID=as.factor(tcellSubDF[,"subjectID"]),
                      parentCount=tcellSubDF[,"parentCount"],
                      childCount=tcellSubDF[,cellPop],
                      sampleName=tcellSubDF[,"sampleName"])
    pv <- tcellSafeMod(mdf)
    if (!is.na(pv)) {
        mePvalVec <- append(mePvalVec,pv[1])
        names(mePvalVec)[length(mePvalVec)] <- paste0(cellPop)
    }
}
qp <- p.adjust(mePvalVec,"bonferroni")
selRTNames <- sort(unique(names(which(qp < 0.1))))
paperLog(sort(qp[selRTNames]),"selected populations with bonferroni p")
methodReportDF <- data.frame(
    phenotypes=names(sort(qp[selRTNames])),
    bonferroni.p=as.numeric(sort(qp[selRTNames]))
)
xtable(methodReportDF,digits=4)
#
#report all phenotypes with FDR < 0.20 for supplement
#
selRTNames <- names(sort(qp[selRTNames]))
qp <- p.adjust(mePvalVec,"BH")
twentyNamesDF <- data.frame(
    Phenotype=names(sort(qp[which(qp < 0.2)])),
    fdr.adjusted.p.value=as.numeric(sort(qp[which(qp < 0.2)])),
    stringsAsFactors=FALSE
)
twentyNamesDF$Phenotype <- gsub("HLA DR","HLADR",twentyNamesDF$Phenotype)
twentyNamesDF$Phenotype <- gsub("CD279 PD1","PD1",twentyNamesDF$Phenotype)
twentyNamesDF$Phenotype <- gsub("CD197 CCR7","CCR7",twentyNamesDF$Phenotype)
twentyNamesDF$Phenotype <- gsub("CD278 ICOS","ICOS",twentyNamesDF$Phenotype)
twentyNamesDF$Phenotype <- gsub("PD1Dim","PD1 Dim ",twentyNamesDF$Phenotype)
twentyNamesDF$Phenotype <- gsub("PD1Bright","PD1 Bright ",twentyNamesDF$Phenotype)
twentyNamesDF$Phenotype <- gsub("CD4Dim","CD4 Dim ",twentyNamesDF$Phenotype)
twentyNamesDF$Phenotype <- gsub("CD4Bright","CD4 Bright ",twentyNamesDF$Phenotype)
twentyNamesDF$Phenotype <- gsub("-","- ",twentyNamesDF$Phenotype)
twentyNamesDF$Phenotype <- gsub("\\+","\\+ ",twentyNamesDF$Phenotype)
print(xtable(twentyNamesDF,digits=3,caption="FAUST Phenotypes associated with outcome in the CITN-09 T cell data at the FDR-adjusted 20\\% level."),include.rownames=FALSE)
twentyNamesDF$Ranking <- seq(nrow(twentyNamesDF))
print(xtable(twentyNamesDF[,c(1,3)],caption="FAUST Phenotypes associated with outcome in the CITN-09 T cell data at the FDR-adjusted 20\\% level."),include.rownames=FALSE)
#
#generate table of effect sizes with 95% CI for supplement
#
qp <- p.adjust(mePvalVec,"bonferroni")
ciList <- list()
for (pName in selRTNames) {
    mdf <- data.frame(
        resType=as.factor(tcellSubDF[,"modelRT"]),
        subjectID=as.factor(tcellSubDF[,"subjectID"]),
        parentCount=tcellSubDF[,"parentCount"],
        childCount=tcellSubDF[,pName],
        sampleName=tcellSubDF[,"sampleName"]
    )
    
    m <- glmer(cbind(childCount,(parentCount-childCount)) ~ resType + (1|subjectID),
               data=mdf,
               family="binomial",
               control=glmerControl(
                   optimizer="bobyqa", optCtrl = list(maxfun = 1e9),
                   check.conv.singular = .makeCC(action = "warning",  tol = 1e-4))
               )
    ciList <- append(ciList,list(cbind(fixef(m),confint(m)[2:3,])))
    names(ciList)[length(ciList)] <- pName
}

print(xtable(
data.frame(
    faustCluster=names(ciList),
    pointEstimate=as.numeric(unlist(lapply(ciList,function(x){x[2,1]}))),
    CI_95_Lower=as.numeric(unlist(lapply(ciList,function(x){x[2,2]}))),
    CI_95_Upper=as.numeric(unlist(lapply(ciList,function(x){x[2,3]})))
),
digits=3))
#
#compute interaction model
#
interactList <- list()
for (pName in selRTNames) {
    mdf <- data.frame(
        resType=as.factor(tcellSubDF[,"modelRT"]),
        viralStatus=as.factor(tcellSubDF[,"tumorViralStatus"]),
        subjectID=as.factor(tcellSubDF[,"subjectID"]),
        parentCount=tcellSubDF[,"parentCount"],
        childCount=tcellSubDF[,pName],
        sampleName=tcellSubDF[,"sampleName"]
    )
    
    m <- glmer(cbind(childCount,(parentCount-childCount)) ~ resType*viralStatus + (1|subjectID),
               data=mdf,
               family="binomial",
               control=glmerControl(
                   optimizer="bobyqa", optCtrl = list(maxfun = 1e9),
                   check.conv.singular = .makeCC(action = "warning",  tol = 1e-4))
               )
    
    interactList <- append(interactList,list(coefficients(summary(m))[4,4]))
    names(interactList)[length(interactList)] <- pName
}
#
#get the baseline data for the populations significant at Bonferroni 0.10 level
#
baselineData <- tcellSubDF[,c("modelRT","OverallResponse","subjectID","sampleName",selRTNames,"parentCount",
                              "tumorViralStatus","PD1_Total","productiveClonality")]
Response <- rep("Non-Responder",nrow(baselineData))
Response[which(baselineData$modelRT==1)] <- "Responder"
baselineData$Response <- as.factor(Response)
baselineData$Response <- relevel(baselineData$Response,"Non-Responder")
#
#plot baseline correlations with external data sources
#      
plotListTitle <- plotList1 <- plotList0 <- list()
set.seed(777)
for (popNum in seq(length(selRTNames))) {
    cellPop <- selRTNames[popNum]
    plotName <- cellPop
    stat.test <- data.frame(
        group1=c("Non-Responder"),
        group2=c("Responder"),
        obsp=round(as.numeric(mePvalVec[cellPop]),6)
    )
    plotName <- gsub("HLA DR","HLADR",plotName)
    plotName <- gsub("CD279 PD1","PD1",plotName)
    plotName <- gsub("CD197 CCR7","CCR7",plotName)
    plotName <- gsub("CD278 ICOS","ICOS",plotName)
    plotName <- gsub("PD1Dim","PD1 Dim ",plotName)
    plotName <- gsub("CD4Dim","CD4 Dim ",plotName)
    plotName <- gsub("PD1Bright","PD1 Bright ",plotName)
    plotName <- gsub("CD4Bright","CD4 Bright ",plotName)
    plotName <- gsub("-","- ",plotName)
    plotName <- gsub("\\+","\\+ ",plotName)
    plotName <- gsub(" PD1","<br>PD1",plotName)
    plotName <- gsub("CD4 Bright","**CD4 Bright**",plotName)
    plotName <- gsub("PD1 Bright","**PD1 Bright**",plotName)
    plotName <- gsub("PD1 Dim","**PD1 Dim**",plotName)
    plotName <- gsub("CD8\\+","**CD8\\+**",plotName)
    plotName <- gsub("HLADR\\+","**HLADR\\+**",plotName)
    plotName <- gsub("CD28\\+","**CD28\\+**",plotName)
    plotName <- gsub("CD3\\+","**CD3\\+**",plotName)
    baselineData$prop <- (baselineData[,cellPop]/baselineData[,"parentCount"])
    baselineData$pct <- (100*baselineData$prop) 
    #
    #viral status plot
    #
    boxDF <- na.omit(baselineData[,c("tumorViralStatus","pct","Response",cellPop,"subjectID","parentCount")])
    colnames(boxDF)[4] <- "Count"
    boxDF$ResVS <- as.factor(paste0(as.character(boxDF[,"Response"]),
                                    "\n",
                                    as.character(boxDF[,"tumorViralStatus"])))
    boxDF$ResVS <- factor(boxDF$ResVS,
                          levels = c("PD/SD\nNegative",
                                     "PD/SD\nPositive",
                                     "CR/PR\nNegative",
                                     "CR/PR\nPositive"))
    p0p <- ggplot(boxDF) +
        geom_boxplot(
            outlier.color=NA,
            aes(
                x=Response,
                y=pct,
                fill=Response
            )
        )+
        geom_point(
            position=position_jitterdodge(),
            aes(
                x=Response,
                y=pct,
                fill=Response,
                shape=tumorViralStatus
            ),
            size=global_plot_point_size ,
            stroke=0.5
        )+
        scale_fill_manual(
            values=c("#932667FF","#F6D645FF"),
            guide=FALSE
        )+
        scale_shape_manual(values=c(2,19))+
        theme_classic(base_size = global_plot_font_size)+
        theme(
            legend.position="none",
            plot.title=element_markdown(
                size = global_plot_font_size,
                hjust=0.5,
                colour = "#000000",
                lineheight=1.2
            )
        )+
        xlab("")+
        guides(shape=guide_legend("Tumor viral status"))+ 
        ggtitle(plotName)+
        scale_size_continuous(range=c(1,10))+
        scale_y_continuous(trans="sqrt",
                           labels=function(x){paste0(x,"%")},
                           breaks = c(0,0.00005,0.0005,0.005,round(seq(0,(max(boxDF$pct) * 1.1),by=(max(boxDF$pct) * 1.1)/7),2)))+
        scale_x_discrete(labels=c("Non-responder","Responder"))+
        coord_cartesian(ylim=c(0,(max(boxDF$pct) * 1.25)))+
        stat_pvalue_manual(stat.test,y.position=max(0.45,(max(boxDF$pct) * 1.01)),label="obsp",label.size=2)+
        ylab("")
    plotList0 <- append(plotList0,list(p0p))
    stat.test.2 <- data.frame(
        group1=0.5,
        group2=2.5,
        obsp=round(as.numeric(mePvalVec[cellPop]),6)
    )
    p1p <- ggplot(boxDF) +
        geom_boxplot(
            outlier.color=NA,
            aes(
                x=tumorViralStatus,
                y=pct,
                fill=interaction(tumorViralStatus,Response)
            )
        )+
        geom_point(
            position=position_jitterdodge(),
            aes(
                x=tumorViralStatus,
                y=pct,
                fill=interaction(tumorViralStatus,Response),
                shape=tumorViralStatus
            ),
            size=global_plot_point_size ,
            stroke=0.5
        )+
        scale_fill_manual(
            values=c("#932667FF","#932667FF","#F6D645FF","#F6D645FF"),
            guide=FALSE
        )+
        scale_shape_manual(values=c(2,19,2,19))+
        theme_classic(base_size = global_plot_font_size)+
        theme(
            legend.position="none",
            plot.title=element_markdown(
                size = global_plot_font_size,
                hjust=0.5,
                colour = "#000000",
                lineheight=1.2
            )
        )+
        xlab("")+
        guides(shape=guide_legend("Tumor viral status"))+ 
        ggtitle(plotName)+
        scale_size_continuous(range=c(1,10))+
        scale_y_continuous(trans="sqrt",
                           labels=function(x){paste0(x,"%")},
                           breaks = c(0,0.00005,0.0005,0.005,round(seq(0,(max(boxDF$pct) * 1.1),by=(max(boxDF$pct) * 1.1)/7),2)))+
        scale_x_discrete(labels=c("Virus negative","Virus positive"))+
        coord_cartesian(ylim=c(0,(max(boxDF$pct) * 1.25)))+
        stat_pvalue_manual(stat.test.2,y.position=max(0.45,(max(boxDF$pct) * 1.01)),label="obsp",size=2)+
        ylab("")
    plotList1 <- append(plotList1,list(p1p))
    pName <- ggdraw()+draw_label(plotName,size=global_plot_font_size)
    plotListTitle <- append(plotListTitle,list(pName))
}
#
#
#compute first CD8 ihc correlations
#
#
baseIHC <- na.omit(baselineData[,c("tumorViralStatus","PD1_Total","CD4-CD8+CD3+CD45RA-HLA DR+CD279 PD1DimCD28+CD127-CD25-CD197 CCR7-","parentCount","modelRT")])
baseIHC$pct <- (100*(baseIHC[,"CD4-CD8+CD3+CD45RA-HLA DR+CD279 PD1DimCD28+CD127-CD25-CD197 CCR7-"]/baseIHC[,"parentCount"]))
ihcCorDF <- baseIHC %>% group_by(tumorViralStatus) %>% summarise(correlation=cor(pct,PD1_Total))
baseIHCSub <- baseIHC[which(baseIHC$tumorViralStatus=="Positive"),]
p3p <- ggplot(baseIHC)+
    geom_point(
        aes(
            x=PD1_Total,
            y=pct,
            color=as.factor(modelRT),
            shape=tumorViralStatus
        ),
        stroke=0.5,
        size=global_plot_point_size 
    )+
    theme_classic(base_size=global_plot_font_size)+
    theme(
        legend.position="none",
        plot.title=element_markdown(
            size=global_plot_font_size,
            hjust=0.5,
            colour = "#000000",
            lineheight = 1.2
        )
    )+
    scale_color_viridis(begin=0.4,end=0.9,option="B",discrete=TRUE)+
    xlab("Total tumor PD-1\n(median number of positive cells/mm^2)")+
    scale_shape_manual(values=c(2,19))+
    guides(shape=guide_legend("Tumor viral status"))+
    ggtitle("")+
    geom_smooth(
        data=baseIHCSub,
        method="lm",
        aes(x=PD1_Total,y=pct),
        color="gray",
        linetype="dashed"
    )+
    scale_size_continuous(range=c(1,2))+
    coord_cartesian(ylim=c(0,(max(baseIHC$pct) * 1.2)))+
    scale_y_continuous(trans="sqrt",labels=function(x){paste0(x,"%")},
                       breaks = c(0.000005,round(seq(0,(max(baseIHC$pct) * 1.1),by=(max(baseIHC$pct) * 1.1)/10),2)))+
    ylab("")
ihcLabelPos <- paste0("Correlation within\nvirus positive subjects:\n",
                      round(as.numeric(ihcCorDF[which(ihcCorDF[,1]=="Positive"),"correlation"]),3))
ihcDim <- ggdraw(p3p+ggtitle(""))+
    draw_label(ihcLabelPos,x=0.4,y=0.7,size=global_plot_font_size)
#
#compute second CD8 ihc correlations
#
baseIHC <- na.omit(baselineData[,c("tumorViralStatus","PD1_Total","CD4-CD8+CD3+CD45RA-HLA DR+CD279 PD1BrightCD28+CD127-CD25-CD197 CCR7-","parentCount","modelRT")])
baseIHC$pct <- (100*(baseIHC[,"CD4-CD8+CD3+CD45RA-HLA DR+CD279 PD1BrightCD28+CD127-CD25-CD197 CCR7-"]/baseIHC[,"parentCount"]))
ihcCorDF <- baseIHC %>% group_by(tumorViralStatus) %>% summarise(correlation=cor(pct,PD1_Total))
baseIHCSub <- baseIHC[which(baseIHC$tumorViralStatus=="Positive"),]
p3p <- ggplot(baseIHC)+
    geom_point(
        aes(
            x=PD1_Total,
            y=pct,
            color=as.factor(modelRT),
            shape=tumorViralStatus
        ),
        size=global_plot_point_size,
        stroke=0.5,
        size=2.5
    )+
    theme_classic(base_size=global_plot_font_size)+
    theme(
        legend.position="none",
        plot.title=element_markdown(
            size=global_plot_font_size,
            hjust=0.5,
            colour = "#000000",
            lineheight = 1.2
        )
    )+
    scale_color_viridis(begin=0.4,end=0.9,option="B",discrete=TRUE)+
    xlab("Total tumor PD-1\n(median number of positive cells/mm^2)")+
    scale_shape_manual(values=c(2,19))+
    guides(shape=guide_legend("Tumor viral status"))+
    ggtitle("")+
    geom_smooth(data=baseIHCSub,
                method="lm",
                aes(x=PD1_Total,
                    y=pct),
                color="gray",
                linetype="dashed")+
    scale_size_continuous(range=c(1,2))+
    coord_cartesian(ylim=c(0,(max(baseIHC$pct) * 1.2)))+
    scale_y_continuous(trans="sqrt",labels=function(x){paste0(x,"%")},
                       breaks = c(0.000005,round(seq(0,(max(baseIHC$pct) * 1.1),by=(max(baseIHC$pct) * 1.1)/10),2)))+
    ylab("")
ihcLabelPos <- paste0("Correlation within\nvirus positive subjects:\n",
                      round(as.numeric(ihcCorDF[which(ihcCorDF[,1]=="Positive"),"correlation"]),3))
ihcBright <- ggdraw(p3p+ggtitle(""))+
    draw_label(ihcLabelPos,x=0.4,y=0.7,size=global_plot_font_size)
#
#compute first CD8 clonality plot
#
baseClone <- na.omit(baselineData[,c("tumorViralStatus","productiveClonality",
                                   "CD4-CD8+CD3+CD45RA-HLA DR+CD279 PD1DimCD28+CD127-CD25-CD197 CCR7-",
                                   "parentCount","modelRT")])
baseClone$pct <- (100*(baseClone[,"CD4-CD8+CD3+CD45RA-HLA DR+CD279 PD1DimCD28+CD127-CD25-CD197 CCR7-"]/baseClone[,"parentCount"]))
cloneCorDF <- baseClone %>% group_by(tumorViralStatus) %>% summarise(correlation=cor(pct,productiveClonality))
baseCloneSub <- baseClone[which(baseClone$tumorViralStatus=="Positive"),]
p4p <- ggplot(baseClone)+
    geom_point(
        aes(
            x=productiveClonality,
            y=pct,
            color=as.factor(modelRT),
            shape=tumorViralStatus
        ),
        size=global_plot_point_size,
        size=2.5,
        stroke=0.5
    )+
    theme_classic(base_size=global_plot_font_size)+
    theme(
        legend.position="none",
        plot.title=element_markdown(
            size=global_plot_font_size,
            hjust=0.5,
            colour = "#000000",
            lineheight = 1.2
        )
    )+
    scale_color_viridis(begin=0.4,end=0.9,option="B",discrete=TRUE)+
    xlab("TCR clonality in the tumor\n(1 - normalized Shannon entropy)")+
    ggtitle("")+
    scale_shape_manual(values=c(2,19))+
    guides(shape=guide_legend("Tumor viral status"))+
    coord_cartesian(ylim=c(0,(max(baseClone$pct) * 1.2)))+
    geom_smooth(data=baseCloneSub,
                method="lm",
                aes(x=productiveClonality,
                    y=pct),
                color="gray",
                linetype="dashed")+
    scale_size_continuous(range=c(1,7))+
    scale_y_continuous(trans="sqrt",labels=function(x){paste0(x,"%")},
                       breaks = c(0.00005,round(seq(0,(max(baseClone$pct) * 1.1),by=(max(baseClone$pct) * 1.1)/10),2)))+
    ylab("")
cloneLabelPos <- paste0("Correlation within\nvirus positive subjects:\n",
                        round(as.numeric(cloneCorDF[which(cloneCorDF[,1]=="Positive"),"correlation"]),3))
cloneDim <- ggdraw(p4p+ggtitle(""))+
    draw_label(cloneLabelPos,x=0.9,y=0.37,size=global_plot_font_size,hjust=1)
#
#compute second CD8 clonality plot
#
baseClone <- na.omit(baselineData[,c("tumorViralStatus","productiveClonality",
                                   "CD4-CD8+CD3+CD45RA-HLA DR+CD279 PD1BrightCD28+CD127-CD25-CD197 CCR7-",
                                   "parentCount","modelRT")])
baseClone$pct <- (100*(baseClone[,"CD4-CD8+CD3+CD45RA-HLA DR+CD279 PD1BrightCD28+CD127-CD25-CD197 CCR7-"]/baseClone[,"parentCount"]))
cloneCorDF <- baseClone %>% group_by(tumorViralStatus) %>% summarise(correlation=cor(pct,productiveClonality))
baseCloneSub <- baseClone[which(baseClone$tumorViralStatus=="Positive"),]
p4p <- ggplot(baseClone)+
    geom_point(
        aes(
            x=productiveClonality,
            y=pct,
            color=as.factor(modelRT),
            shape=tumorViralStatus
        ),
        size=global_plot_point_size,
        stroke=0.5
    )+
    theme_classic(base_size=global_plot_font_size)+
    theme(
        legend.position="none",
        plot.title=element_markdown(
            size=global_plot_font_size,
            hjust=0.5,
            colour = "#000000",
            lineheight = 1.2
        )
    )+
    scale_color_viridis(begin=0.4,end=0.9,option="B",discrete=TRUE)+
    xlab("TCR clonality in the tumor\n(1 - normalized Shannon entropy)")+
    ggtitle("")+
    scale_shape_manual(values=c(2,19))+
    guides(shape=guide_legend("Tumor viral status"))+
    coord_cartesian(ylim=c(0,(max(baseClone$pct) * 1.2)))+
    geom_smooth(data=baseCloneSub,
                method="lm",
                aes(x=productiveClonality,
                    y=pct),
                color="gray",
                linetype="dashed")+
    scale_size_continuous(range=c(1,7))+
    scale_y_continuous(trans="sqrt",labels=function(x){paste0(x,"%")},
                       breaks = c(0.00005,round(seq(0,(max(baseClone$pct) * 1.1),by=(max(baseClone$pct) * 1.1)/10),2)))+
    ylab("")
cloneLabelPos <- paste0("Correlation within\nvirus positive subjects:\n",
                        round(as.numeric(cloneCorDF[which(cloneCorDF[,1]=="Positive"),"correlation"]),3))
cloneBright <- ggdraw(p4p+ggtitle(""))+
    draw_label(cloneLabelPos,x=0.9,y=0.37,size=global_plot_font_size,hjust=1)
#
#create legend for grid plot
#
cellPop <- selRTNames[1]
plotName <- cellPop
plotName <- gsub("HLA DR","HLADR",plotName)
plotName <- gsub("CD279 PD1","PD1",plotName)
plotName <- gsub("CD197 CCR7","CCR7",plotName)
plotName <- gsub("CD278 ICOS","ICOS",plotName)
plotName <- gsub("PD1Dim","PD1 Dim ",plotName)
plotName <- gsub("CD4Dim","CD4 Dim ",plotName)
plotName <- gsub("PD1Bright","PD1 Bright ",plotName)
plotName <- gsub("CD4Bright","CD4 Bright ",plotName)
plotName <- gsub("-","- ",plotName)
plotName <- gsub("\\+","\\+ ",plotName)
plotName <- gsub(" PD1","\\\nPD1",plotName)
baselineData$prop <- baselineData[,cellPop]/baselineData[,"parentCount"]
baselineData$pct <- (100*baselineData$prop) 
boxDF <- na.omit(baselineData[,c("tumorViralStatus","pct","Response",cellPop,"subjectID","parentCount")])
colnames(boxDF)[4] <- "Count"
boxDF$ResVS <- as.factor(paste0(as.character(boxDF[,"Response"]),
                                "\n",
                                as.character(boxDF[,"tumorViralStatus"])))
boxDF$ResVS <- factor(boxDF$ResVS,
                      levels = c("PD/SD\nNegative",
                                 "PD/SD\nPositive",
                                 "CR/PR\nNegative",
                                 "CR/PR\nPositive"))
commonLegend <- ggplot(boxDF) +
    geom_boxplot(
        outlier.color=NA,
        aes(
            x=Response,
            y=pct,
            fill=interaction(Response,tumorViralStatus)
        )
    )+
    geom_point(
        position=position_jitterdodge(),
        aes(
            x=Response,
            y=pct,
            color=Response,
            shape=tumorViralStatus
        ),
        stroke=0.5
    )+
    scale_fill_manual(
        values=c("#932667FF","#F6D645FF","#932667FF","#F6D645FF"),
        guide=FALSE
    )+
    scale_color_manual(values=c("#932667FF","#F6D645FF","#932667FF","#F6D645FF"))+
    scale_shape_manual(values=c(2,19))+
    theme_classic(base_size = global_plot_font_size)+
    theme(legend.position="bottom")+
    xlab("")+
    ylab("FAUST Cluster % of Lymphocytes by Sample")+
    guides(shape=guide_legend("Tumor viral status",override.aes=list(size=global_plot_font_size),ncol=2),
           color=guide_legend("Response",override.aes=list(size=global_plot_font_size),ncol=2))+
    ggtitle("Whole blood flow cytometry FAUST correlate")+
    scale_size_continuous(range=c(1,10))
myLegend <- get_legend(commonLegend)
discTitle <- ggdraw() +
  draw_label(
    "Discovery analysis\nCITN-09 fresh whole blood dataset",
    fontface = 'bold',
    hjust = 0.5
  ) 
#
#
#
#compute targeted frequencies on PBMC dataset
#
#
#
#
#Load the annotation thresholds for the PBMC data
#
set.seed(123)
dataPath <- file.path(normalizePath("../.."),"faustRuns","CITN-09","cryopreserved_pbmcs")
rl <- readRDS(file.path(dataPath,"faustData","gateData","L_resList.rds"))
#
#Target CD8 PD-1+ Counts.
#
firstSample <- TRUE
for (fn in list.files(file.path(dataPath,"faustData","sampleData"))) {
    # 
    #import the expression data for each sample
    #
    fns <- gsub(".fcs_[[:alnum:]]+",".fcs",fn)
    fns <- gsub(".fcs_[[:alnum:]]+",".fcs",fns)
    print(paste0(fns,": ",(fns %in% names(rl[["CD4"]]))))
    exprsMat <- readRDS(file.path(dataPath,"faustData","sampleData",fn,"exprsMat.rds"))
    #  
    #find all cells that match the phenotype of interest
    #
    childLookup <- Reduce(intersect,list(
                                   cd4=which(exprsMat[,"CD4"] < rl[["CD4"]][[fns]][[1]]),
                                   cd8=which(exprsMat[,"CD8"] >= rl[["CD8"]][[fns]]),
                                   cd3=which(exprsMat[,"CD3"] >= rl[["CD3"]][[fns]]),
                                   cd45ra=which(exprsMat[,"CD45RA"] < rl[["CD45RA"]][[fns]]),
                                   dr=which(exprsMat[,"HLA DR"] >= rl[["HLA DR"]][[fns]]),
                                   pd1=which(exprsMat[,"CD279 PD1"] >= rl[["CD279 PD1"]][[fns]]),
                                   cd28=which(exprsMat[,"CD28"] >= rl[["CD28"]][[fns]]),
                                   cd127=which(exprsMat[,"CD127"] < rl[["CD127"]][[fns]]),
                                   cd25=which(exprsMat[,"CD25"] < rl[["CD25"]][[fns]]),
                                   ccr7=which(exprsMat[,"CD197 CCR7"] < rl[["CD197 CCR7"]][[fns]])
                                   ))
    
    parentLookup <- which(exprsMat[,"CD3"] >= rl[["CD3"]][[fns]])
    childCount <- length(childLookup)
    parentCount <- length(parentLookup)
    resDF <- data.frame(
        sampleName=fn,
        childCount=childCount,
        parentCount=parentCount,
        prop=(childCount/parentCount),
        stringsAsFactors=FALSE
    )
    if (firstSample) {
        allRDF <- resDF
        firstSample <- FALSE
    }
    else {
        allRDF <- rbind(allRDF,resDF)
    }
}
allRDF$`subjectID` <- str_extract(allRDF$sampleName,"Subject-\\d\\d")
modelDF <- inner_join(allRDF,meta_data,by=c("subjectID"))
#
#model virus positive subjects
#
resType <- rep(0,nrow(modelDF))
resType[which(modelDF$`OverallResponse` %in% c("CR","PR"))] <- 1
modelDF$resType <- as.factor(resType)
modelDF$Response <- as.factor(resType)
modelDF$pct <- (100*modelDF$prop)
table(modelDF$resType,modelDF$tumorViralStatus)
testDF <- modelDF[which(modelDF$tumorViralStatus=="Positive"),]
pv <- tcellSafeMod(testDF)
stat.test.3 <- data.frame(
    group1=c(1.75),
    group2=c(2.25),
    obsp=c(round(as.numeric(pv),3))
)
#
#plot cd8
#
pbmcCd8Freqs <- ggplot(modelDF)+
        geom_boxplot(
            outlier.color=NA,
            aes(
                x=tumorViralStatus,
                y=pct,
                fill=interaction(tumorViralStatus,Response)
            )
        )+
        geom_point(
            position=position_jitterdodge(),
            aes(
                x=tumorViralStatus,
                y=pct,
                fill=interaction(tumorViralStatus,Response),
                shape=interaction(tumorViralStatus,Response)
            ),
            size=global_plot_point_size ,
            stroke=0.5
        )+
        scale_fill_manual(values=c("#932667FF","#932667FF","#F6D645FF","#F6D645FF"),guide=FALSE)+
        scale_shape_manual(values=c(2,19,2,19))+
           theme_classic(base_size = global_plot_font_size)+
        theme(
          legend.position="none",
          plot.title=element_markdown(
            size=global_plot_font_size,
            hjust=0.5,
            colour = "#000000",
            lineheight=1.2
          )
       )+
    xlab("")+
    ylab("100 • ( Phenotype / Total CD3+ )")+
    ggtitle(paste0("CD4- **CD3+** **CD8+** CD45RA- **PD-1+**<br>**HLA-DR+** **CD28+** CD127- CD25- CCR7-<br>p-value in virus positive subjects: ",round(pv,3)))+
    guides(shape=guide_legend("Tumor viral status"))+ 
    scale_size_continuous(range=c(1,10))+
    scale_y_continuous(trans="sqrt",
                       labels=function(x){paste0(x,"%")})+
    scale_x_discrete(labels=c("Virus Negative","Virus Positive"))+
    stat_pvalue_manual(stat.test.3,y.position=0.375,tip.length=0.01,label="obsp",size=2)+
    coord_cartesian(ylim=c(0,0.25))
#
#
#
#
#compute targeted frequencies on CyTOF dataset
#
#
#
#
#
cytofDataPath <- file.path(normalizePath("../.."),"faustRuns","melanoma_cytof")
set.seed(123)
rl <- readRDS(file.path(cytofDataPath,"faustData","gateData","live_resList.rds"))
#
#gate out the targeted cd8 population by sample.
#
firstSample <- TRUE
for (fn in list.files(file.path(cytofDataPath,"faustData","sampleData"))) {
    # 
    #import the expression data for each sample
    #
    exprsMat <- readRDS(file.path(cytofDataPath,"faustData","sampleData",fn,"exprsMat.rds"))
    #  
    #find all cells that match the phenotype of interest
    #
    lookup <- Reduce(intersect,list(
                                   cd4=which(exprsMat[,"145Nd_CD4"] < rl[["145Nd_CD4"]][[fn]]),
                                   cd8=which(exprsMat[,"146Nd_CD8"] >= rl[["146Nd_CD8"]][[fn]]),
                                   cd3=which(exprsMat[,"154Sm_CD3"] >= rl[["154Sm_CD3"]][[fn]]),
                                   cd45ra=which(exprsMat[,"153Eu_CD45RA"] < rl[["153Eu_CD45RA"]][[fn]]),
                                   dr=which(exprsMat[,"157Gd_HLA-DR"] >= rl[["157Gd_HLA-DR"]][[fn]]),
                                   pd1=which(exprsMat[,"170Er_PD-1"] >= rl[["170Er_PD-1"]][[fn]]),
                                   cd28=which(exprsMat[,"155Gd_CD28"] >= rl[["155Gd_CD28"]][[fn]]),
                                   cd127=which(exprsMat[,"165Ho_CD127"] < rl[["165Ho_CD127"]][[fn]]),
                                   cd25=which(exprsMat[,"141Pr_CD25"] < rl[["141Pr_CD25"]][[fn]]),
                                   ccr7=which(exprsMat[,"169Tm_CCR7"] < rl[["169Tm_CCR7"]][[fn]])
                               ))
    childCount <- length(lookup)
    parentCount <- length(which(exprsMat[,"154Sm_CD3"] >= rl[["154Sm_CD3"]][[fn]]))
    resDF <- data.frame(
        sampleName=fn,
        childCount=childCount,
        parentCount=parentCount,
        prop=(childCount/parentCount),
        stringsAsFactors=FALSE
    )
    if (firstSample) {
        allRDF <- resDF
        firstSample <- FALSE
    }
    else {
        allRDF <- rbind(allRDF,resDF)
    }
}
#
#merge in meta data to identify therapy, responder type.
#
metaDataIn <- read_xlsx(file.path(normalizePath("../.."),
                                  "dataSets",
                                  "publication_gating_sets",
                                  "melanoma_cytof",
                                  "FlowRepository_FR-FCM-ZYKP_files",
                                  "attachments",
                                  "Sample_Mapping_File.xlsx"))
metaData <- as.data.frame(metaDataIn[-which(apply(metaDataIn,1,function(x){length(which(is.na(x)))})==7),])
colnames(metaData)[which(colnames(metaData)=="Sample File name")] <- "sampleName"
colnames(metaData)[which(colnames(metaData)=="Time Point")] <- "timePoint"
colnames(metaData)[which(colnames(metaData)=="Immunotherapy          Type")] <- "immunotherapyType"
countDF <- inner_join(allRDF,metaData,by=c("sampleName"))
#
#subset to pembrolizumab samples, fit model, and plot.
#
plotLookup <- which(countDF[,"immunotherapyType"]=="Pembrolizumab")
propDF <- countDF[plotLookup,]
propDF$pct <- (100*propDF$prop)
m <- glmer(cbind(childCount,(parentCount-childCount)) ~ Response + (1|sampleName),
           data=propDF,family="binomial",
           control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1e9),
                                check.conv.singular = .makeCC(action = "warning",  tol = 1e-4))
           )
pv <- coefficients(summary(m))[2,4]
stat.test.out <- data.frame(
    group1=c("Non-responder"),
    group2=c("Responder"),
    obsp=c(round(as.numeric(pv),3))
)
cytofCd8Freqs <- ggplot(propDF) +
    geom_boxplot(
        outlier.color=NA,
        aes(
            x=Response,
            y=pct,
            fill=Response
        )
    )+
    geom_point(
        position=position_jitterdodge(),
        aes(
            x=Response,
            y=pct,
            fill=Response
        ),
        size=global_plot_point_size ,
        stroke=0.5
    )+
    scale_fill_manual(
        values=c("#932667FF","#F6D645FF"),
        guide=FALSE
    )+
    theme_classic(base_size = global_plot_font_size)+
    theme(
        legend.position="none",
        plot.title=element_markdown(
            size = global_plot_font_size,
            hjust=0.5,
            colour = "#000000",
            lineheight=1.2
        )
    )+
    xlab("")+
    ggtitle(plotName)+
    scale_y_continuous(trans="sqrt")+
    scale_x_discrete(labels=c("Non-responder","Responder"))+
    coord_cartesian(ylim=c(0,0.525))+
    stat_pvalue_manual(stat.test.out,y.position=0.65,label="obsp",size=2)+
    ylab("")
#
#collect all components into a single figure
#
fig_grid <- plot_grid(
    (plotList1[[2]]+ylab("")+ggtitle("PD-1 dim")),(plotList1[[3]] + ylab("") + ggtitle("PD-1 bright")),
    (ihcDim + ylab("")),(ihcBright + ylab("")),
    (cloneDim+ylab("")),(cloneBright+ylab("")),
    (pbmcCd8Freqs+ggtitle("MCC validation<br>PD-1+")+ylab("")),(cytofCd8Freqs+ggtitle("Melanoma validaton<br>PD-1+")+ylab("")),
    nrow=4,
    ncol=2,
    labels=c(
        "A","",
        "B","",
        "C","",
        "D",""
    ),
    label_size=global_plot_label_size
)
#
#add the title
#
axis_info <- ggdraw() +
    draw_label(
        "Frequency of effector memory CD8+ T-cells expressing CD28, HLA-DR, and PD-1 relative to total CD3+",
        fontface = 'bold',
        y = 0.5,
        angle = 90,
        size=global_plot_label_size
    ) 
fig_grid_out <- plot_grid(
    axis_info,
    fig_grid,
    ncol = 2,
    rel_widths = c(0.1, 1)
)
#
#append legend and save
#
figure <- plot_grid(
    fig_grid_out,
    ggdraw(myLegend),
    nrow=2,
    ncol=1,
    rel_heights=c(1,0.05)
)
ggsave(
    filename=file.path(normalizePath("."),"figure_citn09.png"),
    plot=figure,
    width=5,
    height=7.5,
    units="in",
    dpi=300
)
