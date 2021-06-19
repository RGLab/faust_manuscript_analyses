library(flowWorkspace)
library(DT)
library(broom)
library(forcats)
library(scales)
library(data.table)
library(multcomp)
library(xtable)
library(ggplot2)
library(cowplot)
library(lme4)
library(tidyr)
library(readxl)
library(dplyr)
library(stringr)
library(openCyto)
library(viridis)
library(faust)
#
#include following libraries/files to bold titles
#
library(grid)
#library(devtools)
#devtools::install_github("clauswilke/gridtext")
library(gridtext)
library(rlang)
source("./helperFunctions.R")
set.seed(12345)
global_plot_font_size <- 14
global_plot_label_size <- 20
scaleFUN <- function(x) sprintf("%.2f", x)
my_loglog_trans= function() trans_new("my_loglog", function(x){log(log(x))}, function(x){exp(exp(x))})
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
#CITN-09 Myeloid
#
paperLog("Start analysis","CITN-09 Myeloid")
#
#load meta data
#
gs <- load_gs(file.path(normalizePath("../.."),
                        "dataSets",
                        "publication_gating_sets",
                        "CITN-09",
                        "myeloid_panel",
                        "citn09_myeloid_gs"))
meta_data_all <- as.data.frame(pData(gs))
meta_data_all$subjectID <- meta_data_all$publicationID
bioData <- meta_data_all[which(meta_data_all$VISIT=="C01"),]
#
#load faust count matrix. merge metadata and derive variables.
#
citn09CountDF <- as.data.frame(readRDS(file.path(normalizePath("../.."),"faustRuns","CITN-09","myeloid_panel","faustData","faustCountMatrix.rds")))
citn09Pops <- colnames(citn09CountDF)
citn09Pops <- setdiff(citn09Pops,c("0_0_0_0_0"))
paperLog(length(citn09Pops),"Number of populations in CITN-09 Myeloid.")
paperLog(nrow(citn09CountDF),"Number of samples processed in CITN-09 Myeloid.")
paperLog(quantile(1-(citn09CountDF[,"0_0_0_0_0"]/apply(citn09CountDF,1,sum))),"Annotated subsets.")
citn09CountDF$parentCount <- apply(citn09CountDF,1,sum)
citn09CountDF$sampleName <- rownames(citn09CountDF)
citn09CountDF$subjectID <- str_extract(rownames(citn09CountDF),"Subject-\\d\\d")
paperLog(table(citn09CountDF$subjectID),"Samples by subject in CITN-09")
citn09CountDF$visitSTR <- str_extract(rownames(citn09CountDF),"C\\d\\d|EOT")
paperLog(table(citn09CountDF$visitSTR),"Time course in CITN-09")
paperLog(rownames(citn09CountDF)[which(is.na(citn09CountDF$visitSTR))],"check if any visitSTR are undefined")
citn09CountDF$visitSTR[which(is.na(citn09CountDF$visitSTR))] <- c("EOT","C01")
paperLog(rownames(citn09CountDF)[which(citn09CountDF$visitSTR=="C01")],"updated eot,C1D01 visit strings")
citn09AnalysisDF <- left_join(citn09CountDF,bioData,by="subjectID")
citn09AnalysisDF[which(citn09AnalysisDF$OverallResponse == "NE"),"subjectID"]
#remove subject that does not have outcome classified
citn09AnalysisDF <- citn09AnalysisDF[-which(citn09AnalysisDF$OverallResponse == "NE"),]
#
#subset to C01 for modeling and derive CR/PR merged responder type.
#
citn09SubDF <- citn09AnalysisDF[which(citn09AnalysisDF[,c("visitSTR")]=="C01"),]
citn09ResType <- rep(0,nrow(citn09SubDF))
citn09ResType[which(citn09SubDF$OverallResponse %in% c("CR","PR"))] <- 1
citn09SubDF$modelRT <- as.factor(citn09ResType)
paperLog(table(citn09SubDF$subjectID,citn09SubDF$modelRT),"subject responder status in CITN-09")
paperLog(table(citn09SubDF$modelRT),"Response rate in CITN-09")
#
#specify the GLMM
#
citn09SafeMod <- function(dataSet) {
    out <- tryCatch(
    {
        m <- glmer(cbind(childCount,(parentCount-childCount)) ~ modelRT + (1|subjectID),
                   data=dataSet,
                   family="binomial",
                   control=glmerControl(
                       optimizer="bobyqa", optCtrl = list(maxfun = 1e9),
                       check.conv.singular = .makeCC(action = "warning",  tol = 1e-4
                                                     ))
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
#fit GLMM for all populations discovered by faust
#
citn09PvalVec <- c()
for (cellPop in citn09Pops) {
    mdf <- data.frame(modelRT=citn09SubDF[,"modelRT"],
                      subjectID=as.factor(citn09SubDF[,"subjectID"]),
                      parentCount=citn09SubDF[,"parentCount"],
                      childCount=citn09SubDF[,cellPop],
                      sampleName=citn09SubDF[,"sampleName"])
    pv <- citn09SafeMod(mdf)
    if (!is.na(pv)) {
        citn09PvalVec <- append(citn09PvalVec,pv[1])
        names(citn09PvalVec)[length(citn09PvalVec)] <- paste0(cellPop)
    }
}
#
#adjust for bonferonni and select significant populations
#
citn09AdjPvec <- p.adjust(citn09PvalVec,"bonferroni")
#
#generate tables for methods section 
#
methodReportNames <- sort(unique(names(which(citn09AdjPvec < 0.1))))
methodReportDF <- data.frame(
    phenotypes=names(sort(citn09AdjPvec[methodReportNames])),
    bonferroni.p=as.numeric(sort(citn09AdjPvec[methodReportNames]))
)
xtable(methodReportDF,digits=4)
#
#continue with analysis for figure
#
citn09ORStr <-  citn09SubDF[,"OverallResponse"]
citn09ORStr[which(citn09ORStr == "PD")] <- "0_PD/SD"
citn09ORStr[which(citn09ORStr == "SD")] <- "0_PD/SD"
citn09ORStr[which(citn09ORStr == "PR")] <- "1_PR/CR"
citn09ORStr[which(citn09ORStr == "CR")] <- "1_PR/CR"
citn09PlotPop <- "CD33BrightCD16-CD15-CD14+CD3-HLA DRBrightCD20-CD19-CD11B+CD56-CD11C+"
paperLog(citn09AdjPvec[citn09PlotPop],"bonferroni p value of myeloid correlate in CITN-09")
citn09PropDF <- data.frame(
    subjectID=citn09SubDF[,"subjectID"],
    overallResponse=as.factor(citn09ORStr),
    Count= citn09SubDF[,citn09PlotPop],
    Total=  citn09SubDF[,"parentCount"],
    timePoint =  citn09SubDF[,"visitSTR"]
)
citn09PropDF$prop <- citn09PropDF$Count/citn09PropDF$Total
citn09PropDF$pct <- (100*(citn09PropDF$Count/citn09PropDF$Total))
citn09PlotName <- citn09PlotPop
citn09PlotName <- gsub("HLA DR","HLADR",citn09PlotName)
citn09PlotName <- gsub("Medium","Dim ",citn09PlotName)
citn09PlotName <- gsub("HLADRBright","HLADR Bright ",citn09PlotName)
citn09PlotName <- gsub("CD33Bright","CD33 Bright ",citn09PlotName)
citn09PlotName <- gsub("-","- ",citn09PlotName)
citn09PlotName <- gsub("\\+","\\+ ",citn09PlotName)
citn09PlotName <- gsub("CD33 Bright","**CD33 Bright** ",citn09PlotName)
citn09PlotName <- gsub("CD14\\+","**CD14\\+** ",citn09PlotName)
citn09PlotName <- gsub("HLADR Bright","**HLA-DR Bright** ",citn09PlotName)
citn09PlotName <- gsub("CD3-","CD3-<br>",citn09PlotName)
citn09PlotName <- gsub("CD11B\\+","**CD11B\\+** ",citn09PlotName)
citn09PlotName <- gsub("CD11C\\+","**CD11C\\+** ",citn09PlotName)
boxplotCITN09 <- ggplot(citn09PropDF,aes(x=overallResponse,y=pct,fill=overallResponse))+
    geom_boxplot(outlier.color=NA)+
    geom_jitter(aes())+
    scale_fill_viridis(begin=0.4,end=0.9,option="B",discrete=TRUE,guide=FALSE)+
    scale_x_discrete(labels=c("Non-Responder", "Responder"))+
    scale_size(range = c(1, 4))+
    scale_shape_manual(values=c(15,19),guide=FALSE)+
    theme_classic(base_size = 10)+
    theme(
        legend.position="none",
        plot.title=element_markdown(
            size=10,
            hjust=0.5,
            colour = "#000000",
            lineheight=1.2
        )
    )+
    ggtitle(paste0("MCC anti-PD-1 trial<br>",citn09PlotName,"<br>Bonferroni p-value: ",as.numeric(round(citn09AdjPvec[citn09PlotPop],3))))+
    xlab("")+
    ylab("% of CD45+ Cells")+
    guides(size=guide_legend("Number of events"))+
    scale_y_continuous(trans="sqrt",
                       labels=function(x) paste0(x,"%"))
#
#target CD14+ cells for PFDA
#
citn09CD14Pops <- names(citn09PvalVec)
citn09CD14Pops <- citn09CD14Pops[grepl("HLA DRBright",citn09CD14Pops)]
citn09CD14Pops <- citn09CD14Pops[grepl("CD3-",citn09CD14Pops)]
citn09CD14Pops <- citn09CD14Pops[grepl("CD56-",citn09CD14Pops)]
citn09CD14Pops <- citn09CD14Pops[grepl("CD19-",citn09CD14Pops)]
citn09CD14Pops <- citn09CD14Pops[grepl("CD14\\+",citn09CD14Pops)]
citn09CD14Pops <- citn09CD14Pops[grepl("CD16-",citn09CD14Pops)]

#
#fit the PFDA model
#
citn09GlmerDFPrep <- citn09SubDF[,c(citn09CD14Pops,"parentCount","modelRT","subjectID")]
citn09GlmerDF <- gather(citn09GlmerDFPrep,key=cellPop,value=count,-c("parentCount","modelRT","subjectID"))
citn09GlmerDF$popName <- as.factor(paste0("Pop_",as.numeric(as.factor(citn09GlmerDF$cellPop))))
citn09GlmerDF$obsFactor <- as.factor(paste0("Obs_",seq(nrow(citn09GlmerDF))))
citn09GlmerFit <- glmer(
    formula=cbind(count,(parentCount-count))  ~ popName*modelRT + (1|obsFactor),
    data=citn09GlmerDF,
    family="binomial",
    nAGQ=0,
    control=glmerControl(
        optimizer="bobyqa",
        optCtrl = list(maxfun = 1e9),
        boundary.tol = 1e-9,
        tolPwrss=1e-9,
        check.conv.singular = .makeCC(action = "warning",  tol = 1e-12)
    )
)
citn09CoefNames = names(fixef(citn09GlmerFit))
citn09GLHTMat <- matrix(0, nrow = 1, ncol = length(citn09CoefNames))
citn09GLHTMat[,which(citn09CoefNames == "modelRT1")] <- 1
citn09GLHTMat[1,which(citn09CoefNames %like% ":modelRT1")] <- 1/(length(which(citn09CoefNames%like%":modelRT1"))+1)
rownames(citn09GLHTMat) <- c("CD14Monocytes")
citn09Multi <- glht(citn09GlmerFit,citn09GLHTMat, alternative = "greater")
citn09ObsPval <- round(summary(citn09Multi)$test$pvalues[1],9)
citn09CT95 <- confint(citn09Multi,level=0.95)
citn09MultiResults <- as.data.frame(citn09CT95$confint)
citn09MultiResults$upr <- 3
citn09MultiResults$Compartment <- rownames(citn09MultiResults)
citn09MultiResults$ModelType <- "PFDA"
citn09MultiResults$Component <- "PFDA"
#
#get confidence intervals for univariate fits
#
citn09CIListCD14 <- list()
for (cd14pt in citn09CD14Pops) {
    mdf <- data.frame(
        count=citn09SubDF[,cd14pt],
        parentCount=citn09SubDF[,"parentCount"],
        modelRT=citn09SubDF[,"modelRT"],
        subjectID=citn09SubDF[,"subjectID"]
    )
    m <- glmer(
        formula=cbind(count,(parentCount-count))  ~ modelRT + (1|subjectID),
        data = mdf,
        family="binomial",
        control=glmerControl(
            optimizer="bobyqa",
            optCtrl = list(maxfun = 1e9),
            boundary.tol = 1e-9,
            tolPwrss=1e-9,
            check.conv.singular = .makeCC(action = "warning",  tol = 1e-12)
        )
    )
    #90 % CI since testing one-sided hypothesis
    ci1 <- confint(m,level=0.90,devtol=1e-5)
    outRes <- cbind(fixef(m)[1:2],ci1[2:3,])
    colnames(outRes) <- c("Estimate","lwr","upr")
    outRes[,"upr"] <- 3
    citn09CIListCD14 <- append(citn09CIListCD14,list(outRes))
    names(citn09CIListCD14)[length(citn09CIListCD14)] <- cd14pt
}
citn09UniPlotDFPrep <- as.data.frame(Reduce(rbind,
                                       lapply(citn09CIListCD14,
                                              function(x){as.numeric(x[which(rownames(x)=="modelRT1"),,drop=TRUE])})))
colnames(citn09UniPlotDFPrep) <- c("Estimate","lwr","upr")
rownames(citn09UniPlotDFPrep) <- seq(nrow(citn09UniPlotDFPrep))
citn09UniPlotDFPrep$Compartment <- as.character(paste0("Population ",seq(nrow(citn09UniPlotDFPrep))))
citn09UniPlotDFPrep$Method <- "FAUST\npopulations\nused by\nMultivariate"
citn09UniPlotDFPrep[which(citn09UniPlotDFPrep$Compartment=="14+"),"Compartment"] <- "Univariate"
citn09UniResults <- citn09UniPlotDFPrep[order(citn09UniPlotDFPrep[,"lwr"]),]
citn09UniResults$ModelType <- "CD14+ CD16-\nHLA-DR bright\nFAUST Phenotypes\nUnivariate CIs"
citn09UniResults$Component <- paste0("Univariate_",str_pad(seq(nrow(citn09UniResults)),pad="0",width=2))
#
#aggregate over cd14 phenotypes and test
#
citn09_aggregate_meta_data <- data.frame(
    modelRT=citn09SubDF[,"modelRT"],
    subjectID=as.factor(citn09SubDF[,"subjectID"]),
    sampleName=citn09SubDF[,"sampleName"]
)
citn09_aggregate_target_df <- faust::getCountsForTargetMarkers(
    projectPath=file.path(normalizePath("../.."),"faustRuns","CITN-09","myeloid_panel"),
    referencePhenotype="CD33BrightCD16-CD15-CD14+CD3-HLA DRBrightCD20-CD19-CD11B+CD56-CD11C+",
    targetMarkers=c("CD16","CD14","CD3","HLA DR","CD19","CD56"),
    metaDataDF=citn09_aggregate_meta_data
    )

citn09AggregateModel <- glmer(
    cbind(childCount,(parentCount-childCount)) ~ modelRT + (1|subjectID),
    data=citn09_aggregate_target_df,
    family="binomial",
    control=glmerControl(
        optimizer="bobyqa", optCtrl = list(maxfun = 1e9),
        check.conv.singular = .makeCC(action = "warning",  tol = 1e-4)
    )
)
citn09AggregateCI <- confint(citn09AggregateModel,level=0.90,devtol=1e-5)
citn09AggregateOutRes <- cbind(fixef(citn09AggregateModel),citn09AggregateCI[2:3,])
colnames(citn09AggregateOutRes) <- c("Estimate","lwr","upr")
citn09AggregateOutRes[,"upr"] <- 3
citn09AggregateResults <- as.data.frame(citn09AggregateOutRes[which(rownames(citn09AggregateOutRes)=="modelRT1"),,drop=FALSE])
colnames(citn09AggregateResults) <- c("Estimate","lwr","upr")
rownames(citn09AggregateResults) <- c("CD14Monocytes")
citn09AggregateResults$Compartment <- c("CD14Monocytes")
citn09AggregateResults$compartmentSize <- length(citn09CD14Pops)
citn09AggregateResults$ModelType <- "Targeted"
citn09AggregateResults$Component <- "Targeted"

#
#combine and plot
#
citn09PFDAplotDF <- rbind(
    citn09UniResults[,c("Estimate","lwr","upr","ModelType","Component")],
    citn09AggregateResults[,c("Estimate","lwr","upr","ModelType","Component")],
    citn09MultiResults[,c("Estimate","lwr","upr","ModelType","Component")]
)
citn09PFDAplotDF[,c("Estimate","lwr","upr")] <- exp(citn09PFDAplotDF[,c("Estimate","lwr","upr")])
citn09PFDAplotDF$ModelType <- factor(citn09PFDAplotDF$ModelType,levels=c("PFDA",
                                                                         "Targeted",
                                                                         "CD14+ CD16-\nHLA-DR bright\nFAUST Phenotypes\nUnivariate CIs"))
pfdaPlotCITN09 <- ggplot(citn09PFDAplotDF, aes(x=ModelType, y=Estimate))+
    geom_errorbar(
        aes(ymin=lwr, ymax=40, colour=`Component`),
        size=0.05,
        width=0.5,
        position=position_dodge2(width=0.75,reverse=TRUE)
    )+
    geom_point(
        aes(group=`Component`),
        size=0.6,
        shape=15,
        position=position_dodge2(width=0.5,reverse=TRUE)
    )+
    geom_hline(
        yintercept=1,
        size=0.05,
        linetype="solid",
        color="red"
    )+
    coord_flip(xlim=c(1,3),ylim=c(exp(-5),exp(3)))+
    ylab("Point estimate of odds-ratio (square mark)\nwith one-sided 95% confidence interval")+
    scale_color_manual(values=rep("#000000FF",nrow(citn09PFDAplotDF)))+
    theme_classic(base_size=6)+
    xlab("")+
    theme(
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_text(size=9),
        plot.title = element_text(hjust = 0.5)
    )+
    ggtitle(paste0("MCC anti-PD-1 trial, Myeloid Panel"))+
    guides(colour=guide_legend(keywidth = 4, keyheight = 1))+
    theme(legend.position="none")+
    scale_y_continuous(trans="log",labels=scaleFUN,breaks=c(exp(-3),exp(-2),exp(-1),1,exp(1),exp(2)))        
paperLog("End analysis","CITN-09 Myeloid")
#
#
#CITN-07 Myeloid
#
#
paperLog("Start analysis","CITN-07 Myeloid")
#
# Read in FAUST Results
#
citn07CountDF <- as.data.frame(readRDS(file.path(normalizePath("../.."),"faustRuns","CITN-07","faustData","faustCountMatrix.rds")))
citn07CellPops <- colnames(citn07CountDF)
citn07CellPops <- setdiff(citn07CellPops,c("0_0_0_0_0"))
paperLog(length(citn07CellPops),"Number of populations in CITN-07.")
paperLog(nrow(citn07CountDF),"Number of samples processed in CITN-07 Myeloid.")
paperLog(quantile(1-(citn07CountDF[,"0_0_0_0_0"]/apply(citn07CountDF,1,sum))),"Annotated subsets.")
citn07CountDF$parentCount <- apply(citn07CountDF,1,sum)
citn07CountDF$sampleName <- rownames(citn07CountDF)
citn07CountDF$`Subject ID` <- str_extract(rownames(citn07CountDF),"Subject-\\d\\d")
#
#derive visit
#
cytometry_visit <- gsub(".fcs_[[:digit:]]+","",citn07CountDF$sampleName)
cytometry_visit <- gsub("PT_ID_Subject-[[:digit:]][[:digit:]] ","",cytometry_visit)
cytometry_visit <- gsub("PT ID_Subject-[[:digit:]][[:digit:]] ","",cytometry_visit)
cytometry_visit <- gsub("PT_ID_Subject-[[:digit:]][[:digit:]]_","",cytometry_visit)
cytometry_visit <- gsub("PT ID_Subject-[[:digit:]][[:digit:]]_","",cytometry_visit)
cytometry_visit <- gsub("PT_Subject-[[:digit:]][[:digit:]] ","",cytometry_visit)
cytometry_visit <- gsub("C1D-7","C1D-07",cytometry_visit)
cytometry_visit <- gsub("C1D-7 P","C1D-07",cytometry_visit)
cytometry_visit <- gsub("C1D-07 P","C1D-07",cytometry_visit)
cytometry_visit <- gsub("C3D01 P","C3D01",cytometry_visit)
cytometry_visit <- gsub("C3D011","C3D01",cytometry_visit)
cytometry_visit <- gsub("C4D01_[[:digit:]]+","C4D01",cytometry_visit)
cytometry_visit <- gsub("FUW04_[[:digit:]]+","FUW04",cytometry_visit)
cytometry_visit <- gsub("FUW4","FUW04",cytometry_visit)
cytometry_visit <- gsub("FUw04","FUW04",cytometry_visit)
cytometry_visit <- gsub("FUW12_[[:digit:]]+","FUW12",cytometry_visit)
citn07CountDF$`Cytometry_Visit` <- as.factor(cytometry_visit)
#
#Read in demographic data
#
gs <- load_gs(file.path(normalizePath("../.."),
                        "dataSets",
                        "publication_gating_sets",
                        "CITN-07",
                        "citn07_longitudinal_gs"))
citn07DemoDF <- as.data.frame(pData(gs))
citn07DemoDF$`sampleName` <- rownames(citn07DemoDF)
#
# Tidy the count matrix
# Merge the counts with the demographics
#
citn07AnalysisDF <- inner_join(citn07CountDF, citn07DemoDF,by=c("sampleName"))
citn07AnalysisDF[which(citn07AnalysisDF$`Disease Status`=="NA"),"Disease Status"] <- "No RECURRENCE"
citn07AnalysisDF$`Disease Status` <- as.factor(citn07AnalysisDF$`Disease Status`)
citn07AnalysisDF$Cohort <- as.factor(citn07AnalysisDF$Cohort)
citn07AnalysisDF$Cohort <- relevel(factor(citn07AnalysisDF$Cohort),"CH2")
citn07AnalysisDF$Recurrence <- relevel(factor(citn07AnalysisDF$`Disease Status`),"RECURRENCE")
citn07AnalysisDF$`NYESO-1 IHC` <- as.factor(citn07AnalysisDF$`NYESO-1 IHC`)
citn07ModelDF <- citn07AnalysisDF  %>% dplyr::filter((Cytometry_Visit == "C1D-07" & Cohort == "CH1") | (Cytometry_Visit == "C1D01" & Cohort == "CH2"))
#
#specify the GLMM
#
citn07SafeMod <- function(dataSet) {
    out <- tryCatch(
    {
        m <- glmer(
            cbind(Count, Total-Count) ~ Recurrence + Cohort + NYESO_1_IHC + (1|Subject_ID),
            data = dataSet,
            family="binomial",
            control=glmerControl(
                optimizer="bobyqa",
                optCtrl = list(maxfun = 1e9),
                boundary.tol = 1e-9,
                tolPwrss=1e-9,
                check.conv.singular = .makeCC(action = "warning",  tol = 1e-12)
            )
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
#fit GLMM for all populations discovered by faust
#
citn07PvalVec <- c()
for (cellPop in citn07CellPops) {
    mdf <- data.frame(
        Count=citn07ModelDF[,cellPop],
        Total=citn07ModelDF[,"parentCount"],
        Recurrence=citn07ModelDF[,"Recurrence"],
        Cohort=citn07ModelDF[,"Cohort"],
        NYESO_1_IHC=citn07ModelDF[,"NYESO-1 IHC"],
        Subject_ID=citn07ModelDF[,"Subject ID"]
    )
    pv <- citn07SafeMod(mdf)
    if (!is.na(pv)) {
        citn07PvalVec <- append(citn07PvalVec,pv[1])
        names(citn07PvalVec)[length(citn07PvalVec)] <- paste0(cellPop)
    }
}
citn07AdjPvec <- p.adjust(citn07PvalVec,"bonferroni")

#
#generate tables for supplement
#
citn07MethodReportNames <- sort(unique(names(which(citn07AdjPvec < 0.1))))
citn07MethodReportDF <- data.frame(
    phenotypes=names(sort(citn07AdjPvec[citn07MethodReportNames])),
    bonferroni.p=as.numeric(sort(citn07AdjPvec[citn07MethodReportNames]))
)
xtable(citn07MethodReportDF,digits=4)

#
#continue with analysis
#
citn07Plotpop <- "CD8-CD3-HLA DRBrightCD19-CD14+CD11C+CD4-CD123-CD16-CD56-"
citn07BonfPforPlot <- round(citn07AdjPvec[citn07Plotpop],3)
citn07PlotTitle <- "FLT3-ligand + Therapeutic Vx trial<br>CD8- CD3- **HLA-DR Bright** CD19- **CD14+**<br>**CD11C+** CD4- CD123- CD16- CD56-<br>Bonferroni p-value: "
citn07PlotTitle <- paste0(citn07PlotTitle,citn07BonfPforPlot)
citn07InitPlotDF <- data.frame(
    Count=citn07ModelDF[,citn07Plotpop],
    Total=citn07ModelDF[,"parentCount"],
    Recurrence=citn07ModelDF[,"Recurrence"],
    Cohort=citn07ModelDF[,"Cohort"],
    NYESO_1_IHC=citn07ModelDF[,"NYESO-1 IHC"],
    Subject_ID=citn07ModelDF[,"Subject ID"],
    DiseaseStatus=citn07ModelDF[,"Disease Status"]
)
boxplotCITN07 <- ggplot(citn07InitPlotDF,aes(x=DiseaseStatus,y=(100*(Count/Total))))+
    ggtitle(citn07PlotTitle)+
    scale_fill_viridis(
        begin = 0.9,
        end = 0.4,
        option = "B",
        discrete = TRUE,
        guide = FALSE
    ) +
    geom_boxplot(
        aes(fill = DiseaseStatus),
        outlier.color = NA
    )+
    geom_jitter(aes())+
    theme_classic(base_size = 10) + 
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 1),
        plot.title=element_markdown(
            size=10,
            hjust=0.5,
            colour = "#000000",
            lineheight=1.2
        )
    )+
    ylab("% of CD45+ cells")+
    scale_y_continuous(trans="sqrt",labels=function(x) paste0(x,"%"))+
    scale_x_discrete("", labels = c("Recurrence","No Recurrence"), limits = c("RECURRENCE","No RECURRENCE"))
#
#target CD14+ cells for PFDA
#
citn07CD14Pops <- names(citn07PvalVec)
citn07CD14Pops <- citn07CD14Pops[grepl("HLA DRBright",citn07CD14Pops)]
citn07CD14Pops <- citn07CD14Pops[grepl("CD3-",citn07CD14Pops)]
citn07CD14Pops <- citn07CD14Pops[grepl("CD56-",citn07CD14Pops)]
citn07CD14Pops <- citn07CD14Pops[grepl("CD19-",citn07CD14Pops)]
citn07CD14Pops <- citn07CD14Pops[grepl("CD14\\+",citn07CD14Pops)]
citn07CD14Pops <- citn07CD14Pops[grepl("CD16-",citn07CD14Pops)]
#
#fit the PFDA model
#
citn07GlmerDFPrep <- citn07ModelDF[,c(citn07CD14Pops,"parentCount","Recurrence","Cohort","NYESO-1 IHC","Subject ID")]
citn07GlmerDF <- gather(citn07GlmerDFPrep,key=cellPop,value=count,-c("parentCount","Recurrence","Cohort","NYESO-1 IHC","Subject ID"))
citn07GlmerDF$popName <- as.factor(paste0("Pop_",as.numeric(as.factor(citn07GlmerDF$cellPop))))
citn07GlmerDF$obsFactor <- as.factor(paste0("Obs_",seq(nrow(citn07GlmerDF))))
citn07GlmerFit <- glmer(
    formula= cbind(count,(parentCount-count))  ~ popName*Recurrence + Cohort + `NYESO-1 IHC` + (1|obsFactor),
    data=citn07GlmerDF,
    family="binomial",
    control=glmerControl(
        optimizer="bobyqa",
        optCtrl = list(maxfun = 1e9),
        boundary.tol = 1e-9,
        tolPwrss=1e-9,
        check.conv.singular = .makeCC(action = "warning",  tol = 1e-12)
    )
)
citn07CoefNames = names(fixef(citn07GlmerFit))
citn07GLHTMat <- matrix(0, nrow = 1, ncol = length(citn07CoefNames))
citn07GLHTMat[,which(citn07CoefNames == "RecurrenceNo RECURRENCE")] <- 1
citn07GLHTMat[1,which(citn07CoefNames %like% ":RecurrenceNo RECURRENCE")] <- 1/(length(which(citn07CoefNames%like%":RecurrenceNo RECURRENCE"))+1)
rownames(citn07GLHTMat) <- c("CD14Monocytes")
citn07Multi <- glht(citn07GlmerFit,citn07GLHTMat, alternative = "greater")
citn07ObsPval <- round(summary(citn07Multi)$test$pvalues[1],9)
citn07CT95 <- confint(citn07Multi,level=0.95)
citn07MultiResults <- as.data.frame(citn07CT95$confint)
citn07MultiResults$upr <- 3
citn07MultiResults$Compartment <- rownames(citn07MultiResults)
citn07MultiResults$ModelType <- "PFDA"
citn07MultiResults$Component <- "PFDA"
#
#get confidence intervals for univariate fits
#
citn07CIListCD14 <- list()
for (cd14pt in citn07CD14Pops) {
    mdf <- data.frame(
        Count=citn07ModelDF[,cd14pt],
        Total=citn07ModelDF[,"parentCount"],
        Recurrence=citn07ModelDF[,"Recurrence"],
        Cohort=citn07ModelDF[,"Cohort"],
        NYESO_1_IHC=citn07ModelDF[,"NYESO-1 IHC"],
        Subject_ID=citn07ModelDF[,"Subject ID"]
    )
    m <- glmer(
        cbind(Count, Total-Count) ~ Recurrence + Cohort + NYESO_1_IHC + (1 |Subject_ID),
        data = mdf,
        family="binomial",
        control=glmerControl(
            optimizer="bobyqa",
            optCtrl = list(maxfun = 1e9),
            boundary.tol = 1e-9,
            tolPwrss=1e-9,
            check.conv.singular = .makeCC(action = "warning",  tol = 1e-12)
        )
    )
    #90 % CI since testing one-sided hypothesis
    ci1 <- confint(m,level=0.90,devtol=1e-5)
    outRes <- cbind(fixef(m)[1:2],ci1[2:3,])
    colnames(outRes) <- c("Estimate","lwr","upr")
    outRes[,"upr"] <- 3
    citn07CIListCD14 <- append(citn07CIListCD14,list(outRes))
    names(citn07CIListCD14)[length(citn07CIListCD14)] <- cd14pt
}
citn07UniPlotDFPrep <- as.data.frame(Reduce(rbind,
                                       lapply(citn07CIListCD14,
                                              function(x){as.numeric(x[which(rownames(x)=="RecurrenceNo RECURRENCE"),,drop=TRUE])})))
colnames(citn07UniPlotDFPrep) <- c("Estimate","lwr","upr")
rownames(citn07UniPlotDFPrep) <- seq(nrow(citn07UniPlotDFPrep))
citn07UniPlotDFPrep$Compartment <- as.character(paste0("Population ",seq(nrow(citn07UniPlotDFPrep))))
citn07UniPlotDFPrep$Method <- "FAUST\npopulations\nused by\nMultivariate"
citn07UniPlotDFPrep[which(citn07UniPlotDFPrep$Compartment=="14+"),"Compartment"] <- "Aggregate"
citn07UniResults <- citn07UniPlotDFPrep[order(citn07UniPlotDFPrep[,"lwr"]),]
citn07UniResults$ModelType <- "CD14+ CD16-\nHLA-DR bright\nFAUST Phenotypes\nUnivariate CIs"
citn07UniResults$Component <- paste0("Univariate_",str_pad(seq(nrow(citn07UniResults)),pad="0",width=2))
#
#aggregate over cd14 phenotypes and test
#
citn07_aggregate_meta_data <- data.frame(
    sampleName=citn07ModelDF[,"sampleName"],
    Recurrence=citn07ModelDF[,"Recurrence"],
    Cohort=citn07ModelDF[,"Cohort"],
    NYESO_1_IHC=citn07ModelDF[,"NYESO-1 IHC"],
    Subject_ID=citn07ModelDF[,"Subject ID"]
)
citn07_aggregate_target_df <- faust::getCountsForTargetMarkers(
    projectPath=file.path(normalizePath("../.."),"faustRuns","CITN-07"),
    referencePhenotype="CD8-CD3-HLA DRBrightCD19-CD14+CD11C+CD4-CD123-CD16-CD56-",
    targetMarkers=c("HLA DR","CD14","CD16","CD3","CD19","CD56"),
    metaDataDF=citn07_aggregate_meta_data
    )
citn07AggregateModel <- glmer(
    cbind(childCount, parentCount-childCount) ~ Recurrence + Cohort + NYESO_1_IHC + (1 |Subject_ID),
    data = citn07_aggregate_target_df,
    family="binomial",
    control=glmerControl(
        optimizer="bobyqa",
        optCtrl = list(maxfun = 1e9),
        boundary.tol = 1e-9,
        tolPwrss=1e-9,
        check.conv.singular = .makeCC(action = "warning",  tol = 1e-12)
    )
)
citn07AggregateCIs <- confint(citn07AggregateModel,level=0.90,devtol=1e-5)
citn07AggregateOutRes <- cbind(fixef(citn07AggregateModel)[1:2],citn07AggregateCIs[2:3,])
colnames(citn07AggregateOutRes) <- c("Estimate","lwr","upr")
citn07AggregateOutRes[,"upr"] <- 3
citn07AggregateResults <- as.data.frame(citn07AggregateOutRes[which(rownames(citn07AggregateOutRes)=="RecurrenceNo RECURRENCE"),,drop=FALSE])
colnames(citn07AggregateResults) <- c("Estimate","lwr","upr")
rownames(citn07AggregateResults) <- c("CD14Monocytes")
citn07AggregateResults$Compartment <- c("CD14Monocytes")
citn07AggregateResults$compartmentSize <- length(citn07CD14Pops)
citn07AggregateResults$ModelType <- "Targeted"
citn07AggregateResults$Component <- "Targeted"
#
#combine and plot
#
citn07PFDAplotDF <- rbind(
    citn07UniResults[,c("Estimate","lwr","upr","ModelType","Component")],
    citn07MultiResults[,c("Estimate","lwr","upr","ModelType","Component")],
    citn07AggregateResults[,c("Estimate","lwr","upr","ModelType","Component")]
)
citn07PFDAplotDF[,c("Estimate","lwr","upr")] <- exp(citn07PFDAplotDF[,c("Estimate","lwr","upr")])
citn07PFDAplotDF$ModelType <- factor(citn07PFDAplotDF$ModelType,levels=c("PFDA",
                                                                         "Targeted",
                                                                         "CD14+ CD16-\nHLA-DR bright\nFAUST Phenotypes\nUnivariate CIs"))
pfdaPlotCITN07 <- ggplot(citn07PFDAplotDF, aes(x=ModelType, y=Estimate))+
    geom_errorbar(
        aes(ymin=lwr, ymax=40, colour=`Component`),
        size=0.05,
        width=0.5,
        position=position_dodge2(width=0.75,reverse=TRUE)
    )+
    geom_point(
        aes(group=`Component`),
        size=0.6,
        shape=15,
        position=position_dodge2(width=0.5,reverse=TRUE)
    )+
    geom_hline(
        yintercept=1,
        size=0.05,
        linetype="solid",
        color="red"
    )+
    coord_flip(xlim=c(1,3),ylim=c(exp(-5),exp(3)))+
    ylab("Point estimate of odds-ratio (square mark)\nwith one-sided 95% confidence interval")+
    scale_color_manual(values=rep("#000000FF",nrow(citn07PFDAplotDF)))+
    theme_classic(base_size=6)+
    xlab("")+
    theme(
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_text(size=9),
        plot.title = element_text(hjust = 0.5)
    )+
    ggtitle("FLT3-ligand + Therapeutic Vx trial")+
    guides(colour=guide_legend(keywidth = 4, keyheight = 1))+
    theme(legend.position="none")+
    scale_y_continuous(trans="log",labels=scaleFUN,breaks=c(exp(-3),exp(-2),exp(-1),1,exp(1),exp(2)))        
#
#
#
#Krieg FACS dataset
#
#
#
paperLog("Start analysis","Krieg FACS")
kfacsDataPath <- file.path(normalizePath("../.."),"faustRuns","melanoma_facs")
#
#Bring in faust count matrix and derive meta data for modeling.
#
kfacsCountDF <- as.data.frame(readRDS(file.path(kfacsDataPath,"faustData","faustCountMatrix.rds")))
paperLog(nrow(kfacsCountDF),"Number of samples processed in Krieg FACS.")
kfacsCellPops <- colnames(kfacsCountDF)
kfacsCellPops <- setdiff(kfacsCellPops,"0_0_0_0_0")
paperLog(length(kfacsCellPops),"number of cellular populations in Krieg FACS")
paperLog(median(1-kfacsCountDF[,"0_0_0_0_0"]/apply(kfacsCountDF,1,sum)),"Median annotation in Krieg FACS")
kfacsCountDF$parentCount <- apply(kfacsCountDF,1,sum)
kfacsResp <- rep(0,nrow(kfacsCountDF))
kfacsResp[grepl("_R_",rownames(kfacsCountDF))] <- 1
kfacsSubjectType <- rep("1_NonResponder",nrow(kfacsCountDF))
kfacsSubjectType[grepl("_R_",rownames(kfacsCountDF))] <- "2_Responder"
kfacsCountDF$subjectType <- kfacsSubjectType
kfacsCountDF$sampleID <- rownames(kfacsCountDF)
kfacsCountDF$modelRT <- as.factor(kfacsResp)
#
#Specify model, and look for correlates.
#
kfacsSafeMod <- function(dataSet) {
    out <- tryCatch(
    {
        m <- glmer(cbind(Count,(Total-Count)) ~ modelRT + (1|subjectID),
                   data=dataSet,family="binomial",
                   control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1e9),
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

kfacsPvalVec <- c()
for (cellPop in kfacsCellPops) {
    mdf <- data.frame(
        modelRT=kfacsCountDF[,"modelRT"],
        subjectID=as.factor(kfacsCountDF[,"sampleID"]),
        Total=kfacsCountDF[,"parentCount"],
        Count=kfacsCountDF[,cellPop]
    )
    pv <- kfacsSafeMod(mdf)
    if (!is.na(pv)) {
        kfacsPvalVec <- append(kfacsPvalVec,pv[1])
        names(kfacsPvalVec)[length(kfacsPvalVec)] <- paste0(cellPop)
    }
}
kfacsAdjPvec <- p.adjust(kfacsPvalVec,"bonferroni")
#
#generate tables for supplement
#
kfacsMethodReportNames <- sort(unique(names(which(kfacsAdjPvec < 0.1))))
kfacsMethodReportDF <- data.frame(
    phenotypes=names(sort(kfacsAdjPvec[kfacsMethodReportNames])),
    bonferroni.p=as.numeric(sort(kfacsAdjPvec[kfacsMethodReportNames]))
)
xtable(kfacsMethodReportDF,digits=4)
#
#continue with analysis
#
kfacsPlotpop <- "CD3-CD4+HLA-DR+CD14+CD19-CD11b+CD16-CD56-CD45RO+"
kfacsPlotTitle <- "Melanoma anti-PD-1 trial, FACS<br>CD3- **CD4+** **HLA-DR+** **CD14+**<br>CD19- **CD11B+** CD16- CD56- **CD45RO+**<br>Bonferroni p-value: "
kfacsPlotTitle <- paste0(kfacsPlotTitle,as.numeric(kfacsAdjPvec[kfacsPlotpop]))
kfacsPropDF <- data.frame(
    subjectType=as.factor(kfacsCountDF[,"subjectType"]),
    Count=kfacsCountDF[,kfacsPlotpop],
    Total=kfacsCountDF[,"parentCount"],
    stringsAsFactors=FALSE
)
kfacsPropDF$prop <- kfacsPropDF$Count/kfacsPropDF$Total
kfacsPropDF$pct <- (100*(kfacsPropDF$Count/kfacsPropDF$Total))
    
boxplotKFACS <- ggplot(kfacsPropDF,aes(x=subjectType,y=pct,fill=subjectType))+
    geom_boxplot(outlier.color=NA)+
    geom_jitter() +
    scale_fill_viridis(begin=0.4,end=0.9,option="B",discrete=TRUE,guide=FALSE)+
    scale_x_discrete(labels=c("Non-Responder", "Responder"))+
    scale_size(range = c(1, 4))+
    theme_classic(base_size = 8)+
    theme(
        legend.position="none",
        plot.title=element_markdown(
            size=10,
            hjust=0.5,
            colour = "#000000",
            lineheight=1.2
        )
    )+
    ggtitle(kfacsPlotTitle)+
    xlab("")+
    ylab("% of live cells")+
    guides(size=guide_legend("Number of events"))+
    scale_y_continuous(trans="sqrt",
                       labels=function(x) paste0(x,"%"))
    

#
#target CD14+ cells for PFDA
#
kfacsCD14Pops <- names(kfacsPvalVec)
kfacsCD14Pops <- kfacsCD14Pops[grepl("HLA-DR\\+",kfacsCD14Pops)]
kfacsCD14Pops <- kfacsCD14Pops[grepl("CD3-",kfacsCD14Pops)]
kfacsCD14Pops <- kfacsCD14Pops[grepl("CD56-",kfacsCD14Pops)]
kfacsCD14Pops <- kfacsCD14Pops[grepl("CD19-",kfacsCD14Pops)]
kfacsCD14Pops <- kfacsCD14Pops[grepl("CD14\\+",kfacsCD14Pops)]
kfacsCD14Pops <- kfacsCD14Pops[grepl("CD16-",kfacsCD14Pops)]

#
#fit the PFDA model
#
kfacsGlmerDFPrep <- kfacsCountDF[,c(kfacsCD14Pops,"parentCount","modelRT","sampleID")]
kfacsGlmerDF <- gather(kfacsGlmerDFPrep,key=cellPop,value=count,-c("parentCount","modelRT","sampleID"))
kfacsGlmerDF$popName <- as.factor(paste0("Pop_",as.numeric(as.factor(kfacsGlmerDF$cellPop))))
kfacsGlmerDF$obsFactor <- as.factor(paste0("Obs_",seq(nrow(kfacsGlmerDF))))
kfacsGlmerFit <- glmer(
    formula=cbind(count,(parentCount-count))  ~ popName*modelRT + (1|obsFactor),
    data=kfacsGlmerDF,
    family="binomial",
    control=glmerControl(
        optimizer="bobyqa",
        optCtrl = list(maxfun = 1e9),
        boundary.tol = 1e-9,
        tolPwrss=1e-9,
        check.conv.singular = .makeCC(action = "warning",  tol = 1e-12)
    )
)
kfacsCoefNames = names(fixef(kfacsGlmerFit))
kfacsGLHTMat <- matrix(0, nrow = 1, ncol = length(kfacsCoefNames))
kfacsGLHTMat[,which(kfacsCoefNames == "modelRT1")] <- 1
kfacsGLHTMat[1,which(kfacsCoefNames %like% ":modelRT1")] <- 1/(length(which(kfacsCoefNames%like%":modelRT1"))+1)
rownames(kfacsGLHTMat) <- c("CD14Monocytes")
kfacsMulti <- glht(kfacsGlmerFit,kfacsGLHTMat, alternative = "greater")
kfacsObsPval <- round(summary(kfacsMulti)$test$pvalues[1],9)
kfacsCT95 <- confint(kfacsMulti,level=0.95)
kfacsMultiResults <- as.data.frame(kfacsCT95$confint)
kfacsMultiResults$upr <- 3
kfacsMultiResults$Compartment <- rownames(kfacsMultiResults)
kfacsMultiResults$ModelType <- "PFDA"
kfacsMultiResults$Component <- "PFDA"

#
#get confidence intervals for univariate fits
#
kfacsCIListCD14 <- list()
for (cd14pt in kfacsCD14Pops) {
    mdf <- data.frame(
        count=kfacsCountDF[,cd14pt],
        parentCount=kfacsCountDF[,"parentCount"],
        modelRT=kfacsCountDF[,"modelRT"],
        subjectID=kfacsCountDF[,"sampleID"]
    )
    m <- glmer(
        formula=cbind(count,(parentCount-count))  ~ modelRT + (1|subjectID),
        data = mdf,
        family="binomial",
        control=glmerControl(
            optimizer="bobyqa",
            optCtrl = list(maxfun = 1e9),
            boundary.tol = 1e-9,
            tolPwrss=1e-9,
            check.conv.singular = .makeCC(action = "warning",  tol = 1e-12)
        )
    )
    #90 % CI since testing one-sided hypothesis
    ci1 <- confint(m,level=0.90,devtol=1e-5)
    outRes <- cbind(fixef(m)[1:2],ci1[2:3,])
    colnames(outRes) <- c("Estimate","lwr","upr")
    outRes[,"upr"] <- 3
    kfacsCIListCD14 <- append(kfacsCIListCD14,list(outRes))
    names(kfacsCIListCD14)[length(kfacsCIListCD14)] <- cd14pt
}
kfacsUniPlotDFPrep <- as.data.frame(Reduce(rbind,
                                       lapply(kfacsCIListCD14,
                                              function(x){as.numeric(x[which(rownames(x)=="modelRT1"),,drop=TRUE])})))
colnames(kfacsUniPlotDFPrep) <- c("Estimate","lwr","upr")
rownames(kfacsUniPlotDFPrep) <- seq(nrow(kfacsUniPlotDFPrep))
kfacsUniPlotDFPrep$Compartment <- as.character(paste0("Population ",seq(nrow(kfacsUniPlotDFPrep))))
kfacsUniPlotDFPrep$Method <- "FAUST\npopulations\nused by\nMultivariate"
kfacsUniPlotDFPrep[which(kfacsUniPlotDFPrep$Compartment=="14+"),"Compartment"] <- "Aggregate"
kfacsUniResults <- kfacsUniPlotDFPrep[order(kfacsUniPlotDFPrep[,"lwr"]),]
kfacsUniResults$ModelType <- "CD14+ CD16-\nHLA-DR+\nFAUST Phenotypes\nUnivariate CIs"
kfacsUniResults$Component <- paste0("Univariate_",str_pad(seq(nrow(kfacsUniResults)),pad="0",width=2))

#
#aggregate over cd14 phenotypes and test
#
kfacs_aggregate_meta_data <- data.frame(
    modelRT=kfacsCountDF[,"modelRT"],
    sampleName=kfacsCountDF[,"sampleID"],
    subjectID=as.factor(kfacsCountDF[,"sampleID"])
)

kfacs_aggregate_target_df <- faust::getCountsForTargetMarkers(
    projectPath=file.path(normalizePath("../.."),"faustRuns","melanoma_facs"),
    referencePhenotype="CD3-CD4+HLA-DR+CD14+CD19-CD11b+CD16-CD56-CD45RO+",
    targetMarkers=c("CD3","HLA-DR","CD14","CD19","CD16","CD56"),
    metaDataDF=kfacs_aggregate_meta_data
    )

kfacsAggregateModel <- glmer(
    formula=cbind(childCount,(parentCount-childCount))  ~ modelRT + (1|subjectID),
    data = kfacs_aggregate_target_df,
    family="binomial",
    control=glmerControl(
        optimizer="bobyqa",
        optCtrl = list(maxfun = 1e9),
        boundary.tol = 1e-9,
        tolPwrss=1e-9,
        check.conv.singular = .makeCC(action = "warning",  tol = 1e-12)
    )
)

kfacsAggregateCIs <- confint(kfacsAggregateModel,level=0.90,devtol=1e-5)
kfacsAggregateOutRes <- cbind(fixef(kfacsAggregateModel)[1:2],kfacsAggregateCIs[2:3,])
colnames(kfacsAggregateOutRes) <- c("Estimate","lwr","upr")
kfacsAggregateOutRes[,"upr"] <- 3
kfacsAggregateResults <- as.data.frame(kfacsAggregateOutRes[which(rownames(kfacsAggregateOutRes)=="modelRT1"),,drop=FALSE])
colnames(kfacsAggregateResults) <- c("Estimate","lwr","upr")
rownames(kfacsAggregateResults) <- c("CD14Monocytes")
kfacsAggregateResults$Compartment <- c("CD14Monocytes")
kfacsAggregateResults$compartmentSize <- length(kfacsCD14Pops)
kfacsAggregateResults$ModelType <- "Targeted"
kfacsAggregateResults$Component <- "Targeted"
#
#combine and plot
#
kfacsPFDAplotDF <- rbind(
    kfacsUniResults[,c("Estimate","lwr","upr","ModelType","Component")],
    kfacsMultiResults[,c("Estimate","lwr","upr","ModelType","Component")],
    kfacsAggregateResults[,c("Estimate","lwr","upr","ModelType","Component")]
)
kfacsPFDAplotDF[,c("Estimate","lwr","upr")] <- exp(kfacsPFDAplotDF[,c("Estimate","lwr","upr")])
kfacsPFDAplotDF$ModelType <- factor(kfacsPFDAplotDF$ModelType,levels=c("PFDA",
                                                                       "Targeted",
                                                                       "CD14+ CD16-\nHLA-DR+\nFAUST Phenotypes\nUnivariate CIs"))
pfdaPlotKFACS <- ggplot(kfacsPFDAplotDF, aes(x=ModelType, y=Estimate))+
    geom_errorbar(
        aes(ymin=lwr, ymax=40, colour=`Component`),
        size=0.05,
        width=0.5,
        position=position_dodge2(width=0.75,reverse=TRUE)
    )+
    geom_point(
        aes(group=`Component`),
        size=0.6,
        shape=15,
        position=position_dodge2(width=0.5,reverse=TRUE)
    )+
    geom_hline(
        yintercept=1,
        size=0.05,
        linetype="solid",
        color="red"
    )+
    coord_flip(xlim=c(1,3),ylim=c(exp(-5),exp(3)))+
    ylab("Point estimate of odds-ratio (square mark)\nwith one-sided 95% confidence interval")+
    scale_color_manual(values=rep("#000000FF",nrow(kfacsPFDAplotDF)))+
    theme_classic(base_size=6)+
    xlab("")+
    theme(
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_text(size=9),
        plot.title = element_text(hjust = 0.5)
    )+
    ggtitle(paste0("Melanoma anti-PD-1 trial, FACS"))+
    guides(colour=guide_legend(keywidth = 4, keyheight = 1))+
    theme(legend.position="none")+
    scale_y_continuous(trans="log",labels=scaleFUN,breaks=c(exp(-3),exp(-2),exp(-1),1,exp(1),exp(2)))        
########################################################################
#
#
#
#
#T cell panels from CITN-09
#
#
#
#
########################################################################
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
citn09DataPath <- file.path(normalizePath("../.."),"faustRuns","CITN-09","whole_blood_analysis")
allTcellCountDF <- as.data.frame(readRDS(file.path(citn09DataPath,"faustData","faustCountMatrix.rds")))
allTcellPanelPops <- colnames(allTcellCountDF)
tcellPops <- names(which(sapply(allTcellPanelPops,function(x){grepl("CD3\\+",x)})))
dpt <- Reduce(intersect,list(
                            cd4=tcellPops[grepl("CD4Bright",tcellPops)],
                            cd8=tcellPops[grepl("CD8\\+",tcellPops)]
                            ))
dpt2 <- Reduce(intersect,list(
                            cd4=tcellPops[grepl("CD4Dim",tcellPops)],
                            cd8=tcellPops[grepl("CD8\\+",tcellPops)]
                            ))
dnt <- Reduce(intersect,list(
                            cd4=tcellPops[grepl("CD4-",tcellPops)],
                            cd8=tcellPops[grepl("CD8-",tcellPops)]
                            ))
tcellModelPops <- setdiff(tcellPops,c(dpt,dpt2,dnt))
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
    referencePhenotype=tcellModelPops[1],
    targetMarkers=c("CD3","CD4","CD8"),
    metaDataDF=sel_meta_data,
    markersInParentPhenotype=c("CD3")
)
cd3ParentDF <- cd3CountDF[,c("sampleName","parentCount")]
allTcellCountDF <- left_join(allTcellCountDF,cd3ParentDF,by=c("sampleName"))
#
#subset to C01 for modeling and derive CR/PR merged responder type.
#
tcellCountDF <- allTcellCountDF[which(allTcellCountDF[,c("visitSTR")]=="C01"),]
tcellAnalysisDF <- left_join(tcellCountDF,meta_data,by="subjectID")
tcellSubDF <- tcellAnalysisDF[which(tcellAnalysisDF[,c("visitSTR")]=="C01"),]
tcellResType <- rep(0,nrow(tcellSubDF))
tcellResType[which(tcellSubDF$OverallResponse %in% c("CR","PR"))] <- 1
tcellSubDF$modelRT <- tcellResType
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
#fit GLMM for all CD4+ and CD8+ t cell populations discovered by faust
#
mePvalVec <- c()
for (cellPop in tcellModelPops) {
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
#
#fit the PFDA model
#
cd4CITN09ModelPops <- names(mePvalVec)
cd4CITN09ModelPops <- cd4CITN09ModelPops[grepl("CD3\\+",cd4CITN09ModelPops)]
cd4CITN09ModelPops <- cd4CITN09ModelPops[grepl("CD4Bright",cd4CITN09ModelPops)]
cd4CITN09ModelPops <- cd4CITN09ModelPops[grepl("CD8-",cd4CITN09ModelPops)]
cd4CITN09ModelPops <- setdiff(cd4CITN09ModelPops,cd4CITN09ModelPops[grepl("CD279 PD1-",cd4CITN09ModelPops)])
cd4CITN09GlmerDFPrep <- tcellSubDF[,c(cd4CITN09ModelPops,"parentCount","modelRT","subjectID")]
cd4CITN09GlmerDF <- gather(cd4CITN09GlmerDFPrep,key=cellPop,value=count,-c("parentCount","modelRT","subjectID"))
cd4CITN09GlmerDF$popName <- as.factor(paste0("Pop_",as.numeric(as.factor(cd4CITN09GlmerDF$cellPop))))
cd4CITN09GlmerDF$obsFactor <- as.factor(paste0("Obs_",seq(nrow(cd4CITN09GlmerDF))))
cd4CITN09GlmerDF$Response <- as.factor(cd4CITN09GlmerDF$modelRT)
cd4CITN09GlmerFit <- glmer(
    formula=cbind(count,(parentCount-count))  ~ popName*Response + (1|obsFactor),
    data=cd4CITN09GlmerDF,
    family="binomial",
    nAGQ=0,
    control=glmerControl(
        optimizer="bobyqa",
        optCtrl = list(maxfun = 1e9),
        boundary.tol = 1e-9,
        tolPwrss=1e-9,
        check.conv.singular = .makeCC(action = "warning",  tol = 1e-12)
    )
)
cd4CITN09CoefNames <- names(fixef(cd4CITN09GlmerFit))
cd4CITN09GLHTMat <- matrix(0, nrow = 1, ncol = length(cd4CITN09CoefNames))
cd4CITN09GLHTMat[,which(cd4CITN09CoefNames == "Response1")] <- 1
cd4CITN09GLHTMat[1,which(cd4CITN09CoefNames %like% ":Response1")] <- 1/(length(which(cd4CITN09CoefNames%like%":Response1"))+1)
rownames(cd4CITN09GLHTMat) <- c("CD4Tcells")
cd4CITN09Multi <- glht(cd4CITN09GlmerFit,cd4CITN09GLHTMat, alternative = "greater")
cd4CITN09ObsPval <- round(summary(cd4CITN09Multi)$test$pvalues[1],9)
cd4CITN09CT95 <- confint(cd4CITN09Multi,level=0.95)
cd4CITN09MultiResults <- as.data.frame(cd4CITN09CT95$confint)
cd4CITN09MultiResults$upr <- 3
cd4CITN09MultiResults$Compartment <- rownames(cd4CITN09MultiResults)
cd4CITN09MultiResults$ModelType <- "PFDA"
cd4CITN09MultiResults$Component <- "PFDA"
#
#get confidence intervals for univariate fits
#
cd4CITN09CIListCD14 <- list()
for (cd4pt in cd4CITN09ModelPops) {
    mdf <- data.frame(
        childCount=tcellSubDF[,cd4pt],
        parentCount=tcellSubDF[,"parentCount"],
        modelRT=as.factor(tcellSubDF[,"modelRT"]),
        subjectID=tcellSubDF[,"subjectID"]
    )
    m <- glmer(
        formula=cbind(childCount,(parentCount-childCount))  ~ modelRT + (1|subjectID),
        data=mdf,
        family="binomial",
        control=glmerControl(
            optimizer="bobyqa", optCtrl = list(maxfun = 1e9),
            check.conv.singular = .makeCC(action = "warning",  tol = 1e-4))
    )
    #90 % CI since testing one-sided hypothesis
    try({ci1 <- confint(m,level=0.90,devtol=1e-12)
        outRes <- cbind(fixef(m)[1:2],ci1[2:3,])
        colnames(outRes) <- c("Estimate","lwr","upr")
        outRes[,"upr"] <- 3
        cd4CITN09CIListCD14 <- append(cd4CITN09CIListCD14,list(outRes))
        names(cd4CITN09CIListCD14)[length(cd4CITN09CIListCD14)] <- cd14pt
        })
}
cd4CITN09UniPlotDFPrep <- as.data.frame(Reduce(rbind,
                                       lapply(cd4CITN09CIListCD14,
                                              function(x){as.numeric(x[which(rownames(x)=="modelRT1"),,drop=TRUE])})))
colnames(cd4CITN09UniPlotDFPrep) <- c("Estimate","lwr","upr")
rownames(cd4CITN09UniPlotDFPrep) <- seq(nrow(cd4CITN09UniPlotDFPrep))
cd4CITN09UniPlotDFPrep$Compartment <- "Univariate"
cd4CITN09UniPlotDFPrep$Method <- "FAUST\npopulations\nused by\nMultivariate"
cd4CITN09UniResults <- cd4CITN09UniPlotDFPrep[order(cd4CITN09UniPlotDFPrep[,"lwr"]),]
cd4CITN09UniResults$ModelType <- "CD4 Bright CD8- CD3+\n PD-1 Dim or Bright\nFAUST Phenotypes\nUnivariate CIs"
cd4CITN09UniResults$Component <- paste0("Univariate_",str_pad(seq(nrow(cd4CITN09UniResults)),pad="0",width=2))
#
#Aggregate
#
cd4citn09_aggregate_meta_data <- data.frame(
    modelRT=tcellSubDF[,"modelRT"],
    subjectID=as.factor(tcellSubDF[,"subjectID"]),
    sampleName=tcellSubDF[,"sampleName"]
)

cd4citn09_pd1dim_aggregate_target_df <- faust::getCountsForTargetMarkers(
    projectPath=file.path(normalizePath("../.."),"faustRuns","CITN-09","whole_blood_analysis"),
    referencePhenotype="CD4BrightCD8-CD3+CD45RA-HLA DR-CD279 PD1DimCD28+CD127-CD25-CD197 CCR7-",
    targetMarkers=c("CD4","CD8","CD3","CD279 PD1"),
    metaDataDF=cd4citn09_aggregate_meta_data,
    markersInParentPhenotype = c("CD3")
    )

cd4citn09_pd1bright_aggregate_target_df <- faust::getCountsForTargetMarkers(
    projectPath=file.path(normalizePath("../.."),"faustRuns","CITN-09","whole_blood_analysis"),
    referencePhenotype="CD4BrightCD8-CD3+CD45RA+HLA DR+CD279 PD1BrightCD28+CD127+CD25-CD197 CCR7+",
    targetMarkers=c("CD4","CD8","CD3","CD279 PD1"),
    metaDataDF=cd4citn09_aggregate_meta_data,
    markersInParentPhenotype = c("CD3")
    )

cd4citn09CN1 <- colnames(cd4citn09_pd1dim_aggregate_target_df)
cd4citn09CN1[which(cd4citn09CN1=="childCount")] <- "pd1dim_childCount"
colnames(cd4citn09_pd1dim_aggregate_target_df) <- cd4citn09CN1
cd4brightcitn09CN1 <- colnames(cd4citn09_pd1bright_aggregate_target_df)
cd4brightcitn09CN1[which(cd4brightcitn09CN1=="childCount")] <- "pd1bright_childCount"
colnames(cd4citn09_pd1bright_aggregate_target_df) <- cd4brightcitn09CN1
cd4citn09_aggregate_target_df <- inner_join(cd4citn09_pd1dim_aggregate_target_df,cd4citn09_pd1bright_aggregate_target_df,
                                            by=c("sampleName","parentCount","subjectID","modelRT"))
cd4citn09_aggregate_target_df$childCount <- (cd4citn09_aggregate_target_df$`pd1dim_childCount`+cd4citn09_aggregate_target_df$`pd1bright_childCount`)

cd4citn09AggregateModel <- glmer(
    cbind(childCount,(parentCount-childCount)) ~ modelRT + (1|subjectID),
    data=cd4citn09_aggregate_target_df,
    family="binomial",
    control=glmerControl(
        optimizer="bobyqa", optCtrl = list(maxfun = 1e9),
        check.conv.singular = .makeCC(action = "warning",  tol = 1e-4)
    )
)
cd4citn09AggregateCI <- confint(cd4citn09AggregateModel,level=0.90,devtol=1e-5)
cd4citn09AggregateOutRes <- cbind(fixef(cd4citn09AggregateModel),cd4citn09AggregateCI[2:3,])
colnames(cd4citn09AggregateOutRes) <- c("Estimate","lwr","upr")
cd4citn09AggregateOutRes[,"upr"] <- 3
cd4citn09AggregateResults <- as.data.frame(cd4citn09AggregateOutRes[which(rownames(cd4citn09AggregateOutRes)=="modelRT"),,drop=FALSE])
colnames(cd4citn09AggregateResults) <- c("Estimate","lwr","upr")
rownames(cd4citn09AggregateResults) <- c("CD4 Tcells")
cd4citn09AggregateResults$Compartment <- c("CD4 Tcells")
cd4citn09AggregateResults$compartmentSize <- length(cd4CITN09ModelPops) 
cd4citn09AggregateResults$ModelType <- "Targeted"
cd4citn09AggregateResults$Component <- "Targeted"
#
#combine and plot
#
cd4citn09PFDAplotDF <- rbind(
    cd4CITN09UniResults[,c("Estimate","lwr","upr","ModelType","Component")],
    cd4CITN09MultiResults[,c("Estimate","lwr","upr","ModelType","Component")],
    cd4citn09AggregateResults[,c("Estimate","lwr","upr","ModelType","Component")]
)
cd4citn09PFDAplotDF[,c("Estimate","lwr","upr")] <- exp(cd4citn09PFDAplotDF[,c("Estimate","lwr","upr")])
cd4citn09PFDAplotDF$ModelType <- factor(cd4citn09PFDAplotDF$ModelType,levels=c("PFDA",
                                                                       "Targeted",
                                                                       "CD4 Bright CD8- CD3+\n PD-1 Dim or Bright\nFAUST Phenotypes\nUnivariate CIs"))

pfdaPlotCD4CITN09 <- ggplot(cd4citn09PFDAplotDF, aes(x=ModelType, y=Estimate))+
    geom_errorbar(
        aes(ymin=lwr, ymax=40, colour=`Component`),
        size=0.05,
        width=0.5,
        position=position_dodge2(width=0.75,reverse=TRUE)
    )+
    geom_point(
        aes(group=`Component`),
        size=0.6,
        shape=15,
        position=position_dodge2(width=0.5,reverse=TRUE)
    )+
    geom_hline(
        yintercept=1,
        size=0.05,
        linetype="solid",
        color="red"
    )+
    coord_flip(xlim=c(1,3),ylim=c(exp(-5),exp(3)))+
    ylab("Point estimate of odds-ratio (square mark)\nwith one-sided 95% confidence interval")+
    scale_color_manual(values=rep("#000000FF",nrow(cd4citn09PFDAplotDF)))+
    theme_classic(base_size=6)+
    xlab("")+
    theme(
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_text(size=9),
        plot.title = element_text(hjust = 0.5)
    )+
    ggtitle(paste0("MCC anti-PD-1 trial, T cell Panel, CD4"))+
    guides(colour=guide_legend(keywidth = 4, keyheight = 1))+
    theme(legend.position="none")+
    scale_y_continuous(trans="log",labels=scaleFUN,breaks=c(exp(-3),exp(-2),exp(-1),1,exp(1),exp(2)))        

################################################################################################################################################
#
#
#
#
#
#CD8 pfda
#
#
#
#
#
################################################################################################################################################
#
#fit the PFDA model
#
cd8CITN09ModelPops <- names(mePvalVec)
cd8CITN09ModelPops <- cd8CITN09ModelPops[grepl("CD3\\+",cd8CITN09ModelPops)]
cd8CITN09ModelPops <- cd8CITN09ModelPops[grepl("CD4-",cd8CITN09ModelPops)]
cd8CITN09ModelPops <- cd8CITN09ModelPops[grepl("CD8\\+",cd8CITN09ModelPops)]
cd8CITN09ModelPops <- setdiff(cd8CITN09ModelPops,cd8CITN09ModelPops[grepl("CD279 PD1-",cd8CITN09ModelPops)])
cd8CITN09GlmerDFPrep <- tcellSubDF[,c(cd8CITN09ModelPops,"parentCount","modelRT","subjectID")]
cd8CITN09GlmerDF <- gather(cd8CITN09GlmerDFPrep,key=cellPop,value=count,-c("parentCount","modelRT","subjectID"))
cd8CITN09GlmerDF$popName <- as.factor(paste0("Pop_",as.numeric(as.factor(cd8CITN09GlmerDF$cellPop))))
cd8CITN09GlmerDF$obsFactor <- as.factor(paste0("Obs_",seq(nrow(cd8CITN09GlmerDF))))
cd8CITN09GlmerDF$Response <- as.factor(cd8CITN09GlmerDF$modelRT)
cd8CITN09GlmerFit <- glmer(
    formula=cbind(count,(parentCount-count))  ~ popName*Response + (1|obsFactor),
    data=cd8CITN09GlmerDF,
    family="binomial",
    nAGQ=0,
    control=glmerControl(
        optimizer="bobyqa",
        optCtrl = list(maxfun = 1e9),
        boundary.tol = 1e-9,
        tolPwrss=1e-9,
        check.conv.singular = .makeCC(action = "warning",  tol = 1e-12)
    )
)
cd8CITN09CoefNames <- names(fixef(cd8CITN09GlmerFit))
cd8CITN09GLHTMat <- matrix(0, nrow = 1, ncol = length(cd8CITN09CoefNames))
cd8CITN09GLHTMat[,which(cd8CITN09CoefNames == "Response1")] <- 1
cd8CITN09GLHTMat[1,which(cd8CITN09CoefNames %like% ":Response1")] <- 1/(length(which(cd8CITN09CoefNames%like%":Response1"))+1)
rownames(cd8CITN09GLHTMat) <- c("CD8Tcells")
cd8CITN09Multi <- glht(cd8CITN09GlmerFit,cd8CITN09GLHTMat, alternative = "greater")
cd8CITN09ObsPval <- round(summary(cd8CITN09Multi)$test$pvalues[1],9)
cd8CITN09CT95 <- confint(cd8CITN09Multi,level=0.95)
cd8CITN09MultiResults <- as.data.frame(cd8CITN09CT95$confint)
cd8CITN09MultiResults$upr <- 3
cd8CITN09MultiResults$Compartment <- rownames(cd8CITN09MultiResults)
cd8CITN09MultiResults$ModelType <- "PFDA"
cd8CITN09MultiResults$Component <- "PFDA"
#
#get confidence intervals for univariate fits
#
cd8CITN09CIListCD14 <- list()
for (cd8pt in cd8CITN09ModelPops) {
    mdf <- data.frame(
        childCount=tcellSubDF[,cd8pt],
        parentCount=tcellSubDF[,"parentCount"],
        modelRT=as.factor(tcellSubDF[,"modelRT"]),
        subjectID=tcellSubDF[,"subjectID"]
    )
    m <- glmer(
        formula=cbind(childCount,(parentCount-childCount))  ~ modelRT + (1|subjectID),
        data=mdf,
        family="binomial",
        control=glmerControl(
            optimizer="bobyqa", optCtrl = list(maxfun = 1e9),
            check.conv.singular = .makeCC(action = "warning",  tol = 1e-4))
    )
    #90 % CI since testing one-sided hypothesis
    try({ci1 <- confint(m,level=0.90,devtol=1e-12)
        outRes <- cbind(fixef(m)[1:2],ci1[2:3,])
        colnames(outRes) <- c("Estimate","lwr","upr")
        outRes[,"upr"] <- 3
        cd8CITN09CIListCD14 <- append(cd8CITN09CIListCD14,list(outRes))
        names(cd8CITN09CIListCD14)[length(cd8CITN09CIListCD14)] <- cd14pt
        })
}
cd8CITN09UniPlotDFPrep <- as.data.frame(Reduce(rbind,
                                       lapply(cd8CITN09CIListCD14,
                                              function(x){as.numeric(x[which(rownames(x)=="modelRT1"),,drop=TRUE])})))
colnames(cd8CITN09UniPlotDFPrep) <- c("Estimate","lwr","upr")
rownames(cd8CITN09UniPlotDFPrep) <- seq(nrow(cd8CITN09UniPlotDFPrep))
cd8CITN09UniPlotDFPrep$Compartment <- "Univariate"
cd8CITN09UniPlotDFPrep$Method <- "FAUST\npopulations\nused by\nMultivariate"
cd8CITN09UniResults <- cd8CITN09UniPlotDFPrep[order(cd8CITN09UniPlotDFPrep[,"lwr"]),]
cd8CITN09UniResults$ModelType <- "CD8+ CD4- CD3+\n PD-1 Dim or Bright\nFAUST Phenotypes\nUnivariate CIs"
cd8CITN09UniResults$Component <- paste0("Univariate_",str_pad(seq(nrow(cd8CITN09UniResults)),pad="0",width=2))
#
#Aggregate
#
cd8citn09_aggregate_meta_data <- data.frame(
    modelRT=tcellSubDF[,"modelRT"],
    subjectID=as.factor(tcellSubDF[,"subjectID"]),
    sampleName=tcellSubDF[,"sampleName"]
)

cd8citn09_pd1dim_aggregate_target_df <- faust::getCountsForTargetMarkers(
    projectPath=file.path(normalizePath("../.."),"faustRuns","CITN-09","whole_blood_analysis"),
    referencePhenotype="CD4-CD8+CD3+CD45RA-HLA DR+CD279 PD1DimCD28+CD127-CD25-CD197 CCR7-",
    targetMarkers=c("CD4","CD8","CD3","CD279 PD1"),
    metaDataDF=cd8citn09_aggregate_meta_data,
    markersInParentPhenotype = c("CD3")
    )

cd8citn09_pd1bright_aggregate_target_df <- faust::getCountsForTargetMarkers(
    projectPath=file.path(normalizePath("../.."),"faustRuns","CITN-09","whole_blood_analysis"),
    referencePhenotype="CD4-CD8+CD3+CD45RA-HLA DR+CD279 PD1BrightCD28+CD127-CD25-CD197 CCR7-",
    targetMarkers=c("CD4","CD8","CD3","CD279 PD1"),
    metaDataDF=cd8citn09_aggregate_meta_data,
    markersInParentPhenotype = c("CD3")
    )
cd8citn09CN1 <- colnames(cd8citn09_pd1dim_aggregate_target_df)
cd8citn09CN1[which(cd8citn09CN1=="childCount")] <- "pd1dim_childCount"
colnames(cd8citn09_pd1dim_aggregate_target_df) <- cd8citn09CN1
cd8brightcitn09CN1 <- colnames(cd8citn09_pd1bright_aggregate_target_df)
cd8brightcitn09CN1[which(cd8brightcitn09CN1=="childCount")] <- "pd1bright_childCount"
colnames(cd8citn09_pd1bright_aggregate_target_df) <- cd8brightcitn09CN1
cd8citn09_aggregate_target_df <- inner_join(cd8citn09_pd1dim_aggregate_target_df,cd8citn09_pd1bright_aggregate_target_df,
                                            by=c("sampleName","parentCount","subjectID","modelRT"))
cd8citn09_aggregate_target_df$childCount <- (cd8citn09_aggregate_target_df$`pd1dim_childCount`+cd8citn09_aggregate_target_df$`pd1bright_childCount`)
cd8citn09AggregateModel <- glmer(
    cbind(childCount,(parentCount-childCount)) ~ modelRT + (1|subjectID),
    data=cd8citn09_aggregate_target_df,
    family="binomial",
    control=glmerControl(
        optimizer="bobyqa", optCtrl = list(maxfun = 1e9),
        check.conv.singular = .makeCC(action = "warning",  tol = 1e-4)
    )
)
cd8citn09AggregateCI <- confint(cd8citn09AggregateModel,level=0.90,devtol=1e-5)
cd8citn09AggregateOutRes <- cbind(fixef(cd8citn09AggregateModel),cd8citn09AggregateCI[2:3,])
colnames(cd8citn09AggregateOutRes) <- c("Estimate","lwr","upr")
cd8citn09AggregateOutRes[,"upr"] <- 3
cd8citn09AggregateResults <- as.data.frame(cd8citn09AggregateOutRes[which(rownames(cd8citn09AggregateOutRes)=="modelRT"),,drop=FALSE])
colnames(cd8citn09AggregateResults) <- c("Estimate","lwr","upr")
rownames(cd8citn09AggregateResults) <- c("CD8 Tcells")
cd8citn09AggregateResults$Compartment <- c("CD8 Tcells")
cd8citn09AggregateResults$compartmentSize <- length(cd8CITN09ModelPops) 
cd8citn09AggregateResults$ModelType <- "Targeted"
cd8citn09AggregateResults$Component <- "Targeted"
#
#combine and plot
#
cd8citn09PFDAplotDF <- rbind(
    cd8CITN09UniResults[,c("Estimate","lwr","upr","ModelType","Component")],
    cd8CITN09MultiResults[,c("Estimate","lwr","upr","ModelType","Component")],
    cd8citn09AggregateResults[,c("Estimate","lwr","upr","ModelType","Component")]
)
cd8citn09PFDAplotDF[,c("Estimate","lwr","upr")] <- exp(cd8citn09PFDAplotDF[,c("Estimate","lwr","upr")])
cd8citn09PFDAplotDF$ModelType <- factor(cd8citn09PFDAplotDF$ModelType,levels=c("PFDA",
                                                                               "Targeted",
                                                                               "CD8+ CD4- CD3+\n PD-1 Dim or Bright\nFAUST Phenotypes\nUnivariate CIs"))
pfdaPlotCD8CITN09 <- ggplot(cd8citn09PFDAplotDF, aes(x=ModelType, y=Estimate))+
    geom_errorbar(
        aes(ymin=lwr, ymax=40, colour=`Component`),
        size=0.05,
        width=0.5,
        position=position_dodge2(width=0.75,reverse=TRUE)
    )+
    geom_point(
        aes(group=`Component`),
        size=0.6,
        shape=15,
        position=position_dodge2(width=0.5,reverse=TRUE)
    )+
    geom_hline(
        yintercept=1,
        size=0.05,
        linetype="solid",
        color="red"
    )+
    coord_flip(xlim=c(1,3),ylim=c(exp(-5),exp(3)))+
    ylab("Point estimate of odds-ratio (square mark)\nwith one-sided 95% confidence interval")+
    scale_color_manual(values=rep("#000000FF",nrow(cd8citn09PFDAplotDF)))+
    theme_classic(base_size=6)+
    xlab("")+
    theme(
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_text(size=9),
        plot.title = element_text(hjust = 0.5)
    )+
    ggtitle(paste0("MCC anti-PD-1 trial, T cell Panel, CD8"))+
    guides(colour=guide_legend(keywidth = 4, keyheight = 1))+
    theme(legend.position="none")+
    scale_y_continuous(trans="log",labels=scaleFUN,breaks=c(exp(-3),exp(-2),exp(-1),1,exp(1),exp(2)))        
row1 <- plot_grid(pfdaPlotCD4CITN09,pfdaPlotCD8CITN09,nrow=1)
row2 <- plot_grid(pfdaPlotCITN09, pfdaPlotCITN07, pfdaPlotKFACS,nrow=1)
plots <- align_plots(row1,row2,align='v',axis='l')                    
pOut <- plot_grid(plots[[1]],plots[[2]],ncol=1,labels=c("A","B"),label_size=20)
######################################################
#
#
#updated lay out of the figure
#
#
######################################################
citn09DF <- citn09PFDAplotDF
citn07DF <- citn07PFDAplotDF
kfDF <- kfacsPFDAplotDF
citn09DF$mType <- paste0("MCC anti-PD-1 trial\n",
                         gsub("FAUST Phenotypes\n","",gsub("CD14\\+ CD16-\nHLA-DR bright\n","",as.character(citn09DF[,"ModelType"]))))
citn09DF$comp <- paste0("MCC anti-PD-1 trial\n",as.character(citn09DF[,"Component"]))
citn07DF$mType <- paste0("FLT3-ligand +\nTherapeutic Vx trial\n",
                         gsub("FAUST Phenotypes\n","",gsub("CD14\\+ CD16-\nHLA-DR bright\n","",as.character(citn07DF[,"ModelType"]))))
citn07DF$comp <- paste0("FLT3-ligand +\nTherapeutic Vx trial\n",as.character(citn07DF[,"Component"]))
kfDF$mType <- paste0("Melanoma\nanti-PD-1 trial\n",
                     gsub("FAUST Phenotypes\n","",gsub("CD14\\+ CD16-\nHLA-DR\\+\n","",as.character(kfDF[,"ModelType"]))))
kfDF$comp <- paste0("Melanoma\nanti-PD-1 trial\n",as.character(kfDF[,"Component"]))
crossDF <- rbind(rbind(citn09DF,citn07DF),kfDF)
mcc_p_len <- length(intersect(which(grepl("MCC",crossDF$comp)),which(grepl("Univariate",crossDF$comp))))
flt3_p_len <- length(intersect(which(grepl("FLT3",crossDF$comp)),which(grepl("Univariate",crossDF$comp))))
mel_p_len <- length(intersect(which(grepl("Melanoma",crossDF$comp)),which(grepl("Univariate",crossDF$comp))))
pfdaPlotCross <- ggplot(crossDF, aes(x=mType, y=Estimate))+
    geom_errorbar(
        aes(ymin=lwr, ymax=40, colour=`comp`,
            width=c(rep(0.5,mcc_p_len),0.15,0.15,rep(0.5,flt3_p_len),0.15,0.15,rep(0.5,mel_p_len),0.15,0.15)),
        size=0.05,
        position=position_dodge2(width=0.75,reverse=TRUE)
    )+
    geom_point(
        aes(group=`comp`),
        size=0.6,
        shape=15,
        position=position_dodge2(width=0.5,reverse=TRUE)
    )+
    geom_hline(
        yintercept=1,
        size=0.05,
        linetype="solid",
        color="red"
    )+
    coord_flip(xlim=c(1,9),ylim=c(exp(-2.75),exp(3.15)))+
    ylab("Point estimate of odds-ratio (square mark)\nwith one-sided 95% confidence interval")+
    scale_color_manual(values=rep("#000000FF",nrow(crossDF)))+
    theme_classic(base_size=8)+
    xlab("")+
    theme(
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5)
    )+
    ggtitle("CD14+ CD16- HLA-DR+/bright\ncompartment, multiple trials")+
    guides(colour=guide_legend(keywidth = 4, keyheight = 1))+
    theme(legend.position="none")+
    scale_y_continuous(trans="log",labels=scaleFUN,breaks=c(exp(-2),exp(-1),1,exp(1),exp(2)))        
cd4DF <- cd4citn09PFDAplotDF
cd8DF <- cd8citn09PFDAplotDF
cd4DF$mType <- as.character(cd4DF[,"ModelType"])
cd4DF$mType <- gsub("FAUST Phenotypes\n","",gsub("CD4 Bright CD8- CD3\\+\\n PD-1 Dim or Bright\\n","",cd4DF$mType))
cd4DF$mType <- paste0("CD4 compartment\n",cd4DF$mType)
cd4DF$comp <- paste0("CD4 compartment\n",as.character(cd4DF[,"Component"]))
cd8DF$mType <- as.character(cd8DF[,"ModelType"])
cd8DF$mType <- gsub("FAUST Phenotypes\n","",gsub("CD8\\+ CD4- CD3\\+\\n PD-1 Dim or Bright\\n","",cd8DF$mType))
cd8DF$mType <- paste0("CD8 compartment\n",cd8DF$mType)
cd8DF$comp <- paste0("CD8 compartment\n",as.character(cd8DF[,"Component"]))
tcellDF <- rbind(cd4DF,cd8DF)
cd8_p_len <- length(intersect(which(grepl("CD8",tcellDF$comp)),which(grepl("Univariate",tcellDF$comp))))
cd4_p_len <- length(intersect(which(grepl("CD4",tcellDF$comp)),which(grepl("Univariate",tcellDF$comp))))
pfdaPlotTcell <- ggplot(tcellDF, aes(x=as.factor(mType), y=Estimate))+
    geom_errorbar(
        aes(ymin=lwr, ymax=40, colour=`comp`,width=c(rep(1,cd4_p_len),0.1,0.1,rep(1,cd8_p_len),0.1,0.1)),
        size=0.05,
        position=position_dodge2(width=1,reverse=TRUE)
    )+
    geom_point(
        aes(group=`comp`),
        size=0.6,
        shape=15,
        position=position_dodge2(width=1,reverse=TRUE)
    )+
    geom_hline(
        yintercept=1,
        size=0.05,
        linetype="solid",
        color="red"
    )+
    coord_flip(xlim=c(1,6),ylim=c(exp(-2.1),exp(2)))+
    ylab("Point estimate of odds-ratio (square mark)\nwith one-sided 95% confidence interval")+
    scale_color_manual(values=rep("#000000FF",nrow(tcellDF)))+
    theme_classic(base_size=8)+
    xlab("")+
    theme(
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5)
    )+
    ggtitle("T cell compartments\nMCC anti-PD-1 trial")+
    guides(colour=guide_legend(keywidth = 4, keyheight = 1))+
    theme(legend.position="none")+
    scale_y_continuous(trans="log",labels=scaleFUN,breaks=c(exp(-2),exp(-1),1,exp(1),exp(2)))
fig_grid <- plot_grid(
    pfdaPlotTcell,
    pfdaPlotCross,
    nrow=1,
    ncol=2,
    labels=c("A","B"),
    label_size=12
)
ggsave(
    filename=file.path(normalizePath("."),"figure_cross_dataset.png"),
    plot=fig_grid,
    width=6,
    height=7,
    units="in",
    dpi=300
)
#
#generate supplement plot
#
supplement_plot <- plot_grid(
    boxplotCITN09,boxplotCITN07,boxplotKFACS,
    nrow=2,
    labels=c("A","B","C"),
    label_size=12
)

ggsave(
    filename=file.path(normalizePath(".."),"supplementary_figures","supplement_myeloid_plot.png"),
    plot=supplement_plot,
    width=8,
    height=8,
    units="in",
    dpi=300
)



