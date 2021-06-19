library(RColorBrewer)
library(viridis)
library(ggplot2)
library(cowplot)
library(Rphenograph)
library(uwot)
library(faust)
#
#construct the annotation embedding
#
set.seed(1234) 
sampleName <- "PT ID_Subject-45 C01_025.fcs_828269"
projectPath <- file.path(normalizePath("../.."),"faustRuns","CITN-09","whole_blood_analysis")
umapEmbed <- makeAnnotationEmbedding(
    projectPath=projectPath,
    sampleNameVec=sampleName
)
#
#embed the raw expression matrix (collected by the makeAnnotationEmbedding function) for comparison
#
exprsMat <- umapEmbed[,c("CD278 ICOS","CD3", "CD127", "CD197 CCR7", "CD279 PD1", "CD8", "CD4", "CD28", "CD25", "HLA DR", "CD45RA"),
                      drop=FALSE]
set.seed(12345)
sampleUmapEmbed <- umap(exprsMat,n_neighbors=15,metric="euclidean",min_dist=0.2)
umapEmbed$ex <- sampleUmapEmbed[,1]
umapEmbed$ey <- sampleUmapEmbed[,2]
#
#relabel columns and get indices for plotting, set colors, etc.
#
inputCN <- colnames(umapEmbed)
updateCN <- gsub("CD279 PD1","PD1",inputCN)
updateCN <- gsub("CD197 CCR7","CCR7",updateCN)
updateCN <- gsub("CD278 ICOS","ICOS",updateCN)
updateCN <- gsub("HLA DR","HLADR",updateCN)
colnames(umapEmbed) <- updateCN
correlatesLookup <- rep(0,nrow(umapEmbed))
correlatesLookup[which(umapEmbed$faustLabels=="CD4-CD8+CD3+CD45RA-HLA DR+CD279 PD1DimCD28+CD127-CD25-CD197 CCR7-")] <- 1
correlatesLookup[which(umapEmbed$faustLabels=="CD4-CD8+CD3+CD45RA-HLA DR+CD279 PD1BrightCD28+CD127-CD25-CD197 CCR7-")] <- 2
correlatesLookup[which(umapEmbed$faustLabels=="CD4BrightCD8-CD3+CD45RA-HLA DR-CD279 PD1DimCD28+CD127-CD25-CD197 CCR7-")] <- 3
correlatesLookup[which(umapEmbed$faustLabels=="CD4BrightCD8-CD3+CD45RA-HLA DR+CD279 PD1DimCD28+CD127-CD25-CD197 CCR7-")] <- 4
umapEmbed$correlates <- as.factor(correlatesLookup)
nullLookup <- which(umapEmbed$correlates==0)
corPlotLookup <- c(nullLookup,setdiff(seq(nrow(umapEmbed)),nullLookup))
#
#fix seed for color call outs
#
set.seed(123) 
umapEmbed$faustLabels <- as.factor(umapEmbed$faustLabels)
nfClust <- length(table(umapEmbed$faustLabels))
myFaustColors <- colorRampPalette(brewer.pal(12, "Paired"))(nfClust)
myFaustColors <- myFaustColors[sample(seq(length(myFaustColors)),length(myFaustColors))]
names(myFaustColors) <- levels(umapEmbed$faustLabels)
myFaustColors[which(names(myFaustColors)=="0_0_0_0_0")] <- "#ccccccFF"
myFaustColors[which(names(myFaustColors)=="CD4-CD8+CD3+CD45RA-HLA DR+CD279 PD1DimCD28+CD127-CD25-CD197 CCR7-")] <- "#d7191cFF"
myFaustColors[which(names(myFaustColors)=="CD4-CD8+CD3+CD45RA-HLA DR+CD279 PD1BrightCD28+CD127-CD25-CD197 CCR7-")] <- "#abd9e9FF"
myFaustColors[which(names(myFaustColors)=="CD4BrightCD8-CD3+CD45RA-HLA DR-CD279 PD1DimCD28+CD127-CD25-CD197 CCR7-")] <- "#fdae61FF"
myFaustColors[which(names(myFaustColors)=="CD4BrightCD8-CD3+CD45RA-HLA DR+CD279 PD1DimCD28+CD127-CD25-CD197 CCR7-")] <- "#2c7bb6FF"
faustPlotIndex <- seq(nrow(umapEmbed))
faustClusterIndices <- which(umapEmbed$faustLabels=="0_0_0_0_0")
faustPlotIndex <- setdiff(faustPlotIndex,faustClusterIndices)
faustPlotIndex <- append(faustClusterIndices,faustPlotIndex)
#
#extract legends
#
pLeg1 <- ggplot(umapEmbed,aes_string(x="umapX",y="umapY",color="PD1_Windsorized"))+
    geom_point()+
    theme(legend.position="bottom")+
    scale_color_viridis()+
    labs(color="Standardized\nexpression")
legend1 <- get_legend(pLeg1)
ggdraw(legend1)
pd1Leg <- umapEmbed$`PD1_faust_annotation`
pd1Leg[which(pd1Leg=="-")] <- " -"
pd1Leg[which(pd1Leg=="Dim")] <- " Dim"
pd1Leg[which(pd1Leg=="Bright")] <- "+/Bright"
umapEmbed$pd1Leg <- pd1Leg
pLeg2 <- ggplot(umapEmbed,aes_string(x="umapX",y="umapY",color="pd1Leg"))+
    geom_point(size=0.1)+
    theme_classic()+
    theme(legend.position="bottom")+
    scale_color_manual(
        values=c("#ece7f2ff","#a6bddbff","#2b8cbeff"),
        guide=guide_legend(
            title="FAUST Annotation",
            override.aes=list(
                size=8,
                shape=15
            ),
            label.theme = element_text(size=15,face="bold")
        )
    )
legend2 <- get_legend(pLeg2)
pLeg3 <- ggplot(umapEmbed,aes_string(x="umapX",y="umapY",color="CD3_faust_annotation"))+
    geom_point(size=0.1)+
    theme_classic()+
    theme(legend.position="bottom")+
    scale_color_manual(
        values=c("#ece7f2ff","#2b8cbeff"),
        guide=guide_legend(
            title="FAUST Annotation",
            override.aes=list(
                size=8,
                shape=15
            ),
            label.theme = element_text(size=15,face="bold")
        )
    )
legend3 <- get_legend(pLeg3)
ggdraw(legend3)
tLegend <- ggdraw(legend1)
bLegend <- ggdraw(legend2)
all_legends <- plot_grid(ggdraw(legend1),ggdraw(legend2),nrow=1)
#
#plot parameters
#
ptSize <- 0.001

#
#annotation embedding
#
anEmCD4 <- ggplot(umapEmbed,aes_string(x="umapX",y="umapY",color="CD4_Windsorized"))+
    geom_point(size=ptSize)+
    theme_classic()+
    scale_color_viridis(option="D")+
    theme(
        legend.position="none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )+
    xlab("annotation embedding X")+
    ylab("annotation embedding Y")+
    ggtitle("CD4 expression level")

anEmCD3 <- ggplot(umapEmbed,aes_string(x="umapX",y="umapY",color="CD3_Windsorized"))+
    geom_point(size=ptSize)+
    theme_classic()+
    scale_color_viridis(option="D")+
    theme(
        legend.position="none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )+
    xlab("annotation embedding X")+
    ylab("annotation embedding Y")+
    ggtitle("CD3 expression level")

anEmPD1 <- ggplot(umapEmbed,aes_string(x="umapX",y="umapY",color="PD1_Windsorized"))+
    geom_point(size=ptSize)+
    theme_classic()+
    scale_color_viridis(option="D")+
    theme(
        legend.position="none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )+
    xlab("annotation embedding X")+
    ylab("annotation embedding Y")+
    ggtitle("PD-1 expression level")

anEmHLADR <- ggplot(umapEmbed,aes_string(x="umapX",y="umapY",color="HLADR_Windsorized"))+
    geom_point(size=ptSize)+
    theme_classic()+
    scale_color_viridis(option="D")+
    theme(
        legend.position="none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )+
    xlab("annotation embedding X")+
    ylab("annotation embedding Y")+
    ggtitle("HLA-DR expression level")

#
#expression embedding
#
exEmCD4 <- ggplot(umapEmbed,aes_string(x="ex",y="ey",color="CD4_Windsorized"))+
    geom_point(size=ptSize)+
    theme_classic()+
    scale_color_viridis(option="D")+
    theme(
        legend.position="none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )+
    xlab("expression embedding X")+
    ylab("expression embedding Y")+
    ggtitle("CD4 expression level")

exEmCD3 <- ggplot(umapEmbed,aes_string(x="ex",y="ey",color="CD3_Windsorized"))+
    geom_point(size=ptSize)+
    theme_classic()+
    scale_color_viridis(option="D")+
    theme(
        legend.position="none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )+
    xlab("expression embedding X")+
    ylab("expression embedding Y")+
    ggtitle("CD3 expression level")

exEmPD1 <- ggplot(umapEmbed,aes_string(x="ex",y="ey",color="PD1_Windsorized"))+
    geom_point(size=ptSize)+
    theme_classic()+
    scale_color_viridis(option="D")+
    theme(
        legend.position="none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )+
    xlab("expression embedding X")+
    ylab("expression embedding Y")+
    ggtitle("PD-1 expression level")

exEmHLADR <- ggplot(umapEmbed,aes_string(x="ex",y="ey",color="HLADR_Windsorized"))+
    geom_point(size=ptSize)+
    theme_classic()+
    scale_color_viridis(option="D")+
    theme(
        legend.position="none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )+
    xlab("expression embedding X")+
    ylab("expression embedding Y")+
    ggtitle("HLADR expression level")

#
#the associated faust annotations, annotation embedding
#
anEmFACD4<- ggplot(umapEmbed,aes_string(x="umapX",y="umapY",color="CD4_faust_annotation"))+
        geom_point(size=ptSize)+
        theme_classic()+
        scale_color_manual(values=c("#ece7f2ff","#2b8cbeff","#a6bddbff"))+
        theme(
            legend.position="none",
            plot.title = element_text(hjust = 0.5),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
        )+
    xlab("annotation embedding X")+
    ylab("annotation embedding Y")+
    ggtitle("CD4 FAUST annotation")

anEmFACD3 <- ggplot(umapEmbed,aes_string(x="umapX",y="umapY",color="CD3_faust_annotation"))+
        geom_point(size=ptSize)+
        theme_classic()+
        scale_color_manual(values=c("#ece7f2ff","#2b8cbeff"))+
        theme(
            legend.position="none",
            plot.title = element_text(hjust = 0.5),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
        )+
    xlab("annotation embedding X")+
    ylab("annotation embedding Y")+
    ggtitle("CD3 FAUST annotation")

anEmFAPD1 <- ggplot(umapEmbed,aes_string(x="umapX",y="umapY",color="PD1_faust_annotation"))+
        geom_point(size=ptSize)+
        theme_classic()+
        scale_color_manual(values=c("#ece7f2ff","#2b8cbeff","#a6bddbff"))+
        theme(
            legend.position="none",
            plot.title = element_text(hjust = 0.5),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
        )+
    xlab("annotation embedding X")+
    ylab("annotation embedding Y")+
    ggtitle("PD-1 FAUST annotation")

anEmFAHLADR <- ggplot(umapEmbed,aes_string(x="umapX",y="umapY",color="HLADR_faust_annotation"))+
        geom_point(size=ptSize)+
        theme_classic()+
        scale_color_manual(values=c("#ece7f2ff","#2b8cbeff"))+
        theme(
            legend.position="none",
            plot.title = element_text(hjust = 0.5),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
        )+
    xlab("annotation embedding X")+
    ylab("annotation embedding Y")+
    ggtitle("HLA-DR FAUST annotation")


#
#the associated faust annotations, expression embedding
#
exEmFACD4 <- ggplot(umapEmbed,aes_string(x="ex",y="ey",color="CD4_faust_annotation"))+
        geom_point(size=ptSize)+
        theme_classic()+
        scale_color_manual(values=c("#ece7f2ff","#2b8cbeff","#a6bddbff"))+
        theme(
            legend.position="none",
            plot.title = element_text(hjust = 0.5),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
        )+
    xlab("expression embedding X")+
    ylab("expression embedding Y")+
    ggtitle("CD4 FAUST annotation")


exEmFACD3 <- ggplot(umapEmbed,aes_string(x="ex",y="ey",color="CD3_faust_annotation"))+
        geom_point(size=ptSize)+
        theme_classic()+
        scale_color_manual(values=c("#ece7f2ff","#2b8cbeff"))+
        theme(
            legend.position="none",
            plot.title = element_text(hjust = 0.5),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
        )+
    xlab("expression embedding X")+
    ylab("expression embedding Y")+
    ggtitle("CD3 FAUST annotation")

exEmFAPD1 <- ggplot(umapEmbed,aes_string(x="ex",y="ey",color="PD1_faust_annotation"))+
        geom_point(size=ptSize)+
        theme_classic()+
        scale_color_manual(values=c("#ece7f2ff","#2b8cbeff","#a6bddbff"))+
        theme(
            legend.position="none",
            plot.title = element_text(hjust = 0.5),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
        )+
    xlab("expression embedding X")+
    ylab("expression embedding Y")+
    ggtitle("PD-1 FAUST annotation")

exEmFAHLADR <- ggplot(umapEmbed,aes_string(x="ex",y="ey",color="HLADR_faust_annotation"))+
        geom_point(size=ptSize)+
        theme_classic()+
        scale_color_manual(values=c("#ece7f2ff","#2b8cbeff"))+
        theme(
            legend.position="none",
            plot.title = element_text(hjust = 0.5),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
        )+
    xlab("expression embedding X")+
    ylab("expression embedding Y")+
    ggtitle("HLA-DR FAUST annotation")

#
#show all FAUST clusters and the embedding
#
anEmAllFaust <- ggplot(umapEmbed[faustPlotIndex,],aes(x=umapX,y=umapY,color=faustLabels))+
    geom_jitter(size=ptSize)+
    theme_classic()+
    theme(
        legend.position='none',
        plot.title = element_text(hjust = 0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )+
    ggtitle("Selected FAUST phenotypes")+
    scale_color_manual(values=myFaustColors)+
    xlab("annotation embedding X")+
    ylab("annotation embedding Y")

anEmCorrelates <- ggplot(umapEmbed[corPlotLookup,],aes(x=umapX,y=umapY,color=correlates))+
    geom_jitter(size=ptSize)+
    theme_classic()+
    theme(
        legend.position='none',
        plot.title = element_text(hjust = 0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )+
    ggtitle("Significant FAUST phenotypes")+
    scale_color_manual(values=c("#ccccccFF","#d7191cFF","#abd9e9FF","#fdae61FF","#2c7bb6FF"),
                       guide=FALSE)+
    xlab("annotation embedding X")+
    ylab("annotation embedding Y")

anEmCorrelates + geom_point(size=1) + theme_classic() + theme(legend.position="bottom")

corSubsetDF <- umapEmbed[corPlotLookup,]
corxLookup <- intersect(which(corSubsetDF$umapX < 0),which(corSubsetDF$umapX > -10))
coryLookup <- intersect(which(corSubsetDF$umapY < 5),which(corSubsetDF$umapY > -5))
corPlotSubLookup <- intersect(corxLookup,coryLookup)

subAnEmCorrelates <- ggplot(corSubsetDF[corPlotSubLookup,],aes(x=umapX,y=umapY,color=correlates))+
    geom_jitter(size=ptSize)+
    theme_classic()+
    theme(
        legend.position='none',
        plot.title = element_text(hjust = 0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )+
    ggtitle("Significant FAUST phenotypes (magnified)")+
    scale_color_manual(values=c("#ccccccFF","#d7191cFF","#abd9e9FF","#fdae61FF","#2c7bb6FF"),
                       guide=FALSE)+
    xlab("annotation embedding X")+
    ylab("annotation embedding Y")

anCorPlotFinalOut <- subAnEmCorrelates + geom_point(size=1) +
    annotation_custom(
        ggplotGrob(anEmCorrelates + ggtitle("Entire embedding plot") + xlab("") + ylab("")),
        xmin = -4.5, xmax = -0.5, ymin = 0.5, ymax = 4.5
    )

#
#expression all faust, correlates
#
exEmAllFaust <- ggplot(umapEmbed[faustPlotIndex,],aes(x=ex,y=ey,color=faustLabels))+
    geom_jitter(size=ptSize)+
    theme_classic()+
    theme(
        legend.position='none',
        plot.title = element_text(hjust = 0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )+
    ggtitle("Selected FAUST phenotypes")+
    scale_color_manual(values=myFaustColors)+
    xlab("expression embedding X")+
    ylab("expression embedding Y")

exEmCorrelates <- ggplot(umapEmbed[corPlotLookup,],aes(x=ex,y=ey,color=correlates))+
    geom_jitter(size=ptSize)+
    theme_classic()+
    theme(
        legend.position='none',
        plot.title = element_text(hjust = 0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )+
    ggtitle("Significant FAUST phenotypes")+
    scale_color_manual(values=c("#ccccccFF","#d7191cFF","#abd9e9FF","#fdae61FF","#2c7bb6FF"),
                       guide=FALSE)+
    xlab("expression embedding X")+
    ylab("expression embedding Y")


#
#combine and report
#
exprsGrid1 <- plot_grid(
    exEmCD4,exEmCD3,exEmPD1,exEmHLADR,
    nrow=1)

exprsGrid2 <- plot_grid(
    anEmCD4,anEmCD3,anEmPD1,anEmHLADR,
    nrow=1)

exprsTitle <- ggdraw() +
  draw_label(
    "Embeddings colored by standardized expression",
    fontface = 'bold',
    hjust = 0.5
  ) 

exprsOut1 <- plot_grid(
    exprsTitle,exprsGrid1, tLegend,
    ncol = 1,
    rel_heights = c(0.1, 1,0.05)
)

exprsOut2 <- plot_grid(
    exprsTitle,exprsGrid2, tLegend,
    ncol = 1,
    rel_heights = c(0.1, 1,0.05)
)


annGrid1 <- plot_grid(
    exEmFACD4,exEmFACD3,exEmFAPD1,exEmFAHLADR,
    nrow=1)

annGrid2 <- plot_grid(
    anEmFACD4,anEmFACD3,anEmFAPD1,anEmFAHLADR,
    nrow=1)

annTitle <- ggdraw() +
  draw_label(
    "Embeddings colored by FAUST annotations",
    fontface = 'bold',
    hjust = 0.5
  ) 

annOut1 <- plot_grid(
    annTitle, annGrid1, bLegend,
    ncol = 1,
    rel_heights = c(0.1, 1,0.05)
)

annOut2 <- plot_grid(
    annTitle, annGrid2, bLegend,
    ncol = 1,
    rel_heights = c(0.1, 1,0.05)
)


corGrid1 <- plot_grid(exEmAllFaust,exEmCorrelates,labels=c("C","D"),nrow=1,label_size=20)

corGrid2 <- plot_grid(anEmAllFaust,anCorPlotFinalOut,labels=c("G","H"),nrow=1,label_size=20)

plots <- align_plots(exprsGrid1,annGrid1,corGrid1,
                     exprsGrid2,annGrid2,corGrid2,
                     align='v',axis='l')                    

pOutPrep <- plot_grid(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]],ncol=1,labels=c("A","B","","E","F",""),
                      rel_heights=c(0.5,0.5,1,0.5,0.5,1),
                      label_size=18)
legendRow <- plot_grid(tLegend,bLegend,nrow=2)
pOut <- plot_grid(pOutPrep,
                  legendRow,
                  ncol = 1,
                  rel_heights = c(1,0.05))

cowplot::save_plot(
             filename="./figure_embedding.png",
             plot=pOut,
             base_width=12,
             base_height=18
         )
