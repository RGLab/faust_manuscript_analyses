library(flowWorkspace)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(viridis)
#
#Read in FAUST Results 
#
countDF <- as.data.frame(readRDS(file.path(normalizePath("../.."),"faustRuns","CITN-07","faustData","faustCountMatrix.rds")))
cellPops <- colnames(countDF)
dcPops <- Reduce(intersect,list(
                               cd14=cellPops[grepl("CD14-",cellPops)],
                               cd3=cellPops[grepl("CD3-",cellPops)],
                               cd19=cellPops[grepl("CD19-",cellPops)],
                               cd16=cellPops[grepl("CD16-",cellPops)],
                               cd56=cellPops[grepl("CD56-",cellPops)],
                               dr=cellPops[!grepl("HLA DR-",cellPops)]
                           ))

dc_count_df <- data.frame(
    Count=as.numeric(apply(countDF[,dcPops],1,sum)),
    Total=as.numeric(apply(countDF,1,sum)),
    merge_name=rownames(countDF),
    stringsAsFactors=FALSE
)
#
#Read in meta data
#
gs <- load_gs(file.path(normalizePath("../.."),
                        "dataSets",
                        "publication_gating_sets",
                        "CITN-07",
                        "citn07_longitudinal_gs"))

meta_data <- as.data.frame(pData(gs))
meta_data$`merge_name` <- rownames(meta_data)
#
#combine and derive cytometry visit
#
plot_data_wide <- inner_join(dc_count_df,meta_data,by=c("merge_name"))
cytometry_visit <- gsub(".fcs_[[:digit:]]+","",plot_data_wide$merge_name)
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
plot_data_wide$`cytometry_visit` <- as.factor(cytometry_visit)
plot_data <-plot_data_wide[,c("Publication_ID","cytometry_visit","Cohort","Count","Total")]
chtLabs <- c("Cohort 1", "Cohort 2")
names(chtLabs) <- c("CH1", "CH2")
#
#draw the plot and save it
#
p1 <- ggplot(plot_data,aes(y = Count/Total, x = cytometry_visit)) +
    geom_jitter(aes(size = Count), width = 0.2, height = 0, alpha = 0.1) + 
    geom_line(aes(group = `Publication_ID`), alpha = 0.1) +
    facet_grid(Cohort~.,labeller = labeller(Cohort=chtLabs)) +
    scale_y_continuous("Total DC compartment as % of CD45+ cells", labels = scales::percent) +
    scale_x_discrete("Visit") +
    theme_classic(base_size = 6) +
    scale_color_viridis(
        begin = 0.9,
        end = 0.4,
        option = "B",
        discrete = TRUE,
        guide = FALSE
    ) +
    scale_size(range = c(1,4)) +
    theme(legend.position = "bottom", axis.text.x  = element_text(angle=45, hjust = 1))+
    geom_line(data = (plot_data %>% dplyr::group_by(cytometry_visit, Cohort) %>% dplyr::summarize(proportion = median(Count/Total), uci = proportion + qnorm(0.975)*sd(Count/Total), lci = proportion + qnorm(0.025)*sd(Count/Total)) %>% ungroup), 
              aes(x = cytometry_visit, y = proportion, group = Cohort), size = 1) +
    geom_errorbar(data = (plot_data %>%  dplyr::group_by(cytometry_visit, Cohort) %>% dplyr::summarize(proportion = median(Count/Total), uci = proportion - mad(Count/Total), lci = proportion + mad(Count/Total)) %>% ungroup),
                  aes(y = proportion, ymin = lci, ymax = uci), width = 0.2)

pOut <- p1+theme(legend.position="none")
saveRDS(pOut,file.path(normalizePath("."),"artifacts","longitudinal_figure_panel_b.rds"))






