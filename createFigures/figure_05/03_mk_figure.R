#
#need to run this function on a computer that has alpha enabled to properly shade the lines from CITN-07
#
library(ggplot2)
library(cowplot)
library(Cairo)
p1 <- readRDS(file.path(normalizePath("."),"artifacts","longitudinal_figure_panel_a.rds"))+
    theme_classic(base_size = 8)+
    theme(
        legend.position="none",
        plot.title=element_text(size=10,hjust=0.5)
    )+
    ggtitle("MCC anti-PD-1 trial")

p2 <- readRDS(file.path(normalizePath("."),"artifacts","longitudinal_figure_panel_b.rds"))+
    theme_classic(base_size = 8)+
    theme(
        legend.position="none",
        plot.title=element_text(size=10,hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )+
    ggtitle("FLT3-ligand + Therapeutic Vx trial")
p2a <- ggdraw(p2)+
    draw_label("Cohort stimulated with\nFLT3-ligand",x=0.75,y=0.8,size=7) +
    draw_label("Cohort not stimulated\nby FLT3-ligand",x=0.75,y=0.4,size=7)

pOut <- plot_grid(
    p1,p2a,
    nrow=1,
    ncol=2,
    labels=c("A","B"),
    label_size=12
)


ggsave(
    filename="./figure_longitudinal.png",
    plot=pOut,
    width=8,
    height=4,
    units="in",
    dpi=300,
    type="cairo"
)






