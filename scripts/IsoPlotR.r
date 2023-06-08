library(tibble)
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)
library(lattice)
library(cowplot)
library(ggpubr)
library(purrr)

#   seqnames    start       end         width   strand  isoform_id          gene_id             gene_name
#   chr9        32551676    32551807    132     +       ENST00000450093     ENSG00000235453     NA
#   chr9        32552234    32552479    246     +       ENST00000450093     ENSG00000235453     NA
#   chr9        32552801    32553007    207     +       ENST00000450093     ENSG00000235453     NA
#   chr9        32551144    32551693    550     +       ENST00000453396     ENSG00000235453     NA
#   chr9        32552234    32552479    246     +       ENST00000453396     ENSG00000235453     NA
#   chr9        32552801    32553007    207     +       ENST00000453396     ENSG00000235453     NA
#   chr9        32552327    32552479    153     +       ENST00000692500     ENSG00000235453     NA
#   chr9        32552801    32553007    207     +       ENST00000692500     ENSG00000235453     NA

plot_iso <- function(
    df, 
    start = NULL,
    end = NULL,
    output,
    color=colours
    ) 
    {
    colours = c("#AF93B4", "#FFA9A3", "#B9E6FF", "#ACE4AA", "#FFD770", "#FAA00F", "#9893DA", "#FB9DAE", "#157A3F", "A44A3F")
    colnames(df) <-  c("seqnames", "start", "end", "width", "strand", "isoform_id", "gene_id", "gene_name")
    count = 1
    plots <- list()
    options(repr.plot.width=20, repr.plot.height=12)
    for (iso in unique(tab$isoform_id)) {
        df.iso <- df[df$isoform_id == iso,]
        df.iso <- df.iso[order(df.iso$start),] 
        chr_start = min(df.iso$start)
        chr_end = max(df.iso$end)
        chr = df.iso$seqnames[1]
        strand = df.iso$strand[1]
        gene_id = df.iso$gene_id[1]
        gene_name = df.iso$gene_name[1]
        exons_s = c(df.iso$start)
        exons_e = c(df.iso$end)
        exons_size = c(df.iso$width)
        if (strand == "+"){
            arrow = "> > >"
        }
        if (strand == "-"){
            arrow = "< < <"
        }
        exons_coord = sort(c(exons_s, exons_e), decreasing = FALSE)
        cords = c()
        for (i in seq(2, length(exons_coord)-1, by=2)){
            here = (exons_coord[i+1] + exons_coord[i]) / 2 
            cords <- c(cords, here)
        }
        df2 = data.frame(cords, y=rep(1.66, length(cords)))
        plot <- 
            df.iso %>%
            ggplot() +
            scale_x_discrete(limits=exons_coord) +
            # ligne du milieu
            geom_text(
                data=df2,
                aes(x = cords ,
                    y =y,
                    label = arrow),
                size = 6) +
            geom_hline(
                yintercept = 1.65,
                size = 1) +
            geom_rect(
                data = df.iso,
                mapping = aes(
                    xmin = start,
                    xmax = end,
                    ymin = 1.55,
                    ymax = 1.75),
                color = "black",
                fill = colours[count]) +
            geom_text(
                data=df.iso,
                aes(
                    x = (start+(end - start)/2) ,
                    y = 1.9,
                    label = c(1:length(isoform_id))),
                vjust="inward",hjust="inward",
                size = 4)+
            geom_text(
                data = df.iso,
                aes(
                    x = (start+(end - start)/2) ,
                    y = 1.65,
                    label = width),
                size = 4,
                color = "white",
                fontface = "bold") +
            # etiquettes noms exons
            xlab("") +
            ylab("") +
            theme_bw() +
            theme(
                axis.text.x = element_text(angle = 0, size=10),
                axis.text.y = element_blank(),
                axis.ticks = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.border = element_blank(),
                plot.title = element_text(face = "bold", hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5),
                plot.margin = unit(c(1,1,1,1),"cm"), # unit(c(top, right, bottom, left)
                aspect.ratio=1/30) +
            labs(title = iso, subtitle = chr)

    plots[[count]] <- plot
    count = count + 1
    }
    return(plots)
}
#tab = read.csv("/groups/dog/stage/enora/isoformswitch/exons_ENSG00000235453.csv")
#plot_iso(tab, output="/groups/dog/stage/enora/isoformswitch/iso")
