
library(ggthemes)
library(tidyverse) # for plotting and wrangling data
library("qiime2R")
library("phyloseq")
require(dplyr)
library(ggpubr)
library(cowplot)

MO_16S_Merged_MO_unrarified=readRDS("./MO_16S_Merged_MO_unrarified.rds")

sample_names(MO_16S_Merged_MO_unrarified) <- sub("^MO_", "", sample_names(MO_16S_Merged_MO_unrarified))
samp_order=rev(sample_names(MO_16S_Merged_MO_unrarified))

MO_16S_Merged_MO_unrarified_rel  = transform_sample_counts(MO_16S_Merged_MO_unrarified, function(x) x / sum(x) )

Oleispira = subset_taxa(MO_16S_Merged_MO_unrarified_rel, Genus %in% "Oleispira")
Oleispira_glom = tax_glom(Oleispira, "Genus")
# plot_bar(Oleispira_glom, "DOY", "Abundance", "Class", title="Oleispira")

Olei_taxa=Oleispira_glom %>% tax_table() %>% as.data.frame() %>% rownames_to_column( var="ASV")
Olei_table=Oleispira_glom %>% otu_table() %>%  t() %>% as.data.frame()%>%rownames_to_column( var="ASV")
colSums(Olei_table[,2:22]) #Check that it worked


Olei_Rel_abun=merge(Olei_taxa,Olei_table)

Olei_Rel_abun_long=  Olei_Rel_abun %>%
pivot_longer( cols=-c(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species),
names_to = "DOY", values_to = "Rel_abun") %>%
mutate(DOY=as.numeric(DOY), Rel_abun=as.numeric(Rel_abun)*100)


Oleispira_plot=
ggplot(data=Olei_Rel_abun_long) +
    annotate("rect", xmin = 110, xmax = 144, ymin = -Inf, ymax = Inf,
             alpha = .1,fill = "lightgreen")+
    geom_vline(xintercept=c(144), linetype=2, size=0.75)+
    geom_line(aes(x = DOY, y = Rel_abun), size = 0.75)+
    theme_bw()+
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, 100)
    ) +
    theme(
        axis.title.x= element_text(size=10, face="bold"),
        axis.title.y= element_text(size=10, face="bold"),
      axis.text = element_text(face="bold", size=8, colour="black"),
      # axis.text.x = element_blank(),
      axis.ticks=element_line(size=0.5, color="black"),
      axis.ticks.length = unit(0.25, "cm"),
      legend.position = c(0.38, 0.7),
      legend.background=element_rect(linetype="solid", colour = 'black'),
      legend.title=element_text(face="bold", size=8),
      legend.title.align=0.5,
      legend.text=element_text( size=8),
      legend.text.align=0,
      legend.key.height = unit(0.5, "line"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(
          color="white", fill="white",linetype="solid"),
      strip.text.x = element_text(
          size = 13, color = "black", face = "bold.italic"),
      panel.background=element_rect(size=1.5, color="black"),
      legend.key = element_blank(),
      plot.tag = element_text(size=8, face="bold"))+
# geom_hline(yintercept=0.05, linetype=2, size=0.75)+
scale_x_continuous(expand = c(0, 0), minor_breaks= seq(110, 162, 5))+
labs(tag = "a", x = "Day of year", y=expression(bolditalic(Oleispira)~~bold("relative abundance (%)")), color="Sequence variants")+
         theme(plot.margin=unit(c(0.1,0.5,0.2,0.2),"cm"))


Chl_avg=readRDS("./Full_Chlorophyll_avgd.RDS")
Chl_mgm2= Chl_avg%>%
filter(Location=="TM" & Chl_Flag != 0)

Chl_vol= Chl_avg%>%
    filter(Location=="TM" & !(Chl_Flag %in% c(0,7,8)))
Chl_brine= Chl_avg%>%
    filter(Location=="TM" & !(Chl_Flag %in% c(0,7,8)))%>% na.omit()

chl_oleispira=merge(Chl_mgm2, Olei_Rel_abun_long)

fit1=lm(data=chl_oleispira, Rel_abun~chl_mgm2)
summary(fit1)


inlay_plot2=
    ggplot(data = chl_oleispira, aes(x = chl_mgm2, y = Rel_abun)) +
    geom_point(size=0.5) +
    # geom_smooth(method = "lm", se = FALSE, color = "blue") +
    theme_classic() +
    stat_smooth(method = "lm", se = FALSE, color = "black", size=0.5) +
    # scale_x_continuous(expand = c(0, 0), minor_breaks= seq(110, 162, 5))+
    labs(tag="b",y = "Rel. abun. (%)", x = expression(bold("Chlorophyll") ~bolditalic("a")  ~bold("(mg m"^"-2"*")"))) +
    #stat_regline_equation(aes(label = paste(after_stat(adj.rr.label), sep = "~~~"))) +
    stat_cor(aes(label = ..rr.label..), size=2)+
    stat_cor(aes(label = ..p.label..), p.accuracy=0.001, vjust = 3.2, size=2)+
    theme(
        axis.title.x = element_text(size = 7, face = "bold"),
        axis.title.y = element_text(size = 7, face = "bold"),
        axis.text = element_text(size = 4.5, color = "black"),
        axis.ticks = element_line(size = 0.5, color = "black"),
        plot.margin = unit(c(0, 0.1, 0.1, 0.1), "cm"),
        plot.background = element_rect(colour = "black", fill="white", size=1),
        plot.tag = element_text(size=7, face="bold")
    )



final_plot <- ggdraw() +
    draw_plot(Oleispira_plot) +
    draw_plot(inlay_plot2, x = 0.52, y = 0.61, width = 0.4, height = 0.33)

print(final_plot)
