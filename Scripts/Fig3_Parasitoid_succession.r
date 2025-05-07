library("qiime2R")
library("readr")
library("phyloseq")
library(tidyverse)
require(reshape)
require(vegan)
require(indicspecies)
require(cluster)
require(plyr)
library(microViz)
library(patchwork)
library(cowplot)

trace("parse_taxonomy", edit=TRUE)

# function (taxonomy, tax_sep, trim_extra)
# {
#     if (missing(taxonomy)) {
#         stop("Taxonomy Table not supplied.")
#     }
#     if (missing(trim_extra)) {
#         trim_extra = TRUE
#     }
#     if (missing(tax_sep)) {
#         tax_sep = "; |;"
#     }
#     if (sum(colnames(taxonomy) %in% c("Feature.ID", "Taxon")) !=
#         2) {
#         stop("Table does not match expected format. ie does not have columns Feature.ID and Taxon.")
#     }
#     taxonomy <- taxonomy[, c("Feature.ID", "Taxon")]
#     if (trim_extra) {
#         taxonomy$Taxon <- gsub("[kpcofgs]__", "", taxonomy$Taxon)
#         taxonomy$Taxon <- gsub("D_\\d__", "", taxonomy$Taxon)
#     }
#     taxonomy <- suppressWarnings(taxonomy %>% separate(Taxon,
#         c("Domain", "Supergroup", "Division", "Subdivision",
#             "Class", "Order", "Family", "Genus", "Species"),
#         sep = tax_sep))
#     taxonomy <- apply(taxonomy, 2, function(x) if_else(x == "",
#         NA_character_, x))
#     taxonomy <- as.data.frame(taxonomy)
#     rownames(taxonomy) <- taxonomy$Feature.ID
#     taxonomy$Feature.ID <- NULL
#     return(taxonomy)
# }
#
# MO_18S_Merged_MO_rare=readRDS("G:/My Drive/01_MOSAIC/Data/Phyloseq_Objects/MO_18S_Merged_rare.rds")
MO_18S_Merged_MO_unrarified=readRDS("./MO_18S_Merged_MO_unrarified.rds")
sample_names(MO_18S_Merged_MO_unrarified) <- sub("^MO_", "", sample_names(MO_18S_Merged_MO_unrarified))
samp_order=rev(sample_names(MO_18S_Merged_MO_unrarified))
##################################################
#Parasitoid Succession
##################################################

MO_18S_Merged_MO_unrarified_rel  = transform_sample_counts(MO_18S_Merged_MO_unrarified, function(x) (x / sum(x))*100)

Chytrids = subset_taxa(MO_18S_Merged_MO_unrarified_rel, Class %in% "Chytridiomycota")
Chytrid_glom = tax_glom(Chytrids, "Class")
plot_bar(Chytrid_glom, "DOY", "Abundance", "Class", title="Chytridiomycota", facet_grid = "Class")


Cryomonads = subset_taxa(MO_18S_Merged_MO_unrarified_rel, Order %in% "Cryomonadida")
Cryomonad_glom = tax_glom(Cryomonads, "Family")
plot_bar(Cryomonad_glom, "DOY", "Abundance", "Class", title="Cryomonadida", facet_grid = "Family")


Peronosporomycetes = subset_taxa(MO_18S_Merged_MO_unrarified_rel, Class %in% "Peronosporomycetes")
Perono_glom = tax_glom(Peronosporomycetes, "Class")
plot_bar(Perono_glom, "DOY", "Abundance", "Class", title="Peronosporomycetes", facet_grid = "Class")

Syndiniales = subset_taxa(MO_18S_Merged_MO_unrarified_rel, Class %in% "Syndiniales")
Syndiniales_glom = tax_glom(Syndiniales, "Class")
plot_bar(Syndiniales_glom, "DOY", "Abundance", "Class", title="Syndiniales", facet_grid = "Class")

Laby = subset_taxa(MO_18S_Merged_MO_unrarified_rel, Order %in% "Labyrinthulomycetes")
Laby_glom = tax_glom(Laby, "Order")
plot_bar(Laby_glom, "DOY", "Abundance", "Class", title="Labyrinthulomycetes", facet_grid = "Class")


merged_phylo=merge_phyloseq(Chytrid_glom, Cryomonad_glom, Laby_glom, Perono_glom, Syndiniales_glom)
plot_bar(merged_phylo, "DOY", "Abundance", "Class", title="Syndiniales", facet_grid = "Class")

merged_phylo_fixed=merged_phylo %>%
 tax_fix(
  min_length = 4,
  unknowns = c("Cryomonadida_X", "Cryomonadida Order"),
  sep = " ", anon_unique = TRUE,
  suffix_rank = "classified"
 )

 merged_phylo_fixed_glom = tax_glom(merged_phylo_fixed, "Family")
 plot_bar(merged_phylo_fixed_glom, "DOY", "Abundance", "Order", title="merged_phylo_fixed_glom", facet_grid = "Family")


merged_phylo_otu=merged_phylo_fixed_glom %>% otu_table() %>% t() %>% as.data.frame() %>% rownames_to_column( var="ASV")
merged_phylo_taxa=merged_phylo_fixed_glom %>% tax_table() %>% as.data.frame() %>% rownames_to_column( var="ASV")

columns_to_mutate =c("ASV", "Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")


Parasitoid_Taxa_Merger_DF=merged_phylo_taxa %>% inner_join(merged_phylo_otu, by="ASV") %>%
    mutate_at(vars(columns_to_mutate), as.factor) %>%
    pivot_longer(cols = -c(ASV, Domain, Supergroup, Division, Subdivision, Class, Order, Family, Genus, Species), names_to = "DOY", values_to = "Count")%>%
    mutate(Class = fct_recode(Class,  "Labyrinthulea (Sagenista)"="Sagenista")) %>%
    mutate(Class = fct_recode(Class,  "Syndiniales (Dinophyceae)"="Syndiniales")) %>%
    mutate(Class = fct_recode(Class,  "Oomycetes"="Peronosporomycetes")) #%>%
    # mutate(Class = fct_recode(Class,  "Oomycetes Class"="Peronosporomycetes Class"))
Parasitoid_Taxa_Merger_DF$DOY=as.numeric(Parasitoid_Taxa_Merger_DF$DOY)


Parasitoid_Taxa_Merger_DF <- Parasitoid_Taxa_Merger_DF %>%
    dplyr::mutate(
        Class = as.character(Class),  # Convert to character for case_when
        Class = dplyr::case_when(
            Class == "Chytridiomycota" ~ "Chytridiomycetes",
            Class == "Filosa-Thecofilosea" ~ "Thecofilosea",
            Class == "Labyrinthulea (Sagenista)" ~ "Labyrinthulomycetes (Bigyra)",
            TRUE ~ Class  # Keep other values unchanged
        ),
        Class = factor(Class)  # Convert back to factor
    )



Chytridiomycota_panel_white=
ggplot(Parasitoid_Taxa_Merger_DF %>% filter(Class=="Chytridiomycetes"))+
    annotate("rect", xmin = 110, xmax = 144, ymin = -Inf, ymax = Inf,
             alpha = .1,fill = "lightgreen")+
    geom_vline(xintercept=c(144), linetype=2, size=0.5)+
    geom_line(aes(DOY, Count), size=0.75)+
    facet_wrap(~Class, ncol=1)+
    theme_bw()+
    theme(
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text = element_text(face="bold", size=8, colour="black"),
        axis.text.x = element_blank(),
        axis.ticks=element_line(size=0.5, color="black"),
        axis.ticks.length = unit(0.25, "cm"),
        legend.title=element_text(face="bold", size=10),
        legend.title.align=0.5,
        legend.text=element_text( size=8),
        legend.text.align=0,
        legend.key.width = unit(1, "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(
            color="white", fill="white",linetype="solid"),
        strip.text.x = element_text(
            size = 8, color = "black", face = "bold"),
        panel.background=element_rect(size=0.75, color="black"),
        legend.key = element_blank(),
    )+
    scale_y_continuous(expand = c(0, 0))+
    scale_x_continuous(expand = c(0, 0), minor_breaks= seq(110, 160, 5))+
    coord_cartesian(xlim = c(110, 160), ylim=c(-0.1, 30.1))+
annotate("text", x = 113, y = 27, label = "a", size = 5)



    Thecofilosea_panel_white=
    ggplot(Parasitoid_Taxa_Merger_DF %>% filter(Class=="Thecofilosea"))+
        annotate("rect", xmin = 110, xmax = 144, ymin = -Inf, ymax = Inf,
                 alpha = .1,fill = "lightgreen")+
        geom_vline(xintercept=c(144), linetype=2, size=0.5)+
        geom_line(aes(DOY, Count, group=Family, color=Family), size=0.75)+
        scale_color_manual(values = c("#13d4df", "#2e8eb5", "#1a3662"))+
        facet_wrap(~Class, ncol=1)+
        theme_bw()+
        theme(
            axis.title.y= element_blank(),
            axis.title.x= element_blank(),
            axis.text = element_text(face="bold", size=8, colour="black"),
            axis.text.x = element_blank(),
            axis.ticks=element_line(size=0.5, color="black"),
            axis.ticks.length = unit(0.25, "cm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(
                color="white", fill="white",linetype="solid"),
            strip.text.x = element_text(
                size = 8, color = "black", face = "bold"),
            panel.background=element_rect(size=0.75, color="black"),
            legend.position = c(0.27, 0.6),
            legend.background=element_blank(),
            legend.title=element_blank(),
            legend.title.align=0.5,
            legend.text=element_text( size=6),
            legend.text.align=0,
            legend.key.height = unit(0.5, "line"),
            legend.key = element_blank()
        )+
        scale_y_continuous(expand = c(0, 0))+
        scale_x_continuous(expand = c(0, 0), minor_breaks= seq(110, 160, 5))+
        coord_cartesian(xlim = c(110, 160), ylim=c(-0.1, 30.1))+
annotate("text", x = 113, y = 27, label = "b", size = 5)


        Labyrinthulea_panel_white=
        ggplot(Parasitoid_Taxa_Merger_DF %>% filter(Class=="Labyrinthulomycetes (Bigyra)"))+
            annotate("rect", xmin = 110, xmax = 144, ymin = -Inf, ymax = Inf,
                     alpha = .1,fill = "lightgreen")+
                     geom_vline(xintercept=c(144), linetype=2, size=0.5)+
                     geom_line(aes(DOY, Count), size=0.75)+
                     facet_wrap(~Class, ncol=1)+
                     theme_bw()+
                     theme(
                         axis.title.y= element_blank(),
                         axis.title.x= element_blank(),
                         axis.text = element_text(face="bold", size=8, colour="black"),
                         axis.text.x = element_blank(),
                         axis.ticks=element_line(size=0.5, color="black"),
                         axis.ticks.length = unit(0.25, "cm"),
                         legend.title=element_text(face="bold", size=10),
                         legend.title.align=0.5,
                         legend.text=element_text( size=8),
                         legend.text.align=0,
                         legend.key.width = unit(1, "line"),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         strip.background = element_rect(
                             color="white", fill="white",linetype="solid"),
                         strip.text.x = element_text(
                             size = 8, color = "black", face = "bold"),
                         panel.background=element_rect(size=0.75, color="black"),
                         legend.key = element_blank(),
                     )+
                     scale_y_continuous(expand = c(0, 0))+
                     scale_x_continuous(expand = c(0, 0), minor_breaks= seq(110, 160, 5))+
                     coord_cartesian(xlim = c(110, 160), ylim=c(-0.1, 30.1))+
                 # geom_text(aes(label = "a"), x=111, y=27, size=7)
                 annotate("text", x = 113, y = 27, label = "c", size = 5)

      Oomycetes_panel_white=
      ggplot(Parasitoid_Taxa_Merger_DF %>% filter(Class=="Oomycetes"))+
            annotate("rect", xmin = 110, xmax = 144, ymin = -Inf, ymax = Inf,
                     alpha = .1,fill = "lightgreen")+
                     geom_vline(xintercept=c(144), linetype=2, size=0.5)+
                     geom_line(aes(DOY, Count), size=0.75)+
                     facet_wrap(~Class, ncol=1)+
                     theme_bw()+
                     theme(
                         axis.title.y= element_blank(),
                         axis.title.x= element_blank(),
                         axis.text = element_text(face="bold", size=8, colour="black"),
                         axis.text.x = element_blank(),
                         axis.ticks=element_line(size=0.5, color="black"),
                         axis.ticks.length = unit(0.25, "cm"),
                         legend.title=element_text(face="bold", size=10),
                         legend.title.align=0.5,
                         legend.text=element_text( size=8),
                         legend.text.align=0,
                         legend.key.width = unit(1, "line"),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         strip.background = element_rect(
                             color="white", fill="white",linetype="solid"),
                         strip.text.x = element_text(
                             size = 8, color = "black", face = "bold"),
                         panel.background=element_rect(size=0.75, color="black"),
                         legend.key = element_blank(),
                     )+
                     scale_y_continuous(expand = c(0, 0))+
                     scale_x_continuous(expand = c(0, 0), minor_breaks= seq(110, 160, 5))+
                     coord_cartesian(xlim = c(110, 160), ylim=c(-0.1, 30.1))+
                 annotate("text", x = 113, y = 27, label = "d", size = 5)


    Syndiniales_panel_white=
    ggplot(Parasitoid_Taxa_Merger_DF %>% filter(Class=="Syndiniales (Dinophyceae)"))+
          annotate("rect", xmin = 110, xmax = 144, ymin = -Inf, ymax = Inf,
                   alpha = .1,fill = "lightgreen")+
                   geom_vline(xintercept=c(144), linetype=2, size=0.5)+
                   geom_line(aes(DOY, Count), size=0.75)+
                   facet_wrap(~Class, ncol=1)+
                   theme_bw()+
                   theme(
                       axis.title.y= element_blank(),
                       axis.title.x= element_text(size=10, face="bold"),
                       axis.text = element_text(face="bold", size=8, colour="black"),
                       axis.ticks=element_line(size=0.5, color="black"),
                       axis.ticks.length = unit(0.25, "cm"),
                       legend.title=element_text(face="bold", size=10),
                       legend.title.align=0.5,
                       legend.text=element_text( size=8),
                       legend.text.align=0,
                       legend.key.width = unit(1, "line"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       strip.background = element_rect(
                           color="white", fill="white",linetype="solid"),
                       strip.text.x = element_text(
                           size = 8, color = "black", face = "bold"),
                       panel.background=element_rect(size=0.75, color="black"),
                       legend.key = element_blank(),
                   )+
                   scale_y_continuous(expand = c(0, 0))+
                   scale_x_continuous(expand = c(0, 0), minor_breaks= seq(110, 160, 5))+
                   coord_cartesian(xlim = c(110, 160), ylim=c(-0.1, 30.1))+
                   labs(x="Day of year")+
               annotate("text", x = 113, y = 27, label = "e", size = 5)

combined_plot <- Chytridiomycota_panel_white /
  Thecofilosea_panel_white/
        Labyrinthulea_panel_white   / Oomycetes_panel_white /
        Syndiniales_panel_white
y_axis_title <- ggdraw() +
    draw_label("Relative abundance (%)", size = 10, angle = 90, vjust = 1, fontface ="bold")

Para_succ_panel_white=y_axis_title + combined_plot+ plot_layout(widths = c(1,14))&
    theme(plot.margin = margin(t = 1, r = 4, b = 1, l = 2))
