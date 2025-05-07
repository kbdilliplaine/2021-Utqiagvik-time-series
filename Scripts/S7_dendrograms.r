require(reshape)
require(vegan)
require(indicspecies)
require(cluster)
require(plyr)
library(tidyverse)
library(readr)
library(ggthemes)
library(ggalt)
library(patchwork)
library(dendextend)
library(factoextra)

######################################
# Import Data and modify as needed
######################################
#These data are rarified
MO_18S_nometazoa=readRDS(file='./MO_18S_wometazoa.RDS')
#Change rownames to be just the DOY
rownames(MO_18S_nometazoa) <- sub("^MO_", "", rownames(MO_18S_nometazoa))
MO_18S_nometazoa=MO_18S_nometazoa[,c(1:740)]%>% as.matrix()


taxa_18S_nometa=readRDS(file='./MO_18S_Merged_MO_rare_nometa_tax_table.RDS')

################################################################################
#Make a cluster dendrogram using all 21 samples
################################################################################

set.seed()
MO_18S_Full_dist=vegdist(MO_18S_nometazoa, method = "bray")

Full_clust <- hclust(MO_18S_Full_dist, "ward.D2")

fviz_nbclust(MO_18S_nometazoa, pam, method = "silhouette", k.max=5)+ theme_classic()
plot(silhouette(cutree(Full_clust,2), MO_18S_Full_dist))

dend <- as.dendrogram(Full_clust)

# Color branches by cluster and add rectangles around them
dend <- color_branches(dend, k = 2)    # Color branches for clusters
dend <- set(dend, "branches_lwd", 2)   # Set branch line width (increase value to make it thicker)
dend <- set(dend, "labels_cex", 1)

dend_ggplot <- as.ggdend(dend)

Euks=
ggplot(dend_ggplot) +
    geom_rect(aes(xmin = 0, xmax = 13.5, ymin = -Inf, ymax = 1.05),
              color = "black", fill = NA, linetype = "solid", linewidth = 1) +
    geom_rect(aes(xmin = 13.5, xmax = Inf, ymin = -Inf, ymax = 1.05),
              color = "black", fill = NA, linetype = "solid", linewidth = 1) +
    annotate("text", x = 5, y = 1.1, label = "E1", size = 5, fontface = "bold") +
    annotate("text", x = 17, y = 1.1, label = "E2", size = 5, fontface = "bold") +
    theme_few() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank()
    ) +
    coord_cartesian(ylim = c(-0.2, NA))  # Adjust ylim to extend plot area



MO_16S_Merged_MO_rare_otu=readRDS(file='./MO_16S_Merged_MO_rare_otu.RDS')
rownames(MO_16S_Merged_MO_rare_otu) <- sub("^MO_", "", rownames(MO_16S_Merged_MO_rare_otu))

taxa_16S=readRDS(file='./MO_16S_Merged_MO_rare_tax_table.RDS')
################################################################################
#Make a cluster dendrogram using all 21 samples
################################################################################


set.seed()
MO_16S_Full_dist=vegdist(MO_16S_Merged_MO_rare_otu, method = "bray")

Full_clust_16S <- hclust(MO_16S_Full_dist, "ward.D2")
plot(Full_clust_16S)
rect.hclust(Full_clust_16S, k = 3)

dend <- as.dendrogram(Full_clust_16S)

# Color branches by cluster and add rectangles around them
dend <- color_branches(dend, k = 3)    # Color branches for clusters
dend <- set(dend, "branches_lwd", 2)   # Set branch line width (increase value to make it thicker)
dend <- set(dend, "labels_cex", 1)

dend_ggplot <- as.ggdend(dend)

Proks=
ggplot(dend_ggplot) +
    geom_rect(aes(xmin = 0, xmax = 6.5, ymin = -Inf, ymax = 1.1),
              color = "black", fill = NA, linetype = "solid", linewidth = 1) +
    geom_rect(aes(xmin = 6.5, xmax = 16.5, ymin = -Inf, ymax = 1.1),
              color = "black", fill = NA, linetype = "solid", linewidth = 1) +
    geom_rect(aes(xmin = 16.5, xmax = Inf, ymin = -Inf, ymax = 1.1),
              color = "black", fill = NA, linetype = "solid", linewidth = 1) +
    annotate("text", x = 3, y = 1.2, label = "P3", size = 5, fontface = "bold") +
    annotate("text", x = 11, y = 1.2, label = "P2", size = 5, fontface = "bold") +
    annotate("text", x = 19, y = 1.2, label = "P1", size = 5, fontface = "bold") +
    theme_few() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank()
    ) +
    coord_cartesian(ylim = c(-0.2, NA))  # Adjust ylim to extend plot area

DENDROS=Proks+Euks+
  plot_layout(widths = c(1, 1),
  design="1
          2")+
  plot_annotation(tag_levels = 'a')
