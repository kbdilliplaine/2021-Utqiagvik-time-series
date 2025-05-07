
#Prep Qiime2 data for R
#https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121
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
library(forcats)
library(paletteer)
library(readr)
library(ggthemes)
library(ggalt)

MO_16S_Merged_MO_unrarified=readRDS("./MO_16S_Merged_MO_unrarified.rds")
sample_names(MO_16S_Merged_MO_unrarified) <- sub("^MO_", "", sample_names(MO_16S_Merged_MO_unrarified))
samp_order=rev(sample_names(MO_16S_Merged_MO_unrarified))

###############################################################################
#Area plots
###############################################################################
# Step 1: Transform to relative abundance at Genus level
Genus_Rel_16 = MO_16S_Merged_MO_unrarified %>%
    tax_transform("compositional", rank = "Genus")

# Step 2: Extract and reshape OTU table
otu_table_16 = Genus_Rel_16 %>%
    otu_table() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "Genus")

# Step 3: Calculate total abundance for each Genus
taxa_sums_16 = rowSums(otu_table_16[, 2:22])  # assuming samples are in columns 2 to 22
taxa_counts_16 = data.frame(Genus = otu_table_16$Genus, Total_Count = taxa_sums_16)

# Step 4: Identify top 10 genera
top_n_taxa_16 = slice_max(taxa_counts_16, n = 10, order_by = Total_Count) %>%
    pull(Genus)

# Step 5: Melt the full data (no pruning yet)
abundance_melted_16 = melt(otu_table_16, id.vars = "Genus") %>%
    dplyr::rename(DOY = variable) %>%
    mutate(
        DOY = as.numeric(as.character(DOY)),
        value = value * 100,  # Scale to percentages
        Genus = ifelse(Genus %in% top_n_taxa_16, Genus, "Other")  # Lump non-top10 into "Other"
    )

# Step 6: Group by DOY and Genus and sum values
abundance_melted_16 = abundance_melted_16 %>%
    group_by(DOY, Genus) %>%
    dplyr::summarize(value = sum(value), .groups = "drop")

# Step 7: Reorder factor levels
abundance_melted_16$Genus = fct_relevel(
    abundance_melted_16$Genus,
    rev(top_n_taxa_16),  # Top 10 in order
    "Other"           # "Other" last
)

# Step 8: Define color palette
colors_for_plot = c(
    paletteer::paletteer_d("colorBlindness::paletteMartin")[1:10],  # 10 colors
    "gray"  # Gray for "Other"
)

    Area_16S_Genus=
    ggplot(abundance_melted_16, aes(x = DOY, y = value, fill = Genus)) +
        geom_area(linetype = 1, linewidth =0.25 ,colour="black" )+
        # scale_fill_paletteer_d("colorBlindness::paletteMartin")+
         scale_fill_manual(values = colors_for_plot) +
        geom_vline(xintercept=c(144), linetype=2)+
        scale_x_continuous(expand=c(0,0))+
        scale_y_continuous(expand=c(0,0))+
        theme(
            plot.title = element_text(hjust = 0.5, face="bold", size=10),
            strip.text.x = element_text(size=8, face="bold", hjust=0),
            #strip.text.x = element_blank(),
            strip.text.y = element_text(size=8, face="bold.italic"),
            legend.text = element_text( size=8),
            legend.title = element_text( size=10, face="bold"),
            legend.key.size = unit(0.4, 'cm'),
            axis.title = element_text(size=10, face="bold"),
            axis.title.x = element_text(vjust=-1),
            axis.text.x = element_text(size=8, vjust = -1),
            axis.text = element_text(size=8, face="bold",  colour = "black"),
            axis.ticks =element_line(color = "black"),
            panel.spacing = unit(2.5, "lines"),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            # plot.margin = margin(8,18,0,2)
        )+
        labs(tag="c", x = "Day of year", y = "% Relative abundance", fill = "Genus", title="Most abundant prokaryotic taxa (16S rRNA)")


#################################
#16S nMDS
#################################
MO_16S_Merged_MO_rare_otu=readRDS(file='./MO_16S_Merged_MO_rare_otu.RDS')
rownames(MO_16S_Merged_MO_rare_otu) <- sub("^MO_", "", rownames(MO_16S_Merged_MO_rare_otu))
taxa_16S=readRDS(file='./MO_16S_Merged_MO_rare_tax_table.RDS')
Enviro_Chl=readRDS(file='./Enviro_Chl.RDS')


#############################
# Minimum code required to generate 16S nMDS
#############################

set.seed()
MO_16S_Full_dist=vegdist(MO_16S_Merged_MO_rare_otu, method = "bray")

Full_clust <- hclust(MO_16S_Full_dist, "ward.D2")

clust_2<-cutree(Full_clust, k = 2)
clust_3<-cutree(Full_clust, k = 3)

#Square root transformation
# Wisconsin double standardization
MO_16S_Wisc <- wisconsin(MO_16S_Merged_MO_rare_otu)
MDS_Wisc=metaMDS(MO_16S_Wisc, distance='bray', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
MDS_Wisc_3=metaMDS(MO_16S_Wisc, distance='bray', k=3, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)

clust_2_df=as.data.frame(clust_2)
clust_3_df=as.data.frame(clust_3)

clust_2_df <- tibble::rownames_to_column(clust_2_df, var = "DOY")
clust_2_df$DOY <- as.numeric(clust_2_df$DOY)
clust_3_df <- tibble::rownames_to_column(clust_3_df, var = "DOY")
clust_3_df$DOY <- as.numeric(clust_3_df$DOY)

Enviro_Chl_clust <- Enviro_Chl %>%
  full_join(clust_2_df, by = "DOY") %>%
  full_join(clust_3_df, by = "DOY") %>%
  mutate(
    clust_2 = as.factor(clust_2),
    clust_3 = as.factor(clust_3)
  )


point_locations_MDS_Wisc_3 <-
scores(MDS_Wisc_3, display = "sites") %>%
as.data.frame() %>%
rownames_to_column("DOY") %>%
mutate(DOY=as.numeric(DOY))%>%
left_join(Enviro_Chl_clust, by="DOY")

breaks=c(83, 127, 143, 170)
labels=c("Early", "Late", "Post-")
point_locations_MDS_Wisc_3 <- point_locations_MDS_Wisc_3 %>%
    mutate(Phase = cut(DOY, breaks=breaks, labels=labels))

breaks2=c(83, 143, 170)
labels2=c("Bloom", "Post-")
point_locations_MDS_Wisc_3 <- point_locations_MDS_Wisc_3 %>%
    mutate(Phase_2 = cut(DOY, breaks=breaks2, labels=labels2))

point_locations_MDS_Wisc_3 <- point_locations_MDS_Wisc_3 %>%
  mutate(Period = ifelse(DOY > 132, "Melt", "Growth"))


fill_palette_16S <- c("lightgreen", "lightblue", "black")
# point_palette_16S <- c("darkgreen", "blue", "black")


centroids <- point_locations_MDS_Wisc_3 %>%
    group_by(clust_3) %>%
    dplyr::summarize(
        centroid_NMDS1 = mean(NMDS1, na.rm = TRUE),
        centroid_NMDS2 = mean(NMDS2, na.rm = TRUE)
    ) %>%
    mutate(groups = case_when(
        clust_3 == 1 ~ "Early-",
        clust_3 == 2 ~ "Mid-",
        clust_3 == 3 ~ "Post-",
        TRUE ~ NA_character_
    ))%>%
    mutate(
        centroid_NMDS2 = case_when(
            clust_3 == "3" ~ 0.25,  # Adjust value for clust_2 == 1
            TRUE ~ centroid_NMDS2     # Retain existing values otherwise
        ))

  centroids2 <- point_locations_MDS_Wisc_3 %>%
      group_by(clust_3) %>%
      dplyr::summarize(
          centroid_NMDS1 = mean(NMDS1, na.rm = TRUE),
          centroid_NMDS2 = mean(NMDS2, na.rm = TRUE)
      ) %>%
      mutate(groups = case_when(
          clust_3 == 1 ~ "P1",
          clust_3 == 2 ~ "P2",
          clust_3 == 3 ~ "P3",
          TRUE ~ NA_character_
      ))%>%
      mutate(
          centroid_NMDS2 = case_when(
              clust_3 == "3" ~ 0.25,  # Adjust value for clust_2 == 1
              TRUE ~ centroid_NMDS2     # Retain existing values otherwise
          ))


#Poster sizing
  New_MDS_Wisc_16S_clust_3=
  ggplot(data = point_locations_MDS_Wisc_3,
         aes(x =NMDS1, y =  NMDS2, group = clust_3,
             colour = DOY,
             fill = clust_3,
             label = clust_3,
             # shape = Phase_2
         ))+
      geom_encircle( fill=NA, s_shape = 1, expand = 0,
                     alpha = 1, color = "black", size=1, show.legend = FALSE)+
      geom_encircle(aes(fill = clust_3), s_shape = 1, expand = 0,
                    alpha = 0.2, color = "black", size=1, show.legend = FALSE)+
      geom_point(size = 2)+
      scale_colour_gradient(low = "lightblue", high = "darkblue", na.value = NA,
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", frame.linewidth =0.75))+
      annotate("text", x = 0.8, y = -0.7,
               label = paste0("Stress = ",
                              round(MDS_Wisc_3[['stress']], digits = 3)
                            ),
                          size=3)+
                          geom_text(data = centroids2, aes(x = centroid_NMDS1, y = centroid_NMDS2, label = groups),
              size = 4, fontface = "bold", color = "black") +
      theme_few()+
      theme(
        plot.title = element_text(hjust = 0.5, face="bold", size=10),
        strip.text.x = element_text(size=8, face="bold", hjust=0),
        #strip.text.x = element_blank(),
        strip.text.y = element_text(size=8, face="bold.italic"),
        legend.text = element_text( size=8),
        legend.title = element_text( size=10, face="bold"),
        legend.key.size = unit(0.4, 'cm'),
        axis.title = element_text(size=10, face="bold"),
        axis.title.x = element_text(vjust=-1),
        axis.text.x = element_text(size=8, vjust = -1),
        axis.text = element_text(size=8, face="bold",  colour = "black"),
        axis.ticks =element_line(color = "black"),
        panel.spacing = unit(2.5, "lines"),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.margin = margin(10,10,10,10)
      )+
      labs(tag="a", x = "nMDS 1", y = "nMDS 2", color = "Day of year",
       title="Prokaryotic nMDS (16S rRNA)")+
        scale_fill_manual(values=fill_palette_16S)+
        guides(fill="none")

############################
#18S
############################
library(openxlsx)

#NOTE: Run trace and replace the code with this to allow for the PR2 taxonomy.

trace("parse_taxonomy", edit=TRUE)
#
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



MO_18S_Merged_MO_rare=readRDS("./MO_18S_Merged_rare.rds")
sample_names(MO_18S_Merged_MO_rare) <- sub("^MO_", "", sample_names(MO_18S_Merged_MO_rare))
samp_order=rev(sample_names(MO_18S_Merged_MO_rare))


MO_18S_Merged_MO_unrarified=readRDS("./MO_18S_Merged_MO_unrarified.rds")
sample_names(MO_18S_Merged_MO_unrarified) <- sub("^MO_", "", sample_names(MO_18S_Merged_MO_unrarified))
samp_order=rev(sample_names(MO_18S_Merged_MO_unrarified))

###############################################################################
#Load the updated taxonomy
###############################################################################
Taxa_Update <- read.xlsx("./ASV_Table_18S_unrare_Modified_Incorporated.xlsx", sheet = 1)
Taxa_Update=Taxa_Update[,c(1:10)] %>% column_to_rownames("ASV")
tax_table(MO_18S_Merged_MO_unrarified) <- as.matrix(Taxa_Update)
###############################################################################
#Area plots
###############################################################################
# Step 1: Transform to relative abundance at Genus level
Genus_Rel = MO_18S_Merged_MO_unrarified %>%
    tax_transform("compositional", rank = "Genus")

# Step 2: Extract and reshape OTU table
otu_table = Genus_Rel %>%
    otu_table() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "Genus")

# Step 3: Calculate total abundance for each Genus
taxa_sums = rowSums(otu_table[, 2:22])  # assuming samples are in columns 2 to 22
taxa_counts = data.frame(Genus = otu_table$Genus, Total_Count = taxa_sums)

# Step 4: Identify top 10 genera
top_n_taxa = slice_max(taxa_counts, n = 10, order_by = Total_Count) %>%
    pull(Genus)

# Step 5: Melt the full data (no pruning yet)
abundance_melted = melt(otu_table, id.vars = "Genus") %>%
    dplyr::rename(DOY = variable) %>%
    mutate(
        DOY = as.numeric(as.character(DOY)),
        value = value * 100,  # Scale to percentages
        Genus = ifelse(Genus %in% top_n_taxa, Genus, "Other")  # Lump non-top10 into "Other"
    )

# Step 6: Group by DOY and Genus and sum values
abundance_melted = abundance_melted %>%
    group_by(DOY, Genus) %>%
    dplyr::summarize(value = sum(value), .groups = "drop")

# Step 7: Reorder factor levels
abundance_melted$Genus = fct_relevel(
    abundance_melted$Genus,
    rev(top_n_taxa),  # Top 10 in order
    "Other"           # "Other" last
)

# Step 8: Define color palette
colors_for_plot = c(
    paletteer::paletteer_d("colorBlindness::paletteMartin")[1:10],  # 10 colors
    "gray"  # Gray for "Other"
)


    Area_18S_Genus=
    ggplot(abundance_melted, aes(x = DOY, y = value, fill = Genus)) +
        geom_area(linetype = 1, linewidth =0.25 ,colour="black" )+
            scale_fill_manual(values = colors_for_plot) +
        geom_vline(xintercept=c(144), linetype=2)+
        scale_x_continuous(expand=c(0,0))+
        scale_y_continuous(expand=c(0,0))+
        theme(
          plot.title = element_text(hjust = 0.5, face="bold", size=10),
          strip.text.x = element_text(size=8, face="bold", hjust=0),
          #strip.text.x = element_blank(),
          strip.text.y = element_text(size=8, face="bold.italic"),
          legend.text = element_text( size=8),
          legend.title = element_text( size=10, face="bold"),
          legend.key.size = unit(0.4, 'cm'),
          axis.title = element_text(size=10, face="bold"),
          axis.title.x = element_text(vjust=-1),
          axis.text.x = element_text(size=8, vjust = -1),
          axis.text = element_text(size=8, face="bold",  colour = "black"),
          axis.ticks =element_line(color = "black"),
          panel.spacing = unit(2.5, "lines"),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
        )+
        labs(tag="d", x = "Day of year", y = "% Relative abundance", fill = "Genus", title="Most abundant eukaryotic genera (18S rRNA)")

#########################
#nmds
#########################

MO_18S_nometazoa=readRDS(file='./MO_18S_wometazoa.RDS')
rownames(MO_18S_nometazoa) <- sub("^MO_", "", rownames(MO_18S_nometazoa))

MO_18S_nometazoa=MO_18S_nometazoa[,c(1:740)]%>% as.matrix()
rowSums(MO_18S_nometazoa)


taxa_18S_nometa=readRDS(file='./MO_18S_Merged_MO_rare_nometa_tax_table.RDS')
################################################################################
#Make a cluster dendrogram using all 21 samples
################################################################################
set.seed(06041990)
MO_18S_Full_dist=vegdist(MO_18S_nometazoa, method = "bray")

Full_clust <- hclust(MO_18S_Full_dist, "ward.D2")

clust_2<-cutree(Full_clust, k = 2)
clust_3<-cutree(Full_clust, k = 3)

Full_clust_bdisp_2=betadisper(MO_18S_Full_dist, clust_2)
anova(Full_clust_bdisp_2)

MO_18S_nometazoa_Wisc <- wisconsin(MO_18S_nometazoa)
MDS_Wisc=metaMDS(MO_18S_nometazoa_Wisc, distance='bray', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)

clust_2_df=as.data.frame(clust_2)
clust_3_df=as.data.frame(clust_3)

clust_2_df <- tibble::rownames_to_column(clust_2_df, var = "DOY")
clust_2_df$DOY <- as.numeric(clust_2_df$DOY)
clust_3_df <- tibble::rownames_to_column(clust_3_df, var = "DOY")
clust_3_df$DOY <- as.numeric(clust_3_df$DOY)

Enviro_Chl_clust <- Enviro_Chl %>%
  full_join(clust_2_df, by = "DOY") %>%
  full_join(clust_3_df, by = "DOY") %>%
  mutate(
    clust_2 = as.factor(clust_2),
    clust_3 = as.factor(clust_3)
  )

point_locations_MDS_Wisc <- scores(MDS_Wisc, display = "sites") %>%
as.data.frame() %>%
rownames_to_column("DOY") %>%
mutate(DOY=as.numeric(DOY))%>%
left_join(Enviro_Chl_clust, by="DOY")

point_locations_MDS_Wisc
levels(point_locations_MDS_Wisc$clust_3)[levels(point_locations_MDS_Wisc$clust_3)=='2'] <- NA


fill_palette_18S <- c("orange", "purple", "purple")

centroids_18S <- point_locations_MDS_Wisc %>%
    group_by(clust_2) %>%
    dplyr::summarize(
        centroid_NMDS1 = mean(NMDS1, na.rm = TRUE),
        centroid_NMDS2 = mean(NMDS2, na.rm = TRUE)
    ) %>%
    mutate(groups = case_when(
        clust_2 == 1 ~ "Bloom",
        clust_2 == 2 ~ "Post-",
        TRUE ~ NA_character_
    )) %>%
    mutate(
        centroid_NMDS2 = case_when(
            clust_2 == "1" ~ 0.15,  # Adjust value for clust_2 == 1
            TRUE ~ centroid_NMDS2     # Retain existing values otherwise
        ),
        centroid_NMDS1 = case_when(
            clust_2 == "2" ~ 0.09,    # Adjust value for clust_2 == 2
            TRUE ~ centroid_NMDS1     # Retain existing values otherwise
        )
    )


centroids_18S_2 <- point_locations_MDS_Wisc %>%
    group_by(clust_2) %>%
    dplyr::summarize(
        centroid_NMDS1 = mean(NMDS1, na.rm = TRUE),
        centroid_NMDS2 = mean(NMDS2, na.rm = TRUE)
    ) %>%
    mutate(groups = case_when(
        clust_2 == 1 ~ "E1",
        clust_2 == 2 ~ "E2",
        TRUE ~ NA_character_
    )) %>%
    mutate(
        centroid_NMDS2 = case_when(
            clust_2 == "1" ~ 0.15,  # Adjust value for clust_2 == 1
            TRUE ~ centroid_NMDS2     # Retain existing values otherwise
        ),
        centroid_NMDS1 = case_when(
            clust_2 == "2" ~ 0.09,    # Adjust value for clust_2 == 2
            TRUE ~ centroid_NMDS1     # Retain existing values otherwise
        )
    )


New_MDS_Wisc_18S=
ggplot(data = point_locations_MDS_Wisc,
       aes(x =NMDS1, y =  NMDS2, group = clust_2,
           colour = DOY,
           fill = clust_3,
           label = clust_3,
           # shape = Phase_2
       ))+
    geom_encircle( data = point_locations_MDS_Wisc[!is.na(point_locations_MDS_Wisc$clust_3),],fill=NA, s_shape = 1, expand = 0,
                   alpha = 1, color = "black", size=1, show.legend = FALSE, na.rm=T)+
    geom_encircle(data = point_locations_MDS_Wisc[!is.na(point_locations_MDS_Wisc$clust_3),], aes(fill = clust_3), s_shape = 1, expand = 0,
                  alpha = 0.2, color = "black", size=1, show.legend = FALSE, na.rm=T)+
    geom_point(size = 2)+
    scale_colour_gradient(low = "lightblue", high = "darkblue", na.value = NA,
  guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", frame.linewidth =0.75))+
annotate("text", x = -1.25, y = 0.4,
         label = paste0("Stress = ",
                        round(MDS_Wisc[['stress']], digits = 3)
         ),
       size=3)+
      geom_text(data = centroids_18S_2, aes(x = centroid_NMDS1, y = centroid_NMDS2, label = groups, fill=NA),
      size = 4, fontface = "bold", color = "black") +
    theme_few()+
    theme(
      plot.title = element_text(hjust = 0.5, face="bold", size=10),
      strip.text.x = element_text(size=8, face="bold", hjust=0),
      #strip.text.x = element_blank(),
      strip.text.y = element_text(size=8, face="bold.italic"),
      legend.text = element_text( size=8),
      legend.title = element_text( size=10, face="bold"),
      legend.key.size = unit(0.4, 'cm'),
      axis.title = element_text(size=10, face="bold"),
      axis.title.x = element_text(vjust=-1),
      axis.text.x = element_text(size=8, vjust = -1),
      axis.text = element_text(size=8, face="bold",  colour = "black"),
      axis.ticks =element_line(color = "black"),
      panel.spacing = unit(2.5, "lines"),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      plot.margin = margin(10,10,10,10)
    )+
    labs(tag="c", x = "nMDS 1", y = "nMDS 2", color = "Day of year",
     title="Eukaryotic nMDS (18S rRNA)")+
      scale_fill_manual(values=fill_palette_18S)+
      guides(fill="none")

combined=New_MDS_Wisc_16S_clust_3+Area_16S_Genus+New_MDS_Wisc_18S+Area_18S_Genus+
    plot_layout(widths = c(1, 2),
    design="12
            34")+
    plot_annotation(tag_levels = 'a')
