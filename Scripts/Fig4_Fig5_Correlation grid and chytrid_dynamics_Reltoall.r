library(tidyverse) # for plotting and wrangling data
library(SpiecEasi) # Has sparcc and also does clr transforms
library(otuSummary)
library(reshape2) # has the melt function, which I use to wrangle data
library(psych) # for calculating regular correlations with p values
pass <- function(x){x}


#Tutorial https://biovcnet.github.io/_pages/NetworkScience_SparCC.nb.html

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

# MO_18S_Merged_MO_rare=readRDS("G:/My Drive/01_MOSAIC/Data/Phyloseq_Objects/MO_18S_Merged_rare.rds")
MO_18S_Merged_MO_unrarified=readRDS("G:/My Drive/01_MOSAIC/Data/Phyloseq_Objects/MO_18S_Merged_MO_unrarified.rds")
sample_names(MO_18S_Merged_MO_unrarified) <- sub("^MO_", "", sample_names(MO_18S_Merged_MO_unrarified))

samp_order=rev(sample_names(MO_18S_Merged_MO_unrarified))

MO_18S_Merged_MO_unrarified_species=tax_glom(MO_18S_Merged_MO_unrarified, taxrank=rank_names(MO_18S_Merged_MO_unrarified)[9])

#Filters based on at least 0.1% in 50% of samples.
zwh0=genefilter_sample(MO_18S_Merged_MO_unrarified, function(x) x > 0.001 * sum(x), A = 0.5 * nsamples(MO_18S_Merged_MO_unrarified))
zwh0_sp=genefilter_sample(MO_18S_Merged_MO_unrarified_species, function(x) x > 0.001 * sum(x), A = 0.5 * nsamples(MO_18S_Merged_MO_unrarified_species))

zFiltered_18S = prune_taxa(zwh0, MO_18S_Merged_MO_unrarified) ##Remove the selected OTUs
zFiltered_18S_sp = prune_taxa(zwh0_sp, MO_18S_Merged_MO_unrarified_species) ##Remove the selected OTUs

sort(sample_sums(zFiltered_18S))
sort(sample_sums(zFiltered_18S_sp))

zFiltered_18S_Rel <- transform_sample_counts(zFiltered_18S, function(x) x / sum(x))
zFiltered_18S_Rel_sp <- transform_sample_counts(MO_18S_Merged_MO_unrarified_species, function(x) x / sum(x))

zMO_otu=zFiltered_18S_Rel %>% otu_table() %>% as.data.frame() #%>% rownames_to_column( var="ASV")
rowSums(zMO_otu)
zMO_otu_sp=zFiltered_18S_Rel_sp %>% otu_table() %>% as.data.frame() #%>% rownames_to_column( var="ASV")
rowSums(zMO_otu_sp)



# MO_taxa=Filtered_18S %>% tax_table() %>% as.data.frame() %>% rownames_to_column( var="ASV")
zMO_taxa=zFiltered_18S %>% tax_table() %>% as.data.frame() %>% rownames_to_column( var="ASV")
#ASV here is the actual character string designation, so I cant use species as a rowname or it will be duplicated.
zMO_taxa_sp=zFiltered_18S_sp %>% tax_table() %>% as.data.frame() %>% rownames_to_column( var="ASV")


    zMO_taxa=zMO_taxa %>%
        # build grouping by combination of variables
        dplyr::group_by(Genus, Species) %>%
        # add row number which works per group due to prior grouping
        dplyr::mutate(duplicateID = dplyr::row_number()) %>%
         unite(Identifier, c("Species", "duplicateID"), remove=F) %>%
        # ungroup to prevent unexpected behaviour down stream
        dplyr::ungroup()


zparasitoid_otu=otu_table(subset_taxa(zFiltered_18S_Rel, Class =="Chytridiomycota" | Family=="Protaspa-lineage"| Family=="Cryothecomonas-lineage"| Family=="Cryomonadida_X"| Class=="Peronosporomycetes" | Class =="Syndiniales" | Order =="Labyrinthulomycetes")) %>%
    as.data.frame() %>% rownames_to_column( var="ASV")
zparasitoid_tax=tax_table(subset_taxa(zFiltered_18S_Rel, Class =="Chytridiomycota" | Family=="Protaspa-lineage"| Family=="Cryothecomonas-lineage"| Family=="Cryomonadida_X"| Class=="Peronosporomycetes" | Class =="Syndiniales" | Order =="Labyrinthulomycetes")) %>%
    as.data.frame() %>% rownames_to_column( var="ASV")
rownames(zparasitoid_otu) <- zparasitoid_otu[,1]
zparasitoid_otu=zparasitoid_otu[,-1]

#Species
zparasitoid_otu_sp=otu_table(subset_taxa(zFiltered_18S_Rel_sp, Class =="Chytridiomycota" | Family=="Protaspa-lineage"| Family=="Cryothecomonas-lineage"| Family=="Cryomonadida_X"| Class=="Peronosporomycetes" | Class =="Syndiniales" | Order =="Labyrinthulomycetes")) %>%
    as.data.frame() %>% rownames_to_column( var="ASV")
zparasitoid_tax_sp=tax_table(subset_taxa(zFiltered_18S_Rel_sp, Class =="Chytridiomycota" | Family=="Protaspa-lineage"| Family=="Cryothecomonas-lineage"| Family=="Cryomonadida_X"| Class=="Peronosporomycetes" | Class =="Syndiniales" | Order =="Labyrinthulomycetes")) %>%
    as.data.frame() %>% rownames_to_column( var="ASV")
rownames(zparasitoid_otu_sp) <- zparasitoid_otu_sp[,1]
zparasitoid_otu_sp=zparasitoid_otu_sp[,-1]


zPhyto_otu=otu_table(subset_taxa(zFiltered_18S_Rel, Class =="Bacillariophyceae" |Class =="Mediophyceae"|Class =="Coscinodiscophyceae"| Class=="Dinophyceae")) %>%
    as.data.frame() %>% rownames_to_column( var="ASV") #Dinophyceae exludes Syndinialtes and keep all the others
    zPhyto_tax=tax_table(subset_taxa(zFiltered_18S_Rel, Class =="Bacillariophyceae" |Class =="Mediophyceae"|Class =="Coscinodiscophyceae"| Class=="Dinophyceae")) %>%
        as.data.frame() %>% rownames_to_column( var="ASV") #Dinophyceae exludes Syndinialtes and keep all the others

rownames(zPhyto_otu) <- zPhyto_otu[,1]
zPhyto_otu=zPhyto_otu[,-1]

#species level
zPhyto_otu_sp=otu_table(subset_taxa(zFiltered_18S_Rel_sp, Class =="Bacillariophyceae" |Class =="Mediophyceae"|Class =="Coscinodiscophyceae"| Class=="Dinophyceae")) %>%
    as.data.frame() %>% rownames_to_column( var="ASV") #Dinophyceae exludes Syndinialtes and keep all the others
    zPhyto_tax=tax_table(subset_taxa(zFiltered_18S_Rel_sp, Class =="Bacillariophyceae" |Class =="Mediophyceae"|Class =="Coscinodiscophyceae"| Class=="Dinophyceae")) %>%
        as.data.frame() %>% rownames_to_column( var="ASV") #Dinophyceae exludes Syndinialtes and keep all the others

rownames(zPhyto_otu_sp) <- zPhyto_otu_sp[,1]
zPhyto_otu_sp=zPhyto_otu_sp[,-1]



zPhytoPara_otu=otu_table(subset_taxa(zFiltered_18S_Rel, Class =="Bacillariophyceae" |Class =="Mediophyceae"|Class =="Coscinodiscophyceae"| Class=="Dinophyceae" | Class =="Chytridiomycota" | Family=="Protaspa-lineage"| Family=="Cryothecomonas-lineage"| Family=="Cryomonadida_X"| Class=="Peronosporomycetes" | Class =="Syndiniales" | Order =="Labyrinthulomycetes")) %>%
    as.data.frame() %>% rownames_to_column( var="ASV") #Dinophyceae exludes Syndinialtes and keep all the others
    zPhytoPara_otu_tax=tax_table(subset_taxa(zFiltered_18S_Rel, Class =="Bacillariophyceae" |Class =="Mediophyceae"|Class =="Coscinodiscophyceae"| Class=="Dinophyceae" | Class =="Chytridiomycota" | Family=="Protaspa-lineage"| Family=="Cryothecomonas-lineage"| Family=="Cryomonadida_X"| Class=="Peronosporomycetes" | Class =="Syndiniales" | Order =="Labyrinthulomycetes")) %>%
        as.data.frame() %>% rownames_to_column( var="ASV") #Dinophyceae exludes Syndinialtes and keep all the others

rownames(zPhytoPara_otu) <- zPhytoPara_otu[,1]
zPhytoPara_otu=zPhytoPara_otu[,-1]

#species level
zPhytoPara_otu_sp=otu_table(subset_taxa(zFiltered_18S_Rel_sp, Class =="Bacillariophyceae" |Class =="Mediophyceae"|Class =="Coscinodiscophyceae"| Class=="Dinophyceae" | Class =="Chytridiomycota" | Family=="Protaspa-lineage"| Family=="Cryothecomonas-lineage"| Family=="Cryomonadida_X"| Class=="Peronosporomycetes" | Class =="Syndiniales" | Order =="Labyrinthulomycetes")) %>%
    as.data.frame() %>% rownames_to_column( var="ASV") #Dinophyceae exludes Syndinialtes and keep all the others
    zPhytoPara_otu_tax_sp=tax_table(subset_taxa(zFiltered_18S_Rel_sp, Class =="Bacillariophyceae" |Class =="Mediophyceae"|Class =="Coscinodiscophyceae"| Class=="Dinophyceae" | Class =="Chytridiomycota" | Family=="Protaspa-lineage"| Family=="Cryothecomonas-lineage"| Family=="Cryomonadida_X"| Class=="Peronosporomycetes" | Class =="Syndiniales" | Order =="Labyrinthulomycetes")) %>%
        as.data.frame() %>% rownames_to_column( var="ASV") #Dinophyceae exludes Syndinialtes and keep all the others

rownames(zPhytoPara_otu_sp) <- zPhytoPara_otu_sp[,1]
zPhytoPara_otu_sp=zPhytoPara_otu_sp[,-1]



############################################################################
#Directed correlations
############################################################################


zPhyto_Clr <- clr(zPhyto_otu)
zPara_Clr <- clr(zparasitoid_otu)



zspearCorTestClr <- corr.test(zPhyto_Clr, zPara_Clr, method = "spearman", adjust = "none")
zspearCorClr <- zspearCorTestClr$r
zspearCorClr[is.na(zspearCorClr)] = 0
zspearPClr <- zspearCorTestClr$p
zspearPClr[is.na(zspearPClr)] = 0

zPhyto_Clr_sp <- clr(zPhyto_otu_sp)
zPara_Clr_sp <- clr(zparasitoid_otu_sp)

zspearCorTestClr_sp <- corr.test(zPhyto_Clr_sp, zPara_Clr_sp, method = "spearman", adjust = "none")
zspearCorClr_sp <- zspearCorTestClr_sp$r
zspearCorClr_sp[is.na(zspearCorClr_sp)] = 0
zspearPClr_sp <- zspearCorTestClr_sp$p
zspearPClr_sp[is.na(zspearPClr_sp)] = 0

zspearCorClr_processed <- zspearCorClr  %>% reshape2::melt() %>% na.omit() %>% dplyr::rename(rho = value)
zspearPClr_processed <- zspearPClr  %>% reshape2::melt() %>% na.omit() %>% dplyr::rename(p = value)
zspearRhoP_Clr <- left_join(zspearCorClr_processed, zspearPClr_processed, by = c("Var1", "Var2")) %>%
# calculate the false discovery rate to adjust for multiple p values
mutate(fdr = p.adjust(p, method = "BH"))

zresultx_Clr <- zMO_taxa %>%
    inner_join(zspearRhoP_Clr, by = c("ASV" = "Var1"))
zresulty_Clr <- zMO_taxa %>%
    inner_join(zresultx_Clr, by = c("ASV" = "Var2"), suffix = c(".var1", ".var2"))

fdrThreshy <- 0.05 # fdr threshold
zspearOkPy_Clr <- zresulty_Clr%>% filter(fdr < fdrThreshy)


zspearCorClr_processed_sp <- zspearCorClr_sp  %>% reshape2::melt() %>% na.omit() %>% dplyr::rename(rho = value)
zspearPClr_processed_sp <- zspearPClr_sp  %>% reshape2::melt() %>% na.omit() %>% dplyr::rename(p = value)
zspearRhoP_Clr_sp <- left_join(zspearCorClr_processed_sp, zspearPClr_processed_sp, by = c("Var1", "Var2")) %>%
# calculate the false discovery rate to adjust for multiple p values
mutate(fdr = p.adjust(p, method = "BH"))

zresultx_Clr_sp <- zMO_taxa_sp %>%
    inner_join(zspearRhoP_Clr_sp, by = c("ASV" = "Var1"))
zresulty_Clr_sp <- zMO_taxa %>%
    inner_join(zresultx_Clr_sp, by = c("ASV" = "Var2"), suffix = c(".var1", ".var2"))

fdrThreshy <- 0.05 # fdr threshold
zspearOkPy_Clr_sp <- zresulty_Clr_sp%>% filter(fdr < fdrThreshy)


combo_sig_subset_full <- zresulty_Clr %>%
    filter(Identifier.var1 %in% zspearOkPy_Clr$Identifier.var1,
           Identifier.var2 %in% zspearOkPy_Clr$Identifier.var2)

combo_sig_subset_full_ordered <- as.data.frame(combo_sig_subset_full)


combo_sig_subset_full_ordered=combo_sig_subset_full_ordered %>%
mutate(Identifier.var1 = recode(Identifier.var1,
'Algomyces_stechlinensis_1'= "Lobulomycetales Order_1",
'Cryomonadida_XX_sp._1'="Cryomonadida Order_1",
'Cryomonadida_XX_sp._2'="Cryomonadida Order_2",
'Cryothecomonas-lineage_X_sp._1'="Cryothecomonas-lineage_1",
'Cryothecomonas-lineage_X_sp._2'="Cryothecomonas-lineage_2",
'Cryothecomonas-lineage_X_sp._3'="Cryothecomonas-lineage_3",
'Cryothecomonas-lineage_X_sp._4'="Cryothecomonas-lineage_4",
'Protaspa-lineage_X_sp._1'="Protaspa-lineage_1",
'Protaspa-lineage_X_sp._2'="Protaspa-lineage_2",
'Protaspa-lineage_X_sp._3'="Protaspa-lineage_3",
'Protaspa-lineage_X_sp._4'="Protaspa-lineage_4",
'Protaspa-lineage_X_sp._5'="Protaspa-lineage_5",
'Peronosporomycetes_XXX_sp._1'="Oomycetes Class_1",
'Fibrophrys_sp._1'="Labyrinthulomycetes Order_1",
'Oblongichytrium_sp._1'="Oblongichytrium Genus_1",
))
zspearOkPy_Clr=zspearOkPy_Clr %>%
mutate(Identifier.var1 = recode(Identifier.var1,
  'Algomyces_stechlinensis_1'= "Lobulomycetales Order_1",
  'Cryomonadida_XX_sp._1'="Cryomonadida Order_1",
  'Cryomonadida_XX_sp._2'="Cryomonadida Order_2",
  'Cryothecomonas-lineage_X_sp._1'="Cryothecomonas-lineage_1",
  'Cryothecomonas-lineage_X_sp._2'="Cryothecomonas-lineage_2",
  'Cryothecomonas-lineage_X_sp._3'="Cryothecomonas-lineage_3",
  'Cryothecomonas-lineage_X_sp._4'="Cryothecomonas-lineage_4",
  'Protaspa-lineage_X_sp._1'="Protaspa-lineage_1",
  'Protaspa-lineage_X_sp._2'="Protaspa-lineage_2",
  'Protaspa-lineage_X_sp._3'="Protaspa-lineage_3",
  'Protaspa-lineage_X_sp._4'="Protaspa-lineage_4",
  'Protaspa-lineage_X_sp._5'="Protaspa-lineage_5",
  'Peronosporomycetes_XXX_sp._1'="Oomycetes Class_1",
  'Fibrophrys_sp._1'="Labyrinthulomycetes Order_1",
  'Oblongichytrium_sp._1'="Oblongichytrium Genus_1",
  ))
combo_sig_subset_full_ordered=combo_sig_subset_full_ordered %>%
mutate(Identifier.var2 = recode(Identifier.var2,
'Amphiprora_alata_1'="Entomoneis Genus_1",
'Amphiprora_alata_2'="Entomoneis Genus_2",
'Attheya_septentrionalis_1'= "Attheya Genus_1",
'Cymbellales Order_1'="Cymbellales Order_1",
'Fragilariopsis_sp._1'="Fragilariopsis Genus_1",
'Navicula Genus_1'="Naviculaceae Family_2",
'Navicula_amphiceropsis_1'="Navicula Genus_1",
'Navicula_amphiceropsis_2'="Navicula Genus_2",
'Navicula_amphiceropsis_3'="Navicula Genus_3",
'Naviculaceae_X_sp._1'="Naviculaceae Family_1",
'Naviculales Order_1'="Naviculales Order_1",
'Nitzschia Genus_1'="Nitzschia Genus_2",
'Synedropsis_sp._1'="Synedropsis Genus_1",
'Thalassiosiraceae Family_1'="Thalassiosiraceae Family_1",
'Tryblionella Genus_1'="Nitzschia Genus_3",
'Tryblionella Genus_3'="Nitzschia Genus_6",
'Tryblionella Genus_4'="Nitzschia Genus_5",
'Tryblionella Genus_5'="Nitzschia Genus_4",
'Tryblionella_levidensis_1'="Nitzschia Genus_1",
# 'Tryblionella_levidensis_2'="Tryblionella_levidensis_sv2",
'Apocalathium_euryceps_1'="Thoracosphaeraceae Family_1",
# 'Dinophyceae Class_1'="Dinophyceae Class_sv1",
'Gymnodiniaceae Family_1'="Gymnodiniaceae Family_1",
'Gymnodinium_sp._1'="Gymnodinium Genus_1",
# 'Gymnodinium_sp._2'="Gymnodinium_sv2",
'Heterocapsa Genus_1'="Heterocapsa Genus_1",
'Suessiaceae Family_2'="Suessiaceae Family_2",
'Suessiaceae Family_3'="Suessiaceae Family_3"
))

zspearOkPy_Clr=zspearOkPy_Clr %>%
mutate(Identifier.var2 = recode(Identifier.var2,
  'Amphiprora_alata_1'="Entomoneis Genus_1",
  'Amphiprora_alata_2'="Entomoneis Genus_2",
  'Cymbellales Order_1'="Cymbellales Order_1",
  'Fragilariopsis_sp._1'="Fragilariopsis Genus_1",
  'Navicula Genus_1'="Naviculaceae Family_2",
  'Navicula_amphiceropsis_1'="Navicula Genus_1",
  'Navicula_amphiceropsis_2'="Navicula Genus_2",
  'Navicula_amphiceropsis_3'="Navicula Genus_3",
  'Naviculaceae_X_sp._1'="Naviculaceae Family_1",
  'Naviculales Order_1'="Naviculales Order_1",
  'Nitzschia Genus_1'="Nitzschia Genus_2",
  'Synedropsis_sp._1'="Synedropsis Genus_1",
  'Tryblionella Genus_1'="Nitzschia Genus_3",
  'Tryblionella Genus_3'="Nitzschia Genus_6",
  'Tryblionella Genus_4'="Nitzschia Genus_5",
  'Tryblionella Genus_5'="Nitzschia Genus_4",
  'Tryblionella_levidensis_1'="Nitzschia Genus_1",
  # 'Tryblionella_levidensis_2'="Tryblionella_levidensis_sv2",
  'Apocalathium_euryceps_1'="Thoracosphaeraceae Family_1",
  # 'Dinophyceae Class_1'="Dinophyceae Class_sv1",
  'Gymnodiniaceae Family_1'="Gymnodiniaceae Family_1",
  'Gymnodinium_sp._1'="Gymnodinium Genus_1",
  # 'Gymnodinium_sp._2'="Gymnodinium_sv2",
  'Heterocapsa Genus_1'="Heterocapsa Genus_1",
  'Suessiaceae Family_2'="Suessiaceae Family_2",
  'Suessiaceae Family_3'="Suessiaceae Family_3",
  'Attheya_septentrionalis_1'= "Attheya Genus_1",
  'Thalassiosiraceae Family_1'="Thalassiosiraceae Family_1"
  ))










combo_sig_subset_full_ordered$Class.var2 <- factor(combo_sig_subset_full_ordered$Class.var2,
                                                   levels = c("Bacillariophyceae", "Mediophyceae", "Dinophyceae"))

combo_sig_subset_full_ordered <- combo_sig_subset_full_ordered[order(combo_sig_subset_full_ordered$Class.var2,
                                                                     combo_sig_subset_full_ordered$Identifier.var2),]

combo_sig_subset_full_ordered$Identifier.var2 <- factor(combo_sig_subset_full_ordered$Identifier.var2,
                                                        levels = unique(combo_sig_subset_full_ordered$Identifier.var2),
                                                        ordered = TRUE)

combo_sig_subset_full_ordered$Class.var1 <- as.factor(combo_sig_subset_full_ordered$Class.var1)
combo_sig_subset_full_ordered <- combo_sig_subset_full_ordered[order(combo_sig_subset_full_ordered$Class.var1,
                                                                     combo_sig_subset_full_ordered$Identifier.var1),]

combo_sig_subset_full_ordered$Identifier.var1 <- factor(combo_sig_subset_full_ordered$Identifier.var1,
                                                        levels = unique(combo_sig_subset_full_ordered$Identifier.var1),
                                                        ordered = TRUE)


#species level
combo_sig_subset_full_sp=zresulty_Clr_sp %>%
      filter(Species.var1 %in% zspearOkPy_Clr_sp$Species.var1,
      Species.var2 %in% zspearOkPy_Clr_sp$Species.var2)

      combo_sig_subset_full_ordered_sp=as.data.frame(combo_sig_subset_full_sp)
      combo_sig_subset_full_ordered_sp$Class.var2 <- factor(combo_sig_subset_full_ordered_sp$Class.var2,
                                                         levels = c("Bacillariophyceae", "Mediophyceae", "Dinophyceae"))

      # zspearOkPy_Clr_ordered$Species.var2=as.factor(zspearOkPy_Clr_ordered$Species.var2)
      combo_sig_subset_full_ordered_sp <- combo_sig_subset_full_ordered_sp[order(combo_sig_subset_full_ordered_sp$Class.var2,combo_sig_subset_full_ordered_sp$Species.var2),]
      factor(combo_sig_subset_full_ordered_sp$Species.var2,levels = unique(combo_sig_subset_full_ordered_sp$Species.var2),ordered = T)
      combo_sig_subset_full_ordered_sp$Species.var2 <- factor(combo_sig_subset_full_ordered_sp$Species.var2,levels = unique(combo_sig_subset_full_ordered_sp$Species.var2),ordered = T)

      combo_sig_subset_full_ordered_sp$Class.var1=as.factor(combo_sig_subset_full_ordered_sp$Class.var1)
      combo_sig_subset_full_ordered_sp <- combo_sig_subset_full_ordered_sp[order(combo_sig_subset_full_ordered_sp$Class.var1,combo_sig_subset_full_ordered_sp$Species.var1),]
      factor(combo_sig_subset_full_ordered_sp$Species.var1,levels = unique(combo_sig_subset_full_ordered_sp$Species.var1),ordered = T)
      combo_sig_subset_full_ordered_sp$Species.var1 <- factor(combo_sig_subset_full_ordered_sp$Species.var1,levels = unique(combo_sig_subset_full_ordered_sp$Species.var1),ordered = T)


x_axis_manually_curated_pal=c("darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","forestgreen","forestgreen",
"black","black","black","black","black","black","black","black")

y_axis_manually_curated_pal=c("darkred", "darkblue", "darkblue","darkblue","darkblue", "darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue",
"darkorange", "purple", "purple")

#Italicize only if the
italic_labeller <- function(labels) {
  sapply(labels, function(label) {
    if (grepl("Genus", label)) {
      # Find the starting position of "Genus"
      pos <- regexpr("Genus", label)
      pre <- substr(label, 1, pos - 1)
      post <- substr(label, pos, nchar(label))
      # If there is text before "Genus", format it as bolditalic and the rest as bold.
      if (nchar(pre) > 0) {
        expr <- paste0("bolditalic('", pre, "') * bold('", post, "')")
      } else {
        expr <- paste0("bold('", label, "')")
      }
    } else {
      # If "Genus" is not found, render the entire label as bold.
      expr <- paste0("bold('", label, "')")
    }
    expr
  }, USE.NAMES = FALSE)
}

#Fig. 4
  Parasitoid_ordered_subbedw_sigmarker=
  combo_sig_subset_full_ordered %>% filter(rho >= 0.5 | rho <= -0.5)%>%
  ggplot(aes(x = Identifier.var2, y = Identifier.var1, fill = rho)) +
  geom_tile() +
  #use scale discrete to remove the small grey buffer
  geom_tile(color="black") +
  scale_x_discrete( expand=c(0,0), labels = function(x) parse(text = italic_labeller(x))) +
  scale_y_discrete(expand=c(0,0), labels = function(x) parse(text = italic_labeller(x)))+
  geom_point(data = zspearOkPy_Clr, shape = 16)+
  scale_fill_gradient2(limits = c(-1, 1),
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", frame.linewidth =0.5))+
  labs(x="Putative hosts", y="Putative parasitoids", title="Parasitoid-Host interactions", fill=expression(rho))+
  theme(
      # plot.title = element_text(hjust = 0.5, face="bold", size=25),
      plot.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, color=x_axis_manually_curated_pal),
      axis.text.y = element_text(color=y_axis_manually_curated_pal),
      axis.text = element_text(face = "bold", size=8),
      axis.title = element_text(size=10, face="bold"),
      axis.ticks =element_line(color = "black"),
      legend.text = element_text( size=8),
      legend.title = element_text(face="bold", size=12),
      legend.key.size = unit(1, 'lines'),
      legend.text.align = 0,
      panel.spacing = unit(3.5, "lines"),
      strip.background = element_blank(),
      panel.grid = element_blank(),               # Remove grid lines
      panel.border = element_rect(colour = "black", fill=NA, size=0.5),
      plot.margin = margin(2,2,2,2)
  )



####################################
# Sort descending the sum of all rows for "numerically important hosts" to investigate
####################################
zCounts=merge(zMO_taxa %>% select(ASV, Identifier), zMO_otu %>% t()%>%as.data.frame()%>%rownames_to_column(var="ASV"))


zCounts %>%
    mutate(row_sum = rowSums(select(., starts_with("1")))) %>%
    select(Identifier, row_sum) %>% arrange(desc(row_sum))


zCounts %>%
    filter(str_detect(Identifier, "Tryblionella"))

Target="Tryblionella_levidensis_1"
vec=zspearOkPy_Clr %>%
filter(Identifier.var2 ==Target) %>%
pull(Identifier.var1) #%>%


zCounts %>% filter(Identifier== Target | Identifier %in% vec)%>% select(2:23) %>%
    pivot_longer( cols=-c(Identifier), names_to = "DOY", values_to = "Count") %>%
    ggplot()+geom_point(aes(DOY, Count, color=Identifier))+
    geom_line(aes(DOY, Count, group=Identifier, color=Identifier))



zCounts=zCounts %>%
mutate(Identifier = recode(Identifier,
  'Algomyces_stechlinensis_1'= "Lobulomycetales Order_1",
  'Cryomonadida_XX_sp._1'="Cryomonadida Order_1",
  'Cryomonadida_XX_sp._2'="Cryomonadida Order_2",
  'Cryothecomonas-lineage_X_sp._1'="Cryothecomonas Genus_1",
  'Cryothecomonas-lineage_X_sp._2'="Cryothecomonas Genus_2",
  'Cryothecomonas-lineage_X_sp._3'="Cryothecomonas Genus_3",
  'Cryothecomonas-lineage_X_sp._4'="Cryothecomonas Genus_4",
  'Protaspa-lineage_X_sp._1'="Protaspa Genus_1",
  'Protaspa-lineage_X_sp._2'="Protaspa Genus_2",
  'Protaspa-lineage_X_sp._3'="Protaspa Genus_3",
  'Protaspa-lineage_X_sp._4'="Protaspa Genus_4",
  'Protaspa-lineage_X_sp._5'="Protaspa Genus_5",
  'Peronosporomycetes_XXX_sp._1'="Oomycota Class_1",
  'Fibrophrys_sp._1'="Labyrinthulomycetes Order_1",
  'Oblongichytrium_sp._1'="Oblongichytrium Genus_1",
'Amphiprora_alata_1'="Entomoneis Genus_1",
'Amphiprora_alata_2'="Entomoneis Genus_2",
'Attheya_septentrionalis_1'= "Attheya Genus_1",
'Cymbellales Order_1'="Cymbellales Order_1",
'Fragilariopsis_sp._1'="Fragilariopsis Genus_1",
'Navicula Genus_1'="Naviculaceae Family_2",
'Navicula_amphiceropsis_1'="Navicula Genus_1",
'Navicula_amphiceropsis_2'="Navicula Genus_2",
'Navicula_amphiceropsis_3'="Navicula Genus_3",
'Naviculaceae_X_sp._1'="Naviculaceae Family_1",
'Naviculales Order_1'="Naviculales Order_1",
'Nitzschia Genus_1'="Nitzschia Genus_2",
'Synedropsis_sp._1'="Synedropsis Genus_1",
'Thalassiosiraceae Family_1'="Thalassiosiraceae Family_1",
'Tryblionella Genus_1'="Nitzschia Genus_3",
'Tryblionella Genus_3'="Nitzschia Genus_6",
'Tryblionella Genus_4'="Nitzschia Genus_5",
'Tryblionella Genus_5'="Nitzschia Genus_4",
'Tryblionella_levidensis_1'="Nitzschia Genus_1",
# 'Tryblionella_levidensis_2'="Tryblionella_levidensis_2",
'Apocalathium_euryceps_1'="Thoracosphaeraceae Family_1",
# 'Dinophyceae Class_1'="Dinophyceae Class_1",
'Gymnodiniaceae Family_1'="Gymnodiniaceae Family_1",
'Gymnodinium_sp._1'="Gymnodinium Genus_1",
'Gymnodinium_sp._2'="Gymnodinium Genus_2",
'Heterocapsa Genus_1'="Heterocapsa Genus_1",
'Suessiaceae Family_2'="Suessiaceae Family_2",
'Suessiaceae Family_3'="Suessiaceae Family_3"
))

####################
#Chytrid Nitzschia relationship Fig. 5
###################
zCounts_long <- zCounts %>%
  select(2:23) %>%
  pivot_longer(cols = -c(Identifier), names_to = "DOY", values_to = "Prop") %>%
  mutate(DOY=as.numeric(DOY))

Chytrid_Nitzschias=
ggplot() +
    annotate("rect", xmin = 110, xmax = 144, ymin = -Inf, ymax = Inf,
             alpha = .1,fill = "lightgreen")+
    geom_vline(xintercept=c(144), linetype=2, size=0.5)+
    geom_line(data = zCounts_long %>% filter(Identifier == "Nitzschia Genus_4"|Identifier == "Lobulomycetales Order_1"),
              aes(x = DOY, y = Prop*100, group = Identifier, color=Identifier), size = 0.75)+
              geom_line(data = zCounts_long %>% filter(Identifier == "Nitzschia Genus_1"|Identifier == "Nitzschia Genus_5"|Identifier == "Nitzschia Genus_6"),
                            aes(x = DOY, y = Prop*100, group = Identifier, color=Identifier), size = 0.5, linetype="dashed")+
              scale_color_manual( values=c("darkred", "grey40", "black", "orange", "purple"))+

      # geom_line(data = zCounts_long %>% filter(Identifier == "Algomyces_stechlinensis_sv1"),
      #           aes(x = DOY, y = Prop, group = Identifier), size = 1, color = "darkred") +
      theme_few()+
theme(
      axis.title.y= element_text(size=10, face="bold"),
      # axis.title.y=  element_blank(),
      axis.title.x= element_text(size=10, face="bold"),
      # axis.title.x= element_blank(),
      axis.text = element_text(face="bold", size=8, colour="black"),
      # axis.text.x = element_blank(),
      axis.ticks=element_line(size=0.5, color="black"),
      axis.ticks.length = unit(0.25, "cm"),
      legend.position = c(0.65, 0.8),
      # legend.background=element_rect(linetype="solid", colour = 'black'),
      strip.background = element_rect(
          color="white", fill="white",linetype="solid"),
      # legend.title=element_text(face="bold", size=8),
      legend.title=element_blank(),
      legend.title.align=0.5,
      legend.text=element_text( size=6),
      legend.text.align=0,
       legend.key.height = unit(0.5, "line"),
      legend.key.width = unit(1.4, "line"),
      # legend.background=element_blank(),
      legend.background=element_rect(linetype="solid", colour = 'black', size=0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(fill = "transparent", colour = NA),
      strip.text.x = element_text(
          size = 10, color = "black", face = "bold.italic"),
          panel.background = element_rect(colour = "black", size=0.5),
      )+
      scale_y_continuous(expand = c(0, 0),
      limits=c(0,40),
      breaks = seq(0, 40, by = 10)
      )+
      scale_x_continuous(expand = c(0, 0), minor_breaks= seq(110, 162, 5))+
      labs(x = "Day of year", y="Relative abundance (%)", color="Sequence variants")+
    theme(plot.margin=unit(c(0.2,0.5,0.2,0.2),"cm"))
