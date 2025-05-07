
library(tidyverse)
library(ggthemes)
library(patchwork)
library(cowplot)
# MO_16S_Merged_MO_unrarified=readRDS("/Users/kyle/Library/CloudStorage/GoogleDrive-kbdilliplaine@alaska.edu/My Drive/01_MOSAIC/Data/Phyloseq_Objects/MO_16S_Merged_MO_unrarified.rds")

Chl_avg=readRDS(file = 'G:/My Drive/01_MOSAIC/Enviro_New_Bottom10s_2024.RDS')
Enviro_Chl=readRDS(file='G:/My Drive/01_MOSAIC/Github Scripts/Enviro_Chl.RDS')

#####################################################################################################
# Seawater only
#####################################################################################################
sw_nuts_short <- readRDS("G:/My Drive/01_MOSAIC/Data/sw_nuts_short.rds")

Disc=c("#26E8EB", "#3A4EC2", "#040291", "darkgray", "darkgreen", "black")

cols    <- c( "c1" = "#26E8EB", "c2" = "#3A4EC2", "c3"="#040291", "c4"="black" )
lines  <- c("l1" = 1, "l2" = 2, 'l3'=3, 'l4'=4)

nuts_sw_mulitipanel_notation_nut_axes =
    ggplot(sw_nuts_short) +
    theme_few() +
    annotate("rect", xmin = 110, xmax = 144, ymin = -Inf, ymax = Inf, alpha = .1, fill = "lightgreen") +

    stat_smooth(aes(DOY, Ammonia, color = "c1", linetype = "l1"), linewidth = 0.75, se = FALSE, span = 0.5) +
    stat_smooth(aes(DOY, Nitrogen, color = "c2", linetype = "l2"), linewidth = 0.75, se = FALSE, span = 0.5) +
    stat_smooth(aes(DOY, Phosphate, color = "c3", linetype = "l3"), linewidth = 0.75, se = FALSE, span = 0.5) +
    stat_smooth(aes(DOY, Silicate / 3, color = "c4", linetype = "l4"), linewidth = 0.75, se = FALSE, span = 0.5) +

    scale_y_continuous(
        name = expression(
            bold(atop("NH"[4]^"+" * "," ~ "NO"[2]^"-" * "+" * "NO"[3]^"-" * "," ~ "PO"[4]^"3-",
                      "(" * mu * "M)"))
        ),
        sec.axis = sec_axis(
            trans = ~ . * 3,
            name = expression(bold(atop("Si(OH)"[4], "(" * mu * "M)")))
        ),
        limits = c(0, 10),
        expand = c(0, 0)
    ) +

    geom_vline(xintercept = 144, linetype = 2, linewidth = 0.75) +

    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_blank(),
        axis.text = element_text(face = "bold", size = 8, colour = "black"),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        legend.background = element_rect(linetype = "solid", colour = 'black'),
        legend.title = element_text(face = "bold", size = 8),
        legend.title.align = 0.5,
        legend.text = element_text(size = 8),
        legend.text.align = 0,
        legend.key.width = unit(2, "line"),
        legend.key.height = unit(0.5, "line"),
        legend.spacing.y = unit(0.1, "cm"),
        panel.border = element_rect(size = 1, color = "black")
    ) +

    scale_x_continuous(expand = c(0, 0), minor_breaks = seq(110, 162, 5)) +

    scale_linetype_manual(
        name = "Seawater Nutrients",
        breaks = c("l1", "l2", "l3", "l4"),
        values = lines,
        labels = c(
            expression("NH"[4]^"+"),
            expression("NO"[2]^"-" ~ "+" ~ "NO"[3]^"-"),
            expression("PO"[4]^"3-"),
            expression("Si(OH)"[4])
        )
    ) +

    scale_color_manual(
        name = "Seawater Nutrients",
        breaks = c("c1", "c2", "c3", "c4"),
        values = cols,
        labels = c(
            expression("NH"[4]^"+"),
            expression("NO"[2]^"-" ~ "+" ~ "NO"[3]^"-"),
            expression("PO"[4]^"3-"),
            expression("Si(OH)"[4])
        )
    ) +

    annotate("text", x = 113, y = 9, label = "b", size = 5)


#####################################################################################################
# Chl/Phaeo plot panel
#####################################################################################################
#####################################################
#algae
#####################################################

dweibull_density <- function(x) {
   1000 * dweibull(x, shape = 16.1979, scale = 128.2009)
 }

#This route allows for a legend when pulling from diff columns etc.
cols    <- c( "c1" = "#2C9436", "c2" = "#705E0E" )
shapes  <- c("s1" = 16, "s2" = 17)



algae=
  ggplot(Chl_avg)+ theme_few()+
  annotate("rect", xmin = 110, xmax = 144, ymin = -Inf, ymax = Inf,
           alpha = .1,fill = "lightgreen")+
      geom_vline(xintercept=c(144), linetype=2, size=0.75)+
      geom_point(aes(DOY, chl_mgm2, color="c1", shape="s1"), size=1)+
      geom_point(aes(DOY, phaeo_prop*50, color="c2", shape="s2"), size=1)+
      #stat_smooth(aes(DOY, phaeo_prop*3600), se=F, span=2,, color="#705E0E")+
      geom_line(stat = "function", fun = dweibull_density, color = "#2C9436", size = 0.75, linetype=5) +
      scale_y_continuous(name = expression(bold(atop("Chlorophyll" ~bolditalic("a"), "(mg m"^"-2"*")"))),
        sec.axis = sec_axis( trans=~.*0.02,
        name=expression(bold("Proportion Phaeophytin"))),
        breaks = seq(0, 50, 10),
                   limits = c(-0.2,50),
                   expand = c(0,0))+
      #trans value is just the  1/maxchl y axis, i.e., 1/50
      theme(#legend.position = c(0.8, 0.7),
            axis.title.x = element_blank(),
            axis.title.y= element_text(size=10, face="bold"),
            axis.text.x = element_blank(),
            axis.text = element_text(face="bold", size=8, colour="black"),
            axis.ticks=element_line(size=0.5, color="black"),
            axis.ticks.length = unit(0.25, "cm"),
            legend.background=element_rect(linetype="solid", colour = 'black'),
            legend.title=element_text(face="bold", size=8),
            legend.title.align=0.5,
            legend.text=element_text( size=8),
            legend.text.align=0,
            legend.key.width = unit(2, "line"),
            legend.key.height = unit(0.75, "line"),
            panel.border=element_rect(size=1, color="black"))+
      scale_x_continuous(expand = c(0, 0), minor_breaks= seq(110, 162, 5))+
      # coord_cartesian(ylim=c(0, 3.600), xlim=c(110,NA)) + #must use cartesian coords to set the limits so that it doesn't eliminate points that change the loess line.
      labs(x = "Day of Year")+
      scale_shape_manual(name = "Pigment",
                         breaks = c("s1", "s2"),
                         values = shapes,
                         labels = c(expression("Chl"~italic("a")), "Phaeo"))+
      scale_color_manual(name = "Pigment",
                         breaks = c("c1", "c2"),
                         values = cols,
                         labels = c(expression("Chl"~italic("a")), "Phaeo"))+
                         # geom_text(aes(label = "a"), x=111, y=47.5, size=9)
                         annotate("text", x = 113, y = 47.5, label = "a", size = 5)


#fix the tickmark overhang
g=ggplotGrob(algae)
panel_id <- g$layout$name == "panel"
border <- tail(g$grobs[panel_id][[1]][["children"]], 1)
g <- gtable::gtable_add_grob(g, border,
                             t = min(g$layout$t[panel_id]),
                             l = min(g$layout$l[panel_id]),
                             name = "border",
                             clip = "off")
algae2=ggdraw(g)
algae2



##################################################
#Third panel, physical
#################################################

#Missing snow and ice thicknesses for some.

phys_sub=Enviro_Chl %>%
select(DOY, PAR_B, Snow_Depth, Temperature, Salinity_Brine)%>%
pivot_longer( cols=-c(DOY), names_to = "Variable", values_to = "Value")

Disc=c("#26E8EB", "#3A4EC2", "#040291", "darkgray", "black")


    cols    <- c( "c1" = "Black", "c2" = "Black", "c3"="Black" )
    lines  <- c("l1" = 1, "l2" = 2, 'l3'=3)
PAR_Snow=
    ggplot(Enviro_Chl)+ theme_few()+
    annotate("rect", xmin = 110, xmax = 144, ymin = -Inf, ymax = Inf,
             alpha = .1,fill = "lightgreen")+
        geom_vline(xintercept=c(144), linetype=2, size=0.75)+
        stat_smooth(aes(DOY, PAR_B, color="c1", linetype="l1"), linewidth=0.75, se=F, span=0.65)+
        stat_smooth(aes(DOY, Snow_Depth*13.24503, color="c2", linetype="l2"), linewidth=0.75, se=F, span=0.65)+
        # geom_point(aes(DOY, PAR_B, color="c1", shape="l1"), size=2.5)+
        # geom_point(aes(DOY, Snow_Depth*13.24503, color="c2", shape="l2"), size=2.5)+
        #stat_smooth(aes(DOY, phaeo_prop*3600), se=F, span=2,, color="#705E0E")+
        scale_y_continuous(
          name = expression(
            atop(
              bold("Irradiance"),
              bold("(" * paste(mu, "mols") ~ "photons" ~ m^-2 ~ s^-1 * ")")
            )
          ),
          sec.axis = sec_axis(
            trans = ~ . / 13.24503,
            name = expression(
              atop(
                bold("Snow Depth"),
                bold("(cm)")
              )
            )
          ),
          expand = c(0, 0),
          limits = c(0, 80)
        ) +
        theme(#legend.position = c(0.8, 0.8),
              # axis.title.x = element_blank(),
              axis.title.y= element_text(size=10, face="bold"),
              # axis.text.x = element_blank(),
              axis.text = element_text(face="bold", size=8, colour="black"),
              axis.ticks=element_line(size=0.5, color="black"),
              axis.ticks.length = unit(0.25, "cm"),
              axis.title = element_text(size=10, face="bold"),
              strip.text.x = element_text(size = 10),
              legend.background=element_rect(linetype="solid", colour = 'black'),
              legend.title=element_text(face="bold", size=8),
              legend.title.align=0.5,
              legend.text=element_text( size=8),
              legend.text.align=0,
              legend.key.width = unit(2, "line"),
              panel.border=element_rect(size=1, color="black"))+
        scale_x_continuous(expand = c(0, 0), minor_breaks= seq(110, 162, 5))+
        # coord_cartesian(ylim=c(0, 3.600), xlim=c(110,NA)) + #must use cartesian coords to set the limits so that it doesn't eliminate points that change the loess line.
        labs(x = "Day of Year")+
        scale_linetype_manual(name = "Variable",
                           breaks = c("l1", "l2"),
                           values = lines,
                           labels = c(expression("PAR"),expression("Snow")))+
        scale_color_manual(name = "Variable",
                           breaks = c("c1", "c2"),
                           values = cols,
                           labels = c(expression("PAR"),expression("Snow")))+
                           # geom_text(aes(label = "c"), x=111, y=76, size=9)
                           annotate("text", x = 113, y = 74, label = "c", size = 5)

combo_enviro_sw=algae/nuts_sw_mulitipanel_notation_nut_axes/PAR_Snow
