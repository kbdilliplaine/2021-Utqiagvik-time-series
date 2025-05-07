require(reshape)
require(vegan)
require(indicspecies)
require(cluster)
require(plyr)
library(tidyverse)
library(readr)
library(ggthemes)
library(stringr)
library(patchwork)

#Load in the RDS file containing Enviro_long Enviro_Chl, Chl_avg, sw_nuts_short and sw_nuts_trimmed, as processed in 04_Enviro_Plots_24.r
load("./Enivro_nuts_all.RDS")

Brine_nuts=
    Enviro_long %>%
    filter(Variable %in% c("Ammonia_Brine", "Nitrate_Brine", "Nitrite_Brine", "Phosphate_Brine", "Silicate_Brine", "Nitrogen_Brine"))%>%
    pivot_wider(
      names_from = Variable,
      values_from = Value
    )

##################################################
#Brine nut panel
#################################################

    cols    <- c( "c1" = "#26E8EB", "c2" = "#3A4EC2", "c3"="#040291", "c4"="black" )
    lines  <- c("l1" = 1, "l2" = 2, 'l3'=3, 'l4'=4)


nuts_Brine_mulitipanel_notation_nut_axes =
ggplot(Brine_nuts) +
  theme_few() +
  annotate("rect", xmin = 110, xmax = 144, ymin = -Inf, ymax = Inf, alpha = .1, fill = "lightgreen") +
  stat_smooth(aes(DOY, Ammonia_Brine, color = "c1", linetype = "l1"), linewidth = 1.5, se = F, span = 0.5) +
  stat_smooth(aes(DOY, Nitrogen_Brine, color = "c2", linetype = "l2"), linewidth = 1.5, se = F, span = 0.5) +
  stat_smooth(aes(DOY, Phosphate_Brine, color = "c3", linetype = "l3"), linewidth = 1.5, se = F, span = 0.5) +
  stat_smooth(aes(DOY, Silicate_Brine / 2, color = "c4", linetype = "l4"), linewidth = 1.5, se = F, span = 0.5) +
  scale_y_continuous(
    name = expression("NH"[4]^"+" ~ "," ~ "NO"[2]^"-" ~ "+" ~ "NO"[3]^"-" ~ "," ~ "PO"[4]^"3-" ~ "," ~ "(" * mu * "M)"),
    sec.axis = sec_axis(~ .*2, name = expression("Si(OH)"[4] ~ "(" * mu * "M)"))
  ) +
  coord_cartesian(ylim = c(0, 45), xlim = c(110, NA)) +
  geom_vline(xintercept = c(144), linetype = 2, size = 0.75) +
  theme(
    legend.position = c(0.8, 0.7),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 17, face = "bold"),
    axis.text.x = element_blank(),
    axis.text = element_text(face = "bold", size = 15, colour = "black"),
    axis.ticks = element_line(size = 1, color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    legend.background = element_rect(linetype = "solid", colour = 'black'),
    legend.title = element_text(face = "bold", size = 12),
    legend.title.align = 0.5,
    legend.text = element_text(size = 12),
    legend.text.align = 0,
    legend.key.width = unit(3, "line"),
    panel.border = element_rect(size = 1, color = "black")
  ) +
  scale_x_continuous(expand = c(0, 0), minor_breaks = seq(110, 162, 5)) +
  scale_linetype_manual(
    name = "Brine Nutrients",
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
    name = "Brine Nutrients",
    breaks = c("c1", "c2", "c3", "c4"),
    values = cols,
    labels = c(
      expression("NH"[4]^"+"),
      expression("NO"[2]^"-" ~ "+" ~ "NO"[3]^"-"),
      expression("PO"[4]^"3-"),
      expression("Si(OH)"[4])
    )
  ) +
  annotate("text", label = "b", x = 111, y = 42.75, size = 9)



#####################################################
#algae
#####################################################

mu <-  4.839e+00
sigma <- 5.564e-02
log_normal_density <- function(x) {
# 6.080e+04  * dlnorm(x, meanlog = mu, sdlog = sigma)} #this for ug/L below for mg/L
60.79755* dlnorm(x, meanlog = mu, sdlog = sigma)}
cols    <- c( "c1" = "#2C9436", "c2" = "#705E0E" )
shapes  <- c("s1" = 16, "s2" = 17)

##################################################
#Panel 1 chl.
#################################################
# mg L Brine
algae=
ggplot(Enviro_Chl)+ theme_few()+
annotate("rect", xmin = 110, xmax = 144, ymin = -Inf, ymax = Inf,
         alpha = .1,fill = "lightgreen")+
    geom_vline(xintercept=c(144), linetype=2, size=0.75)+
    geom_point(aes(DOY, chl_Lbrine/1000, color="c1", shape="s1"), size=3)+
    geom_point(aes(DOY, phaeo_prop*3.600, color="c2", shape="s2"), size=3)+
    #stat_smooth(aes(DOY, phaeo_prop*3600), se=F, span=2,, color="#705E0E")+
    geom_line(stat = "function", fun = log_normal_density, color = "#2C9436", size = 1.5, linetype=5) +
    scale_y_continuous(name = expression("Chlorophyll" ~italic("a")  ~"(mg L"^"-1"*"Brine)"), sec.axis = sec_axis( trans=~.*0.2777778, name=expression("Proportion Phaeophytin")), expand = c(0, 0))+
    theme(legend.position = c(0.8, 0.7),
          axis.title.x = element_blank(),
          axis.title.y= element_text(size=17, face="bold"),
          axis.text.x = element_blank(),
          axis.text = element_text(face="bold", size=15, colour="black"),
          axis.ticks=element_line(size=1, color="black"),
          axis.ticks.length = unit(0.25, "cm"),
          legend.background=element_rect(linetype="solid", colour = 'black'),
          legend.title=element_text(face="bold", size=12),
          legend.title.align=0.5,
          legend.text=element_text( size=12),
          legend.text.align=0,
          legend.key.width = unit(3, "line"),
          legend.key.height = unit(0.75, "line"),
          panel.border=element_rect(size=1, color="black"))+
    scale_x_continuous(expand = c(0, 0), minor_breaks= seq(110, 162, 5))+
    coord_cartesian(ylim=c(0, 3.600), xlim=c(110,NA)) + #must use cartesian coords to set the limits so that it doesn't eliminate points that change the loess line.
    labs(x = "Day of Year")+
    scale_shape_manual(name = "Pigment",
                       breaks = c("s1", "s2"),
                       values = shapes,
                       labels = c(expression("Chl"~italic("a")), "Phaeo"))+
    scale_color_manual(name = "Pigment",
                       breaks = c("c1", "c2"),
                       values = cols,
                       labels = c(expression("Chl"~italic("a")), "Phaeo"))+
                       annotate("text", label = "a", x = 111, y = 3.42, size = 9)

 g=ggplotGrob(algae)
 panel_id <- g$layout$name == "panel"
 border <- tail(g$grobs[panel_id][[1]][["children"]], 1)
 g <- gtable::gtable_add_grob(g, border,
                              t = min(g$layout$t[panel_id]),
                              l = min(g$layout$l[panel_id]),
                              name = "border",
                              clip = "off")


##################################################
#Third panel, physical
#################################################
Brine_S=
    Enviro_long %>%
    filter(Variable =="Salinity_Brine")%>%
    pivot_wider(
      names_from = Variable,
      values_from = Value
    )

Brine_Salinity=
    ggplot(Brine_S)+ theme_few()+
    annotate("rect", xmin = 110, xmax = 144, ymin = -Inf, ymax = Inf,
             alpha = .1,fill = "lightgreen")+
        geom_vline(xintercept=c(144), linetype=2, size=0.75)+
        stat_smooth(aes(DOY, Salinity_Brine), color="Black", linetype=1, linewidth=1.5, se=F, span=0.65)+
        geom_point(aes(DOY, Salinity_Brine), color="Black", size=2.5)+
        # geom_point(aes(DOY, Snow_Depth*13.24503, color="c2", shape="l2"), size=2.5)+
        #stat_smooth(aes(DOY, phaeo_prop*3600), se=F, span=2,, color="#705E0E")+
         scale_y_continuous(name = "Salinity")+
        theme(legend.position = c(0.8, 0.8),
        axis.title.y= element_text(size=17, face="bold"),
        # axis.text.x = element_blank(),
        axis.text = element_text(face="bold", size=15, colour="black"),
        axis.ticks=element_line(size=1, color="black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.title.x = element_text(size=20, face="bold"),
        strip.text.x = element_text(size = 15),
        legend.background=element_rect(linetype="solid", colour = 'black'),
        legend.title=element_text(face="bold", size=12),
        legend.title.align=0.5,
        legend.text=element_text( size=12),
        legend.text.align=0,
        legend.key.width = unit(3, "line"),
        panel.border=element_rect(size=1, color="black"))+
        scale_x_continuous(expand = c(0, 0), minor_breaks= seq(110, 162, 5))+
        coord_cartesian(ylim = c(0, 55), xlim = c(110, NA)) +
        # coord_cartesian(ylim=c(0, 3.600), xlim=c(110,NA)) + #must use cartesian coords to set the limits so that it doesn't eliminate points that change the loess line.
        labs(x = "Day of Year")+
        annotate("text", label = "c", x = 111, y = 52.25, size = 9)





#Brine nuts with brine chl and brine salinity
combo_enviro_Brine=algae/nuts_Brine_mulitipanel_notation_nut_axes/Brine_Salinity
