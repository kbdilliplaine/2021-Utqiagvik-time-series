library(ggplot2)
library(dplyr)
library(openxlsx)
library(tidyverse)
library(tidyr)
library(ggpubr)

chlnut=readRDS(file = './Enviro_New_Bottom10s_2024.RDS')

subnut=chlnut %>% select(DOY, chl_mgm2, Nitrate_Bulk, Nitrite_Bulk, Ammonia_Bulk, Nitrogen_Bulk, Phosphate_Bulk, Silicate_Bulk) %>%
  column_to_rownames('DOY')


# Get the column names of the dataframe except for "chl_mgm2"
predictors <- setdiff(names(subnut), "chl_mgm2")

# Run linear regression for each predictor against "chl_mgm2"
regressions <- lapply(predictors, function(predictor) {
  formula <- formula(paste(predictor, " ~ chl_mgm2"))
  lm(formula, data = subnut)
})

  summaries <- lapply(regressions, summary)
  print(summaries)

    adj_r2=chlnut %>%
      select(DOY, chl_mgm2, Nitrate_Bulk, Nitrite_Bulk, Ammonia_Bulk, Nitrogen_Bulk, Phosphate_Bulk, Silicate_Bulk) %>%
      pivot_longer(-c(DOY, chl_mgm2), names_to = "Variable", values_to = "measure") %>%
      group_by(Variable) %>%
      mutate(adj_r2 = summary(lm(measure ~ chl_mgm2))$adj.r.squared) %>% filter(DOY==111)%>%
      mutate(Variable_Label = factor(Variable, levels = c(
       "Ammonia_Bulk", "Nitrite_Bulk", "Nitrate_Bulk", "Phosphate_Bulk", "Silicate_Bulk", "Nitrogen_Bulk"),
       labels = c(
           expression("NH"[4]^"+"),
           expression("NO"[2]^"-"),
           expression("NO"[3]^"-"),
           expression("PO"[4]^"3-"),
           expression("Si(OH)"[4]),
           expression("NO"[2]^"-"~"+"~"NO"[3]^"-")  # Plain text for Nitrogen_Bulk
       )
   ))

      chlnut_with_labels <- chlnut %>%
select(DOY, chl_mgm2, Nitrate_Bulk, Nitrite_Bulk, Ammonia_Bulk, Nitrogen_Bulk, Phosphate_Bulk, Silicate_Bulk) %>%
        pivot_longer(-c(DOY, chl_mgm2), names_to = "Variable", values_to = "measure") %>%
        mutate(Variable_Label = factor(Variable, levels = c(
          "Ammonia_Bulk", "Nitrite_Bulk", "Nitrate_Bulk", "Phosphate_Bulk", "Silicate_Bulk", "Nitrogen_Bulk"),
          labels = c(
            expression("NH"[4]^"+"),
            expression("NO"[2]^"-"),
            expression("NO"[3]^"-"),
            expression("PO"[4]^"3-"),
            expression("Si(OH)"[4]),
            expression("NO"[2]^"-"~"+"~"NO"[3]^"-")  # Plain text for Nitrogen_Bulk
          )))

    linreg=
        ggplot(chlnut_with_labels, aes(x = chl_mgm2, y = measure)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE) +
        theme_few() +
        facet_wrap(~Variable_Label, scales = "free_y", labeller = label_parsed) +  # Apply parsed labels
        labs(x = expression("Chlorophyll" ~italic("a")  ~"(mg m"^"-2"*")"), y = expression(paste("Nutrients ("*mu, "M)"))) +
        geom_text(data=adj_r2, aes(label = paste0("adj. r² = ", round(adj_r2, 2))), color = "black", hjust = "inward", vjust = "top")+
        theme(axis.title = element_text(size = 17, face = "bold"),
              axis.text=element_text(color="black"))


##########
#seawater
##########
seawater_nuts <- read.csv("./Seawater_nutrients.csv")

####################
# N:P and Si:N Redfield
####################

NP=
ggplot(seawater_nuts) +
    geom_point(aes(DOY, NP_Ratio)) +
    geom_line(aes(DOY, NP_Ratio)) +
    geom_hline(yintercept = 16, color = "red")+
        labs(x = "Day of year", y = "Nitrate:phosphate") +
        guides(color = guide_legend(title = "Source"))+
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        guides(color = guide_legend(title = "Source"))+
        theme(axis.title = element_text(face="bold", size=12),
              plot.title = element_text(hjust = 0.5, size=15),
              strip.text.x =element_text(size = 10, face="bold"),
              #strip.text.y = element_text(size = 15, face="bold.italic"),
              axis.text = element_text(face="bold", size=10, colour = "black"),
              legend.text = element_text( face = "bold", size = 10),
           legend.title = element_blank(),
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              # legend.position="bottom",
              # panel.spacing = unit(2, "lines"),
              # strip.text = element_text(hjust = -0),
              # strip.text.x = element_text(size = 20, face="bold", hjust=0),
              panel.background = element_blank(),
              panel.grid = element_blank(),
              # panel.border = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              # axis.title.x = element_text(vjust=-0.5),
              # plot.margin = margin(2,2,2,2)
)

SiN=
ggplot(seawater_nuts) +
    geom_point(aes(DOY, SiN_Ratio)) +
    geom_line(aes(DOY, SiN_Ratio)) +
    geom_hline(yintercept = 0.9375, color = "red")+
            labs( x = "Day of year", y = "Silicate:nitrate") +
            guides(color = guide_legend(title = "Source"))+
            guides(color = guide_legend(title = "Source"))+
            theme(axis.title = element_text(face="bold", size=12),
                  plot.title = element_text(hjust = 0.5, size=15),
                  strip.text.x =element_text(size = 10, face="bold"),
                  #strip.text.y = element_text(size = 15, face="bold.italic"),
                  axis.text = element_text(face="bold", size=10, colour = "black"),
                  legend.text = element_text( face = "bold", size = 10),
               legend.title = element_blank(),
                  # legend.position="bottom",
                  # panel.spacing = unit(2, "lines"),
                  # strip.text = element_text(hjust = -0),
                  # strip.text.x = element_text(size = 20, face="bold", hjust=0),
                  panel.background = element_blank(),
                  panel.grid = element_blank(),
                  # panel.border = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=1),
                  # axis.title.x = element_text(vjust=-0.5),
                  # plot.margin = margin(2,2,2,2)
            )




Combo_Redfield=ggarrange(NP,SiN, ncol=1, labels=c("a", "b"), align="v")
