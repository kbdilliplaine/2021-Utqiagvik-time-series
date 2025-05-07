library(tidyverse)
library(lubridate)
library(ggpubr)

# File paths
SIMB_fp <- './SIMB3_278860.csv'
PABR_fp <- './PABRweather/PABR2.csv'

# Start and end dates
date_start <- ymd("2021-04-21")
date_end <- ymd("2021-06-11")

# Read data
#Set and convert excel date.
data_SIMB <- read_csv(SIMB_fp) %>%
    mutate(datetime=as.POSIXct(time_stamp * (60*60*24)
               , origin="1899-12-30"
               , tz="America/Anchorage"))%>%
  filter(datetime >= date_start, datetime <= date_end)

data_PABR <- read_csv(PABR_fp) %>%
  mutate(datetime = ymd_hms(DATE, tz="America/Anchorage")) %>%
  filter(datetime >= date_start, datetime <= date_end, !is.na(HourlyDryBulbTemperature)) %>%
  mutate(air_temp = (HourlyDryBulbTemperature - 32) * 5.0/9.0)


AirT=
ggplot() +
    geom_line(data = data_SIMB, aes(x = decimal_day, y = air_temp, color = "SIMB3 Buoy"), size = 1, linetype = "solid", alpha = 0.7) +
    geom_line(data = data_PABR, aes(x = decimal_day, y = air_temp, color = "Airport"), size = 1, linetype = "solid", alpha = 0.7) +
    geom_hline(yintercept = 0, color = "black", linetype="dashed") +
    geom_hline(yintercept = -1.8, color = "red")+
    labs( x = "Date", y = "Air temperature (°C)") +
    # scale_x_datetime(date_breaks = "1 week", date_labels = "%b %d") +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_manual(name = "Source", values = c("SIMB3 Buoy" = "blue", "Airport" = "red")) +
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
          legend.position=c(.8,.3),
          # axis.title.y=element_blank(),
          # legend.position="bottom",
          # panel.spacing = unit(2, "lines"),
          # strip.text = element_text(hjust = -0),
          # strip.text.x = element_text(size = 20, face="bold", hjust=0),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          # panel.border = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          # axis.title.x = element_text(vjust=-0.5),
          plot.margin = margin(2, 5, 0, 20),  # Add more space on the left side
          panel.spacing = unit(0.5, "lines"), # Reduce the spacing between panels
          plot.tag.position = c(0, 0.95)   # Move tag label slightly more left
    )

temp=readRDS( file = './utq_temperature.RDS')

Icesnow=
temp %>% group_by(DOY) %>%
  filter(Type=="Ice-snow")


  icesnowT=ggplot(data = Icesnow) +
    geom_point(aes(DOY, T)) +
    geom_line(aes(DOY, T))+
    geom_hline(yintercept = -1.8, color = "red")+
    geom_hline(yintercept = 0, color = "black", linetype="dashed")+
      labs( x = "Day of year", y = "Ice-snow temperature (°C)") +
      theme(axis.title = element_text(face="bold", size=12),
            plot.title = element_text(hjust = 0.5, size=15),
            strip.text.x =element_text(size = 10, face="bold"),
            axis.text = element_text(face="bold", size=10, colour = "black"),
            legend.text = element_text( face = "bold", size = 10),
         legend.title = element_blank(),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            plot.margin = margin(0, 5, 0, 20),  # Add more space on the left side
            panel.spacing = unit(0.5, "lines"), # Reduce the spacing between panels
            plot.tag.position = c(0, 0.95)   # Move tag label slightly more left
      )


  Snowair=
  temp %>% group_by(DOY) %>%
    filter(Type=="Snow-Air")

    SnowairT=
    ggplot(data = Snowair) +
        geom_point(aes(DOY, T)) +
        geom_line(aes(DOY, T))+
        geom_hline(yintercept = -1.8, color = "red")+
        geom_hline(yintercept = 0, color = "black", linetype="dashed")+
        labs( x = "Day of year", y = "Snow-air temperature (°C)") +
        guides(color = guide_legend(title = "Source"))+
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        guides(color = guide_legend(title = "Source"))+
        theme(axis.title = element_text(face="bold", size=12),
              plot.title = element_text(hjust = 0.5, size=15),
              strip.text.x =element_text(size = 10, face="bold"),
              axis.text = element_text(face="bold", size=10, colour = "black"),
              legend.text = element_text( face = "bold", size = 10),
           legend.title = element_blank(),
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              panel.background = element_blank(),
              panel.grid = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              plot.margin = margin(0, 5, 0, 20),  # Add more space on the left side
              panel.spacing = unit(0.5, "lines"), # Reduce the spacing between panels
              plot.tag.position = c(0, 0.95)   # Move tag label slightly more left
        )


botts=temp %>% group_by(DOY) %>%
          filter(Type=="Ice" & Depth=="0.025")

bott2.5=botts %>% group_by(DOY) %>%
      dplyr::summarize(temperature.avg=mean(T))



    Icebott=
        ggplot(data = bott2.5) +
          geom_point(aes(DOY, temperature.avg)) +
          geom_line(aes(DOY, temperature.avg))+
          geom_hline(yintercept = -1.8, color = "red")+
          labs( x = "Day of year", y = "Bottom-ice temperature (°C)") +
          guides(color = guide_legend(title = "Source"))+
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          guides(color = guide_legend(title = "Source"))+
          theme(axis.title = element_text(face="bold", size=12),
                plot.title = element_text(hjust = 0.5, size=15),
                strip.text.x =element_text(size = 10, face="bold"),
                axis.text = element_text(face="bold", size=10, colour = "black"),
                legend.text = element_text( face = "bold", size = 10),
                legend.title = element_blank(),
                panel.background = element_blank(),
                panel.grid = element_blank(),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                plot.margin = margin(0, 5, 2, 20),  # Add more space on the left side
                panel.spacing = unit(0.5, "lines"), # Reduce the spacing between panels
                plot.tag.position = c(0, 0.95)   # Move tag label slightly more left
          )


ComboTemp=ggarrange(AirT,SnowairT,icesnowT,Icebott, ncol=1, align="hv", labels=c("a", "b", "c", "d"))
