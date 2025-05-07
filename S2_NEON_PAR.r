#Neonutilities tutorial https://www.neonscience.org/resources/learning-hub/tutorials/download-explore-neon-data
# load packages
library(neonUtilities)
library(raster)
library(readxl)
library(dplyr)
library(tidyverse)
library(openxlsx)
library(lubridate)
library(ggthemes)
library(egg)


# Set global option to NOT convert all character variables to factors
options(stringsAsFactors=F)

#################################################33
# After doing this once, dont need to do it again.
# The stackByTable() function will unzip and join the files in the downloaded zip file.
#stackByTable("C:/Users/Dilliplaine/Desktop/MOSAiC/Data/NEON_par")
#####################################################33
#Read in the stacked files
par30 <- readTableNEON(
  dataFile="./PARPAR_30min.csv",
  varFile="./variables_00024.csv")
head(par30)

#The variables file shows you the definition and units for each column of data.
parvar <- read.csv("./variables_00024.csv")
View(parvar)

plot(PARMean~startDateTime,
     data=par30[which(par30$verticalPosition=="010"),],
     type="l")

     #Cleaning data and subsetting https://www.neonscience.org/resources/learning-hub/tutorials/dc-subset-data-no-data-values-r
     library(lubridate)
     library(ggplot2)


     # convert to POSIX date time class
par30$datetime <- as.POSIXct(par30$startDateTime,
                                format = "%Y-%m-%dT%H:%M",
                                tz = "America/Anchorage")
# subset data - 2009-2011
 par30.2021 <- subset(par30,
      datetime >= as.POSIXct('2021-04-21 00:00') &
      datetime <= as.POSIXct('2021-06-11 23:30'))

  # plot(PARMean~startDateTime,
  #            data=par30.2021%>%filter(verticalPosition=="030"),
  #            type="l")
#################
#Convert to decimal date; add DOY
hrs <- hours(12)
par30.2021$adjTime=par30.2021$startDateTime-hrs

par30.2021$date=decimal_date(par30.2021$adjTime)
par30.2021$DOY=lubridate::yday(par30.2021$adjTime)
plot(PARMean~date,
  data=par30.2021%>%filter(verticalPosition=="030"),
  type="l")
########


#subset and extract just the vertposition 030 PAR data for this specific segment of 1 year with decimal date for the phytotools
par30.2021_2=par30.2021 %>% filter(verticalPosition== "030", PARMean >0) %>% dplyr::select(startDateTime, adjTime, endDateTime, date, DOY, PARMean )
par30.2021_2=na.omit(par30.2021_2)

###########
#Load ice and snow thickness
###########
DF=read_excel("./MO_Snow_Icethick_Bot10S_2.xlsx", sheet = 1)
DF$Date=as_datetime(DF$Date, tz = "UTC")
DF$DOY=yday(DF$Date) #Date is off by one.

DF$DOY=as.numeric(DF$DOY)
Snow_Thickness=DF %>%
filter(Variable=="Snow") %>%
group_by(DOY) %>%
summarize(avg_snowthick=mean(Value))

interpolated_values <- approx(x = Snow_Thickness$DOY, y = Snow_Thickness$avg_snowthick, xout = par30.2021_2$DOY)$y

par30.2021_2$snow_depth=interpolated_values

# ggplot(par30.2021_2)+geom_line(aes(x=startDateTime, y=PARMean))+geom_line(aes(x=startDateTime, y=snow_depth*100), col="blue")


Ice_Thickness=DF %>%
filter(Variable=="Ice" ) %>%
group_by(DOY) %>%
summarize(avg_icethick=mean(Value))

ice_depth_interp <- approx(x = Ice_Thickness$DOY, y = Ice_Thickness$avg_icethick, xout = par30.2021_2$DOY)$y
par30.2021_2$ice_depth=ice_depth_interp

#Underice PAR based on attenuation coefficients derived from A functional regression model for predicting optical depth and estimating attenuation coefficients in sea-ice covers near Resolute Passage, Canada

#Set attenuation coefficients
ui= -2.6#-3.5*.7 #Matthes 2019 change this to -2.6
us=-11
uchl=-0.06 #Also from sea-ice covers near Resolute Passage, Canada
#The paper itself states 11 for snow and 3.6 for ice.
#units should be in meters
par30.2021_2$UI_PAR=round(par30.2021_2$PARMean*exp(((par30.2021_2$ice_depth/100)*ui)+(par30.2021_2$snow_depth/100*us)),2)


##################################################
#Load in field PAR data for pseudotransmissivity
##################################################
Light=read_excel("./PAR2.xlsx", sheet = 1)
library(lubridate)
Light$Date=as_datetime(Light$Date, tz = "UTC")
Light$DOY=yday(Light$Date)+1 #Date is off by one.
Light$PseudoTrans=(Light$PAR_B/Light$PAR_A)*100
Light2=Light %>% select(-c("hi","mean snow", "std. snow"))

ggplot(Light2%>%select(c("DOY","PAR_A","PAR_B")))+
geom_point(aes(x=DOY, y=PAR_A))+
geom_point(aes(x=DOY, y=PAR_B))



Pseudotrans_interp <- approx(x = Light2$DOY, y = Light2$PseudoTrans, xout = par30.2021_2$DOY)$y
par30.2021_2$PseudoTrans=Pseudotrans_interp/100
par30.2021_2$PseudoTrans_PAR=round(par30.2021_2$PARMean*par30.2021_2$PseudoTrans,2)


# ggplot(par30.2021_2)+
# geom_line(aes(x=startDateTime, y=PARMean))+
# geom_line(aes(x=startDateTime, y=snow_depth*100), col="blue")+
# geom_line(aes(x=startDateTime, y=ice_depth*10), col="red")+
# geom_line(aes(x=startDateTime, y=PseudoTrans_PAR*10), col="orange")+
# geom_line(aes(x=startDateTime, y=UI_PAR*10), col="purple")


Chl_avg=readRDS(file = './Full_Chlorophyll_avgd.RDS')
Chl_mgm2= Chl_avg%>%
filter(Location=="TM" & Chl_Flag != 0)

chl_interp <- approx(x = Chl_mgm2$DOY, y = Chl_mgm2$chl_mgm2, xout = par30.2021_2$DOY)$y
par30.2021_2$chl=chl_interp

par30.2021_2$UI_PAR_chl=round(par30.2021_2$PARMean*exp(((par30.2021_2$ice_depth/100)*ui)+((par30.2021_2$snow_depth/100)*us)+(par30.2021_2$chl*uchl)),2)


# ggplot(par30.2021_2)+
#     geom_line(aes(x=startDateTime, y=PseudoTrans_PAR), col="black", size=1)+
#     geom_line(aes(x=startDateTime, y=UI_PAR_chl), col="blue", size=1)


#Daily Light Integral Calculation
#n()=measure per day (30 min intervals); 30=mins and 60=seconds. the 1Million is the conversion of umols to mols.
DLI_sum=par30.2021_2 %>% group_by(DOY) %>% summarize(DLI_pseudo=sum((30*60*PseudoTrans_PAR)/1000000), DLI_atten=sum((30*60*UI_PAR_chl)/1000000))%>% filter(DOY !=110)
DLI_sum <- DLI_sum %>%
    rename(  "DLI_PT" = "DLI_pseudo",
             "DLI_Atten" = "DLI_atten")
DLI_sum$DOY=as.numeric(DLI_sum$DOY)

 par30.2021_2=par30.2021_2 %>% group_by(DOY) %>% mutate(DLI_PT=sum((30*60*PseudoTrans_PAR)/1000000), DLI_Atten=sum((30*60*UI_PAR_chl)/1000000))

# write.csv(DLI_sum, "C:/Users/Dilliplaine/Desktop/MOSAiC/mosaic/DLI.csv")
#End result is Mols/m2/d

PAR30=par30.2021_2 %>% select(startDateTime, endDateTime, adjTime, DOY, PARMean, snow_depth, ice_depth, PseudoTrans, PseudoTrans_PAR, chl, UI_PAR_chl, UI_PAR) %>% filter(DOY !=110)
PAR30 <- PAR30 %>%
        rename("startDateTime_UTC"="startDateTime",
              "endDateTime_UTC"="endDateTime",
              "adjTime_AK"="adjTime",
               "Incident PAR"="PARMean",
              "Snow depth"="snow_depth",
              "Ice depth"="ice_depth",
              "PT"="PseudoTrans",
              "PT PAR"="PseudoTrans_PAR",
              "Chl a" = "chl",
            "Atten. PAR chl"="UI_PAR_chl",
          "Atten. PAR"="UI_PAR")


PAR30_LONG=PAR30 %>%pivot_longer(!c("startDateTime_UTC", "endDateTime_UTC", "adjTime_AK", "DOY"), names_to = "Variable", values_to ="measure" )
PAR30_LONG$Variable <- factor(PAR30_LONG$Variable,
                               levels = c("Incident PAR",
                                         "Snow depth",
                                         "Ice depth",
                                         "PT",
                                         "PT PAR",
                                         "Chl a",
                                       "Atten. PAR chl",
                                     "Atten. PAR"))

########################################
#Supplemental Irradiance figure.
########################################


    # Define the sequence of breaks for 2021
    breaks <- seq(from = as.POSIXct("2021-04-20"),  # April 20 is day 110 in 2021
                  to = as.POSIXct("2021-06-09"),  # June 9 is day 160 in 2021
                  by = "10 days")
    # y_labels <- c("Chl_a_interp" = "8888",
    #               "Var2" = "Label for Var2",
    #               "Var3" = "Label for Var3",
    #             "Var4"="Label for Var3")
                library(ggplot2)
                library(lubridate)
                library(ggthemes)

                # Ensure that 'Tags' dataframe has the right structure
variable_levels <- c("Incident PAR",
          "Snow depth",
          "Ice depth",
          "PT",
          "PT PAR",
          "Chl a",
        "Atten. PAR chl",
      "Atten. PAR")
# c("Incident_PAR_Mean", "ice_depth_interp",
#       "snow_depth_interp", "Chl_a_interp",
#       "PT_interp", "PT_PAR", "Atten_PAR",
#       "Atten_PAR_wochl")
Tags <- data.frame(
  Variable = factor(variable_levels, levels = variable_levels),
  label = c("a", "b", "c", "d", "e", "f", "g", "h"),
  x = as.POSIXct("2021-04-20 20:00:00")  # Position for the text
)


#Custom second y axis units
y_units=c(expression(paste(mu, "mols")~"photons"~m^"-2"~s^"-1"*")"),
"cm",
"cm",
"Transmissivity",
expression(paste(mu, "mols")~"photons"~m^"-2"~s^"-1"*")"),
expression("mg m"^"-2"),
expression(paste(mu, "mols")~"photons"~m^"-2"~s^"-1"*")"),
expression(paste(mu, "mols")~"photons"~m^"-2"~s^"-1"*")"))


  # Define breaks
  breaks <- seq(from = as.POSIXct("2021-04-20"), to = as.POSIXct("2021-06-09"), by = "10 days")
  lim=c(as.POSIXct("2021-04-20", "2021-06-09"))

  # Define the limits for the x-axis: 110 to 170 days of the year
  start_date <- as.POSIXct("2021-04-20 00:00:00")  # Day 110
  end_date <- as.POSIXct("2021-06-14 23:59:59")    # Day 165

  # Updated plot with x-axis limits and no extra space
PAR=
  ggplot(PAR30_LONG) +
    geom_line(aes(x = adjTime_AK, y = measure)) +
    theme_few() +
    labs(x = "Day of year", y = "") +
    theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
          axis.text = element_text(colour = "black", size = 15),
          axis.title = element_text(size = 20),
          panel.spacing = unit(1, "lines"),  # Adjust space between facets
          plot.margin = margin(10, 10, 10, 10)) +  # Control plot margins
    scale_x_datetime(labels = function(x) yday(x),
                     breaks = breaks,  # Use custom breaks
                     limits = c(start_date, end_date),
                     expand = c(0, 0)) +  # Remove extra space on x-axis
    facet_grid(rows = vars(Variable), scales = "free") +
    geom_text(data = Tags, aes(x = x, y = Inf, label = label), vjust = 1.5)  # Adjust text
