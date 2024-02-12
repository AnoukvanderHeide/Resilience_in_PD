### GENERAL PLOTS ###

# This script loads cleaned data files, then makes histograms and time plots for participant characteristics and
# changes in variables over time during the survey period.

rm(list = ls())

# Load wide and long version of data file including covid survey data combined with variables from POM visits
library(readr)

Covid <- read_csv("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/Covid_averages.csv")
SR_data <- read_csv("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models/SR_classes_mean.csv",  # Res_classes_mean is made with Latent Class Trajectory script
                     col_types = cols(ID = col_skip(),
                                      SR_spl2 = col_factor(levels = c("1","2")), SR_spl2ns = col_factor(levels = c("1","2")), 
                                      SR_bet2 = col_factor(levels = c("1","2")), SR_bet2ns = col_factor(levels = c("1","2")), 
                                      SR_qua2 = col_factor(levels = c("1","2"))))
names(SR_data)[names(SR_data) == "SR_spl2"] <- "SR_class"
POM_wide <- read_csv("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/POM_precovid_visit.csv",
                     col_types = cols(Sex = col_factor(levels = c("Male", "Female")), 
                                      # 1=male, 2=female
                                      LivingSituat = col_factor(levels = c("1", "2", "3", "4", "5")), 
                                      # 1=alone, 2=with partner, 3=with family, 4=nursing home, 5=other
                                      DailyActivity = col_factor(levels = c("1", "2", "3", "4", "5", "6", "7", "8"))))
                                      # 1=paid job, 2=household, 3=retired, 4=student, 5=no paid job due to illness, 
                                      # 6=involuntarily no job, 7=volunteer work, 8=other)
Covid_wide <- cbind(Covid, Res_data[1])

Data_long <-  read_csv("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/Covid_POM_long.csv")

#------------------------------------------------------------------------------------------#
#                                 Basic histograms                                         #
#------------------------------------------------------------------------------------------#

library(ggplot2)
library(extrafont)
library(scales) # (Vaal groen/blauw = #5F9EA0, turquoise = #008B8B)
# Color blind red and green combination: red = #CC3311, green = #009988

# Histograms of longitudinal variables from covid survey
Res_hist <- ggplot(Covid_wide, aes(x=SR)) + 
                   geom_histogram(binwidth=1, color = "black", fill="#008B8B") +
                   labs(x="Stressor-reactivity", y="Count") + 
                   theme_minimal(base_size = 17)

PSS_hist <- ggplot(Covid_wide, aes(x=PSSmean)) + 
                   geom_histogram(binwidth=1, color = "black", fill="#008B8B") +
                   labs(x="Perceived Stress Scale", y="Frequency") + 
                   theme_minimal(base_size = 17)
SL_hist <- ggplot(Data_wide, aes(x=SLmean)) + 
                   geom_histogram(bins=25, color = "black", fill="#008B8B") +
                   labs(x="Stressor load", y="Frequency") + 
                   scale_x_continuous(breaks = seq(0,70,10))+
                   theme_minimal(base_size = 17)

# Histograms of disease variables
UPDRS3_hist <- ggplot(POM_wide, aes(x=UPDRS3_off)) + 
                      geom_histogram(bins=20, color = "black", fill="#008B8B") +
                      labs(x="UPDRS-III (off medication)", y="Frequency") + 
                      theme_minimal(base_size = 17)
Disdur_hist <- ggplot(POM_wide, aes(x=DisDur)) + 
                      geom_histogram(bins=20, color = "black", fill="#008B8B") +
                      labs(x="PD disease duration (years)", y="Frequency") + 
                      scale_x_continuous(breaks = seq(0,8,1))+
                      theme_minimal(base_size = 17)

# Bar chart living situation (should be possible to reorder based on group size but that does not work)
Bar_livsit <- ggplot(POM_wide, aes(x=LivingSituat)) +
                     geom_bar(color = "black", fill="#008B8B") +
                     geom_text(aes(label = scales::percent(..prop..)), stat= "count", vjust = -0.5, size=5) +
                     #scale_x_discrete(limits=c("Alone", "With partner", "With family", "Other")) +
                     labs(x="Living Situation", y="Frequency") + 
                     theme_minimal(base_size = 17)

# Bar chart daily activity (recode first because no participants have score 4 (student) so remove this option)
library(dplyr)
Data_wide$Act <- as.numeric(dplyr::recode(as.character(Data_wide$Daily_act), "1"="1", "2"="2", "3"="3", "5"="4", "6"="5", "7"="6", "8"="7"))
Bar_act <- ggplot(POM_wide, aes(x=DailyActivity)) +
                  geom_bar(color = "black", fill="#5F9EA0") +
                  geom_text(aes(label = scales::percent(..prop..)), stat= "count", vjust = -0.4, size=5) +
                  #scale_x_discrete(limits=c("Paid job", "Household \n duties", "Retired", "No paid job \n due to PD", "Unvoluntarily \n no job", "Volunteer \n work", "Other")) +
                  labs(x="", y="Frequency") + 
                  theme_minimal(base_size = 17)


#-----------------------------------------------------------------------------------------#
#                                    Resilience                                           #
#-----------------------------------------------------------------------------------------#

# Calculate equation of PSS-SL linear relation in y = ax+b format
lm_eqn <- function(Data_long){
    m <- lm(PSS ~ SL, Data_long);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                     list(a = format(unname(coef(m)[1]), digits = 2),
                          b = format(unname(coef(m)[2]), digits = 2),
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
  } 
  
# Histogram of the 2 classes from latent class analysis
library(plyr)
mu <- ddply(Covid_wide, "SR_class", summarise, grp.mean=mean(SR))
colors <- c("#CC3311", "#009988");
Res_groups <- ggplot(Covid_wide, aes(x=SR, fill=SR_class, color=SR_class)) +
                    geom_histogram(position="identity",  alpha=0.4) + 
                    labs(x="Stressor-reactivity", y="Frequency") + 
                    geom_vline(data=mu, aes(xintercept=grp.mean, color=SR_class), linetype="dashed") +
                    theme_minimal(base_size = 17) +
                    scale_x_continuous(breaks = seq(-15,15,5)) +
                    scale_colour_manual(name = "",
                                        values = colors, aesthetics = c("colour", "fill"),
                                        breaks = c("1","2"),
                                        labels = c("Low stressor-reactivity","High stressor-reactivity"),
                                        guide = guide_legend(reverse=T))

# Plot timeline separate for all individuals (spaghetti plots)
SR_all <-    ggplot(data = Data_long, aes(x = Weeknr, y = SR, group = ID, color = factor(ID))) +
                    geom_point() + 
                    geom_line(data = Data_long[which(Data_long$SR !="NA"),]) +
                    xlab("Week number after baseline") + 
                    ylab("Individual resilience score") +
                    theme_minimal(base_size = 16) +
                    scale_x_continuous(breaks = seq(0,27,5)) +
                    scale_y_continuous(breaks = seq(-25,20,5))

# Plot only subject 189 (low), 30 (high), 325 (around 0) as example figure
Example <- Data_long[which(Data_long$Subj_ID == "sub-POMU12FA62414399DD8F" | Data_long$Subj_ID == 	
                             "sub-POMUEEB7307F823DB346" | Data_long$Subj_ID == "sub-POMU889979EA743C6466"),]
Example$ID <- ifelse(Example$Subj_ID=="sub-POMU12FA62414399DD8F", "A", ifelse(Example$Subj_ID=="sub-POMUEEB7307F823DB346", "B", "C"))
SR_example <- ggplot(data = Example, aes(x = Weeknr, y = SR, group = ID, color = ID)) +
                       geom_point(size=3, aes(shape=ID)) + 
                       scale_shape_manual(name = "Example participants",
                                          values=c(15, 16, 17),
                                          labels = c("Low SR","Average SR","High SR"))+                        
                       geom_line(linewidth = 0.9, aes(linetype=ID, color = ID)) +
                       scale_linetype_manual(name = "Example participants",
                                             values=c("solid", "dotted", "dashed"),
                                             labels = c("Low SR","Average SR","High SR"))+                     
                       xlab("Week number after baseline") + 
                       ylab("Individual stressor reactivity score") +
                       theme_minimal(base_size = 16) +
                       scale_x_continuous(breaks = seq(0,27,5)) +
                       scale_y_continuous(breaks = seq(-25,20,5)) +
                       scale_color_manual(name = "Example participants",
                                          values = c("#009988", "#EE7733", "#CC3311"),
                                          labels = c("Low SR","Average SR","High SR"))
SR_example <- Res_example + theme(legend.key.size = unit(0.7, "cm"), legend.key.width = unit(1,"cm"),
                    legend.text = element_text(size = 15), legend.title = element_text(size = 15, face = "bold"),
                    legend.margin=margin(0,0,0,0),
                    legend.box.margin=margin(-10,-5,-10,-10))

# Density distribution by ID
Res_density <- ggplot(data=Data_long, aes(x=Res)) + 
                      geom_density(aes(group=factor(Subj_ID), colour=factor(Subj_ID))) +
                      xlab("Resilience") + theme(legend.position = "none")

# Figure displaying way of SR measure calculation
SR_scatter <- ggplot(Data_long, aes(x = SL, y = PSS)) + 
               geom_point(size = 1.5, color = "#595959") +
               geom_smooth(method = lm, size = 1, se = F, color = "#EE7733") + 
               labs(x="Stressor load (SL)", y="Perceived stress scale (PSS)") + 
               scale_x_continuous(breaks = seq(0,75,10)) +
               #geom_text(x = 60, y = 31, size = 5, label = lm_eqn(Data_long), parse = TRUE) +
               #geom_text(x = 58, y = 34, size = 5.5, label = "PSS = 0.18*SL + 5.90") +
               theme_minimal(base_size = 16)
  
ggsave("SR_examples.svg", plot = SR_example, width=6.19, height=4.45, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )
ggsave("Resilience_calculation_scatter.svg", plot = SR_scatter, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )
ggsave("Resilience_time.jpg", plot = Res_time, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )
ggsave("Resilience_spaghettiplot.jpg", plot = Res_all, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )
ggsave("Resilience_boxplot_time_dot.tiff", plot = Res_boxtime, device='tiff', dpi=300, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )
ggsave("Resilience_hist_classes.tiff", plot = Res_groups, device='tiff', dpi=300, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )
