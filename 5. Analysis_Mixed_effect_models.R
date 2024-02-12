# This script models trajectories of longitudinal measures (11 timepoints of COVID-survey) and plots predicted trajectories. 
# The nlme-package is used for mixed-effect models, to explore trajectories of the longitudinal measures stressor reactivity (SR),
# COVID-related stressor exposure, perceived stress (PSS total score), episodic anxiety (PAS) and daily PD symptoms (UPDRSIb+II).
# Mixed-effects models are robust in case of missing data, so models were estimated using all available data, even if timepoints were missing. 
# Weeks since baseline survey (T1) were used as timescale, and natural cubic splines were included for time to allow non-linearity in the model 
# at timepoints of important changes in governmental COVID-19 regulations. Plots are created with predicted trajectories for all variables.

rm(list = ls())

###################     Load and restructure data files     ###################

library(readr)
library(dplyr)

# Long data set needed for mixed-effect models
Data_long <-  read_csv("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/Covid_POM_long.csv")
Covid_long <- Data_long[c("Subj_ID", "Age_cov", "Sex", "LivingSituat", "Weeknr", "timepoint", "Res", "PSS", "SL", "UPDRS", "PAS")] # keep only relevant variables
Classes <- read_csv("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models/Res_classes_long.csv") # load results from latent class trajectory models to have low and high SR classes
Covid_long <- cbind(Covid_long, Classes[5]) # add column with SR class (from latent class analysis)
rm(Data_long, Classes) # remove data file with all variables that are not needed for mixed-models
Covid_long$SR <- -Covid_long$Res 
Covid_long$SR_class <- as.factor(Covid_long$RES_spl2) # 1=low SR (high SR), 2=high SR (low SR)
names(Covid_long)[names(Covid_long) == "Age_cov"] <- "Age"
Covid_long$Age_bin <- as.factor(ifelse(Covid_long$Age<63, 1, 2))
Covid_long$IDnr <- match(Covid_long$Subj_ID, unique(Covid_long$Subj_ID))
Covid_long$Sex <- as.factor(Covid_long$Sex)
Covid_long$LivingSituat <- as.factor(Covid_long$LivingSituat)

###################     Model trajectory of stressor-reactivity (SR)    ###################

library(nlme)
library(lmerTest)
library(splines)

# Explore different random effect structures, using fixed covariates age, sex and living situation. Allow age and sex to interact with time.
SR_base <-  lme(SR ~ Weeknr * (Age + Sex) + LivingSituat,                                       random = ~ 1 | Subj_ID, # Model for linear SR trajectory
                data = Covid_long, na.action = na.omit, keep.data = FALSE) 
SR_spl  <-  lme(SR ~ ns(Weeknr,knots=c(3,6,10,17,25)) * (Age + Sex) + LivingSituat,             random = ~ 1 | Subj_ID, # Model for SR trajectory with splines
                data = Covid_long, na.action = na.omit, keep.data = FALSE) 
SR_fac  <-  lme(SR ~  ns(Weeknr,knots=c(3,6,10,17,25)) * (SR_class + Age + Sex) + LivingSituat, random = ~ 1 | Subj_ID, # Model with SR class as factor added
                data = Covid_long, na.action = na.omit, keep.data = FALSE) 

anova(update(SR_base, . ~ ., method = "ML"), # model with splines on dates of changes in governmental COVID-19 regulations is better (p<.0001)
      update(SR_spl, . ~ ., method = "ML")) 

################### Model trajectories of all longitudinal variables stratified for SR class ###################

# divide whole sample in two groups based on latent SR classes, to explore influence of covariates
Data_SR1 <- Covid_long %>% filter(SR_class ==1)
Data_SR2 <- Covid_long %>% filter(SR_class ==2)

# SR models
SR_high <- lme(SR ~ (Age + Sex) *ns(Weeknr,knots=c(3,6,10,17,24)) + LivingSituat, random = ~ 1 | Subj_ID,
                data = Data_SR1, na.action = na.omit, keep.data = FALSE)
anova(SR_high) # interaction age * time sign.

SR_low  <- lme(SR ~ (Age + Sex) *ns(Weeknr,knots=c(3,6,10,17,24)) + LivingSituat, random = ~ 1 | Subj_ID,
                data = Data_SR2, na.action = na.omit, keep.data = FALSE)
anova(SR_low) # interaction age * time not sign.

# SL models
SL_high <- lme(SL ~ (Age + Sex) *ns(Weeknr,knots=c(3,6,10,17,24)) + LivingSituat, random = ~ 1 | Subj_ID,
               data = Data_SR1, na.action = na.omit, keep.data = FALSE)
anova(SL_high) # interaction age * time sign

SL_low <-  lme(SL ~ (Age + Sex) *ns(Weeknr,knots=c(3,6,10,17,24)) + LivingSituat, random = ~ 1 | Subj_ID,
              data = Data_SR2, na.action = na.omit, keep.data = FALSE)
anova(SL_low) # interaction age*time and sex*time sign

# PSS models
PSS_high <- lme(PSS ~ (Age + Sex) *ns(Weeknr,knots=c(3,6,10,17,24)) + LivingSituat, random = ~ 1 | Subj_ID,
                data = Data_SR1, na.action = na.omit, keep.data = FALSE)
anova(PSS_high) # interaction age * time sign

PSS_low <- lme(PSS ~ (Age + Sex) *ns(Weeknr,knots=c(3,6,10,17,24)) + LivingSituat, random = ~ 1 | Subj_ID,
              data = Data_SR2, na.action = na.omit, keep.data = FALSE)
anova(PSS_low) # interactions not sign.

# UPDRS models
UPD_high <- lme(UPDRS ~ (Age + Sex) *ns(Weeknr,knots=c(3,6,10,17,24)) + LivingSituat, random = ~ 1 | Subj_ID,
                data = Data_SR1, na.action = na.omit, keep.data = FALSE)
anova(UPD_high) # interactions not sign

UPD_low <- lme(UPDRS ~ (Age + Sex) *ns(Weeknr,knots=c(3,6,10,17,24)) + LivingSituat, random = ~ 1 | Subj_ID,
                data = Data_SR2, na.action = na.omit, keep.data = FALSE)
anova(UPD_high) # interactionS not sign

# PAS models
PAS_high <- lme(PAS ~ (Age + Sex) *ns(Weeknr,knots=c(3,6,10,17,24)) + LivingSituat, random = ~ 1 | Subj_ID,
                data = Data_SR1, na.action = na.omit, keep.data = FALSE)
anova(PAS_high) # interactionS not sign

PAS_low <- lme(PAS ~ (Age + Sex) *ns(Weeknr,knots=c(3,6,10,17,24)) + LivingSituat, random = ~ 1 | Subj_ID,
               data = Data_SR2, na.action = na.omit, keep.data = FALSE)
anova(PAS_low) # interaction sex * time sign

################### Model trajectories of with SR class as factor  ###################

# use lm() for this, because ggpredict function does not work for lme or lmer
SR_all <- lm(SR ~ SR_class *ns(Weeknr,knots=c(3,6,10,17,24)) + Age + Sex + LivingSituat + (1|IDnr), data = Covid_long, na.action = na.omit)
SL_all  <- lm(SL ~ (Age + Sex + SR_class) *ns(Weeknr,knots=c(3,6,10,17,24)) + LivingSituat + (1|IDnr), data = Covid_long, na.action = na.omit)
PSS_all <- lm(PSS ~ (Age + Sex + SR_class) *ns(Weeknr,knots=c(3,6,10,17,24)) + LivingSituat + (1|IDnr), data = Covid_long, na.action = na.omit)
PAS_all <- lm(PAS ~ (Age + Sex + SR_class) *ns(Weeknr,knots=c(3,6,10,17,24)) + LivingSituat + (1|IDnr), data = Covid_long, na.action = na.omit)
UPD_all <- lm(UPDRS ~ (Age + Sex + SR_class) *ns(Weeknr,knots=c(3,6,10,17,24)) + LivingSituat + (1|IDnr), data = Covid_long, na.action = na.omit)

anova(SR_all)
anova(SL_all)
anova(PSS_all)
anova(PAS_all)
anova(UPD_all)

################### Plot predicted trajectories of longitudinal variables ###################

library(ggeffects)
library(ggplot2)

# Plot predicted values per SR class
# ggpredict holds non-focal terms constant at mean value or at reference level (for factors). With condition= you can specify other levels to use. 
# Here, we used the most common levels for sex (male) and living situation (living with partner), mean age is used by default.
SR_pred <- ggpredict(SR_all, terms = c("Weeknr", "SR_class"), condition = c(Sex = "Male", LivingSituat = 2))
SL_pred <-  ggpredict(SL_all, terms = c("Weeknr", "SR_class"), condition = c(Sex = "Male", LivingSituat = 2))
PSS_pred <- ggpredict(PSS_all, terms = c("Weeknr", "SR_class"), condition = c(Sex = "Male", LivingSituat = 2))
PAS_pred <- ggpredict(PAS_all, terms = c("Weeknr", "SR_class"), condition = c(Sex = "Male", LivingSituat = 2))
UPD_pred <- ggpredict(UPD_all, terms = c("Weeknr", "SR_class"), condition = c(Sex = "Male", LivingSituat = 2))

# Make plot of change in SR per SR class
SR <- ggplot(SR_pred, aes(x = x, y = predicted, color=group, group = group)) +
             geom_line(linewidth = 1, show.legend = F) +
             geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .15, colour = NA) +
             labs(x = "Weeks since baseline", y = "Stressor reactivity", color = "SR class\n") +
             scale_x_continuous(breaks = seq(0,25,5)) + 
             scale_y_continuous(breaks = seq(-6,6,2)) + 
             scale_color_manual(labels = c("Low SR", "High SR"), values = c("#009988", "#CC3311")) +
             theme_minimal() + 
             theme(legend.key.size = unit(1, "cm"), legend.key.width = unit(1,"cm"),
                  legend.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold"),
                  axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 5)),
                  axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 5)),
                  axis.text = element_text(size = 13))

# Make plot of change in stressor exposure (SL) per SR class
sl <- ggplot(SL_pred, aes(x = x, y = predicted, color=group, group = group)) +
            geom_line(linewidth = 1, show.legend = F) +
            geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .15, colour = NA) +
            labs(x = "Weeks since baseline", y = "Stressor-exposure", color = "SR class\n") +
            scale_x_continuous(breaks = seq(0,25,5)) + 
            scale_y_continuous(breaks = seq(10,28,2)) + 
            scale_color_manual(labels = c("Low SR", "High SR"), values = c("#009988", "#CC3311")) +
            theme_minimal() + 
            theme(legend.key.size = unit(1, "cm"), legend.key.width = unit(1,"cm"),
                  legend.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold"),
                  axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 5)),
                  axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 5)),
                  axis.text = element_text(size = 13))

# Make plot of change in perceived stress (PSS) per SR class
pss <- ggplot(PSS_pred, aes(x = x, y = predicted, color=group, group = group)) +
            geom_line(linewidth = 1, show.legend = F) +
            geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .15, colour = NA) +
            labs(x = "Weeks since baseline", y = "Perceived stress scale", color = "SR class\n") +
            scale_x_continuous(breaks = seq(0,25,5)) + 
            scale_y_continuous(breaks = seq(0,20,2)) + 
            scale_color_manual(labels = c("Low SR", "High SR"), values = c("#009988", "#CC3311")) +
            theme_minimal() + 
            theme(legend.key.size = unit(1, "cm"), legend.key.width = unit(1,"cm"),
                  legend.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold"),
                  axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 5)),
                  axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 5)),
                  axis.text = element_text(size = 13))

# Make plots of change in PAS over survey follow up time per SR class
pas <- ggplot(PAS_pred, aes(x = x, y = predicted, color=group, group = group)) +
             geom_line(linewidth = 1, show.legend = F) +
             geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .15, colour = NA) +
             labs(x = "Weeks since baseline", y = "Parkinson Anxiety Scale", color = "SR class\n") +
             scale_x_continuous(breaks = seq(0,25,5)) + 
             #scale_y_continuous(breaks = seq(0,4,0.5)) + 
             scale_color_manual(labels = c("High SR", "Low SR"), values = c("#009988", "#CC3311")) +
             theme_minimal() + 
             theme(legend.key.size = unit(1, "cm"), legend.key.width = unit(1,"cm"),
                  legend.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold"),
                  axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 5)),
                  axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 5)),
                  axis.text = element_text(size = 13))

# Make plot of change in daily PD symptoms (UPDRSIb+II) over survey follow up time per SR class
upd <- ggplot(UPD_pred, aes(x = x, y = predicted, color=group, group = group)) +
            geom_line(linewidth = 1, show.legend = F) +
            geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .15, colour = NA) +
            labs(x = "Weeks since baseline", y = "Daily PD symptoms", color = "SR class\n") +
            scale_x_continuous(breaks = seq(0,25,5)) + 
            scale_y_continuous(breaks = seq(10,24,2)) + 
            scale_color_manual(labels = c("Low SR", "High SR"), values = c("#009988", "#CC3311")) +
            theme_minimal() + 
            theme(legend.key.size = unit(1, "cm"), legend.key.width = unit(1,"cm"),
                  legend.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold"),
                  axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 5)),
                  axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 5)),
                  axis.text = element_text(size = 13))

# Save trajectory plots (figures 3b-e in manuscript)
ggsave(file = "UPDRS_SRclass_pred.svg", plot = upd, width=3.41, height=3.59, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )
ggsave(file = "PSS_SRclass_pred.svg", plot = pss, width=3.41, height=3.59, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )
ggsave(file = "PAS_SRclass_pred.svg", plot = pas, width=3.41, height=3.59, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )
ggsave(file = "SL_SRclass_pred.svg", plot = sl, width=3.41, height=3.59, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )
ggsave(file = "SR_SRclass_pred.svg", plot = SR, width=3.41, height=3.59, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )

# Create density plot for the two SR-classes (with count instead of density)
SRclass_density <- ggplot(Covid_long, aes(y = SR, after_stat(count), fill = SR_class)) +
                          geom_density(alpha = 0.6, adjust = 0.7, show.legend = F, size = 0.4) + 
                          geom_boxplot(aes(x = -20, y = SR, fill = SR_class), alpha = 1, width = 30, colour = "black", size = 0.4, show.legend = F) + 
                          scale_y_continuous(breaks=c(-25, -20, -15, -10, -5, 0, 5, 10, 15)) + 
                          labs(x = "Count per SR class", y = "Stressor reactivity") +  
                          scale_fill_manual(name = "SR class", 
                                            values = c("#009988", "#CC3311"), 
                                            labels = c("Low SR", "High SR")) + 
                          theme_minimal() + 
                          theme(legend.key.size = unit(1.0, "cm"), legend.key.width = unit(1.0,"cm"),
                                legend.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold"), legend.box.margin=margin(-10,-5,-10,-10),
                                axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 5)), 
                                axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 5)),
                                axis.text = element_text(size = 13))

# Extra: Option with histogram, shows division SR classes (but less smooth); not used in manuscript
ggplot(Covid_long, aes(y = SR, fill = SR_class)) +
        geom_histogram(binwidth = 0.5, position = "identity", alpha = 0.5, show.legend = F) + 
        geom_boxplot(aes(x = -15, y = SR, fill = SR_class), alpha = 0.6, width = 20, colour = "black", show.legend = F) + 
        scale_y_continuous(breaks=c(-25, -20, -15, -10, -5, 0, 5, 10, 15)) + 
        labs(x = "Count per SR class", y = "SR", color = "SR class\n") +  
        scale_fill_manual(values=c("#009988", "#CC3311"), labels = c("High SR", "Low SR")) + 
        theme_minimal() + 
        theme(legend.key.size = unit(1.0, "cm"), legend.key.width = unit(1.0,"cm"),
              legend.text = element_text(size = 16), legend.title = element_text(size = 16, face = "bold"),
              axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),
              axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 10)),
              axis.text = element_text(size = 18))

# Save density plot (figure 3a in manuscript)
ggsave("SRclass_density.svg", plot = SRclass_density, width=3.41, height=3.59, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )


################### Extra: Make measure for univariate intra-individual variability (MSSD) ###################

# Not used in manuscript; used for exploration of variability in SR measure
library(data.table)
library(plyr)

# Function to calculate MSSD 
my.mssd <- function(data)
{  diffToNext <- data[2:length(data)]-data[1:(length(data)-1)] # Computes the difference between each value and the next
diffToNext2 <- diffToNext^2                  # Square the difference
SSdiff <- sum(diffToNext2,na.rm=TRUE)        # Take sum of the squared differences
denominator <- sum(!is.na(diffToNext))       # Compute the number of non-missing elements (denominator)
mssd <- SSdiff/denominator                   # Compute MSSD
return(mssd) }

# Calculate and plot MSSD for stressor reactivity (SR) and anxiety (PAS) measures 
mssd.stats <- ddply(Covid_long, "Subj_ID", summarize, MSSD_SR=my.mssd(SR), MSSD_PAS=my.mssd(PAS))

p <- ggplot(data=mssd.stats, aes(x=MSSD_SR)) +
            geom_histogram(fill="#b3dcdc", color="black") + 
            labs(x = "MSSD SR", y = "Frequency (N)") +
            theme_minimal() 
p + theme(text = element_text(size = 16), axis.text = element_text(size = 16))  
            ggplot(data=mssd.stats, aes(x=MSSD_PAS)) +
            geom_histogram(fill="white", color="black") + 
            labs(x = "MSSD episodic anxiety")

# Calculate correlations of MSSD with UPDRS, BDI and STAI (correlations are very low)
Data_wide <- merge(Data_wide, mssd.stats, by = "Subj_ID")
Data_wide = Data_wide[!Data_wide$Subj_ID =="sub-POMU900F78E54F00A78A",] # Remove patient with extreme high BDI (extreme outlier)
            
            
cor(Data_wide$UPDRSmean, Data_wide$MSSD_SR,  method = "pearson", use = "complete.obs")
cor(Data_wide2$BDI, Data_wide2$MSSD_SR,  method = "pearson", use = "complete.obs")    
cor(Data_wide$STAI_s, Data_wide$MSSD_SR,  method = "pearson", use = "complete.obs") 

p <- ggplot(Data_wide2, aes(x = BDI, y = MSSD_SR)) + 
            geom_point() + geom_smooth(method = lm, size = 1.25) + 
            labs(x="BDI", y="MSSD SR") +
            theme_minimal(base_size = 16)
