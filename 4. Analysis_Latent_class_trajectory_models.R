# We conducted a latent class trajectory analysis using the lcmm-package, to explore whether the sample could be divided into meaningful 
# heterogeneous subgroups (latent classes) with different longitudinal trajectories of stressor-reactivity. The lcmm function estimates models 
# by assuming a parameterized link function for linking the outcome Y(t) with the underlying latent process L(t) it measures.

# lcmm(fixed = y ~ x + cov, # 2-sided linear formula for specifying fixed-effects at latent process level
      #mixture = ~ x,       # 1-sided formula for class-specific fixed effects 
      #random = ~ x,        # formula for random-effects (time: random slope)
      #subject = "ID",      # name of covariate representing grouping structure (ID)
      #classmb,             # formula describing covariates in class-membership model
      #ng = 3,              # number of latent classes
      #link = "linear",     # link function to estimate (linear, beta, default spline, quantile spline)
      #data = Covid_long,   # optional data frame with variables named in fixed, mixture, random, classmb and subject
      #nsim = 100,          # number of points used to plot the estimated link function (default = 100)
      #na.action = 1)}      # how NA's are managed: 1 for na.omit, 2 for na.fail

rm(list = ls())

library(readr)
library(dplyr)
Data_long <-  read_csv("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/Covid_POM_long.csv")
Covid_long <- Data_long[c("Subj_ID", "Age_cov", "Sex", "LivingSituat", "Weeknr", "timepoint", "SR", "PSS", "SL", "UPDRS", "PAS")]
Covid_long <- Covid_long %>% mutate(ID = as.factor(Subj_ID))
Covid_long$ID <- as.numeric(Covid_long$ID)
rm(Data_long)

library(lcmm)
library(splines)

### Step 1: check which residual structure fits best ###
SR_check1 <- lcmm(SR~1+Weeknr, random =~-1, subject = "ID", ng = 1, data = Covid_long)
SR_check2 <- lcmm(SR~1+Weeknr, mixture=~1+Weeknr, random =~-1, subject = "ID", ng = 2, B = SR_check1, data = Covid_long)
#SR_check3 <- lcmm(SR~1+Weeknr+I(Weeknr^2), mixture=~1+Weeknr+I(Weeknr^2), random =~-1, subject = "ID", ng = 2, data = Covid_long)

plot(SR_check1) # Plot residuals
plot(SR_check2) # Plot residuals
plot(SR_check2,which="postprob") # postprob plot

################ Test models to define classes with different trajectories in stress-reactivity ###################  

SR_lin1 <- lcmm(SR~Weeknr, random =~ Weeknr, subject = "ID", ng = 1, link = "linear", data = Covid_long)
SR_lin2 <- lcmm(SR~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 2, link = "linear", B = SR_lin1, data = Covid_long)      # model with 2 classes
SR_lin3 <- lcmm(SR~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "linear", B = SR_lin1, data = Covid_long)      # model with 3 classes
SR_lin4 <- lcmm(SR~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "linear", B = SR_lin1, data = Covid_long)      # model with 4 classes

# models with spline link function
SR_spl1   <- lcmm(SR~Weeknr, random =~ Weeknr, subject = "ID", ng = 1, link = "splines", data = Covid_long)                                    # model with 1 class
SR_spl2   <- lcmm(SR~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 2, link = "splines", B = SR_spl1, data = Covid_long)    # model with 2 classes
SR_spl3   <- lcmm(SR~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "splines", B = SR_spl1, data = Covid_long)    # model with 3 classes
SR_spl4   <- lcmm(SR~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "splines", B = SR_spl1, data = Covid_long)    # model with 4 classes
SR_spl1ns <- lcmm(SR~ ns(Weeknr, knots=c(3,6,10,25)), random = ~ Weeknr, subject = "ID", ng = 1, link = "splines", data = Covid_long)
SR_spl2ns <- lcmm(SR~ ns(Weeknr, knots=c(3,6,10,25)), mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                   random = ~ Weeknr, subject = "ID", ng = 2, link = "splines", B = SR_spl1ns, data = Covid_long)                   # model with 2 classes and manual splines
SR_spl3ns <- lcmm(SR~ ns(Weeknr, knots=c(3,6,10,25)), mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                   random = ~ Weeknr, subject = "ID", ng = 3, link = "splines", data = Covid_long)                                  # model with 3 classes and manual splines
SR_spl4ns <- lcmm(SR~ ns(Weeknr, knots=c(3,6,10,25)), mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                   random = ~ Weeknr, subject = "ID", ng = 4, link = "splines", data = Covid_long)                                  # model with 4 classes and manual splines

#  models with beta link function
SR_bet1   <- lcmm(SR~Weeknr, random =~ Weeknr, subject = "ID", ng = 1, link = "beta", data = Covid_long)                          # model with 1 class
SR_bet2   <- lcmm(SR~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 2, link = "beta", data = Covid_long)       # model with 2 classes
SR_bet3   <- lcmm(SR~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "beta", data = Covid_long)       # model with 3 classes
SR_bet4   <- lcmm(SR~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "beta", data = Covid_long)       # model with 4 classes
SR_bet2ns <- lcmm(SR~ ns(Weeknr, knots=c(3,6,10,17,25)), mixture = ~ ns(Weeknr, knots=c(3,6,10,17,25)),                          
                   random = ~ Weeknr, subject = "ID", ng = 2, link = "beta", data = Covid_long)                                     # model with 2 classes and manual splines
SR_bet3ns <- lcmm(SR~ ns(Weeknr, knots=c(3,6,10,25)), mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                   random = ~ Weeknr, subject = "ID", ng = 3, link = "beta", data = Covid_long)                                     # model with 3 classes and manual splines
SR_bet4ns <- lcmm(SR~ ns(Weeknr, knots=c(3,6,10,25)), mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                   random = ~ Weeknr, subject = "ID", ng = 4, link = "beta", data = Covid_long)                                     # model with 4 classes and manual splines

# models with quantile splines link function
SR_qua2   <- lcmm(SR~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 2, link = "2-quant-splines", data = Covid_long) # model with 2 classes
SR_qua3   <- lcmm(SR~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "2-quant-splines", data = Covid_long) # model with 3 classes
SR_qua4   <- lcmm(SR~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "2-quant-splines", data = Covid_long) # model with 4 classes
SR_qua2ns <- lcmm(SR~ ns(Weeknr, knots=c(3,6,10,25)), mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                   random = ~ Weeknr, subject = "ID", ng = 2, link = "2-quant-splines", data = Covid_long)                               # model with 2 classes and manual splines
SR_qua3ns <- lcmm(SR~ ns(Weeknr, knots=c(3,6,10,25)), mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                 random = ~ Weeknr, subject = "ID", ng = 3, link = "2-quant-splines", data = Covid_long)                                 # model with 3 classes and manual splines
SR_qua4ns <- lcmm(SR~ ns(Weeknr, knots=c(3,6,10,25)), mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                 random = ~ Weeknr, subject = "ID", ng = 4, link = "2-quant-splines", data = Covid_long)                                 # model with 4 classes and manual splines

# to compare model fit results, use summarytable. For class division, use postprob
summarytable(SR_lin1, SR_lin2, SR_lin3, SR_spl1, SR_spl2, which = c("loglik", "AIC", "BIC"))
summarytable(SR_spl1, SR_spl2, which = c("loglik", "AIC", "BIC"))
postprob(SR_spl2)

# Best models based on BIC: SR_spl2ns, SR_bet2ns; but based on class division: SR_qua2, SR_spl2 and SR_bet2
# pull out who is in which class
classes <- as.data.frame(c(SR_spl2$pprob[,1:2],SR_spl2ns$pprob[,1:2], SR_bet2$pprob[,1:2], SR_bet2ns$pprob[,1:2], SR_qua2$pprob[,1:2]))
names(classes)[names(classes) == "class"] <- "SR_spl2"
names(classes)[names(classes) == "class.1"] <- "SR_spl2ns"
names(classes)[names(classes) == "class.2"] <- "SR_bet2"
names(classes)[names(classes) == "class.3"] <- "SR_bet2ns"
names(classes)[names(classes) == "class.4"] <- "SR_qua2"
classes <- classes[c("ID", "SR_spl2", "SR_spl2ns", "SR_bet2", "SR_bet2ns", "SR_qua2")]
classes_long = merge(Covid_long, classes, by = "ID")
classes_long <- subset(classes_long, select=c("ID", "Weeknr", "SR", "SR_spl2", "SR_spl2ns", "SR_bet2", "SR_bet2ns", "SR_qua2"))
classes_long$SR_spl2 <- as.factor(classes_long$SR_spl2)
classes_long$SR_spl2ns <- as.factor(classes_long$SR_spl2ns)
classes_long$SR_bet2 <- as.factor(classes_long$SR_bet2)
classes_long$SR_bet2ns <- as.factor(classes_long$SR_bet2ns)
classes_long$SR_qua2 <- as.factor(classes_long$SR_qua2)

# plot SR_spl2 model
library(ggplot2)
p_mean <- ggplot(classes_long, aes(Weeknr, SR, group=ID, colour=SR_spl2)) + 
          geom_smooth(aes(group=SR_spl2), size=2, se=T)  + 
          labs(x="Week number", y="Stressor reactivity (SR)", colour="Latent Class") + 
          scale_x_continuous(breaks=seq(0,26,5)) +
          scale_y_continuous(breaks=seq(-8,5)) +
          theme_minimal(base_size = 17)
p_all <-  ggplot(classes_long, aes(Weeknr, SR, group=ID, colour=SR_spl2)) + 
          geom_line() + geom_smooth(aes(group=SR_spl2), size=2, se=F)  + 
          labs(x="Week number",y="Stressor reactivity (SR)",colour="Latent Class") +
          scale_x_continuous(breaks=seq(0,26,5)) +
          scale_y_continuous(breaks=seq(-25,15,5)) +
          theme_minimal(base_size = 17)

write.csv(classes_long,'M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models/SR_classes_long.csv')
write.csv(classes,'M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models/SR_classes_mean.csv')
ggsave("SR_2class_all.jpg", plot = p_all, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )
ggsave("SR_2class_mean.jpg", plot = p_mean, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )


#---------------------------------------------------------------------------------------------------------------------#
#                        EXTRA: latent class models for other longitudinal variables                                  #
#---------------------------------------------------------------------------------------------------------------------#


############# Test models to define classes with different trajectories in stressor exposure (SL) ################

SL_bet3_new <- lcmm(SL~ 1+ Weeknr + I(Weeknr^2) + Sex, mixture = ~1+ Weeknr + I(Weeknr^2), random =~ Weeknr, subject = "ID", ng = 3, link = "beta", data = Covid_long)

SL_lin1 <- lcmm(SL~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "linear", data = Covid_long)
SL_lin2 <- lcmm(SL~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "linear", data = Covid_long)
SL_lin3 <- lcmm(SL~Weeknr + Sex, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "linear", data = Covid_long) 
SL_lin4 <- lcmm(SL~Weeknr + Sex, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "linear", data = Covid_long) 
SL_spl1 <- lcmm(SL~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "splines", data = Covid_long) 
SL_spl2 <- lcmm(SL~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "splines", data = Covid_long) 
SL_spl3 <- lcmm(SL~Weeknr + Sex, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "splines", data = Covid_long) 
SL_spl4 <- lcmm(SL~Weeknr + Sex, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "splines", data = Covid_long)
SL_bet1 <- lcmm(SL~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "beta", data = Covid_long)
SL_bet2 <- lcmm(SL~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "beta", data = Covid_long)
SL_bet3 <- lcmm(SL~Weeknr + Sex, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "beta", data = Covid_long)
SL_bet4 <- lcmm(SL~Weeknr + Sex, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "beta", data = Covid_long)
SL_bet5 <- lcmm(SL~ ns(Weeknr, knots=c(3,6,10,25)) + Sex, mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                random = ~ Weeknr, subject = "ID", ng = 3, link = "beta", data = Covid_long)
SL_bet6 <- lcmm(SL~ ns(Weeknr, knots=c(3,6,10,25)) + Sex, mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                random = ~ Weeknr, subject = "ID", ng = 4, link = "beta", data = Covid_long)
SL_bet7 <- lcmm(SL~ ns(Weeknr, knots=c(3,6,10,25)), mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                random = ~ Weeknr, subject = "ID", ng = 3, link = "beta", data = Covid_long)
SL_bet8 <- lcmm(SL~ ns(Weeknr, knots=c(3,6,10,25)), mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                random = ~ Weeknr, subject = "ID", ng = 4, link = "beta", data = Covid_long)
SL_qua1 <- lcmm(SL~Weeknr + Sex, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "2-quant-splines", data = Covid_long)
SL_qua2 <- lcmm(SL~Weeknr + Sex, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "2-quant-splines", data = Covid_long)

summarytable(SL_bet7, SL_bet8, which = c("loglik", "AIC", "BIC"))
postprob(SL_bet5)

################ Test models to define classes with different trajectories in PSS ###################

PSS_lin1 <- lcmm(PSS~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "linear", data = Covid_long)
PSS_lin2 <- lcmm(PSS~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "linear", data = Covid_long)
PSS_lin3 <- lcmm(PSS~Weeknr + Sex, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "linear", data = Covid_long) 
PSS_lin4 <- lcmm(PSS~Weeknr + Sex, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "linear", data = Covid_long) 
PSS_spl1 <- lcmm(PSS~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "splines", data = Covid_long) 
PSS_spl2 <- lcmm(PSS~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "splines", data = Covid_long) 
PSS_spl3 <- lcmm(PSS~Weeknr + Sex, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "splines", data = Covid_long) 
PSS_spl4 <- lcmm(PSS~Weeknr + Sex, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "splines", data = Covid_long)
PSS_spl5 <- lcmm(PSS~ ns(Weeknr, knots=c(3,6,10,25)) + Sex, mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                 random = ~ Weeknr, subject = "ID", ng = 3, link = "splines", data = Covid_long)
PSS_spl6 <- lcmm(PSS~ ns(Weeknr, knots=c(3,6,10,25)) + Sex, mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                 random = ~ Weeknr, subject = "ID", ng = 4, link = "splines", data = Covid_long)
PSS_bet1 <- lcmm(PSS~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "beta", data = Covid_long)
PSS_bet2 <- lcmm(PSS~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "beta", data = Covid_long)
PSS_bet3 <- lcmm(PSS~Weeknr + Sex, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "beta", data = Covid_long)
PSS_bet4 <- lcmm(PSS~Weeknr + Sex, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "beta", data = Covid_long)
PSS_bet5 <- lcmm(PSS~ ns(Weeknr, knots=c(3,6,10,25)) + Sex, mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                 random = ~ Weeknr, subject = "ID", ng = 3, link = "beta", data = Covid_long)
PSS_bet6 <- lcmm(PSS~ ns(Weeknr, knots=c(3,6,10,25)) + Sex, mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                 random = ~ Weeknr, subject = "ID", ng = 4, link = "beta", data = Covid_long)
PSS_bet7 <- lcmm(PSS~ ns(Weeknr, knots=c(3,6,10,25)), mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                 random = ~ Weeknr, subject = "ID", ng = 3, link = "beta", data = Covid_long)
PSS_bet8 <- lcmm(PSS~ ns(Weeknr, knots=c(3,6,10,25)), mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                 random = ~ Weeknr, subject = "ID", ng = 4, link = "beta", data = Covid_long)
PSS_qua1 <- lcmm(PSS~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "2-quant-splines", data = Covid_long)
PSS_qua2 <- lcmm(PSS~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "2-quant-splines", data = Covid_long)

################ Test models to define classes with different trajectories in UPDRS2 ###################

UP_lin1 <- lcmm(UPDRS2~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "linear", data = Covid_long)
UP_lin2 <- lcmm(UPDRS2~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "linear", data = Covid_long)
UP_spl1 <- lcmm(UPDRS2~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "splines", data = Covid_long) 
UP_spl2 <- lcmm(UPDRS2~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "splines", data = Covid_long) 
UP_spl3 <- lcmm(UPDRS2~ ns(Weeknr, knots=c(3,6,10,25)), mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                random = ~ Weeknr, subject = "ID", ng = 3, link = "splines", data = Covid_long)
UP_spl4 <- lcmm(UPDRS2~ ns(Weeknr, knots=c(3,6,10,25)), mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                random = ~ Weeknr, subject = "ID", ng = 4, link = "splines", data = Covid_long)
UP_bet1 <- lcmm(UPDRS2~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "beta", data = Covid_long)
UP_bet2 <- lcmm(UPDRS2~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "beta", data = Covid_long)
UP_bet3 <- lcmm(UPDRS2~ ns(Weeknr, knots=c(3,6,10,25)), mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                random = ~ Weeknr, subject = "ID", ng = 3, link = "beta", data = Covid_long)
UP_bet4 <- lcmm(UPDRS2~ ns(Weeknr, knots=c(3,6,10,25)), mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                random = ~ Weeknr, subject = "ID", ng = 4, link = "beta", data = Covid_long)

################ Test models to define classes with different trajectories in PAS ###################

PAS_lin1 <- lcmm(PAS~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "linear", data = Covid_long)
PAS_lin2 <- lcmm(PAS~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "linear", data = Covid_long)
PAS_lin3 <- lcmm(PAS~Weeknr + Sex, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "linear", data = Covid_long) 
PAS_lin4 <- lcmm(PAS~Weeknr + Sex, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "linear", data = Covid_long) 
PAS_spl1 <- lcmm(PAS~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "splines", data = Covid_long) 
PAS_spl2 <- lcmm(PAS~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "splines", data = Covid_long) 
PAS_spl3 <- lcmm(PAS~Weeknr + Sex, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "splines", data = Covid_long) 
PAS_spl4 <- lcmm(PAS~Weeknr + Sex, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "splines", data = Covid_long)
PAS_spl5 <- lcmm(PAS~ ns(Weeknr, knots=c(3,6,10,25)) + Sex, mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                 random = ~ Weeknr, subject = "ID", ng = 3, link = "splines", data = Covid_long)
PAS_spl6 <- lcmm(PAS~ ns(Weeknr, knots=c(3,6,10,25)) + Sex, mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                 random = ~ Weeknr, subject = "ID", ng = 4, link = "splines", data = Covid_long)
PAS_bet1 <- lcmm(PAS~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "beta", data = Covid_long)
PAS_bet2 <- lcmm(PAS~Weeknr, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "beta", data = Covid_long)
PAS_bet3 <- lcmm(PAS~Weeknr + Sex, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 3, link = "beta", data = Covid_long)
PAS_bet4 <- lcmm(PAS~Weeknr + Sex, mixture =~ Weeknr, random =~ Weeknr, subject = "ID", ng = 4, link = "beta", data = Covid_long)
PAS_bet5 <- lcmm(PAS~ ns(Weeknr, knots=c(3,6,10,25)) + Sex, mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                 random = ~ Weeknr, subject = "ID", ng = 3, link = "beta", data = Covid_long)
PAS_bet6 <- lcmm(PAS~ ns(Weeknr, knots=c(3,6,10,25)) + Sex, mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                 random = ~ Weeknr, subject = "ID", ng = 4, link = "beta", data = Covid_long)
PAS_bet7 <- lcmm(PAS~ ns(Weeknr, knots=c(3,6,10,25)), mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                 random = ~ Weeknr, subject = "ID", ng = 3, link = "beta", data = Covid_long)
PAS_bet8 <- lcmm(PAS~ ns(Weeknr, knots=c(3,6,10,25)), mixture = ~ ns(Weeknr, knots=c(3,6,10,25)), 
                 random = ~ Weeknr, subject = "ID", ng = 4, link = "beta", data = Covid_long)

############### Plot winning models #################

### PAS winning model ###

people1 <- as.data.frame(PAS_bet8$pprob[,1:2]) # pull out who is in which class
d       <- merge(Covid_long, people, by = "ID")
groups  <- subset(d, select=c("ID", "Weeknr", "PAS", "class"))
groups$class <- as.factor(groups$class)

p1_all <-  ggplot(groups, aes(Weeknr, PAS, group=ID, colour=class)) + # plot model for PAS with beta link and manual splines (PAS_bet5 and PAS_bet6)
           geom_line() + geom_smooth(aes(group=class), size=2, se=F)  + 
           labs(x="Time",y="PAS",colour="Latent Class")
p1_mean <- ggplot(groups, aes(Weeknr, PAS, group=ID, colour=class)) + 
           geom_smooth(aes(group=class), size=2, se=T)  + 
           labs(title="PAS beta link incl. government measures splines", x="Week number",y="Episodic anxiety (PAS)",colour="Latent Class") + 
           theme_minimal(base_size = 17)

### SL winning model ###

people2 <- as.data.frame(SL_bet8$pprob[,1:2]) # pull out who is in which class
d       <- merge(Covid_long, people1, by = "ID")
groups  <- subset(d, select=c("ID", "Weeknr", "SL", "class"))
groups$class <- as.factor(groups$class)
 
p2_all <-  ggplot(groups, aes(Weeknr, SL, group=ID, colour=class)) + # plot model for SL with beta link and timepoint spline
           geom_line() + geom_smooth(aes(group=class), size=2, se=F)  + 
           labs(x="Week number",y="Stressor load",colour="Latent Class")
p2_mean <- ggplot(groups, aes(Weeknr, SL, group=ID, colour=class)) + 
           geom_smooth(aes(group=class), size=2, se=T)  + 
           labs(title="SL beta link incl. government measures splines", x="Week number",y="Stressor load",colour="Latent Class")+ 
           theme_minimal(base_size = 17)

### UPDRS2 winning model ###

people3 <- as.data.frame(UP_bet3$pprob[,1:2])
d       <- merge(Covid_long, people4, by = "ID")
groups  <- subset(d, select=c("ID", "Weeknr", "UPDRS2", "class"))
groups$class <- as.factor(groups$class)

p3_mean <- ggplot(groups, aes(Weeknr, UPDRS2, group=ID, colour=class)) + 
           geom_smooth(aes(group=class), size=2, se=T)  + 
           labs(title="UPDRS2 (beta link + government measures spline) ", x="Week number",y="UPDRS2",colour="Latent Class")+ 
           theme_minimal(base_size = 17)
