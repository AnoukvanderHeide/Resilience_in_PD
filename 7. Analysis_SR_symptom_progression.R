
# This script explores whether the COVID-19 pandemic changed the PD symptom progression rate for motor symptoms (UPDRS-III), 
# depressive symptoms (BDI) and anxiety symptoms (STAI), and whether this differs for low vs. high stress-reactive participants.
# We perform the analysis twice for all variables; once for the complete group comparing the last pre-COVID PPP visit with the first
# in-COVID PPP visit, once for the subgroup with 2 pre-COVID and 1 in-COVID PPP visit. Boxplots are created for all analyses.

rm(list = ls())

############## Load POM visit data (made with POM_visit_data.R) and resilience score data ###############

# Load libraries
library(readr)
library(tidyr)
library(dplyr)
library(rstatix)
library(zoo)
library(nlme)
library(ggplot2)
library(ggpubr)
library(emmeans)

POM_data_all <- read_csv("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/POM_allvisits.csv")
Covid <- read_csv("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/Covid_averages.csv")

############## UPDRS3 progression: rm-ANOVA comparing the 3 time points, no SR group comparison ###############

# Make subset of data specific for analysis of UPDRS3 progression
POM_Upd <- POM_data_all[c("Subj_ID", "Visit", "Age_cov", "Sex", "DisDur", "precov", "UPDRS3_off")] # Subset of variables needed for test progression motor sypmtoms

# add some UPDRS total scores manually that were not computed automatically
POM_Upd$UPDRS3_off[POM_Upd$Subj_ID =="sub-POMUDCE36875FFBC2199" & POM_Upd$Visit == 1] <- 33
POM_Upd$UPDRS3_off[POM_Upd$Subj_ID =="sub-POMU78FC9EC1D822238D" & POM_Upd$Visit == 2] <- 25

# Reformat variable names and types suitable for analysis
POM_Upd$Visit <- as.integer(POM_Upd$Visit)
Upd_wide <- pivot_wider(POM_Upd, names_from = Visit, values_from = UPDRS3_off)
names(Upd_wide)[names(Upd_wide) == "1"] <- "UPDRS3.1"
names(Upd_wide)[names(Upd_wide) == "2"] <- "UPDRS3.2"
names(Upd_wide)[names(Upd_wide) == "3"] <- "UPDRS3.3"

# Determine who misses UPDRS3 total scores for each visit
sum(is.na(Upd_wide$UPDRS3.1)) # 4 with missing UPDRS3 for visit 1
sum(is.na(Upd_wide$UPDRS3.2)) # 55 with missing UPDRS3 for visit 2
sum(is.na(Upd_wide$UPDRS3.3)) # 33 with missing UPDRS3 for visit 3

# Combine data and make separate groups for 1 or 2 pre-covid visits
Upd_wide <- cbind(Upd_wide, Covid$SR, Covid$SR_class)
pre1_Upd <- Upd_wide %>% filter(precov == 1) # group with 1 visit before covid
pre2_Upd <- Upd_wide %>% filter(precov == 2) # group with 2 visits before covid

# Put separate groups with 2 pre-covid visits back to long format to have one UPDRS3 variable
pre2_Upd_long <- pivot_longer(pre2_Upd, cols = starts_with("UPDRS3."), names_to = "Visit", values_to = "UPDRS3")
pre2_Upd_long$Visit <- as.factor(pre2_Upd_long$Visit)
levels(pre2_Upd_long$Visit) <- c("1", "2", "3")

# Show mean and SD for all 3 visits (whole group, no comparison resilience yet)
pre2_Upd_long %>% group_by(Visit) %>% get_summary_stats(UPDRS3, type = "mean_sd")

## Before rm-ANOVA, test assumptions: no outliers, normality per group, sphericity. Then: run model and look at results
pre2_Upd_long %>% group_by(Visit) %>% identify_outliers(UPDRS3) # Outliers assumption: Determine if there are extreme outliers
ggqqplot(pre2_Upd_long, "UPDRS3", facet.by = "Visit") # Normality assumption: because sample size >50, the normal QQ plot is used (Shapiro-Wilk test becomes very sensitive even to a minor deviation from normality)
model_UPD_pre2 <- aov(UPDRS3 ~ SR_class*Visit + Error(Subj_ID) + DisDur + Age_cov + Sex, data=pre2_Upd_long)
summary(model_UPD_pre2) 

# Pairwise comparisons
pwc <- pre2_Upd_long %>% pairwise_t_test(UPDRS3 ~ Visit, paired = TRUE, p.adjust.method = "bonferroni")
data.frame(pwc) #Based on the rm-ANOVA, all the tested pairwise differences between time points are statistically significant.
pwc <- pwc %>% add_xy_position(x = "Visit")

# Test if UPDRS delta visit 3-2 is larger than visit 2-1
Upd_wide$delta_upd_3_2 <- Upd_wide$UPDRS3.3 - Upd_wide$UPDRS3.2
Upd_wide$delta_upd_2_1 <- Upd_wide$UPDRS3.2 - Upd_wide$UPDRS3.1
delta <- Upd_wide %>% filter(precov == 2) # group with 2 visits before covid (151 patients)
t.test(delta$delta_upd_3_2, delta$delta_upd_2_1, paired = TRUE, alternative = "greater") # mean diff = 1.50 ; p-value = 0.096

# Plot progression for the groups of participants
Upd_pre2 <- ggplot(data=pre2_Upd_long, aes(x=Visit, y=UPDRS3)) +
                  geom_boxplot(fill="#008B8B", color="black") +
                  geom_jitter(color="black", size=1, alpha=0.4) +
                  xlab("POM Visit") + 
                  ylab("UPDRS-3") +
                  ggtitle("PD progression in group with 2 pre-covid visit") +
                  stat_pvalue_manual(pwc) +
                  labs(subtitle = get_test_label(test, detailed = TRUE),
                       caption = get_pwc_label(pwc)) +
                  theme_minimal(base_size = 16)

ggsave("UPDRS3_change_pre2.tiff", plot = Upd_pre2, device='tiff', dpi=300, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )

############## UPDRS3 progression: rm-ANOVA comparing the 2 stressor-reactivity classes for 2 pre-COVID visits separately ###############

# Summary stats and outlier detection: 
# Of the 151 with 2 pre-COVID visits, 1 person has missing UPDRS3 for visit 2+3, 10 missing UPDRS3 only for visit 3 (unable to do UPDRS3 in off).
# In total, that leaves 150-11=139 complete cases to use in analysis
pre2_Upd_long %>%
  group_by(SR_class,Visit) %>%
  get_summary_stats(UPDRS3, type = "mean_sd") # UPDRS clearly higher for resilience class 2 (1=high resilience, 2=low resilience)

pre2_Upd_long %>%
  group_by(SR_class,Visit) %>%
  identify_outliers(UPDRS3)    # No extreme outliers

# Posthoc tests effect res class at each timepoint: effect resilience class sign. for all visits
pwc1 <- pre2_Upd_long %>% group_by(Visit) %>% pairwise_t_test(UPDRS3 ~ SR_class, p.adjust.method="bonferroni")

# Posthoc tests time for each res class: for both high and low res class sign. effect of visit
pwc2 <- pre2_Upd_long %>% group_by(SR_class) %>% pairwise_t_test(UPDRS3 ~ Visit, paired=TRUE, p.adjust.method="bonferroni")

# Boxplot for all 3 visits divided by resilience class
p <- ggplot(data=pre2_Upd_long, aes(x=SR_class, y=UPDRS3, fill=interaction(SR_class, Visit), alpha=Visit)) + # resilience: 1=high, 2=low
            geom_boxplot() +
            scale_fill_manual(values=c("#009988", "#CC3311", "#009988", "#CC3311", "#009988", "#CC3311")) +
            scale_alpha_manual(values=c(1, 0.5, 0.2)) +
            geom_point(position=position_jitterdodge(),alpha=0.4) +          
            scale_x_discrete(labels=c("\n \n \n Low SR","\n \n \n High SR")) +
            scale_y_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80)) +
            xlab("") +
            ylab("Motor symptoms (MDS-UPDRS III)") +
            theme_minimal()
Upd_res2 <- p + theme(legend.position = "none", axis.text = element_text(size = 17, colour = "black"), axis.title.y = element_text(size = 18, margin = margin(r=10)))
print(Upd_res2)

# Save boxplot (figure 5d in manuscript)
ggsave("UPDRS3_pre2.svg", plot = Upd_res2, width=5, height=5, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )


### UPDRS3 progression: rm-ANOVA comparing the stressor-reactivity classes for all complete cases together ###

# Make a delta LEDD variable to use as covariate
LEDD <- POM_data_all[c("Subj_ID", "Visit", "precov", "LEDD")] 
LEDD$Visit <- as.integer(LEDD$Visit)
LEDD_wide <- pivot_wider(LEDD, names_from = Visit, values_from = LEDD)
names(LEDD_wide)[names(LEDD_wide) == "1"] <- "LEDD.1"
names(LEDD_wide)[names(LEDD_wide) == "2"] <- "LEDD.2"
names(LEDD_wide)[names(LEDD_wide) == "3"] <- "LEDD.3"
LEDD_wide <- LEDD_wide %>% filter(precov == 1 | precov ==2) # group with 1 or visits before covid
LEDD_wide$LEDD_1   <- ifelse(LEDD_wide$precov==1, LEDD_wide$LEDD.1, LEDD_wide$LEDD.2) #For precov==1 use visit 1 as pre, voor precov==2 (all others) use visit 2
LEDD_wide$LEDD_2  <- ifelse(LEDD_wide$precov==1, LEDD_wide$LEDD.2, LEDD_wide$LEDD.3) #For precov==1 use visit 2 as post, voor precov==2 (all others) use visit 3
LEDD_wide$dLEDD <- LEDD_wide$LEDD_2 - LEDD_wide$LEDD_1

# Select complete cases with ALL UPDRS3 scores and with either 1 or 2 pre-covid visits
all_UPD <- Upd_wide %>% filter(precov == 1 | precov ==2) # group with 1 or 2 visits before covid
all_UPD <- cbind(all_UPD, LEDD_wide[8]) # add delta LEDD variable 

# Make new variables for pre-covid UPDRS3, post-covid UPDRS3, and delta UPDRS3(post-pre). Also a simple binary resilience variable
all_UPD$UPDRS_1   <- ifelse(all_UPD$precov==1, all_UPD$UPDRS3.1, all_UPD$UPDRS3.2) #For precov==1 use UPDRS3 of visit 1, voor precov==2 (all others) use UPDRS3 of visit 2
all_UPD$UPDRS_2  <- ifelse(all_UPD$precov==1, all_UPD$UPDRS3.2, all_UPD$UPDRS3.3) #For precov==1 use UPDRS3 of visit 2, voor precov==2 (all others) use UPDRS3 of visit 3
all_UPD = all_UPD[,!(names(all_UPD) %in% c("UPDRS3.1","UPDRS3.2", "UPDRS3.3"))]

# Check if all is complete
all_UPD <- all_UPD %>% filter(!is.na(all_UPD$UPDRS_1), !is.na(all_UPD$UPDRS_2)) # 258 complete cases
sum(is.na(all_UPD$UPDRS_1)) # moet 0 zijn
sum(is.na(all_UPD$UPDRS_2)) # moet 0 zijn
all_UPD$dLEDD <- na.aggregate(all_UPD$dLEDD) # mean imputation of LEDD covariate (5 cases)

# Put back in long format
all_UPD_long <- pivot_longer(all_UPD, cols = starts_with("UPDRS_"), names_to = "covtime", values_to = "UPDRS3")
all_UPD_long$covtime <- as.factor(all_UPD_long$covtime)
all_UPD_long$UPDRS3 <- as.double(all_UPD_long$UPDRS3)

# Test difference between resilience classes
summary <-all_UPD_long %>%
          group_by(SR_class,covtime) %>%
          get_summary_stats(UPDRS3, type = "mean_sd")
data.frame(summary) # UPDRS clearly higher for resilience class 2 (1=high resilience, 2=low resilience), but progression seems slower...

# Test assumptions: 
all_UPD_long %>% group_by(covtime) %>% identify_outliers(UPDRS3)       # Identify outliers (no extreme outliers)   
ggqqplot(all_UPD_long, "UPDRS3") + facet_grid(covtime ~ SR_class)     # All the points fall along reference line, so we assume normality
all_UPD_long %>% group_by(covtime) %>% levene_test(UPDRS3 ~ SR_class) # Homogeneity of variance assumption of between subj. factor
box_m(all_UPD_long[, "UPDRS3", drop = FALSE], all_UPD_long$SR_class)    # Homogeneity of covariances (p=0.69)

# Test 2x2 mixed design ANOVA with factors resilience class and time
all_UPD_long$Subj_ID <- as.factor(all_UPD_long$Subj_ID)
model_UPD <- aov(UPDRS3 ~ SR_class*covtime + Error(Subj_ID/covtime) + DisDur + dLEDD + Age_cov + Sex, data=all_UPD_long)
summary(model_UPD)

# Fit linear mixed-effects model with continuous score for SR
all_UPD_long$IDnr <- match(all_UPD_long$Subj_ID, unique(all_UPD_long$Subj_ID))
model_UPD_con <- lmer(UPDRS3 ~ SR * covtime + dLEDD + DisDur + Age_cov + Sex + (1|IDnr), #Random error comes from the fact that it is a within subjects design
                  data = all_UPD_long, na.action = na.omit)
anova(model_UPD_con)
UPD_pred <- ggpredict(model_con, terms = c("covtime")) #, condition = c(Sex = "Male", LivingSituat = 2)

# Make boxplot for 2x2 design 
p <- ggplot(data=all_UPD_long, aes(x=SR_class, y=UPDRS3, fill=interaction(SR_class, covtime), alpha=covtime)) + # resilience: 1=high, 2=low
            geom_boxplot() +
            scale_fill_manual(values=c("#009988", "#CC3311", "#009988", "#CC3311")) +
            scale_alpha_manual(values=c(1, 0.4, 1, 0.4)) +
            geom_point(position=position_jitterdodge(),alpha=0.4) +
            scale_x_discrete(labels=c("\n \n Low SR","\n \n High SR")) +
            scale_y_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80)) +
            xlab("") +
            ylab("Motor symptoms (MDS-UPDRS III)") +
            theme_minimal()
Upd_resall <- p + theme(legend.position = "none", axis.text = element_text(size = 17, colour = "black"), axis.title.y = element_text(size = 18, margin = margin(r=10)))
print(Upd_resall)

# Save boxplot (figure 5a in manuscript)
ggsave(file = "UPDRS3_prepost.svg", plot = Upd_resall, width=5, height=5, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )


# Extra: Check if number of days between POM visit and start COVID pandemic is related to accelerated UPDRS decrease
Perc <- POM_data_all[c("Subj_ID", "Visit", "Visit_date", "Age_cov", "Sex", "precov", "UPDRS3_off")] 
Perc$Visit <- as.integer(Perc$Visit)
Perc <- pivot_wider(Perc, names_from = Visit, values_from = c("UPDRS3_off","Visit_date"))
Perc$Res <- Covid$Resilience
Perc <- Perc %>% filter(precov == 1 | precov == 2) # 323 cases
Perc <- Perc %>% filter(!is.na(Perc$UPDRS3_off_1), !is.na(Perc$UPDRS3_off_2), !is.na(Perc$UPDRS3_off_3)) # 257 complete cases
Perc$days <- ifelse(Perc$precov==2, (as.Date("2020-04-21") - Perc$Visit_date_2), (as.Date("2020-04-21") - Perc$Visit_date_1))

Perc$UPDRS_1   <- ifelse(Perc$precov==1, Perc$UPDRS3_off_1, Perc$UPDRS3_off_2)  # For precov==1 use UPDRS3 of visit 1, voor precov==2 (all others) use UPDRS3 of visit 2
Perc$UPDRS_2  <- ifelse(Perc$precov==1, Perc$UPDRS3_off_2, Perc$UPDRS3_off_3)   # For precov==1 use UPDRS3 of visit 2, voor precov==2 (all others) use UPDRS3 of visit 3
Perc$UPDRSdelta <- Perc$UPDRS_2 - Perc$UPDRS_1                                  # <0: 65, >0: 165, 0: 7
Perc = Perc[,!(names(Perc) %in% c("UPDRS3_off_1","UPDRS3_off_2", "UPDRS3_off_3"))]

cor(Perc$days, Perc$UPDRSdelta,  method = "pearson") # Hypothesis: the more days between visits overlap with covid-pandemic, the higher delta UPDRS3
cor(Perc$Res, Perc$UPDRSdelta,  method = "pearson")

p <- ggscatter(Perc, x = "UPDRSdelta", y = "UPDRSdelta", 
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "Days of covid between visits", ylab = "delta UPDRS-III (post-covid - pre-covid)") +
              scale_y_continuous(breaks = c(-20, -10, 0, 10, 20, 30))
p + theme(axis.text = element_text(size = 16), axis.title = element_text(size = 17))

############## Relationship STAI progression with stressor-reactivity ###############

POM_STAI <- POM_data_all[c("Subj_ID", "Visit", "Age_cov", "Sex", "DisDur", "precov", "STAI_t", "STAI_s")] #Variables needed for comparison resilience and STAI change
POM_STAI$STAI <- POM_STAI$STAI_t + POM_STAI$STAI_s
POM_STAI = POM_STAI[,!(names(POM_STAI) %in% c("STAI_t","STAI_s"))]

POM_STAI$Visit <- as.integer(POM_STAI$Visit)
STAI_wide <- pivot_wider(POM_STAI, names_from = Visit, values_from = STAI)
names(STAI_wide)[names(STAI_wide) == "1"] <- "STAI.1"
names(STAI_wide)[names(STAI_wide) == "2"] <- "STAI.2"
names(STAI_wide)[names(STAI_wide) == "3"] <- "STAI.3"

# Determine who misses STAI for each visit
sum(is.na(STAI_wide$STAI.1))                           # 10 with missing STAI for visit 1
sum(is.na(STAI_wide$STAI.2))                           # 68 with missing STAI for visit 2
sum(is.na(STAI_wide$STAI.3))                           # 41 with missing STAI for visit 3

# Combine data and make separate groups for 1 or 2 pre-covid visits
STAI_wide <- cbind(STAI_wide, Covid$SR, Covid$SR_class)
pre1_STAI <- STAI_wide %>% filter(precov == 1, !is.na(STAI_wide$STAI.1), !is.na(STAI_wide$STAI.2)) # group with 1 visit before covid
pre2_STAI <- STAI_wide %>% filter(precov == 2, !is.na(STAI_wide$STAI.2), !is.na(STAI_wide$STAI.3)) # group with 2 visits before covid

# Put separate groups with 1 or 2 precovid visits back to long format to have one STAI variable
pre2_STAI_long <- pivot_longer(na.omit(pre2_STAI), cols = starts_with("STAI."), names_to = "Visit", values_to = "STAI")
pre2_STAI_long$Visit <- as.factor(pre2_STAI_long$Visit)
levels(pre2_STAI_long$Visit) <- c("1", "2", "3")

# summary stats
pre2_STAI_long %>%
  group_by(SR_class,Visit) %>%
  get_summary_stats(STAI, type = "mean_sd") # STAI clearly higher for resilience class 2 (1=high resilience, 2=low resilience)

# test for latent classes of resilience
model_STAI_pre2 <- aov(STAI ~ SR_class*Visit + Error(Subj_ID) + DisDur + Age_cov + Sex, data=pre2_STAI_long)
summary(model_STAI_pre2)

# posthoc tests effect res class at each timepoint
one.way1 <- pre2_STAI_long %>%
  group_by(Visit) %>%
  anova_test(dv = STAI, wid = Subj_ID, between = SR_class) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way1 # the simple main effect of resilience class is sign. for all visits

pwc1 <- pre2_STAI_long %>% group_by(Visit) %>% pairwise_t_test(STAI ~ SR_class, p.adjust.method="bonferroni")
pwc1 # Pairwise comparisons show that the mean STAI score was significantly different for all pairs

# posthoc tests time for each res class
one.way2 <- pre2_STAI_long %>%
  group_by(SR_class) %>%
  anova_test(dv = STAI, wid = Subj_ID, within = Visit, covariat = c(Age_cov, Sex, DisDur)) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2 # Effect of visit only sign. for the high resilient group

pwc2 <- pre2_STAI_long %>% group_by(SR_class) %>% pairwise_t_test(STAI ~ Visit, paired=TRUE, p.adjust.method="bonferroni")
pwc2 # The pairwise comparisons visit 1 vs 2 and visit 1 vs 3 were statistically significantly different only for the high resilient group

# Boxplot for all 3 visits divided by resilience class
p <- ggplot(data=pre2_STAI_long, aes(x=SR_class, y=STAI, fill=interaction(SR_class, Visit), alpha=Visit)) + # resilience: 1=high, 2=low
  geom_boxplot() +
  scale_fill_manual(values=c("#009988", "#CC3311", "#009988", "#CC3311", "#009988", "#CC3311")) +
  scale_alpha_manual(values=c(1, 0.5, 0.2)) +
  geom_point(position=position_jitterdodge(),alpha=0.4) +          
  scale_x_discrete(labels=c("\n \n Low SR","\n \n High SR")) +
  scale_y_continuous(breaks = c(40, 60, 80, 100, 120, 140)) +
  xlab("") +
  ylab("State and trait anxiety (STAI)") +
  theme_minimal()
STAI_res2 <- p + theme(legend.position = "none", axis.text = element_text(size = 17, colour = "black"), axis.title.y = element_text(size = 18, margin = margin(r=10)))
print(STAI_res2)

# Save boxplot (figure 5f in manuscript)
ggsave("STAI_pre2.svg", plot = STAI_res2, width=4.5, height=6, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )


### STAI progression: rm-ANOVA comparing the stressor-reactivity classes for all complete cases together ###
all_STAI <- rbind(pre1_STAI, pre2_STAI) # 240 cases

# Make new variables for pre-covid STAI, post-covid STAI, and delta STAI(post-pre). Also a simple binary resilience variable
all_STAI$STAI_1 <- ifelse(all_STAI$precov==1, all_STAI$STAI.1, all_STAI$STAI.2) #For precov==1 use STAI of visit 1, voor precov==2 (all others) use STAI of visit 2
all_STAI$STAI_2 <- ifelse(all_STAI$precov==1, all_STAI$STAI.2, all_STAI$STAI.3) #For precov==1 use STAI of visit 2, voor precov==2 (all others) use STAI of visit 3
all_STAI$STAIdelta <- all_STAI$STAI_2 - all_STAI$STAI_1
all_STAI = all_STAI[,!(names(all_STAI) %in% c("STAI.1","STAI.2", "STAI.3"))]

# Put back in long format
all_STAI_long <- pivot_longer(all_STAI, cols = starts_with("STAI_"), names_to = "covtime", values_to = "STAI")
all_STAI_long$covtime <- as.factor(all_STAI_long$covtime)
all_STAI_long$STAI <- as.double(all_STAI_long$STAI)

# Get summary table to explore difference between SR groups
summary<-all_STAI_long %>%
  group_by(SR_class,covtime) %>%
  get_summary_stats(STAI, type = "mean_sd")
data.frame(summary) # STAI clearly higher for patients with high SR, but decrease from pre- to in-COVID in low SR group

# Test assumptions
out <- all_STAI_long %>% group_by(SR_class,covtime) %>% identify_outliers(STAI)  # No extreme outliers
ggqqplot(all_STAI_long, "STAI") + facet_grid(covtime ~ SR_class)                 # All the points fall along reference line, so we assume normality
all_STAI_long %>% group_by(covtime) %>% levene_test(STAI ~ SR_class)             # Test homoscedasticity: significant
all_STAI_long %>% group_by(covtime) %>% levene_test(log(STAI) ~ SR_class)        

# Test 2x2 mixed design ANOVA with factors SR class and time
all_STAI_long$Subj_ID <- as.factor(all_STAI_long$Subj_ID)
model_STAI <- aov(STAI ~ SR_class*covtime + Error(Subj_ID/covtime) + DisDur + Age_cov + Sex, data=all_STAI_long)
model_STAI_log <- aov(log(STAI) ~ SR_class*covtime + Error(Subj_ID/covtime) + DisDur + Age_cov + Sex, data=all_STAI_long)
summary(model_STAI)

# Fit linear mixed-effects model with continuous score for SR (supplementary analysis)
all_STAI_long$IDnr <- match(all_STAI_long$Subj_ID, unique(all_STAI_long$Subj_ID))
model_STAI_con <- lmer(STAI ~ SR * covtime + DisDur + Age_cov + Sex + (1|IDnr), #Random error comes from the fact that it is a within subjects design
                        data = all_STAI_long, na.action = na.omit)
anova(model_STAI_con)

# Plot with 2x2 boxplots for high versus low SR
p <- ggplot(data=all_STAI_long, aes(x=SR_class, y=STAI, fill=interaction(SR_class, covtime), alpha=covtime)) + # SR: 1=low, 2=high
  geom_boxplot() +
  scale_fill_manual(values=c("#009988", "#CC3311", "#009988", "#CC3311")) +
  scale_alpha_manual(values=c(1, 0.4, 1, 0.4)) +
  geom_point(position=position_jitterdodge(),alpha=0.4) +
  scale_x_discrete(labels=c("\n \n Low SR","\n \n High SR")) +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100, 120, 140)) +
  xlab("") +
  ylab("State and trait anxiety (STAI)") +
  theme_minimal()
STAI_resall <- p + theme(legend.position = "none", axis.text = element_text(size = 17, colour = "black"), axis.title.y = element_text(size = 18, margin = margin(r=5)))
print(STAI_resall)

# Save boxplot (figure 5c in manuscript)
ggsave(file = "STAI_prepost.svg", plot = STAI_resall, width=4.5, height=6, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )


############## Relationship BDI progression with stressor-reactivity ###############

POM_BDI <- POM_data_all[c("Subj_ID", "Visit", "Age_cov", "Sex", "DisDur", "precov", "BDI")] #Variables needed for comparison resilience and BDI change

POM_BDI$Visit <- as.integer(POM_BDI$Visit)
BDI_wide <- pivot_wider(POM_BDI, names_from = Visit, values_from = BDI)
names(BDI_wide)[names(BDI_wide) == "1"] <- "BDI.1"
names(BDI_wide)[names(BDI_wide) == "2"] <- "BDI.2"
names(BDI_wide)[names(BDI_wide) == "3"] <- "BDI.3"

# Determine who misses BDI for each visit
sum(is.na(BDI_wide$BDI.1))                           # 10 with missing BDI for visit 1
sum(is.na(BDI_wide$BDI.2))                           # 68 with missing BDI for visit 2
sum(is.na(BDI_wide$BDI.3))                           # 41 with missing BDI for visit 3

# Combine data and make separate groups for 1 or 2 pre-covid visits
BDI_wide <- cbind(BDI_wide, Covid$SR, Covid$SR_class)
pre1_BDI <- BDI_wide %>% filter(precov == 1, !is.na(BDI_wide$BDI.1), !is.na(BDI_wide$BDI.2)) # group with 1 visit before covid
pre2_BDI <- BDI_wide %>% filter(precov == 2, !is.na(BDI_wide$BDI.2), !is.na(BDI_wide$BDI.3)) # group with 2 visits before covid

# Put separate group with 2 precovid visits back to long format to have one BDI variable
pre2_BDI_long <- pivot_longer(na.omit(pre2_BDI), cols = starts_with("BDI."), names_to = "Visit", values_to = "BDI")
pre2_BDI_long$Visit <- as.factor(pre2_BDI_long$Visit)
levels(pre2_BDI_long$Visit) <- c("1", "2", "3")

# summary stats
pre2_BDI_long %>%
  group_by(SR_class,Visit) %>%
  get_summary_stats(BDI, type = "mean_sd") # BDI clearly higher for resilience class 2 (1=high resilience, 2=low resilience)

# test for latent classes of resilience
pre2_BDI_long = pre2_BDI_long[!pre2_BDI_long$Subj_ID =="sub-POMU900F78E54F00A78A",] # remove extreme outlier
model_BDI_pre2 <- aov(BDI ~ SR_class*Visit + Error(Subj_ID) + DisDur + Age_cov + Sex, data=pre2_BDI_long)
summary(model_BDI_pre2)

# posthoc tests effect res class at each timepoint
one.way1 <- pre2_BDI_long %>%
            group_by(Visit) %>%
            anova_test(dv = BDI, wid = Subj_ID, between = SR_class) %>%
            get_anova_table() %>%
            adjust_pvalue(method = "bonferroni")
one.way1 # the simple main effect of resilience class is sign. for all visits

pwc1 <- pre2_BDI_long %>% group_by(Visit) %>% pairwise_t_test(BDI ~ SR_class, p.adjust.method="bonferroni")
pwc1 # Pairwise comparisons show that the mean BDI score was significantly different for all pairs

# posthoc tests time for each res class
one.way2 <- pre2_BDI_long %>%
            group_by(SR_class) %>%
            anova_test(dv = BDI, wid = Subj_ID, within = Visit) %>%
            get_anova_table() %>%
            adjust_pvalue(method = "bonferroni")
one.way2 # Effect of visit only sign. for the low resilient group

pwc2 <- pre2_BDI_long %>% group_by(SR_class) %>% pairwise_t_test(BDI ~ Visit, paired=TRUE, p.adjust.method="bonferroni")
pwc2 # The pairwise comparisons visit 1 vs 3 and visit 2 vs 3 were statistically significantly different only for the low resilient group

BDI_wide$delta_BDI_3_2 <- BDI_wide$BDI.3 - BDI_wide$BDI.2
BDI_wide$delta_BDI_2_1 <- BDI_wide$BDI.2 - BDI_wide$BDI.1
delta <- BDI_wide %>% filter(precov == 2) # group with 2 visits before covid (150 patients)

# Boxplot for all 3 visits divided by resilience class
p <- ggplot(data=pre2_BDI_long, aes(x=SR_class, y=BDI, fill=interaction(SR_class, Visit), alpha=Visit)) + # resilience: 1=high, 2=low
            geom_boxplot() +
            scale_fill_manual(values=c("#009988", "#CC3311", "#009988", "#CC3311", "#009988", "#CC3311")) +
            scale_alpha_manual(values=c(1, 0.5, 0.2)) +
            geom_point(position=position_jitterdodge(),alpha=0.4) +          
            scale_x_discrete(labels=c("\n \n Low SR","\n \n High SR")) +
            scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30, 35)) +
            xlab("") +
            ylab("Depressive symptoms (BDI)") +
            theme_minimal()
BDI_res2 <- p + theme(legend.position = "none", axis.text = element_text(size = 17, colour = "black"), axis.title.y = element_text(size = 18, margin = margin(r=10)))
print(BDI_res2)

# Save boxplot (figure 5e in manuscript)
ggsave("BDI_pre2.svg", plot = BDI_res2, width=5, height=5, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )


### BDI progression: rm-ANOVA comparing the stressor-reactivity classes for all complete cases together ###
all_BDI <- rbind(pre1_BDI, pre2_BDI) # 239 cases

# Make new variables for pre-covid BDI, post-covid BDI, and delta BDI(post-pre). Also a simple binary resilience variable
all_BDI$BDI_1   <- ifelse(all_BDI$precov==1, all_BDI$BDI.1, all_BDI$BDI.2) #For precov==1 use BDI of visit 1, voor precov==2 (all others) use BDI of visit 2
all_BDI$BDI_2  <- ifelse(all_BDI$precov==1, all_BDI$BDI.2, all_BDI$BDI.3) #For precov==1 use BDI of visit 2, voor precov==2 (all others) use BDI of visit 3
all_BDI = all_BDI[,!(names(all_BDI) %in% c("BDI.1","BDI.2", "BDI.3"))]

# Put back in long format
all_BDI_long <- pivot_longer(all_BDI, cols = starts_with("BDI_"), names_to = "covtime", values_to = "BDI")
all_BDI_long$covtime <- as.factor(all_BDI_long$covtime)
all_BDI_long$BDI <- as.double(all_BDI_long$BDI)

# Summary statistics per group
summary<- all_BDI_long %>%
          group_by(SR_class,covtime) %>%
          get_summary_stats(BDI, type = "mean_sd")
data.frame(summary) # BDI clearly higher for resilience class 2 (1=high resilience, 2=low resilience), but progression seems slower...

# Test assumptions for analysis
all_BDI_long %>% group_by(covtime) %>% identify_outliers(BDI)       # Identify outliers (1 extreme outliers) 
all_BDI_long = all_BDI_long[!all_BDI_long$Subj_ID =="sub-POMU900F78E54F00A78A",] # remove 1 extreme outlier 
ggqqplot(all_BDI_long, "BDI") + facet_grid(covtime ~ SR_class)     # All the points fall along reference line, so we assume normality
all_BDI_long %>% group_by(covtime) %>% levene_test(BDI ~ SR_class) # Homogeneity of variance assumption of between subj. factor: SIGN!!
box_m(all_BDI_long[, "BDI", drop = FALSE], all_BDI_long$SR_class)    # Homogeneity of covariances: SIGN!!

# Test 2x2 mixed design ANOVA with factors resilience class and time
all_BDI_long$Subj_ID <- as.factor(all_BDI_long$Subj_ID)
model_BDI <- aov(BDI ~ SR_class*covtime + Error(Subj_ID/covtime) + DisDur + Age_cov + Sex, data=all_BDI_long)
#model_BDI_log <- aov(log(BDI) ~ SR_class*covtime + Error(Subj_ID/covtime) + DisDur + Age_cov + Sex, data=all_BDI_long)
summary(model_BDI)

# Post hoc test for effect time
posthoc <- all_BDI_long %>%
            group_by(SR_class) %>%
            anova_test(dv = BDI, wid = Subj_ID, within = covtime) %>%
            get_anova_table() %>%
            adjust_pvalue(method = "bonferroni")

pwc <- all_BDI_long %>% 
        group_by(SR_class) %>%
        emmeans_test(BDI ~ covtime, p.adjust.method = "bonferroni") 
pwc <- pwc %>% add_xy_position(x = "covtime")

# Fit linear mixed-effects model with continuous score for SR (supplementary analysis)
all_BDI_long$IDnr <- match(all_BDI_long$Subj_ID, unique(all_BDI_long$Subj_ID))
model_BDI_con <- lmer(BDI ~ SR * covtime + DisDur + Age_cov + Sex + (1|IDnr), #Random error comes from the fact that it is a within subjects design
                      data = all_BDI_long, na.action = na.omit)
anova(model_BDI_con)

posthoc <- emmeans(model_BDI_con, pairwise ~ covtime | SR)

# Plot 2x2 boxplots
p <- ggplot(data=all_BDI_long, aes(x=SR_class, y=BDI, fill=interaction(SR_class, covtime), alpha=covtime)) + # resilience: 1=high, 2=low
  geom_boxplot() +
  scale_fill_manual(values=c("#009988", "#CC3311", "#009988", "#CC3311")) +
  scale_alpha_manual(values=c(1, 0.4, 1, 0.4)) +
  geom_point(position=position_jitterdodge(),alpha=0.4) +
  scale_x_discrete(labels=c("\n \n Low SR","\n \n High SR")) +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30, 35)) +
  xlab("") +
  ylab("Depressive symptoms (BDI)") +
  theme_minimal()

BDI_resall <- p + theme(legend.position = "none", axis.text = element_text(size = 17, colour = "black"), axis.title.y = element_text(size = 18, margin = margin(r=10)))
print(BDI_resall)

# Save boxplot (figure 5b in manuscript)
ggsave(file = "BDI_prepost.svg", plot = BDI_resall, width=5, height=5, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )


############## Extra: Change start vs. end covid survey in PASS, SOZU and RRS ################

All <- read_csv("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/Covid_all_wide.csv")
names(Covid)[names(Covid) == "Age_cov"] <- "Age"
SoZU <- cbind(Covid[1:3], POM_data_all[10], All[46:47], Res_data[2]) # Select columns with Sozu
RRS <- cbind(Covid[1:3], All[48:49], Res_data[2]) # Select columns with RRS
PASS <- cbind(Covid[1:3], All[54:55], Res_data[2]) # Select columns with Sozu

SoZU_long <- pivot_longer(SoZU, cols = starts_with("SoZU_"), names_to = "time", values_to = "SoZU")
RRS_long <- pivot_longer(RRS, cols = starts_with("RRS_"), names_to = "time", values_to = "RRS")
PASS_long <- pivot_longer(PASS, cols = starts_with("PASS_"), names_to = "time", values_to = "PASS")

ano_PASS <- anova_test(PASS_long, dv = PASS, wid = Subj_ID, between = SR_class, within = time, covariate = c(Age, Sex))
ano_PASS # only sign. effect of resilience class
ano_Soc <- anova_test(SoZU_long, dv = SoZU, wid = Subj_ID, between = SR_class, within = time, covariate = c(Age, Sex))
ano_Soc # only sign. effect of resilience class, interaction sex*time
ano_RRS <- anova_test(RRS_long, dv = RRS, wid = Subj_ID, between = SR_class, within = time, covariate = c(Age, Sex))
ano_RRS # only sign. effect of resilience class

# Pairwise comparison social support with sex to explore sex*time interaction
pwc2 <- SoZU_long %>% 
        group_by(Sex) %>% 
        pairwise_t_test(SoZU ~ time, paired=TRUE, p.adjust.method="bonferroni")
pwc2 # for men social support went down and for women it went up, but this is for both sexes not significant.

# Plot change between baseline and final survey in positive appraisal style (PASS)
detach(package:plyr)
sum1 <- PASS_long %>%
  group_by(time, SR_class) %>%
  summarise(sd = sd(PASS, na.rm=TRUE), PASS = mean(PASS, na.rm=TRUE))

p<- ggplot(PASS_long, aes(y=PASS, x=time, color = SR_class)) +
          geom_point(alpha = 0.3, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.3)) +
          geom_line(aes(group = SR_class),data = sum1, size = 1.5, position=position_dodge(width=0.3)) +
          geom_errorbar(aes(ymin = PASS-sd, ymax = PASS+sd), data = sum1, width = 0.5, size = 1.5, position=position_dodge(width=0.3))+
          scale_color_manual(name = "Resilience class", 
                             values = c("#00AFBB", "#E7B800"),
                             breaks = c("1", "2"),
                             labels = c("High resilience","Low resilience")) +
          scale_x_discrete(labels=c("Baseline","Final")) +
          xlab("") + ylab("Positive appraisal style") + theme_minimal()
PASS <- p + theme(legend.key.size = unit(1.5, "cm"), legend.key.width = unit(1.5,"cm"),
                 legend.text = element_text(size = 16), legend.title = element_text(size = 16, face = "bold"),
                 axis.text = element_text(size = 16), axis.title.y = element_text(size = 17, margin = margin(r=10)))
print(PASS)

# Plot change between baseline and final survey in rumination (RRS)
sum2 <- RRS_long %>%
  group_by(time, SR_class) %>%
  summarise(sd = sd(RRS, na.rm=TRUE), RRS = mean(RRS, na.rm=TRUE))
sum2

p<- ggplot(RRS_long, aes(y=RRS, x=time, color = SR_class)) +
          geom_point(alpha = 0.3, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.3)) +
          geom_line(aes(group = SR_class),data = sum2, size = 1.5, position=position_dodge(width=0.3)) +
          geom_errorbar(aes(ymin = RRS-sd, ymax = RRS+sd), data = sum2, width = 0.5, size = 1.5, position=position_dodge(width=0.3))+
          scale_color_manual(name = "Resilience class", 
                             values = c("#00AFBB", "#E7B800"),
                             breaks = c("1", "2"),
                             labels = c("High resilience","Low resilience")) +
          scale_x_discrete(labels=c("Baseline","Final")) +
          xlab("") + ylab("Ruminative response") + theme_minimal()
RRS <- p + theme(legend.key.size = unit(1.5, "cm"), legend.key.width = unit(1.5,"cm"),
                 legend.text = element_text(size = 16), legend.title = element_text(size = 16, face = "bold"),
                 axis.text = element_text(size = 16), axis.title.y = element_text(size = 17, margin = margin(r=10)))
print(RRS)

# Plot change between baseline and final survey in perceived social support (SoZU)
sum3 <- SoZU_long %>%
        group_by(time, SR_class) %>%
        summarise(sd = sd(SoZU, na.rm=TRUE), SoZU = mean(SoZU, na.rm=TRUE))

p<- ggplot(SoZU_long, aes(y=SoZU, x=time, color = SR_class)) +
          geom_point(alpha = 0.3, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.3)) +
          geom_line(aes(group = SR_class),data = sum3, size = 1.5, position=position_dodge(width=0.3)) +
          geom_errorbar(aes(ymin = SoZU-sd, ymax = SoZU+sd), data = sum3, width = 0.5, size = 1.5, position=position_dodge(width=0.3))+
          scale_color_manual(name = "Resilience class", 
                             values = c("#00AFBB", "#E7B800"),
                             breaks = c("1", "2"),
                             labels = c("High resilience","Low resilience")) +
          scale_x_discrete(labels=c("Baseline","Final")) +
          xlab("") + ylab("Social support") + theme_minimal()
Soc <- p + theme(legend.key.size = unit(1.5, "cm"), legend.key.width = unit(1.5,"cm"),
                 legend.text = element_text(size = 16), legend.title = element_text(size = 16, face = "bold"),
                 axis.text = element_text(size = 16), axis.title.y = element_text(size = 17, margin = margin(r=10)))
print(Soc)

ggsave("PASS_change.jpg", plot = PASS, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )
ggsave("SoZU_change.jpg", plot = Soc, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )
ggsave("RRS_change.jpg", plot = RRS, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )

