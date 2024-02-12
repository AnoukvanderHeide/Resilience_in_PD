# This script calculates COVID-related stressor-reactivity from all PSS and SL scores over all time points. 
# The script also produces long and wide data structure of covid survey data, and a data file with average scores 
# to be used in subsequent analysis steps.

rm(list = ls())

library(readr)

# Load all data files: POM_precovid_visit.csv is made in Preprocessing_POM_visit_data.R and Covid_longitudinal.csv is created in Matlab

POM_data <- read_csv("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/POM_precovid_visit.csv",
                     col_types = cols(Sex = col_factor(levels = c("Male", "Female")), 
                                      # 1=male, 2=female
                                      LivingSituat = col_factor(levels = c("1", "2", "3", "4", "5")), 
                                      # 1=alone, 2=with partner, 3=with family, 4=nursing home, 5=other
                                      DailyActivity = col_factor(levels = c("1", "2", "3", "4", "5", "6", "7", "8"))))
                                      # 1=paid job, 2=household, 3=retired, 4=student, 5=no paid job due to illness, 
                                      # 6=involuntarily no job, 7=volunteer work, 8=other)

Covid_data <- read_csv("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/Covid_longitudinal.csv")
names(Covid_data)[names(Covid_data) == "name"] <- "Subj_ID"

is.nan.data.frame <- function(x) # In the Matlab file NaN is used for empty cells, in R this has to be changed to NA (NaN: not a number, NA: not applicable)
    do.call(cbind, lapply(x, is.nan))
Covid_data[is.nan(Covid_data)] <- NA

library(magrittr)
library(tidyr)

Weeknr <- read_csv("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/Covid_weeknrs.csv")
Weeknr[is.nan(Weeknr)] <- NA
Weeknr_long <- Weeknr %>% pivot_longer(cols = starts_with("wk_"), names_to = "timepoint", values_to = "Weeknr", values_drop_na = TRUE) # put weeknrs in long format

#PSS <-  Covid_data[c("PSS_1", "PSS_2", "PSS_3", "PSS_4", "PSS_5", "PSS_6", "PSS_7", "PSS_8", "PSS_9", "PSS_10", "PSS_11")]
#PSS$Reps <- rowSums(!is.na(PSS)) # Count how many surveys were filled in (repetitions)

# Make long format of stressor exposure (SL) and perceived stress (PSS) to calculate normative regression line
data = subset(Covid_data, select = c(1:12,24:45))

data_long <- data %>%
             pivot_longer(cols = starts_with(c("PSS", "UPDRS", "SL")),
                          names_to = c(".value", "timepoint"),
                          names_sep = "_",
                          values_drop_na = TRUE)

# First check whether linear relationship outperforms non-linear relationship
set.seed(123)
training.samples <- data_long$PSS %>% createdata_longPartition(p = 0.8, list = FALSE)
train.data_long  <- data_long[training.samples, ]
test.data_long   <- data_long[-training.samples, ]

model_lin   <- lm(PSS ~ SL, data_long = train.data_long)
predictions <- model %>% predict(test.data_long)
data_long.frame(RMSE = RMSE(predictions, test.data_long$PSS),
           R2 = R2(predictions, test.data_long$PSS))

model_poly  <- lm(PSS ~ SL + I(SL^2), data_long = train.data_long)
predictions <- model_ploy %>% predict(test.data_long)
data_long.frame(RMSE = RMSE(predictions, test.data_long$PSS),
           R2 = R2(predictions, test.data_long$PSS))

model_log   <- lm(log(PSS) ~ SL, data_long = train.data_long)
predictions <- model_log %>% predict(test.data_long)
data_long.frame(RMSE = RMSE(predictions, test.data_long$PSS),
           R2 = R2(predictions, test.data_long$PSS))

library(splines)
knots       <- quantile(train.data_long$lstat, p = c(0.25, 0.5, 0.75))
model_spl   <- lm (PSS ~ bs(SL, knots = knots), data_long = train.data_long)
predictions <- model_spl %>% predict(test.data_long)
data_long.frame(RMSE = RMSE(predictions, test.data_long$PSS),
           R2 = R2(predictions, test.data_long$PSS))

# Calculate regression line
library(ggpubr)
ggplot(data_long, aes(x = SL, y = PSS)) +
       geom_point() +
       stat_smooth()

# Check for and remove extreme outlier(s)
library(rstatix)
out_SL  <- data_long %>% identify_outliers(SL) 
out_PSS <- data_long %>% identify_outliers(PSS) 
out_UPD <- data_long %>% identify_outliers(UPDRS) 
out_SL$Subj_ID[out_SL$is.extreme == T]    # sub-POMUC449B70F434BE59C has extreme value for SL
out_PSS$Subj_ID[out_PSS$is.extreme == T]  # No extreme outliers PSS
out_UPD$Subj_ID[out_UPD$is.extreme == T]  # No extreme outliers UPDRS
outval <- out_SL$SL[out_SL$is.extreme == T]
data_long["SL"][data_long["SL"] == outval] <- NA  # remove outlier SL
rm(out_SL, out_PSS, out_UPD, outval)

cor(data_long$PSS, data_long$SL, use = "pairwise.complete.obs", method = "pearson") # Correlation is 0.34
model <- lm(PSS ~ SL, data = data_long) # y = ax + b: PSS = 0.184*SL + 5.88

b <- model$coefficients[1] 
a <- model$coefficients[2]

# Check assumptions regression
plot(model, which = 2) # Q-Q plot: residuals roughly normally distributed
res <- resid(model)    # Get list of residuals
plot(density(res))     # Create density plot of residuals: bit left skewed

plot(model, which = 1)
library(lmtest)
lmtest::bptest(model) # Breusch-Pagan test: BP = 41.892, df = 1, p-value = 9.646e-11

# Calculate stressor-reactivity across all participants and all time points
Number = 1:11
Create_SR <- function(i){
  SL <- data[paste0("SL_",i)]
  PSS <- data[paste0("PSS_",i)]
  Covid_data[paste0("SR_",i)] <<- -((a*SL+b)-PSS)
  }

lapply(Number,Create_SR)

# Calculate average scores per subject: this should be redone with column names instead of numbers!
Covid_averages <- POM_data %>%
    select(Subj_ID, Sex, Age_cov) %>%
    mutate(PSSmean   = rowMeans(select(Covid_data, starts_with("PSS")), na.rm = TRUE),
           PASmean   = rowMeans(select(Covid_data, starts_with("PAS")), na.rm = TRUE),
           UPDRSmean = rowMeans(select(Covid_data, starts_with("UPDRS")), na.rm = TRUE),
           SLmean    = rowMeans(select(Covid_data, starts_with("SL")), na.rm = TRUE),
           SR        = rowMeans(select(Covid_data, starts_with("SR")), na.rm = TRUE),
           SoZU      = rowMeans(select(Covid_data, starts_with("SoZU")), na.rm = TRUE),
           RRS       = rowMeans(select(Covid_data, starts_with("RRS")), na.rm = TRUE),
           Symp      = rowMeans(select(Covid_data, starts_with("Symp")), na.rm = TRUE),
           BHC       = rowMeans(select(Covid_data, starts_with("BHC")), na.rm = TRUE),
           PASS      = rowMeans(select(Covid_data, starts_with("PASS")), na.rm = TRUE),
           BRS       = rowMeans(select(Covid_data, starts_with("BRS")), na.rm = TRUE),
           Opti      = rowMeans(select(Covid_data, starts_with("Opti")), na.rm = TRUE),
           Neuro     = rowMeans(select(Covid_data, starts_with("Neuro")), na.rm = TRUE),
           LE        = Covid_data$LE)

Alldata_wide <- merge(Covid_averages, subset(POM_data, select = -c(3,6,8)), by = "Subj_ID")
Covid_data <- cbind(Covid_data, Weeknr[2:12])

write.csv(Covid_averages,"M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/Covid_averages.csv", row.names = FALSE)
write.csv(Covid_data,"M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/Covid_all_wide.csv", row.names = FALSE)
write.csv(Alldata_wide,"M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/Covid_POM_wide.csv", row.names = FALSE)

# Make long data file, including pre-covid variables
# First select Subj_ID + all scores that are measured longitudinally.
Sel_cov <- Covid_data %>%
  select(Subj_ID, starts_with("PSS_"), starts_with("PAS_"), starts_with("UPDRS_"), starts_with("SL_"), starts_with("SR_"))

All <- merge(Sel_cov, POM_data, by = "Subj_ID")

# Make long data file using multiple column sets at once and also keep other variables.
data_long2 <- All %>%
              pivot_longer(cols = starts_with(c("PSS_", "PAS_", "UPDRS_", "SL_", "SR_")),
                           names_to = c(".value", "timepoint"),
                           names_sep = "_",
                           values_drop_na = TRUE)

# Check the association between PSS and UPDRS, and between PAS and PSS
ggplot(data = Covid_avg, aes(x = PSSmean, y = UPDRSmean)) + 
        geom_point(color = "#595959") +
        geom_smooth(method = lm, size = 1, se = F, color = "red")

ggplot(data = Covid_long, aes(x = PSS, y = UPDRS)) + 
        geom_point(color = "#595959") +
        geom_smooth(method = lm, size = 1, se = F, color = "red")

cor.test(Covid_long$UPDRS, Covid_long$PSS)
cor.test(Covid_long$PAS, Covid_long$PSS)

# Make one data file with the week numbers since baseline added as column
Covid_long <- cbind(data_long2, Weeknr_long[3]) # add column with the nr of weeks since baseline that surveys are completed
write.csv(Covid_long,"M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/Covid_POM_long.csv", row.names = FALSE)
