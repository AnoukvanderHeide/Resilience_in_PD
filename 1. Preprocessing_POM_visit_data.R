# This script loads POM visit data, restructures it and takes the last precovid-score for each measure of interest.
# The used input file is created with the script that are available at 
# https://github.com/mejoh/Personalized-Parkinson-Project-Motor/blob/master/R/generate_castor_csv.R

rm(list = ls())

############## Load and restructure POM data ###############

library(readr)

POM_all <- read_csv("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/merged_manipulated_2022-08-05.csv")

# Select only the variables that are relevant for analysis of covid data
POM_data <- POM_all[c("pseudonym",         "Timepoint",       "TimepointNr",    "Age",	             "Gender",          "WeeksSinceLastVisit", 
                      "WeeksToFollowUp",   "Race",	          "PayedJob",       "LivingSituat",      "NpsEducYears",	  "DailyActivity",   
                      "MostAffSide",       "ParkinMedUser",   "STAITraitSum",	  "STAIStateSum",      "BDI2Sum",         "MoCASum",	
                      "PDQ39_SingleIndex", "Up1_1to6",        "Up1_7to13",      "Up2Total",          "Up3OfTotal",	    "Up3OnTotal",      
                      "Up3OfHoeYah",       "Up3OnHoeYah",     "DiagParkMonth",	"DiagParkYear",	     "AssessWeekNum",   "AssessMonth",     
                      "AssessYear",        "LEDD",            "Subtype_DiagEx3_DisDurSplit",         "QUIPrsSum",       "BodMasInd",
                      "SmokeCurrent",      "SmokeLastYear")]
rm(POM_all)

# To select only the covid participants, load the list of pseudonyms corresponding to covid participants from the file with covid data
Names <- read_csv("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/Covid_longitudinal.csv")["name"]
colnames(Names)[1] <- "Subj_ID"

# remove PIT visits (select only POM visit 1, 2 and 3) - 1414 observations
POM_data <- POM_data[(POM_data$Timepoint == "ses-POMVisit1" |  POM_data$Timepoint == "ses-POMVisit2" | POM_data$Timepoint == "ses-POMVisit3"),]

# Change some variable names and data class
library(magrittr)
library(dplyr)
POM_data <- POM_data %>% mutate(TimepointNr = dplyr::recode(TimepointNr, '0'='1', '1'='2','2'='3'))
POM_data$TimepointNr <- as.factor(POM_data$TimepointNr)
POM_data$Smoking <- ifelse(!is.na(POM_data$SmokeCurrent), POM_data$SmokeCurrent, POM_data$SmokeLastYear)
rename <- c("pseudonym" = "Subj_ID",
            "TimepointNr" = "Visit",
            "Gender" = "Sex",
            "NpsEducYears" = "Years_edu",
            "BDI2Sum" = "BDI",
            "MoCASum" = "MoCA",
            "PDQ39_SingleIndex" = "PDQ",
            "STAITraitSum" = "STAI_t",
            "STAIStateSum" = "STAI_s",
            "Up1_1to6" = "UPDRS1a",
            "Up1_7to13" = "UPDRS1b",
            "Up2Total" = "UPDRS2",
            "Up3OfTotal" = "UPDRS3_off",
            "Up3OnTotal" = "UPDRS3_on",
            "Up3OfHoeYah" = "HY_off",
            "Up3OnHoeYah" = "HY_on",
            "QUIPrsSum" = "QUIP",
            "Subtype_DiagEx3_DisDurSplit" = "Subtype",
            "BodMasInd" = "BMI")

for (old_name in names(POM_data)) {
  if (old_name %in% names(rename)) {
    new_name <- rename[old_name]
    names(POM_data)[names(POM_data) == old_name] <- new_name
  }
}

POM_data <- POM_data %>% mutate(Subtype = dplyr::recode(Subtype, '1_Mild-Motor'='1', '2_Intermediate'='2','3_Diffuse-Malignant'='3', '4_Undefined'='4'))
POM_data$Subtype <- as.factor(POM_data$Subtype)
POM_data$DiagParkMonth = as.double(POM_data$DiagParkMonth)
POM_data <- POM_data %>% mutate(Weeks = coalesce(WeeksSinceLastVisit, WeeksToFollowUp)) # Merge into one variable with weeks after POM baseline

# Add LEDD missings that WEre separately calculated in Matlab script
POM_data$LEDD[POM_data$Subj_ID=="sub-POMU7080B399FC1F46AF" & POM_data$Visit == "2"] <- 600
POM_data$LEDD[POM_data$Subj_ID=="sub-POMU7080B399FC1F46AF" & POM_data$Visit == "3"] <- 600
POM_data$LEDD[POM_data$Subj_ID=="sub-POMUF425DA29DA955CA4" & POM_data$Visit == "3"] <- 172.5981
POM_data$LEDD[POM_data$Subj_ID=="sub-POMUCC9B68541F5DA80C" & POM_data$Visit == "2"] <- 979.188
POM_data$LEDD[POM_data$Subj_ID=="sub-POMU602ABC476136B813" & POM_data$Visit == "2"] <- 575
POM_data$LEDD[POM_data$Subj_ID=="sub-POMU1A50F30F4A977983" & POM_data$Visit == "3"] <- 1080
POM_data$LEDD[POM_data$Subj_ID=="sub-POMUC6F238EF1A578EA5" & POM_data$Visit == "2"] <- 500
POM_data$LEDD[POM_data$Subj_ID=="sub-POMUC6F238EF1A578EA5" & POM_data$Visit == "3"] <- 500
POM_data$LEDD[POM_data$Subj_ID=="sub-POMUA2B1CD3314B0AEA3" & POM_data$Visit == "3"] <- 650.3759
POM_data$LEDD[POM_data$Subj_ID=="sub-POMU5C5E8FB23B76EE52" & POM_data$Visit == "2"] <- 350
POM_data$LEDD[POM_data$Subj_ID=="sub-POMU70B8157191B09234" & POM_data$Visit == "3"] <- 525.188

# Calculate visit 1 date
POM_covid <- merge(POM_data, Names, by="Subj_ID")
POM_covid$AssessYear[POM_covid$Subj_ID=="sub-POMUD4F57C78CB909834" & POM_covid$Visit == "1"] <- 2018   # For one subject assessyear is missing but it is clearly 2018
POM_covid$vis1 <- as.Date(paste(POM_covid$AssessYear, POM_covid$AssessWeekNum, 1, sep="-"), "%Y-%U-%u")

# Restructure to wide, but add age and visit 1 date later otherwise it does not work
library(tidyr)
Covid_short <- POM_covid[c("Subj_ID", "Visit", "Age", "Race", "Weeks", "DiagParkMonth", "DiagParkYear", "vis1")]
Dates_wide <- pivot_wider(Covid_short, id_cols = c('Subj_ID'), names_from = 'Visit', values_from = c('Weeks', 'DiagParkMonth', 'DiagParkYear', 'Race'))
add <- Covid_short[c("Subj_ID", "Visit", "Age", "vis1")]
add <- add[add$Visit == "1",]
Dates_wide <- merge(Dates_wide, add, by="Subj_ID")

# Check if week variable is present for most people (did restructuring go well?)
sum(is.na(Dates_wide$vis1)) # 0
sum(is.na(Dates_wide$Weeks_2)) # 58 ---> now 65, no idea why. I added 8 manually below from AssessTime in UPDRS3 folder
sum(is.na(Dates_wide$Weeks_3)) # 23

# Calculate visit 2 dates
Dates_wide$vis2 <- Dates_wide$vis1 + Dates_wide$Weeks_2*7 # add number of weeks (days*7) since visit 1

# For several participants, AssessWeek and Year are missing, but AssessTime is available with exact date
Dates_wide$vis2[Dates_wide$Subj_ID == "sub-POMU1CA754739B57B126"] <- "2019-01-21" # add manually because this date is found in Assesstime
Dates_wide$vis2[Dates_wide$Subj_ID == "sub-POMU022823FBC6EBD9D9"] <- "2019-09-02" # add manually because this date is found in Assesstime
Dates_wide$vis2[Dates_wide$Subj_ID == "sub-POMU6AB50AF4C627380A"] <- "2019-01-17" # add manually because this date is found in Assesstime
Dates_wide$vis2[Dates_wide$Subj_ID == "sub-POMU7E7448F5C57F585A"] <- "2020-09-22" # add manually because this date is found in Assesstime
Dates_wide$vis2[Dates_wide$Subj_ID == "sub-POMU8A34CEE519A04260"] <- "2019-04-23" # add manually because this date is found in Assesstime
Dates_wide$vis2[Dates_wide$Subj_ID == "sub-POMU9CC9CDA4417A6629"] <- "2020-11-03" # add manually because this date is found in Assesstime
Dates_wide$vis2[Dates_wide$Subj_ID == "sub-POMUB6FD0417ADB272A7"] <- "2020-09-22" # add manually because this date is found in Assesstime
Dates_wide$vis2[Dates_wide$Subj_ID == "sub-POMUEC079158D61E2B02"] <- "2021-01-12" # add manually because this date is found in Assesstime
Dates_wide$vis2[Dates_wide$Subj_ID == "sub-POMU7EE8061949005CB8"] <- "2021-02-01" # add manually because this date is found in Assesstime
Dates_wide$vis2[Dates_wide$Subj_ID == "sub-POMUF83D2E4BD665A11E"] <- "2020-10-20" # add manually because this date is found in Assesstime
Dates_wide$vis2[Dates_wide$Subj_ID == "sub-POMU60DB808C16A39F2A"] <- "2020-10-01" # add manually because this date is found in Assesstime
Dates_wide$vis2[Dates_wide$Subj_ID == "sub-POMU48DEEC15C708B762"] <- "2020-11-23" # add manually because this date is found in Assesstime
Dates_wide$vis2[Dates_wide$Subj_ID == "sub-POMUD8F51CA936C59315"] <- "2020-10-29" # add manually because this date is found in Assesstime
Dates_wide$vis2[Dates_wide$Subj_ID == "sub-POMU20F386A713FBC96C"] <- "2020-09-08" # add manually because this date is found in Assesstime

# Calculate visit 3 date
Dates_wide$vis3 <- if_else(!is.na(Dates_wide$Weeks_2), # if visit 2 took place and is in week variable
                           Dates_wide$vis2 + Dates_wide$Weeks_3*7, # add weeks since baseline at visit 3 to visit 2 date
                           if_else(!is.na(Dates_wide$vis2), # if visit 2 date is known anyways (due to AssessTime variable)
                                    Dates_wide$vis2 + Dates_wide$Weeks_3*7, # add nr. of weeks to visit 2 date
                                    Dates_wide$vis1 + Dates_wide$Weeks_3*7)) # if no visit 2, add weeks since baseline to visit 1 date 
Dates_wide$vis3[Dates_wide$Subj_ID == "sub-POMU1EA14620C7A20443"] <- "2020-08-04" # add manually because this date is found in Assesstime
#Dates_wide$vis3[Dates_wide$Subj_ID == "sub-POMU022823FBC6EBD9D9"] <- Dates_wide$vis2[Dates_wide$Subj_ID == "sub-POMU022823FBC6EBD9D9"] + Dates_wide$Weeks_3[Dates_wide$Subj_ID == "sub-POMU022823FBC6EBD9D9"]*7

# Make new variables for difference (in days) between visit dates and start COVID pandemic (March 11, 2020 is day that pandemic was declared)
Dates_wide$dif1 <- as.Date("2020-03-11") - Dates_wide$vis1
Dates_wide$dif2 <- as.Date("2020-03-11") - Dates_wide$vis2
Dates_wide$dif3 <- as.Date("2020-03-11") - Dates_wide$vis3

# Make a new variable containing the number of the last POM visit before start covid surveys
Dates_wide$precov <- ifelse(Dates_wide$dif3 < 0 | is.na(Dates_wide$dif3), 
                            (ifelse(Dates_wide$dif2 > 0, 2, 1)), 
                            3)
Dates_wide$precov[is.na(Dates_wide$precov)] = 1 
table(Dates_wide$precov)

# Calculate age at start covid pandemic
Dates_wide$Age_cov <- Dates_wide$Age + (Dates_wide$dif1/365)
Dates_wide$Age_cov <- as.double(Dates_wide$Age_cov)

# Calculate disease duration at start covid pandemic (March 11 2020)
Dates_wide$DiagParkMonth <- ifelse(is.na(Dates_wide$DiagParkMonth_1), 1, Dates_wide$DiagParkMonth_1) # For 4 patients month of diagnosis is missing: use January to calculate disease duration
Dates_wide$Diag_date <- as.Date(paste(Dates_wide$DiagParkYear_1, round(Dates_wide$DiagParkMonth*4.34812141), 1, sep = "-"), "%Y-%U-%u") 
Dates_wide$DisDur <- as.double((as.Date("2020-03-11") - Dates_wide$Diag_date)/365)

# Change variable names and select only useful variables for melt
Dates_wide$Visit_1 <- Dates_wide$vis1
Dates_wide$Visit_2 <- Dates_wide$vis2
Dates_wide$Visit_3 <- Dates_wide$vis3
Dates_wide$Race <- Dates_wide$Race_1
Dates_wide <- Dates_wide[,c("Subj_ID", "Age_cov", "DisDur", "Race", "vis1", "vis2", "vis3", "precov")]

# Put back to long format
Dates_long <- pivot_longer(Dates_wide, cols = starts_with("vis"), names_to = "Visit", values_to = "Visit_date")
Dates_long <- Dates_long %>% mutate(Visit=dplyr::recode(Visit, 'vis1'='1', 'vis2'='2','vis3'='3'))
Dates_long$Visit <- as.factor(Dates_long$Visit)
names(Dates_long)[names(Dates_long) == "value"] <- "Date"
POM_covid <- POM_covid[c("Subj_ID", "Sex", "Visit", "PayedJob", "LivingSituat", "Years_edu", "DailyActivity", "MostAffSide", "LEDD", "STAI_t",
                         "STAI_s", "BDI", "MoCA", "PDQ", "UPDRS1a", "UPDRS1b", "UPDRS2", "UPDRS3_off", "UPDRS3_on", "HY_off", 
                         "HY_on", "ParkinMedUser", "Subtype", "QUIP", "BMI", "Smoking")]

# Select precovid POM scores
Complete <- merge(Dates_long, POM_covid, by=c("Subj_ID", "Visit"))
Complete$LEDD <- if_else(Complete$ParkinMedUser=="No", 0, Complete$LEDD) # For non-med users, NA is automatically filled in, but should be 0
Complete$precov <- as.factor(Complete$precov)
Precov_POM <- Complete[(Complete$Visit == Complete$precov),]
Precov_POM$LEDD <- if_else(Precov_POM$ParkinMedUser=="No", 0, Precov_POM$LEDD) # For non-med users, NA is automatically filled in, but should be 0

# Change the "other" level for living situation based on info in text fields: 1=Alone, 2=With_partner, 3=With_family
Precov_POM$LivingSituat[Precov_POM$Subj_ID=="sub-POMU00094252BA30B84F"] <- 1
Precov_POM$LivingSituat[Precov_POM$Subj_ID=="sub-POMU0AEE0E7E9F195659"] <- 3
Precov_POM$LivingSituat[Precov_POM$Subj_ID=="sub-POMU3C5174C7A51309E9"] <- 1
Precov_POM$LivingSituat[Precov_POM$Subj_ID=="sub-POMU60DB808C16A39F2A"] <- 3
Precov_POM$LivingSituat[Precov_POM$Subj_ID=="sub-POMU7842A231D1D8252C"] <- 3
Precov_POM$LivingSituat[Precov_POM$Subj_ID=="sub-POMU7E7448F5C57F585A"] <- 3
Precov_POM$LivingSituat[Precov_POM$Subj_ID=="sub-POMUAAE0C39E47C0BBA9"] <- 2
Precov_POM$LivingSituat[Precov_POM$Subj_ID=="sub-POMUCE5C9ED9AF4FB994"] <- 3
Precov_POM$LivingSituat[Precov_POM$Subj_ID=="sub-POMUD84A16A5685258CF"] <- 1
Precov_POM$LivingSituat[Precov_POM$Subj_ID=="sub-POMUDCE36875FFBC2199"] <- 2
Precov_POM$LivingSituat[Precov_POM$Subj_ID=="sub-POMUEC079158D61E2B02"] <- 3

# Change the "other" level for daily activity based on text fields
Precov_POM$DailyActivity[Precov_POM$Subj_ID=="sub-POMU1CF55A8FB405D10A"] <- 1
Precov_POM$DailyActivity[Precov_POM$Subj_ID=="sub-POMU1E1A426E03FD63B6"] <- 5
Precov_POM$DailyActivity[Precov_POM$Subj_ID=="sub-POMU1EC01A53575D7104"] <- 7
Precov_POM$DailyActivity[Precov_POM$Subj_ID=="sub-POMU242851EF4B22A600"] <- 3
Precov_POM$DailyActivity[Precov_POM$Subj_ID=="sub-POMU26BDF3AC72A04667"] <- 1
Precov_POM$DailyActivity[Precov_POM$Subj_ID=="sub-POMU4676B34286BA2D6E"] <- 1
Precov_POM$DailyActivity[Precov_POM$Subj_ID=="sub-POMU46D404208FF82848"] <- 5
Precov_POM$DailyActivity[Precov_POM$Subj_ID=="sub-POMU4CC057B13DBB2927"] <- 3
Precov_POM$DailyActivity[Precov_POM$Subj_ID=="sub-POMU6AB50AF4C627380A"] <- 3
Precov_POM$DailyActivity[Precov_POM$Subj_ID=="sub-POMU6E6C9C134160CB8D"] <- 5
Precov_POM$DailyActivity[Precov_POM$Subj_ID=="sub-POMU7C93AFD65302F67C"] <- 7
Precov_POM$DailyActivity[Precov_POM$Subj_ID=="sub-POMU8067BDE54D1B1B4A"] <- 7
Precov_POM$DailyActivity[Precov_POM$Subj_ID=="sub-POMU94E6A93D782CE718"] <- 5
Precov_POM$DailyActivity[Precov_POM$Subj_ID=="sub-POMU98ADD50DEBE46D5D"] <- 1
Precov_POM$DailyActivity[Precov_POM$Subj_ID=="sub-POMUB0B32692E1393605"] <- 7
Precov_POM$DailyActivity[Precov_POM$Subj_ID=="sub-POMUC86C6B41DE5A61DB"] <- 6
Precov_POM$DailyActivity[Precov_POM$Subj_ID=="sub-POMUD8F51CA936C59315"] <- 3
Precov_POM$DailyActivity[Precov_POM$Subj_ID=="sub-POMUDCE36875FFBC2199"] <- 7
Precov_POM$DailyActivity[Precov_POM$Subj_ID=="sub-POMUF69ADAA4CB4CE8EF"] <- 5

# Get characteristics table
library(tableone)
Variables    <- c("Age_cov", "Sex", "DisDur", "Race", "PayedJob", "LivingSituat", "Years_edu", "DailyActivity", 
                  "LEDD", "STAI_t", "STAI_s", "BDI", "MoCA", "PDQ", "UPDRS1a", "UPDRS1b", "UPDRS2", "UPDRS3_off", "UPDRS3_on", "Subtype")
Categorical  <- c("Sex","LivingSituat","DailyActivity", "Race", "Subtype") #Define categorical variables

table1 <- CreateTableOne(vars = Variables, data = Precov_POM, factorVars = Categorical)
table1

write.csv(Precov_POM,"M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/POM_precovid_visit.csv", row.names = FALSE)
write.csv(Complete,"M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/POM_allvisits.csv", row.names = FALSE)


#Precov_POM$group <- rep(1,nrow(Precov_POM))
#Precov_whole$group <- rep(2,nrow(Precov_whole)) # For this comparison, first run the whole script for whole sample (don't make covid subset in line 69)
#POMtot <- rbind(Precov_POM, Precov_whole)
#Sel <- POMtot[!duplicated(POMtot$Subj_ID), ]

#Sel$DailyActivity[Sel$DailyActivity == '8'] <- '7'
#Sel$LivingSituat[Sel$LivingSituat== '5'] <- '3'
#Sel$DailyActivity <- factor(Sel$DailyActivity)
#Sel$LivingSituat <- factor(Sel$LivingSituat)

### Test difference between responders and non-responders ###
#stats::t.test(Sel$Age_cov[ which(Sel$group==2)], Precov_POM$Age_cov, var.equal=TRUE)
#stats::t.test(Sel$Years_edu[ which(Sel$group==2)], Precov_POM$Years_edu, var.equal=TRUE)
#stats::t.test(Sel$BDI[ which(Sel$group==2)], Precov_POM$BDI, var.equal=TRUE)
#stats::t.test(Sel$MoCA[ which(Sel$group==2)], Precov_POM$MoCA, var.equal=TRUE)
#stats::t.test(Sel$UPDRS1a[ which(Sel$group==2)], Precov_POM$UPDRS1a, var.equal=TRUE)
#stats::t.test(Sel$UPDRS1b[ which(Sel$group==2)], Precov_POM$UPDRS1b, var.equal=TRUE)
#stats::t.test(Sel$UPDRS2[ which(Sel$group==2)], Precov_POM$UPDRS2, var.equal=TRUE)
#stats::t.test(Sel$UPDRS3_off[ which(Sel$group==2)], Precov_POM$UPDRS3_off, var.equal=TRUE)
#stats::t.test(Sel$UPDRS3_on[ which(Sel$group==2)], Precov_POM$UPDRS3_on, var.equal=TRUE)
#stats::t.test(Sel$LEDD[ which(Sel$group==2)], Precov_POM$LEDD, var.equal=TRUE)
#stats::t.test(Sel$STAI_s[ which(Sel$group==2)], Precov_POM$STAI_s, var.equal=TRUE)
#stats::t.test(Sel$STAI_t[ which(Sel$group==2)], Precov_POM$STAI_t, var.equal=TRUE)
#stats::t.test(Sel$PDQ[ which(Sel$group==2)], Precov_POM$PDQ, var.equal=TRUE)

#chisq.test(table(Sel$group, Sel$Sex))
#chisq.test(table(Sel$group, Sel$DailyActivity))
#chisq.test(table(Sel$group, Sel$LivingSituat))
#chisq.test(table(Sel$group, Sel$HY_off))
#chisq.test(table(Sel$group, Sel$HY_on))
