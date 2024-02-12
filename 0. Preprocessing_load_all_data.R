# Load COVID-survey data from PEP download 

# Load required libraries
library(readr)
library(dplyr)
library(jsonlite)

# Get a list of all PPP COVID sub-study participants
base_dir <- "P:/3022026.01/pep/ClinVars4" # Define the base directory
subdirectories <- list.dirs(base_dir, full.names = TRUE, recursive = TRUE)
subject_ids <- character(0) # Initialize an empty vector to store subject IDs

# Loop through the sub-directories and extract subject IDs
for (subdir in subdirectories) {
  # Check if the current sub-directory contains "ses-COVIDbasic"
  if (grepl("ses-COVIDbasic$", subdir)) {
    # Extract the subject ID from the path
    subject_id <- dirname(subdir)
    subject_id <- sub(base_dir, "", subject_id) # Remove the base directory
    subject_ids <- c(subject_ids, subject_id)
  }
}

#### LOAD AND SAVE DATA FOR BASELINE AND FINAL SURVEY ####

# Create vector of JSON file names
json_files <- c('PSSPAS.PSSPAS.AnswerSet0.json',
                'PSSQ.PSSQ.AnswerSet0.json',
                'RRS.RRS.AnswerSet0.json',
                'CorImpact.CorImpact.AnswerSet0.json',
                'UPDRS2.Part1UPDRS1b.AnswerSet0.json',
                'BRSCOPE.COPECERQ.AnswerSet0.json',
                'CLEQ.CLEQ.AnswerSet0.json')

# Define survey types
survey_types <- c("basic", "final")

# Initialize an empty lists to store data
data_basic <- list() 
data_final <- list()

# Loop over each subject ID
for (s in subject_ids) {
    
  # Loop over each survey type
  for (survey in survey_types) {
    
    # Define the directory path for "ses-COVIDbasic" or "ses-COVIDfinal" 
    dirc <- file.path(base_dir, s, paste0("ses-COVID", survey))
    
    # Check if the directory exists
    if (dir.exists(dirc)) {
      
      variables <- list()
      
      # Loop over JSON file names and check if they exist
      for (json_file_name in json_files) {
        
        if (survey == "basic") {
          json_cur <- file.path(dirc, paste0('COVID.Castor.CovidQuestionnaires.CovPackBasic.', json_file_name))
        } else if (survey == "final") {
          json_cur <- file.path(dirc, paste0('COVID.Castor.CovidQuestionnaires.CovPackFinal.', json_file_name))
        }
        
        if (file.exists(json_cur) && json_file_name == 'PSSPAS.PSSPAS.AnswerSet0.json') {
          pss_pas <- jsonlite::fromJSON(json_cur)
          
          # Calculate PSS score
          PSS <- sum(as.numeric(pss_pas$crf$PercStreSca01) - 1,
                     as.numeric(pss_pas$crf$PercStreSca02) - 1,
                     as.numeric(pss_pas$crf$PercStreSca03) - 1,
                     4 - (as.numeric(pss_pas$crf$PercStreSca04) - 1),
                     4 - (as.numeric(pss_pas$crf$PercStreSca05) - 1),
                     as.numeric(pss_pas$crf$PercStreSca06) - 1,
                     4 - (as.numeric(pss_pas$crf$PercStreSca07) - 1),
                     4 - (as.numeric(pss_pas$crf$PercStreSca08) - 1),
                     as.numeric(pss_pas$crf$PercStreSca09) - 1,
                     as.numeric(pss_pas$crf$PercStreSca10) - 1)
          
          # Calculate PAS score
          PAS_names <- paste0("ParAnxSca", sprintf("%02d", 1:4)) # List of all variable names to be counted
          PAS <- sum(sapply(PAS_names, function(var) as.numeric(pss_pas$crf[[var]]) -1 ))} # for each value do -1
        
        # Calculate social support score (F-SoZU)
        if (file.exists(json_cur) && json_file_name == 'PSSQ.PSSQ.AnswerSet0.json') {
          pssq <- jsonlite::fromJSON(json_cur)
          PSSQ_names <- paste0("PerSocSup", sprintf("%02d", 1:8)) # List of all variable names to be counted
          PSSQ <- sum(sapply(PSSQ_names, function(var) as.numeric(pssq$crf[[var]])))}
        
        # Calculate ruminative response score (RRS)
        if (file.exists(json_cur) && json_file_name =='RRS.RRS.AnswerSet0.json') {
          rrs <- jsonlite::fromJSON(json_cur)
          RRS_names <- paste0("RumResSca", sprintf("%02d", 1:5)) # List of all variable names to be counted
          RRS <- sum(sapply(RRS_names, function(var) as.numeric(rrs$crf[[var]])))}
        
        # Calculate sum score for stressor exposure    
        if (file.exists(json_cur) && json_file_name == 'CorImpact.CorImpact.AnswerSet0.json') {
          sl <- jsonlite::fromJSON(json_cur)
          SL_names <- paste0("Impact", sprintf("%02d", 3:20)) # List of all variable names to be counted
          SL <- sum(sapply(SL_names, function(var) as.numeric(sl$crf[[var]]))) - 18 }
        
        # Calculate Behavioural coping and Positive appraisal scores
        if (file.exists(json_cur) && json_file_name == 'BRSCOPE.COPECERQ.AnswerSet0.json') {
          cope <- jsonlite::fromJSON(json_cur) # Read and decode the JSON file
          BHC_names <- paste0("CopInvSho", sprintf("%02d", 1:9))
          BHC_names <- BHC_names[BHC_names != "CopInvSho06"] # Exclude "CopInvSho06"
          BHC <- sum(as.numeric(cope$crf[BHC_names]))
          CERQ_names <- paste0("CogEmoReg", sprintf("%02d", 1:12))
          
          # Calculate the sum based on the variable names
          PASS <- sum(4 * as.numeric(cope$crf[CERQ_names]), 
                      5 * as.numeric(cope$crf$CopInvSho06),
                      5 * as.numeric(cope$crf$CopInvSho10)) /14 }
        
        # Calculate sum score for UPDRS-II  
        if (file.exists(json_cur) && json_file_name == 'UPDRS2.Part1UPDRS1b.AnswerSet0.json') {
          UPDRS2_a <- jsonlite::fromJSON(json_cur) # Read and decode the JSON file part 1
          
          if (survey == "basic") {
            file2 <- file.path(dirc, paste0('COVID.Castor.CovidQuestionnaires.CovPackBasic.UPDRS2.Part2MotAsp1.AnswerSet0.json'))
            file3 <- file.path(dirc, paste0('COVID.Castor.CovidQuestionnaires.CovPackBasic.UPDRS2.Part2MotAsp2.AnswerSet0.json'))
          } else if (survey == "final") {
            file2 <- file.path(dirc, paste0('COVID.Castor.CovidQuestionnaires.CovPackFinal.UPDRS2.Part2MotAsp1.AnswerSet0.json'))
            file3 <- file.path(dirc, paste0('COVID.Castor.CovidQuestionnaires.CovPackFinal.UPDRS2.Part2MotAsp2.AnswerSet0.json'))
          }
          
          UPDRS2_b <- jsonlite::fromJSON(file2)  # Read and decode the JSON file part 2
          UPDRS2_c <- jsonlite::fromJSON(file3)  # Read and decode the JSON file part 3 
          
          UPDRS2 <- sum(as.numeric(UPDRS2_a$crf$Updrs2Micturition), as.numeric(UPDRS2_a$crf$Updrs2Dizziness),
                        as.numeric(UPDRS2_a$crf$Updrs2Insomnia),    as.numeric(UPDRS2_a$crf$Updrs2Sensations),
                        as.numeric(UPDRS2_a$crf$Updrs2Sleepy),      as.numeric(UPDRS2_a$crf$Updrs2Fatigue),
                        as.numeric(UPDRS2_a$crf$Updrs2Obstipation), as.numeric(UPDRS2_b$crf$Updrs2Salivation),
                        as.numeric(UPDRS2_b$crf$Updrs2Speech),      as.numeric(UPDRS2_b$crf$Updrs2Swallow),
                        as.numeric(UPDRS2_b$crf$Updrs2Cutlery),     as.numeric(UPDRS2_b$crf$Updrs2Dress),
                        as.numeric(UPDRS2_b$crf$Updrs2Hygiene),     as.numeric(UPDRS2_b$crf$Updrs2Writing),
                        as.numeric(UPDRS2_c$crf$Updrs2Hobbies),     as.numeric(UPDRS2_c$crf$Updrs2BedTurn),
                        as.numeric(UPDRS2_c$crf$Updrs2Tremor),      as.numeric(UPDRS2_c$crf$Updrs2Stepping),
                        as.numeric(UPDRS2_c$crf$Updrs2Balance),     as.numeric(UPDRS2_c$crf$Updrs2Freezing),
                        na.rm = TRUE)}
        
        # Only for the final survey, there is a questionnaire for life events during survey period
        if (file.exists(json_cur) && json_file_name == 'CLEQ.CLEQ.AnswerSet0.json') {
          le <- jsonlite::fromJSON(json_cur)
          LE_names <- paste0("RadicEvent", sprintf("%02da", 1:27)) # List of all variable names to be counted
          LE <- sum(sapply(LE_names, function(var) as.numeric(le$crf[[var]])))}
      } 
      
      # Put the data for this subject in a data frame
      if (tolower(survey) == "basic") {
        prefix <- "1"
      } else if (tolower(survey) == "final") {
        prefix <- "11"
      }
      
      subject_data <- data.frame(Subj_ID = s)
      column_names <- c(paste0("PSS_",    prefix),
                        paste0("PAS_",    prefix),
                        paste0("UPDRS2_", prefix),
                        paste0("SL_",     prefix), 
                        paste0("SoZu_",   prefix),
                        paste0("RRS_",    prefix),
                        paste0("BHC_",    prefix),
                        paste0("PASS_",   prefix))
      col_values <- c(PSS, PAS, UPDRS2, SL, PSSQ, RRS, BHC, PASS)
      
      for (i in seq_along(column_names)) {
        subject_data[[column_names[i]]] <- col_values[i]
      }
      
      # Append the data frame to the list
      if (tolower(survey) == "basic") {
        data_basic <- append(data_basic, list(subject_data))
      } else if (tolower(survey) == "final") {
        data_final <- append(data_final, list(subject_data))
      }
      
      # Remove the decoded JSON objects (only remove variables that exist)
      vars_to_remove <- c("pss_pas", "pssq", "rrs", "sl", "cope", "UPDRS2_a", "UPDRS2_b", "UPDRS2_c", "le")
      for (var_name in vars_to_remove) {if (exists(var_name)) {rm(list = var_name)}}
      
    }
  } 
}
  
# Store the data
baseline_data <- do.call(rbind, data_basic)
final_data <- do.call(rbind, data_final)
all_data <- merge(baseline_data, final_data, by = "Subj_ID", all.x = TRUE) # Merge, but with all the rows in the baseline data
all_data$Subj_ID <- gsub("^/", "", all_data$Subj_ID)

write.csv(all_data, file = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/Covid_base_final_data.csv", row.names = FALSE)



#### LOAD AND SAVE DATA FOR BIWEEKLY SURVEY ####

# Loop through the sub-directories and extract subject IDs for people with biweekly survey responsses
for (subdir in subdirectories) {
  # Check if the current sub-directory contains "ses-COVIDbasic"
  if (grepl("ses-COVIDweek2$", subdir)) {
    # Extract the subject ID from the path
    subject_id <- dirname(subdir)
    subject_id <- sub(base_dir, "", subject_id) # Remove the base directory in name
    subject_ids_bi <- c(subject_ids, subject_id)
  }
}

# Initialize an empty list to store data
data_biweekly <- list() 

# Loop over each subject ID
for (s in subject_ids_bi) {
  
  # Define the directory path for "ses-COVIDbasic" or "ses-COVIDfinal" 
  dirc <- file.path(base_dir, s, paste0("ses-COVIDweek2"))
  
  if (dir.exists(dirc)) {
    
    variables <- list()
    
      # Construct the full path to the JSON file
      json_file_path <- file.path(dirc, json_file_name)
      
      if (file.exists(json_file_path)) {
        
        # Read the JSON file
        json_data <- jsonlite::fromJSON(json_file_path)
        
        # Check the structure of json_data
        str(json_data)
      }
    # Store the 'variables' list in the 'data_biweekly' list for this subject
 #   data_biweekly[[s]] <- variables
#      if (file.exists(json_file_path)) {
 #       sl <- fromJSON(file = json_file_path)
  #      if (!is.null(sl$crf)) {
   #       SL_values <- sapply(paste0("Impact", sprintf("%02d", 3:20)), function(var))}}
      
      
      
      # Longitudinal scores for stressor load
      impact_files <- list.files(path = dirc, pattern = 'COVID.Castor.CovidQuestionnaires.CovPackWeek2.CorImpact.CorImpact.AnswerSet', full.names = TRUE)
      repetitions_impact <- length(impact_files) / 2
      
      impact_sum <- sapply(1:repetitions_impact, function(j) {
        sl <- fromJSON(impact_files[j])
        SL_values <- sapply(paste0("Impact", sprintf("%02d", 3:20)), function(var) as.numeric(sl$crf[[var]]))
        sum(SL_values) - 18})

  }
}
  

  
  week <- sapply(1:repetitions_impact, function(j) {
    Json <- fromJSON(file.path(survey_type_dir, paste0("COVID.Castor.CovidQuestionnaires.CovPackWeek2.CorImpact.CorImpact.AnswerSet", j - 1, ".WeekNumber.json")))
    as.numeric(Json)
  })
  
  # Longitudinal scores for PSS and PAS
  PSSPAS_files <- list.files(path = survey_type_dir, pattern = 'COVID.Castor.CovidQuestionnaires.CovPackWeek2.PSSPAS.PSSPAS.AnswerSet', full.names = TRUE)
  repetitions_PSSPAS <- length(PSSPAS_files) / 2
  
  PSS <- sapply(1:repetitions_PSSPAS, function(k) {
    Json <- fromJSON(PSSPAS_files[k])
    PSS_values <- sapply(paste0("PercStreSca", sprintf("%02d", 1:10)), function(var) as.numeric(Json$crf[[var]]))
    sum(PSS_values) - 1
  })
  
  PAS <- sapply(1:repetitions_PSSPAS, function(k) {
    Json <- fromJSON(PSSPAS_files[k + repetitions_PSSPAS])
    PAS_values <- sapply(paste0("ParAnxSca", sprintf("%02d", 1:4)), function(var) as.numeric(Json$crf[[var]]))
    sum(PAS_values) - 1
  })
  
  # Longitudinal scores for UPDRS2
  UPDRS2_files <- list.files(path = survey_type_dir, pattern = 'COVID.Castor.CovidQuestionnaires.CovPackWeek2.UPDRS2.Part2MotAsp2.AnswerSet', full.names = TRUE)
  repetitions_UPDRS <- length(UPDRS2_files) / 2
  
  UPDRS2 <- sapply(1:repetitions_UPDRS, function(n) {
    Json1 <- fromJSON(UPDRS2_files[n])
    Json2 <- fromJSON(UPDRS2_files[n + repetitions_UPDRS])
    Json3 <- fromJSON(UPDRS2_files[n + 2 * repetitions_UPDRS])
    
    UPDRS2 <- sum(as.numeric(UPDRS2_a$crf$Updrs2Micturition), as.numeric(UPDRS2_a$crf$Updrs2Dizziness),
                  as.numeric(UPDRS2_a$crf$Updrs2Insomnia),    as.numeric(UPDRS2_a$crf$Updrs2Sensations),
                  as.numeric(UPDRS2_a$crf$Updrs2Sleepy),      as.numeric(UPDRS2_a$crf$Updrs2Fatigue),
                  as.numeric(UPDRS2_a$crf$Updrs2Obstipation), as.numeric(UPDRS2_b$crf$Updrs2Salivation),
                  as.numeric(UPDRS2_b$crf$Updrs2Speech),      as.numeric(UPDRS2_b$crf$Updrs2Swallow),
                  as.numeric(UPDRS2_b$crf$Updrs2Cutlery),     as.numeric(UPDRS2_b$crf$Updrs2Dress),
                  as.numeric(UPDRS2_b$crf$Updrs2Hygiene),     as.numeric(UPDRS2_b$crf$Updrs2Writing),
                  as.numeric(UPDRS2_c$crf$Updrs2Hobbies),     as.numeric(UPDRS2_c$crf$Updrs2BedTurn),
                  as.numeric(UPDRS2_c$crf$Updrs2Tremor),      as.numeric(UPDRS2_c$crf$Updrs2Stepping),
                  as.numeric(UPDRS2_c$crf$Updrs2Balance),     as.numeric(UPDRS2_c$crf$Updrs2Freezing),
                  na.rm = TRUE)
    
    sum(UPDRS2_values)
  })
  
  # Fill out all scores for all timepoints
  Covid_longitudinal[i]$PSS <- PSS
  Covid_longitudinal[i]$PAS <- PAS
  Covid_longitudinal[i]$SL <- impact_sum
  Covid_longitudinal[i]$UPDRS2 <- UPDRS2
  Covid_longitudinal[i]$rep <- repetitions_PSSPAS
  
  PSS_all[i,] <- PSS
  week_nr_PSS[i,] <- week
  
  UPDRS_all[i,] <- UPDRS2
  week_nr_up[i,] <- week
  
  impact_all[i,] <- impact_sum
  week_nr_imp[i,] <- week
  
  PAS_all[i,] <- PAS
  week_nr_PAS[i,] <- week


