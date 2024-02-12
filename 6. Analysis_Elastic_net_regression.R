# This script performs elastic-net regression analysis on 17 variables with SR as outcome measure

# The elastic net penalty is controlled by α, and bridges the gap between lasso regression (α=1, the default) and ridge regression (α=0).
# We want a combination of the two methods, since ridge regression deals best with multicollinearity and lasso with variable selection; hence
# we used α=0.5 to use both methods equally. Tuning parameter λ controls strength of the penalty. The ridge penalty shrinks the coefficients 
# of correlated predictors towards each other while lasso tends to pick one of them and discard the others.


rm(list = ls())

###################     Prepare data and settings    ###################

# Load data set and select the variables of interest
library(readr)
POM_data  <- read_csv("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/POM_precovid_visit.csv",
                     col_types = cols(Sex            = col_factor(levels = c("Male", "Female")), 
                                      LivingSituat   = col_factor(levels = c("1", "2", "3", "4", "5")), 
                                      DailyActivity  = col_factor(levels = c("1", "2", "3", "4", "5", "6", "7", "8"))))
Covid_averages <- read_csv("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/Covid_averages.csv")
cov = subset(Covid_averages, select = -c(Age_cov,Sex))
Data_wide <- merge(POM_data, cov, by="Subj_ID")

names(Data_wide)[names(Data_wide) == "Age_cov"]        <- "Age"
names(Data_wide)[names(Data_wide) == "Comorb_present"] <- "Comorbidity"
rm(cov, POM_data, Covid_averages)

levels(Data_wide$LivingSituat) <- c("Alone", "With_partner", "With_family", "Nursinghome", "Other")
levels(Data_wide$DailyActivity) <- c("Paid", "Household", "Retired", "Student", "Not_paid", "Invol_nojob", "Volunteer", "Other", "Not_working")
# original levels daily activity: 1=paid job, 2=household, 3=retired, 4=student, 5=no paid job due to illness, 6=involuntarily no job, 
# 7=volunteer work, 8=other; renamed the levels to remove the levels that contained no participants)

# In elastic net, categorical predictors with many levels or levels with low freqs are not ideal. Therefore, combine some levels (Student and other are not present)
Data_wide$DailyActivity[Data_wide$DailyActivity == 'Invol_nojob'] <- 'Not_working'
Data_wide$DailyActivity[Data_wide$DailyActivity == 'Household'] <- 'Not_working'
Data_wide$DailyActivity[Data_wide$DailyActivity == 'Volunteer'] <- 'Not_paid'

Data_wide$Daily_Act <- factor(Data_wide$DailyActivity)
Data_wide$Liv_Sit <- factor(Data_wide$LivingSituat)
Data_wide$SR_class <- factor(Data_wide$SR_class)
levels(Data_wide$SR_class) <- c("highSR","lowSR") # rename factor levels

# Prepare the data (one selection with continuous resilience measure, one with the latent classes)
library(dplyr)
Data_wide$SR               <- -Data_wide$Res
Data_wide$STAI             <- (Data_wide$STAI_s + Data_wide$STAI_t)/2
Data_wide$UPDRSnm          <- (Data_wide$UPDRS1a + Data_wide$UPDRS1b + Data_wide$UPDRS2)
Data_wide$Quality_of_life  <- -Data_wide$PDQ
Data_wide$Years_education  <- Data_wide$Years_edu
Data_wide$Disease_duration <- Data_wide$DisDur
Data_wide$Comorbidity      <- as.factor(Data_wide$Comorbidity)
levels(Data_wide$Comorbidity) <- c("No", "Yes")

Data_con <- na.omit(subset(Data_wide, select = c(SR, Age, Sex, Disease_duration, Years_education, LEDD, MoCA, Quality_of_life, STAI, BDI, UPDRS3_off, # only select complete cases
                                                 UPDRSnm, SoZU, RRS, PASS, Liv_Sit, Daily_Act, Comorbidity))) 
Data_bin <- na.omit(subset(Data_wide, select = c(SR_class, Age, Sex, Disease_duration, Years_education, LEDD, MoCA, Quality_of_life, STAI, BDI, UPDRS3_off, # only select complete cases
                                                 UPDRSnm, SoZU, RRS, PASS, Liv_Sit, Daily_Act, Comorbidity))) 
rm(Data_wide)

# Load libraries
library(tidymodels)
library(glmnet)
library(future)
library(furrr)
library(viridis)

# Settings for analysis
set.seed(12)
gsLambda <- 101         # Grid size for lambda (glmnet default = 100)
gsAlpha <- 11           # Grid size for alpha (0.1 increments)
nFolds <- 10            # K for K-fold cross validations (standard = 10)
nFoldrepeats <- 100     # Number of repeats for cross-validation (standard = 100)
nPermutations <- 1000   # Number of permutation test (standard = 1000)

###################     Analysis continuous SR measure       ################### 

# Pre-processing pipeline for continuous SR measure
preprocessing_recipe_con <- recipe(Data_con) %>%
        update_role(everything(), new_role = "predictor") %>%   # Set predictors
        update_role(SR, new_role = "outcome") %>%               # Set outcome
        step_naomit(everything(), skip = T) %>%                 # Every row with missing data should be removed
        step_dummy(Sex) %>%                                     # Factors need to become a dummy variable
        step_dummy(Comorbidity) %>%                             # Factors need to become a dummy variable
        step_dummy(Liv_Sit) %>%                                 # Factors need to become a dummy variable
        step_dummy(Daily_Act) %>%                               # Factors need to become a dummy variable
        step_zv(recipes::all_predictors()) %>%                  # Zero Variance Filter: Possible disbalance in some predictors might give a CV fold with eg. only men, that should be removed
        step_YeoJohnson(recipes::all_numeric_predictors()) %>%  # Makes the predictors more Gaussian
        step_normalize(recipes::all_numeric_predictors())       # Z-score everything

preprocessing_recipe_prep_con <- prep(preprocessing_recipe_con, training = Data_con)
DataSet_prep_con              <- bake(preprocessing_recipe_prep_con, Data_con) # Apply pre-processing steps to data

# Run analysis looking for both optimal lambda and alpha
hyperparam_grid <- grid_regular(penalty(range=c(0,3), trans=NULL), mixture(), levels = c(gsLambda, gsAlpha), original = T) # Set grid for grid search alpha and lambda
Data_cv <- rsample::vfold_cv(data = Data_con, v = nFolds, repeats = nFoldrepeats)                                          # Set cross validation for data
model_engine <- linear_reg() %>%  set_args(penalty = tune(), mixture = tune()) %>% set_engine("glmnet")                    # Set model type
modelling_workflow <- workflow() %>% add_recipe(preprocessing_recipe_con) %>% add_model(model_engine)                      # Set workflow 
model_fit_con <- tune::tune_grid(modelling_workflow, resamples = Data_cv, grid = hyperparam_grid, metrics = metric_set(rmse, rsq_trad, rsq)) # Run model fit

# Run analysis looking for optimal lambda with alpha fixed at 0.5 (this is what was used for manuscript)
lambda_grid <- grid_regular(penalty(range=c(0,3), trans=NULL), levels = gsLambda, original = T)                            # Set grid for grid search lambda
Data_cv <- rsample::vfold_cv(data = Data_con, v = nFolds, repeats = nFoldrepeats)                                          # Set cross validation for data
model_engine <- linear_reg() %>%  set_args(penalty = tune(), mixture = 0.5) %>% set_engine("glmnet")                       # Set model type
modelling_workflow <- workflow() %>% add_recipe(preprocessing_recipe_con) %>% add_model(model_engine)                      # Set workflow 
model_fit_con <- tune::tune_grid(modelling_workflow, resamples = Data_cv, grid = lambda_grid, metrics = metric_set(rmse, rsq_trad, rsq)) # Run model fit

# Save model fit
save(model_fit_con, file = file.path("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/", "model_fit_con_new2.RData"))

# Permutation testing settings
best_hyperparams <- tune::select_best(model_fit_con, metric = "rsq_trad")  # Get best hyperparameters & update workflow
final_model      <- finalize_model(model_engine, best_hyperparams)
final_workflow   <- modelling_workflow %>% update_model(final_model)

Data_perm_con    <- rsample::nested_cv(data = Data_con, # Make a set of permuted data and for each set a cross validation (CV)
                                       outside = permutations(permute = SR, times = nPermutations),
                                       inside = vfold_cv(v = nFolds, repeats = 1)) # To make sure computing does not explode, use 1 repeat for permutation

fit_perm_resamples <- function(inner_resamples, wflow) {tune::fit_resamples(wflow, resamples = inner_resamples, 
                                                                            metrics = yardstick::metric_set(rsq_trad, rmse))} # Apply the workflow per mutation per CV

collect_mean_performance <- function(performance_metrics) {performance_metrics %>% # Output the correct metrics per mutation (and mean over the CV)
    dplyr::select(.metric, mean) %>%
    tidyr::pivot_wider(names_from = ".metric", values_from = "mean")}

# Run permutation tests
model_fit_con_perm <- Data_perm_con %>%
  dplyr::mutate(
    fits = furrr::future_map( # furrr::future_map is like purrr::map but with parallel processing
      .x = inner_resamples,   # .x is the list to which the function should be applied for each element
      .f = ~ fit_perm_resamples(inner_resamples = .x, wflow = final_workflow), # .f is the function which will be applied to each element of .x
      .options = furrr::furrr_options(stdout = T, seed = NULL, packages = c("tidyverse", "tidymodels", "glmnet", "vctrs", "furrr", "rlang")),
      .progress = F)) %>%
  dplyr::mutate(performance = purrr::map(fits, tune::collect_metrics)) %>% 
  dplyr::mutate(mean_performance = purrr::map(performance, collect_mean_performance)) %>% # Summarize the performance
  tidyr::unnest(mean_performance)

# Save model fit
save(model_fit_con_perm, file = file.path("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/", "model_fit_con_perm_new2.RData"))

###################       Analysis binary SR-class measure       ################### 

# Pre-processing pipeline for binary SR-class measure
preprocessing_recipe_bin <- recipe(Data_bin) %>%
        update_role(everything(), new_role = "predictor") %>%   # Set predictors
        update_role(SR_class, new_role = "outcome") %>%         # Set outcome
        step_naomit(everything(), skip = T) %>%                 # Every row with missing data should be removed
        step_dummy(Sex) %>%                                     # Factors need to become a dummy variable
        step_dummy(Comorbidity) %>%                             # Factors need to become a dummy variable
        step_dummy(Liv_Sit) %>%                                 # Factors need to become a dummy variable
        step_dummy(Daily_Act) %>%                               # Factors need to become a dummy variable
        step_zv(recipes::all_predictors()) %>%                  # Zero Variance Filter: Possible disbalance in some predictors might give a CV fold with eg. only men, that should be removed
        step_YeoJohnson(recipes::all_numeric_predictors()) %>%  # Makes the predictors more Gaussian
        step_normalize(recipes::all_numeric_predictors())       # Z-score everything

preprocessing_recipe_prep_bin <- prep(preprocessing_recipe_bin, training = Data_bin)
DataSet_prep_bin <- bake(preprocessing_recipe_prep_bin, Data_bin)

# Run analysis looking for both optimal lambda and alpha
hyperparam_grid <- grid_regular(penalty(range=c(0,2), trans=NULL), mixture(), levels = c(gsLambda, gsAlpha), original = T) # Set grid for gidsearch alpha and lambda
Data_cv <- rsample::vfold_cv(data = Data_bin, v = nFolds, repeats = nFoldrepeats)                           # Set cross validation for data
model_engine <- logistic_reg() %>% set_args(penalty = tune(), mixture = tune()) %>% set_engine("glmnet")    # Set model type
modelling_workflow <- workflow() %>% add_recipe(preprocessing_recipe_bin) %>% add_model(model_engine)       # Set workflow 
model_fit_bin <- tune::tune_grid(modelling_workflow, resamples = Data_cv, grid = hyperparam_grid, metrics = metric_set(roc_auc)) # Run new model fit

# Run analysis looking for optimal lambda with alpha fixed at 0.5
lambda_grid <- grid_regular(penalty(range=c(0,2), trans=NULL), levels = gsLambda, original = T)                            # Set grid for grid search lambda
Data_cv <- rsample::vfold_cv(data = Data_bin, v = nFolds, repeats = nFoldrepeats)                                          # Set cross validation for data
model_engine <- logistic_reg() %>% set_args(penalty = tune(), mixture = 0.5) %>% set_engine("glmnet")                      # Set model type
modelling_workflow <- workflow() %>% add_recipe(preprocessing_recipe_bin) %>% add_model(model_engine)                      # Set workflow
model_fit_bin <- tune::tune_grid(modelling_workflow, resamples = Data_cv, grid = lambda_grid, metrics = metric_set(roc_auc)) # Run model fit

# Save model fit
save(model_fit_bin, file = file.path("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/", "model_fit_bin_new.RData"))

# Permutation testing settings
best_hyperparams <- tune::select_best(model_fit_bin, metric = "roc_auc")  #Get best hyperparameters & update workflow
final_model <- finalize_model(model_engine, best_hyperparams)
final_workflow <- modelling_workflow %>% update_model(final_model)

Data_perm_bin <- rsample::nested_cv(data = Data_bin, #Make a set of permutated data + for each set a CV
                                outside = permutations(permute = SR_class, times = nPermutations),
                                inside = vfold_cv(v = nFolds, repeats = 1)) #To make sure the computing does not explode, we always use 1 repeat for permutation

fit_perm_resamples <- function(inner_resamples, wflow) {tune::fit_resamples(wflow, resamples = inner_resamples, 
                                                                            metrics = yardstick::metric_set(roc_auc))} # Apply the workflow per mutaiton per CV

collect_mean_performance <- function(performance_metrics) {performance_metrics %>% # Output the correct metrics per mutation (and mean over the CV)
                                                           dplyr::select(.metric, mean) %>%
                                                           tidyr::pivot_wider(names_from = ".metric", values_from = "mean")}

# Run permutation tests
model_fit_bin_perm <- Data_perm_bin %>%
  dplyr::mutate(
    fits = furrr::future_map( # furrr::future_map is like purrr::map but with parallel processing
      .x = inner_resamples, # .x is the list to which the function should be applied for each element
      .f = ~ fit_perm_resamples(inner_resamples = .x, wflow = final_workflow), # .f is the function which will be applied to each element of .x
      .options = furrr::furrr_options(stdout = T, seed = NULL, packages = c("tidyverse", "tidymodels", "glmnet", "vctrs", "furrr", "rlang")),
      .progress = F)) %>%
  dplyr::mutate(performance = purrr::map(fits, tune::collect_metrics)) %>% 
  dplyr::mutate(mean_performance = purrr::map(performance, collect_mean_performance)) %>% # Summarise the performance
  tidyr::unnest(mean_performance)

# Save model fit
save(model_fit_bin_perm, file = file.path("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/", "model_fit_bin_perm_new.RData"))

###################     Looking at the results         ################### 

# Load results from previous steps, running them again takes a long time
load("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/model_fit_con.RData")
load("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/model_fit_bin.RData")
load("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/model_fit_con_perm.RData")
load("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/model_fit_bin_perm.RData")

## Continuous model ##

model_fit_result_con <- estimate_tune_results(model_fit_con)
best_hyperparams <- tune::select_best(model_fit_con, metric = "rsq_trad") # penalty 0.54
show_best(model_fit_con, metric = "rmse")
best_hyperparams

# Plot heatmap of all parameter solutions (only when alpha is also determined with CV)
model_fit_result_con %>% dplyr::filter(.metric == 'rmse') %>%
        ggplot(aes(x = as.factor(mixture), y = as.factor(penalty))) +
        geom_tile(aes(fill = mean)) +
        scale_fill_viridis(option = 'A') +
        labs(x = 'elastic net penalty', y = 'penalty strength')

model_fit_result_con %>% dplyr::filter(.metric == 'rsq_trad') %>%
        ggplot(aes(x = as.factor(mixture), y = as.factor(penalty))) +
        geom_tile(aes(fill = mean)) +
        scale_fill_viridis(option = 'A') +
        labs(x = 'elastic net penalty', y = 'penalty strength')

# Run model 
PredVars <- names(DataSet_prep_con) %in% c("SR")
bestFIT <- glmnet(as.matrix(DataSet_prep_con[,!PredVars]), # x = all predictors
                  as.matrix(DataSet_prep_con[,PredVars]),  # y = response vector
                  type.measure = "rsq_trad",  
                  alpha = 0.5, # with preset mixture (=alpha), which was optimal after tuning step in binary analysis
                  #alpha = best_hyperparams$mixture, # this takes the mixture as was determined optimal by tuning
                  family = "gaussian")

CI <- boot.glmnet(as.matrix(DataSet_prep_con[,!PredVars]), # x = all predictors
                  as.matrix(DataSet_prep_con[,PredVars]))
covered <- DataSet_prep_con$beta >= CI[,1] & DataSet_prep_con$beta <= CI[,2]
table(covered)

plot(x=bestFIT)                              # Plot the coefficients (shows path of coefficients against ℓ1-norm of the whole coefficient vector as λ varies)
print(bestFIT)                               # summary of glmnet path at each step (df=nr of nonzero coefficients, %Dev= % (of null) deviance explained)
plot(bestFIT, xvar = "dev", label = TRUE)    # Plot variance explained (measure of complexity of the model)
plot(bestFIT, xvar = "lambda", label = TRUE) # Plot coefficients over lambda

# Show coefficients
cCOEF <- as.data.frame(as.matrix(coef(bestFIT, s = best_hyperparams$penalty))) # s is lambda of the best fit (the "optimal" penalty)
cCOEF_named <- cbind(Predictor=rownames(cCOEF), b=cCOEF)
cCOEF_named <- cCOEF_named[-1,]
cCOEF_named$Predictors <- ifelse(cCOEF_named$s1 < 0, "Resilience factors", "Risk factors")

p <- ggplot(data = cCOEF_named, mapping = aes(x = reorder(Predictor,abs(s1)), y = s1)) + 
            geom_bar(stat = "identity", aes(fill = Predictors)) +
            scale_fill_manual(values=c("Resilience factors"="#009988", "Risk factors"="#CC3311")) +         
            coord_flip() + 
            ggtitle("Coefficient contribution to optimal fit") +
            ylab("Coefficient") + xlab("Variable Name") + 
            scale_y_continuous(breaks = c(-1.0, -0.5, 0, 0.5, 1.0, 1.5)) +
            theme_minimal()
Coefs_plot <- p + theme(legend.text = element_text(size = 12, colour = "black"), 
                        legend.title = element_text(size = 14),
                        axis.text = element_text(size = 12, colour = "black"), 
                        axis.title = element_text(size = 14, margin = margin(r=10)))

# Get R-squared of the best model
true_performance <- model_fit_result_con %>% dplyr::filter(.metric == 'rsq_trad') %>% dplyr::pull(mean) %>% max(na.rm = TRUE)

# Plot null distribution
{hist(model_fit_con_perm$rsq_trad, n = 100, xlim = c(-0.2, 0.4)) 
  abline(v=true_performance, col="red")}  # add the predictive performance that was achieved with the real data

# Calculate p-value
Pvalue <- mean(model_fit_con_perm$rsq_trad >= true_performance)
print(Pvalue)

# Save plot with coefficients of elastic net regression (figure 4 in manuscript)
ggsave("Coefs_continuousSR.svg", plot = Coefs_plot, path = "M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Figures and models" )


## Binary SR model ##

# Output results 
model_fit_result_bin <- estimate_tune_results(model_fit_bin)
best_hyperparams <- tune::select_best(model_fit_bin, metric = "roc_auc")
show_best(model_fit_bin, metric = "roc_auc")
best_hyperparams

# Plot heatmap of all parameter solutions
model_fit_result_bin %>% dplyr::filter(.metric == 'roc_auc') %>%
        ggplot(aes(x = as.factor(mixture), y = as.factor(penalty))) +
        geom_tile(aes(fill = mean)) +
        scale_fill_viridis(option = 'A') +
        labs(x = 'elastic net penalty', y = 'penalty strength')

# Run model 
PredVars <- names(DataSet_prep_bin) %in% c("SR_class")
bestFIT <- glmnet(as.matrix(DataSet_prep_bin[,!PredVars]), # x = all predictors
                  as.matrix(DataSet_prep_bin[,PredVars]),  # y = response vector
                  type.measure="roc_auc",  
                  alpha = 0.3,
                  #alpha=best_hyperparams$mixture, 
                  family="binomial")

plot(x=bestFIT)                              # Plot the coefficients (shows path of coefficients against ℓ1-norm of the whole coefficient vector as λ varies)
print(bestFIT)                               # summary of glmnet path at each step (df=nr of nonzero coefficients, %Dev= % (of null) deviance explained)
plot(bestFIT, xvar = "dev", label = TRUE)    # Plot variance explained (measure of complexity of the model)
plot(bestFIT, xvar = "lambda", label = TRUE) # Plot coefficient over lambda

# Show coefficients
cCOEF <- as.data.frame(as.matrix(coef(bestFIT, s = best_hyperparams$penalty))) # s is lambda of the best fit (the "optimal" penalty)
cCOEF_named <- cbind(Predictor=rownames(cCOEF), b=cCOEF)
cCOEF_named <- cCOEF_named[-1,]
cCOEF_named$Predictors <- ifelse(cCOEF_named$s1 < 0, "Resilience factors", "Risk factors")

p2 <- ggplot(data = cCOEF_named, mapping = aes(x = reorder(Predictor,abs(s1)), y = s1)) + 
             geom_bar(stat = "identity", aes(fill = Predictors)) +
             scale_fill_manual(values=c("Resilience factors"="#009988", "Risk factors"="#CC3311")) +         
             coord_flip() + 
             ggtitle("Coefficient contribution to optimal fit") +
             ylab("Coefficient") + xlab("Variable Name") + 
             scale_y_continuous(breaks = c(-1.0, -0.5, 0, 0.5, 1.0, 1.5)) +
             theme_minimal()
Coefs_plot_bin <- p2 + theme(legend.text = element_text(size = 12, colour = "black"), 
                            legend.title = element_text(size = 14),
                            axis.text = element_text(size = 12, colour = "black"), 
                            axis.title = element_text(size = 14, margin = margin(r=10)))

# Get R-squared of the best model
true_performance <- model_fit_result_bin %>% dplyr::filter(.metric == 'roc_auc') %>% dplyr::pull(mean) %>% max(na.rm = TRUE)

# Plot null distribution and add the predictive performance of the data
{hist(model_fit_bin_perm$roc_auc, n = 100, xlim = c(0, 1)) 
  abline(v=true_performance, col="red")}  # add the predictive performance that was achieved with the real data

# Calculate p-value
Pvalue <- mean(model_fit_bin_perm$roc_auc >= true_performance)
print(Pvalue)


###################     Mediation analysis      ################### 

# Mediation analysis to determine relationship between SR, perceived social support and positive appraisal style
 
hist(Data_con$SR)
hist(Data_con$PASS)
hist(Data_con$SoZU)
hist(Data_con$Years_edu)

# Relation PSS -> PASS -> SR
summary(med_rep1_a <- lm(scale(PASS)~scale(SoZU) + Age + Sex + Years_edu,
                         data=Data_con))

summary(med_rep1_b<- lm(SR ~ scale(PASS) + Age + Sex + Years_edu,
                        data=Data_con))

indirect_effect_socsupp <- medci(mu.y = summary(med_rep1_b)$coef[2,1], se.y=summary(med_rep1_b)$coef[2,2],  mu.x = summary(med_rep1_a)$coef[2,1], se.x = summary(med_rep1_a)$coef[2,2], alpha = 0.05, type="dop")

# Calculate sign. level at alpha = .001
medci(mu.y = summary(med_rep1_b)$coef[2,1], se.y=summary(med_rep1_b)$coef[2,2],  mu.x = summary(med_rep1_a)$coef[2,1], se.x = summary(med_rep1_a)$coef[2,2], alpha = 0.001, type="dop")

indirect_effect_socsupp # interval excludes 0, i.e. sign. mediation effect

