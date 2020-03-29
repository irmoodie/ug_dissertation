# Thesis Script - Analysis
# IR Moodie
# Last Modified = 29/03/2020

# ---- Start ----

rm(list=ls()) # clear environment

# ---- Packages ----

library(survival) # survival objects
library(lme4) # mixed models
library(MuMIn) # AICc
library(sjPlot) # diagnostics and table outputs
library(tidyverse) # the good stuff

# ---- Data Import ----

exp1 <- read_csv("exp1.csv", na = "na") # Experiment 1 dataset
exp2 <- read_csv("exp2.csv", na = "na") # Experiment 2 dataset
exp3 <- read_csv("exp2_effective_f.csv", na = "na") # Modified Exp 2 dataset so day is now equal to feeding day so a day where no food was give is not counted (gave up trying to do this in R so done manually)

# ---- Data tidying ----

mass1 <- exp1 %>%
  select(-survival, -censored, -label) # make mass dataset

mass1$mass_diff <- mass1$mass_f-mass1$mass_i # get change in mass
mass1$f_total <- mass1$f0+mass1$f1+mass1$f2+mass1$f3+mass1$f4 # get total mass of food consumed

# --- Survival Exp 1 ----


# ---- Survival Exp 2 ----


# ---- Survival Exp 3 ----


# ---- Consumption Exp 1 ----


# ---- Consumption Exp 2 ----


# ---- Consumption Exp 3 ----


# ---- Mass Experiment 1 ----

mass1_lm1 <- lm(f_total~mass_diff*fungus, data = mass1) # plot full model
summary(mass1_lm1) # doesnt seem to be effect of fungus or interaction
plot_model(mass1_lm1, type = "diag") # diag seems ok

mass1_lm2 <- update(mass1_lm1,~. -mass_diff:fungus) # remove interaction
summary(mass1_lm2) # std, error of fungus is huge
plot_model(mass1_lm2, type = "diag") # diag seems better

mass1_lm3 <- update(mass1_lm2,~. -fungus) # remove fungus
summary(mass1_lm3) # all sig with decent effect size of mass_diff - most likely best model
plot_model(mass1_lm3, type = "diag") # all fine

AIC(mass1_lm1,mass1_lm2,mass1_lm3) # confirm with AIC (AIC confirms lm3)

mass_model <- mass1_lm3 # id model for ease
rm(mass1_lm1,mass1_lm2,mass1_lm3) # remove candidate models

tab_model(mass_model,
          p.val = "kr",
          file = "model_summaries/Mass_Experiment_1.html")
