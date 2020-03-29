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

consump1 <- exp1 %>%
  gather(key = "day", value = "consumption", f0:f4) %>%
  separate(day,into=c("junk", "day"),-1) %>% 
  mutate(day = as.numeric(day)) %>%
  select(-mass_i, -mass_f, -junk, -survival, -censored) # transform consumption data from wide to long format

consump2 <- gather(exp2, key = "day", value = "consumption", f0:f4) %>%
  separate(day,into=c("junk", "day"),-1) %>% 
  mutate(day = as.numeric(day)) %>%
  select(-junk, -survival, -censored) %>%
  filter(feeding == "constant",
         day != 4) # wide to long for exp 2, day 4 removed as no values in 1 treatment

mass1 <- exp1 %>%
  select(-survival, -censored, -label) # make mass dataset

mass1$mass_diff <- mass1$mass_f-mass1$mass_i # get change in mass
mass1$f_total <- mass1$f0+mass1$f1+mass1$f2+mass1$f3+mass1$f4 # get total mass of food consumed

# --- Survival Exp 1 ----


# ---- Survival Exp 2 ----


# ---- Survival Exp 3 ----


# ---- Consumption Exp 1 ----

hist(consump1$consumption) # skewed
hist(log(consump1$consumption)) # log most promising transformation

consump1_lme1 <- lmer(log(consumption)~factor(day)*scale(spore)*fungus+(1|id), data = consump1) # maximum model
summary(consump1_lme1)
plot_model(consump1_lme1, type = "diag")

consump1_lme2 <- update(consump1_lme1,~. -factor(day):fungus)
summary(consump1_lme2)
plot_model(consump1_lme2, type = "diag")

consump1_lme3 <- update(consump1_lme2,~. -factor(day):scale(spore))
summary(consump1_lme3)
plot_model(consump1_lme3, type = "diag")

consump1_lme4 <- update(consump1_lme3,~. -scale(spore):fungus)
summary(consump1_lme4)
plot_model(consump1_lme4, type = "diag")

consump1_lme5 <- update(consump1_lme4,~. -fungus)
summary(consump1_lme5)
plot_model(consump1_lme5, type = "diag")

consump1_lme6 <- update(consump1_lme5,~. -scale(spore))
summary(consump1_lme6)
plot_model(consump1_lme6, type = "diag")

consump1_lme7 <- update(consump1_lme6,~. -factor(day):scale(spore):fungus)
summary(consump1_lme7)
plot_model(consump1_lme7, type = "diag") # really not fantasitc, but the data convinces me that there is no effect of fungus or spore anayway

AICc(consump1_lme1,consump1_lme2,consump1_lme3,consump1_lme4,consump1_lme5,consump1_lme6,consump1_lme7)

consump1_lme <- consump1_lme7 # save for ease
rm(consump1_lme1,consump1_lme2,consump1_lme3,consump1_lme4,consump1_lme5,consump1_lme6,consump1_lme7) # clean up

tab_model(consump1_lme,
          p.val = "kr",
          file = "model_summaries/Consumption_Exp_1.html") # save model

# ---- Consumption Exp 2 ----

hist(consump2$consumption)
hist(log(consump2$consumption)) # bit better

consump2_lme1 <- lmer(log(consumption)~factor(day)*method*fungus+(1|id), data = consump2)
summary(consump2_lme1)
plot_model(consump2_lme1, type = "diag")

consump2_lme2 <- update(consump2_lme1,~. -factor(day):method:fungus)
summary(consump2_lme2)
plot_model(consump2_lme2, type = "diag")

consump2_lme3 <- update(consump2_lme2,~. -factor(day):fungus)
summary(consump2_lme3)
plot_model(consump2_lme3, type = "diag")

consump2_lme4 <- update(consump2_lme3,~. -factor(day):method)
summary(consump2_lme4)
plot_model(consump2_lme4, type = "diag")

consump2_lme5 <- update(consump2_lme4,~. -fungus:method)
summary(consump2_lme5)
plot_model(consump2_lme5, type = "diag")

consump2_lme6 <- update(consump2_lme5,~. -fungus)
summary(consump2_lme6)
plot_model(consump2_lme6, type = "diag")

consump2_lme7 <- update(consump2_lme6,~. -method)
summary(consump2_lme7)
plot_model(consump2_lme7, type = "diag")

consump2_lme8 <- update(consump2_lme7,~. -factor(day))
summary(consump2_lme8)

AICc(consump2_lme1,consump2_lme2,consump2_lme3,consump2_lme4,consump2_lme5,consump2_lme6,consump2_lme7,consump2_lme8)

consump2_lme <- consump2_lme1
rm(consump2_lme1,consump2_lme2,consump2_lme3,consump2_lme4,consump2_lme5,consump2_lme6,consump2_lme7,consump2_lme8)

tab_model(consump2_lme,
          p.val = "kr",
          file = "model_summaries/Consumption_Exp_2.html")

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
          file = "model_summaries/Mass_Experiment_1.html") # output as html table
