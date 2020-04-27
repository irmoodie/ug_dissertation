# Thesis Script - Analysis
# IR Moodie
# Last Modified = 20/04/2020

# ---- Start ----

rm(list=ls()) # clear environment

# ---- Packages ----

library(survival) # survival objects
library(survMisc) # psuedo r2
library(lme4) # mixed models
library(MuMIn) # AICc
library(sjPlot) # diagnostics and table outputs
library(tidyverse) # the good stuff

# ---- Data Import ----

exp1 <- read_csv("exp1.csv", na = "na") # Experiment 1 dataset
exp2 <- read_csv("exp2.csv", na = "na") # Experiment 2 dataset
exp3 <- read_csv("exp2_effective_f.csv", na = "na") # Modified Exp 2 dataset so day is now equal to feeding day so a day where no food was give is not counted (gave up trying to do this in R so done manually)

# ---- Data tidying ----

survival1 <- exp1 %>%
  select(-mass_i, -mass_f, -f0, -f1, -f2, -f3, -f4) # survival data of interest in exp1

survival1$fungus <- factor(survival1$fungus, levels = c("dead", "alive")) # match order with previous graphs

survival2 <- exp2 %>%
  filter(feeding == "constant") %>%
  select(-f0, -f1, -f2, -f3, -f4) # survival data of interest in exp2

survival2$fungus <- factor(survival2$fungus, levels = c("control", "alive")) # match order with previous graphs

survival3 <- exp2 %>%
  filter(method == "spray") %>%
  select(-f0, -f1, -f2, -f3, -f4) # survival data of interest in exp3

survival3$fungus <- factor(survival3$fungus, levels = c("control", "alive")) # match order with previous graphs

consump1 <- exp1 %>%
  gather(key = "day", value = "consumption", f0:f4) %>%
  separate(day,into=c("junk", "day"),-1) %>% 
  mutate(day = as.numeric(day)) %>%
  select(-mass_i, -mass_f, -junk, -survival, -censored) # transform consumption data from wide to long format

consump1$fungus <- factor(consump1$fungus, levels = c("dead", "alive")) # match order with previous graphs

consump2 <- gather(exp2, key = "day", value = "consumption", f0:f4) %>%
  separate(day,into=c("junk", "day"),-1) %>% 
  mutate(day = as.numeric(day)) %>%
  select(-junk, -censored, -survival) %>%
  filter(feeding == "constant",
         day != 4) # wide to long for exp 2, day 4 removed as no values in 1 treatment

consump2$fungus <- factor(consump2$fungus, levels = c("control", "alive")) # match order with previous graphs

consump3 <- gather(exp3, key = "day", value = "consumption", f0:f4) %>%
  separate(day,into=c("junk", "day"),-1) %>% 
  mutate(day = as.numeric(day)) %>%
  select(-junk, -survival, -censored) %>%
  filter(method == "spray", day != 3, day != 4)

consump3$fungus <- factor(consump3$fungus, levels = c("control", "alive")) # match order with previous graphs

mass1 <- exp1 %>%
  select(-survival, -censored, -label) # make mass dataset

mass1$mass_diff <- mass1$mass_f-mass1$mass_i # get change in mass
mass1$f_total <- mass1$f0+mass1$f1+mass1$f2+mass1$f3+mass1$f4 # get total mass of food consumed

mass1$fungus <- factor(mass1$fungus, levels = c("dead", "alive")) # match order with previous graphs

# --- Survival Exp 1 ----

cox1 <- coxph(Surv(survival, censored) ~ fungus*factor(spore), data = survival1) # build max model
summary(cox1) # nothing sig

cox1.1 <- update(cox1,~. -fungus:factor(spore))
summary(cox1.1)

cox1.2 <- update(cox1.1,~. -factor(spore))
summary(cox1.2)

cox1.3 <- update(cox1.1,~. -fungus)
summary(cox1.3)

cox1.4 <- update(cox1.2,~. -fungus)
summary(cox1.4)

cox1modsum <- model.sel(cox1,cox1.1,cox1.2,cox1.3,cox1.4)
cox1modavg <- model.avg(cox1modsum, subset = delta < 2)
summary(cox1modavg) #full = zero method

exp(0.004531)
exp(-0.415969)
exp(0.406460)

# 95% CI
exp(0.004531-(0.473222*1.96))
exp(0.004531+(0.473222*1.96))

exp(-0.415969-(0.644456*1.96))
exp(-0.415969+(0.644456*1.96))

exp(0.406460-(0.548247*1.96))
exp(0.406460+(0.548247*1.96))

rsq(cox1.3)

# ---- Survival Exp 2 ----

cox2 <- coxph(Surv(survival, censored) ~ fungus*method, data = survival2) # build max model
summary(cox2)

cox2.1 <- update(cox2,~. -fungus:method) 
summary(cox2.1) 

rsq(cox2.1)

cox2.2 <- update(cox2.1,~. -method) 
summary(cox2.2)

cox2.3 <- update(cox2.1,~. -fungus) 
summary(cox2.3) 

cox2modsum <- model.sel(cox2,cox2.1,cox2.2,cox2.3)
cox2modavg <- model.avg(cox2modsum, subset = delta < 2)
summary(cox2modavg) #full = zero method

exp(0.6589)
exp(-2.6599)
exp(-0.4843)

# 95% CI
exp(0.6589+c(-1.96,1.96)*0.4799)
exp(-2.6599+c(-1.96,1.96)*0.7577)
exp(-0.4843+c(-1.96,1.96)*1.0178)

rsq(cox2)
rsq(cox2.1)
rsq(cox2.3)


# ---- Survival Exp 3 ----

cox3 <- coxph(Surv(survival, censored) ~ fungus*feeding, data = survival3) # max model
summary(cox3)

cox3.1 <- update(cox3,~. - fungus:feeding) 
summary(cox3.1)

cox3.2 <- update(cox3.1,~. -fungus) 
summary(cox3.2) =

cox3.3 <- update(cox3.1,~. -feeding)
summary(cox3.3)

cox3.4 <- update(cox3.3,~. -fungus)
summary(cox3.4)

cox3modsum <- model.sel(cox3,cox3.1,cox3.2,cox3.3,cox3.4)
cox3modavg <- model.avg(cox3modsum, subset = delta < 2)
summary(cox3modavg) #full = zero method

rsq(cox3.2)
rsq(cox3.4)

exp(0.29)
exp(0.29)-(1.96*0.5743)
exp(0.29)+(1.96*0.5743)

# ---- Consumption Exp 1 ----

consump1_lme1 <- lmer(consumption~factor(day)*factor(spore)*fungus+(1|id), data = consump1) # maximum model
plot_model(consump1_lme1, type = "diag")
consump1_lme1log <- lmer(log(consumption)~factor(day)*factor(spore)*fungus+(1|id), data = consump1) # maximum model
plot_model(consump1_lme1log, type = "diag")

# log better

c11 <- lmer(log(consumption)~factor(day)*factor(spore)*fungus+(1|id), data = consump1)

c12 <- update(c11,~. -factor(day):factor(spore):fungus)

c13 <- update(c12,~. -factor(day):factor(spore))

c14 <- update(c12,~. -factor(day):fungus)

c15 <- update(c12,~. -factor(spore):fungus)

c16 <- update(c14,~. -factor(day):factor(spore))

c17 <- update(c13,~. -factor(spore):fungus)

c18 <- update(c14,~. -factor(day):factor(spore))

c19 <- update(c18,~. -factor(spore):fungus)

c20 <- update(c19,~. -fungus)

c21 <- update(c19,~. -factor(spore))

c22 <- update(c19,~. -factor(day))

c23 <- update(c20,~. -factor(spore))

c24 <- update(c20,~. -factor(day))

c25 <- update(c21,~. -factor(day))

c26 <- update(c25,~. -fungus)


consump1modsum <- model.sel(c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26)
consump1modavg <- model.avg(consump1modsum, subset = delta < 2)
AICc(c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26)
summary(c23) 
plot_model(c23, type = "diag")

tab_model(c23,
          p.val = "kr",
          file = "model_summaries/Consumption_Exp_1.html") # save model

# ---- Consumption Exp 2 ----

consump2_lme1 <- lmer(consumption~factor(day)*method*fungus+(1|id), data = consump2)
consump2_lme1log <- lmer(log(consumption)~factor(day)*method*fungus+(1|id), data = consump2)
plot_model(consump2_lme1, type = "diag")
plot_model(consump2_lme1log, type = "diag")

# log better

d11 <- lmer(log(consumption)~factor(day)*method*fungus+(1|id), data = consump2)

d12 <- update(d11,~. -factor(day):method:fungus)

d13 <- update(d12,~. -factor(day):method)

d14 <- update(d12,~. -factor(day):fungus)

d15 <- update(d12,~. -method:fungus)

d16 <- update(d14,~. -factor(day):method)

d17 <- update(d13,~. -method:fungus)

d18 <- update(d14,~. -factor(day):method)

d19 <- update(d18,~. -method:fungus)

d20 <- update(d19,~. -fungus)

d21 <- update(d19,~. -method)

d22 <- update(d19,~. -factor(day))

d23 <- update(d20,~. -method)

d24 <- update(d20,~. -factor(day))

d25 <- update(d21,~. -factor(day))

d26 <- update(d25,~. -fungus)

consump2modsum <- model.sel(d11,d12,d13,d14,d15,d16,d17,d18,d19,d20,d21,d22,d23,d24,d25,d26)
consump2modavg <- model.avg(consump2modsum, subset = delta < 2)
AICc(d11,d12,d13,d14,d15,d16,d17,d18,d19,d20,d21,d22,d23,d24,d25,d26)
summary(d11)

tab_model(d11,
          p.val = "kr",
          file = "model_summaries/Consumption_Exp_2.html") # save model

# ---- Consumption Exp 3 ----

consump3_lme1 <- lmer(consumption~factor(day)*feeding*fungus+(1|id), data = consump3) # max model
consump3_lme1log <- lmer(log(consumption)~factor(day)*feeding*fungus+(1|id), data = consump3) # max model
plot_model(consump3_lme1, type = "diag") 
plot_model(consump3_lme1log, type = "diag")

# log worse

e11 <- lmer(consumption~factor(day)*feeding*fungus+(1|id), data = consump3)

e12 <- update(e11,~. -factor(day):feeding:fungus)

e13 <- update(e12,~. -factor(day):feeding)

e14 <- update(e12,~. -factor(day):fungus)

e15 <- update(e12,~. -feeding:fungus)

e16 <- update(e14,~. -factor(day):feeding)

e17 <- update(e13,~. -feeding:fungus)

e18 <- update(e14,~. -factor(day):feeding)

e19 <- update(e18,~. -feeding:fungus)

e20 <- update(e19,~. -fungus)

e21 <- update(e19,~. -feeding)

e22 <- update(e19,~. -factor(day))

e23 <- update(e20,~. -feeding)

e24 <- update(e20,~. -factor(day))

e25 <- update(e21,~. -factor(day))

e26 <- update(e25,~. -fungus)

consump3modsum <- model.sel(e11,e12,e13,e14,e15,e16,e17,e18,e19,e20,e21,e22,e23,e24,e25,e26)
consump3modavg <- model.avg(consump3modsum, subset = delta < 2)
AICc(e11,e12,e13,e14,e15,e16,e17,e18,e19,e20,e21,e22,e23,e24,e25,e26)
summary(e11) 

tab_model(e11,
          p.val = "kr",
          file = "model_summaries/Consumption_Exp_3.html")


# ---- Mass Experiment 1 ----

mass1_lm1 <- lm(mass_diff~f_total*fungus*factor(spore), data = mass1) # plot full model
summary(mass1_lm1) 
plot_model(mass1_lm1, type = "diag") 

mass1_lm2 <- update(mass1_lm1,~. -f_total:fungus:factor(spore)) # remove interaction
summary(mass1_lm2) # std, error of fungus is huge
plot_model(mass1_lm2, type = "diag") # diag seems better

mass1_lm3 <- update(mass1_lm2,~. -fungus:factor(spore)) # remove fungus
summary(mass1_lm3) # all sig with decent effect size of mass_diff - most likely best model
plot_model(mass1_lm3, type = "diag") # all fine

mass1_lm4 <- update(mass1_lm3,~. -f_total:factor(spore))
summary(mass1_lm4)

mass1_lm5 <- update(mass1_lm4,~. -f_total:fungus)
summary(mass1_lm5)

mass1_lm6 <- update(mass1_lm5,~. -factor(spore))
summary(mass1_lm6)

mass1_lm7 <- update(mass1_lm6,~. -fungus)
summary(mass1_lm7)

AICc(mass1_lm1,mass1_lm2,mass1_lm3,mass1_lm4,mass1_lm5,mass1_lm6,mass1_lm7) # confirm with AIC

massmodsum <- model.sel(mass1_lm1,mass1_lm2,mass1_lm3,mass1_lm4,mass1_lm5,mass1_lm6,mass1_lm7)
massmodavg <- model.avg(massmodsum, subset = delta <2, fit = TRUE)
summary(massmodavg)

mass1$pred <- 1

mass1 <- mass1 %>%
  drop_na()
  
masspred1 <- predict(massmodavg, newdata = mass1, se = TRUE)
mass1$pred <- masspred1$fit
mass1$se <- masspred1$se.fit

mass_model <- mass1_lm7 # id model for ease
rm(mass1_lm1,mass1_lm2,mass1_lm3,mass1_lm4,mass1_lm5,mass1_lm6,mass1_lm7) # remove candidate models

# ---- Agar Experiment ----

agar <- prop.test(x = c(0, 6), n = c(9, 9))
prop.test(x = c(0, 0), n = c(6, 9))

# ---- End ----