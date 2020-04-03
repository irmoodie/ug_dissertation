# Thesis Script - Analysis
# IR Moodie
# Last Modified = 30/03/2020

# ---- Start ----

rm(list=ls()) # clear environment

# ---- Packages ----

library(survival) # survival objects
library(lme4) # mixed models
library(MuMIn) # AICc
library(sjPlot) # diagnostics and table outputs
library(glmmTMB)
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

survival1result <- AICc(cox1,cox1.1,cox1.2,cox1.3,cox1.4)
survival1result <- tibble::rownames_to_column(survival1result, "Model") # add model names
survival1result$Model <- fct_recode(survival1result$Model,
                          "fungus*spore" = "cox1",
                          "fungus+spore" = "cox1.1",
                          "fungus" = "cox1.2",
                          "spore" = "cox1.3",
                          "none" = "cox1.4") # add helpful names
tab_df(survival1result,
       file = "model_summaries/Suvival_Experiment_1.html")

rm(cox1,cox1.1,cox1.2,cox1.3,cox1.4)


# ---- Survival Exp 2 ----

cox2 <- coxph(Surv(survival, censored) ~ fungus*method, data = survival2) # build max model
summary(cox2) # interaction not significant

cox2.1 <- update(cox2,~. -fungus:method) # remove interaction term
summary(cox2.1) # suspect that this is best model due to high significance and effect sizes

cox2.2 <- update(cox2.1,~. -method) # try remove method
summary(cox2.2) # much worse

cox2.3 <- update(cox2.1,~. -fungus) # try remove fungus
summary(cox2.3) # not much worse but still think 2.1 is best

AICc(cox2,cox2.1,cox2.2,cox2.3) # check with AICc

survival2result <- cox2.1 # save for ease
rm(cox2,cox2.1,cox2.2,cox2.3) # remove

tab_model(survival2result,
          p.val = "kr",
          file = "model_summaries/Survival_Experiment_2.html") # output model summary

# ---- Survival Exp 3 ----

cox3 <- coxph(Surv(survival, censored) ~ fungus*feeding, data = survival3) # max model
summary(cox3) # nothing sig

cox3.1 <- update(cox3,~. - fungus:feeding) # remove interaction
summary(cox3.1)

cox3.2 <- update(cox3.1,~. -fungus) # remove fungus
summary(cox3.2) # nothing sig

cox3.3 <- update(cox3.1,~. -feeding)
summary(cox3.3)

cox3.4 <- update(cox3.3,~. -fungus)
summary(cox3.4)

survival3result<-AICc(cox3,cox3.1,cox3.2,cox3.3,cox3.4) # no effect, as best model is 3.3
survival3result <- tibble::rownames_to_column(survival3result, "Model") # add model names
survival3result$Model <- fct_recode(survival3result$Model,
                                    "fungus*feeding" = "cox3",
                                    "fungus+feeding" = "cox3.1",
                                    "feeding" = "cox3.2",
                                    "fungus" = "cox3.3",
                                    "none" = "cox3.4") # helpful names
tab_df(survival3result,
       file = "model_summaries/Suvival_Experiement_3.html") # output AICc values

rm(cox3,cox3.1,cox3.2,cox3.3,cox3.4)

# ---- Consumption Exp 1 ----

hist(consump1$consumption) # skewed
hist(log(consump1$consumption)) # log most promising transformation

# log vs norm
consump1_lme1 <- lmer(log(consumption)~factor(day)*factor(spore)*fungus+(1|id), data = consump1) # maximum model
summary(consump1_lme1)
tab_model(consump1_lme1)
p <- plot_model(consump1_lme1, type = "diag")
p[[2]] <- p[[2]]$id
pplot <- plot_grid(p)

ggsave("logConsumption_Experiment_1_diag.png",
             plot = pplot,
             path = "model_summaries",
             dpi = "retina", 
             width = 20,
             height = 16,
             type = "cairo")

consump1_lme2 <- update(consump1_lme1,~. -factor(day):fungus)
summary(consump1_lme2)
plot_model(consump1_lme2, type = "diag")

consump1_lme3 <- update(consump1_lme2,~. -factor(day):factor(spore))
summary(consump1_lme3)
plot_model(consump1_lme3, type = "diag")

consump1_lme4 <- update(consump1_lme3,~. -factor(spore):fungus)
summary(consump1_lme4)
plot_model(consump1_lme4, type = "diag")

consump1_lme5 <- update(consump1_lme4,~. -fungus)
summary(consump1_lme5)
plot_model(consump1_lme5, type = "diag")

consump1_lme6 <- update(consump1_lme5,~. -factor(spore))
summary(consump1_lme6)
plot_model(consump1_lme6, type = "diag")

consump1_lme7 <- update(consump1_lme6,~. -factor(day):factor(spore):fungus)
summary(consump1_lme7)
plot_model(consump1_lme7, type = "diag") # really not fantasitc, but the data convinces me that there is no effect of fungus or spore anayway

AICc(consump1_lme1,consump1_lme2,consump1_lme3,consump1_lme4,consump1_lme5,consump1_lme6,consump1_lme7)

consump1_lme <- consump1_lme7 # save for ease
rm(consump1_lme1,consump1_lme2,consump1_lme3,consump1_lme4,consump1_lme5,consump1_lme6,consump1_lme7) # clean up

tab_model(consump1_lme,
          p.val = "kr",
          file = "model_summaries/Consumption_Exp_1.html") # save model

# ---- Consumption Exp 2 ----

hist(consump2$consumption, breaks = 20)
hist(log(consump2$consumption), breaks = 20) # bit better

consump2_lme1 <- lmer(log(consumption)~factor(day)*method*fungus+(1|id), data = consump2)
summary(consump2_lme1)
q <- plot_model(consump2_lme1, type = "diag")
q[[2]] <- q[[2]]$id
qplot <- plot_grid(q)

ggsave("Consumption_Experiment_2_diag.png",
       plot = qplot,
       path = "model_summaries",
       dpi = "retina", 
       width = 20,
       height = 16,
       type = "cairo")

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

# ---- 2 GLMM ----
glmm1 <- glmmTMB(consumption~factor(day)*method*fungus+(1|id), data = consump2, family = Gamma(link = "inverse"))
summary(glmm1)

numcols <- grep("^c\\.",names(glmm1))
glmm1s <- glmm1
glmm1s <- scale(glmm1s)
m1_sc <- update(m1,data=dfs)

tab_model(glmm1)

library("DHARMa")
check_gamma_model <- simulateResiduals(fittedModel = glmm1, n = 500)
plot(check_gamma_model)

sim_glmm1 <- simulate(glmm1, nsim = 100)
sim_glmm1 <- do.call(cbind, sim_glmm1)
head(sim_glmm1)
sim_res_glmm1 <- createDHARMa(simulatedResponse = sim_glmm1,
                              observedResponse = consump2$consumption,
                              fittedPredictedResponse = predict(glmm1),
                              integerResponse = TRUE)


# ---- Consumption Exp 3 ----

hist(consump3$consumption)
hist(log(consump3$consumption))

consump3_lme1 <- lmer(consumption~factor(day)*feeding*fungus+(1|id), data = consump3) # max model
summary(consump3_lme1) # will try removing 3 way interaction
s <- plot_model(consump3_lme1, type = "diag") # not amazing but ok
s[[2]] <- s[[2]]$id
splot <- plot_grid(s)

ggsave("Consumption_Experiment_3_diag.png",
       plot = splot,
       path = "model_summaries",
       dpi = "retina", 
       width = 20,
       height = 16,
       type = "cairo")

consump3_lme2 <- update(consump3_lme1,~. -factor(day):feeding:fungus)
summary(consump3_lme2)
plot_model(consump3_lme2, type = "diag") # not amazing but ok

consump3_lme3 <- update(consump3_lme2,~. -feeding:fungus)
summary(consump3_lme3)
plot_model(consump3_lme3, type = "diag")

consump3_lme4 <- update(consump3_lme3,~. -factor(day):fungus)
summary(consump3_lme4)
plot_model(consump3_lme4, type = "diag")

consump3_lme5 <- update(consump3_lme4,~. -factor(day):feeding)
summary(consump3_lme5)
plot_model(consump3_lme5, type = "diag")

consump3_lme6 <- update(consump3_lme5,~. -feeding)
summary(consump3_lme6)
plot_model(consump3_lme6, type = "diag")

consump3_lme7 <- update(consump3_lme6,~. -fungus)
summary(consump3_lme7)
plot_model(consump3_lme7, type = "diag")

consump3_lme8 <- update(consump3_lme7,~. -factor(day))
summary(consump3_lme8)

AICc(consump3_lme1,consump3_lme2,consump3_lme3,consump3_lme4,consump3_lme5,consump3_lme6,consump3_lme7,consump3_lme8)
consump3_lme <-  lmer(log(consumption)~factor(day)+(1|id), data = consump3) # something was weird in the update
rm(consump3_lme1,consump3_lme2,consump3_lme3,consump3_lme4,consump3_lme5,consump3_lme6,consump3_lme7,consump3_lme8)

summary(consump3_lme)
tab_model(consump3_lme,
          p.val = "kr",
          file = "model_summaries/Consumption_Exp_3.html") # no random effect

# ---- Mass Experiment 1 ----

mass1_lm1 <- lm(mass_diff~f_total*fungus*spore, data = mass1) # plot full model
summary(mass1_lm1) # doesnt seem to be effect of fungus or interaction
plot_model(mass1_lm1, type = "diag") # diag seems ok

mass1_lm2 <- update(mass1_lm1,~. -f_total:fungus:spore) # remove interaction
summary(mass1_lm2) # std, error of fungus is huge
plot_model(mass1_lm2, type = "diag") # diag seems better

mass1_lm3 <- update(mass1_lm2,~. -fungus:spore) # remove fungus
summary(mass1_lm3) # all sig with decent effect size of mass_diff - most likely best model
plot_model(mass1_lm3, type = "diag") # all fine

mass1_lm4 <- update(mass1_lm3,~. -f_total:spore)
summary(mass1_lm4)

mass1_lm5 <- update(mass1_lm4,~. -f_total:fungus)
summary(mass1_lm5)

mass1_lm6 <- update(mass1_lm5,~. -spore)
summary(mass1_lm6)

mass1_lm7 <- update(mass1_lm6,~. -fungus)
summary(mass1_lm7)



AIC(mass1_lm1,mass1_lm2,mass1_lm3,mass1_lm4,mass1_lm5,mass1_lm6,mass1_lm7) # confirm with AIC

mass_model <- mass1_lm7 # id model for ease
rm(mass1_lm1,mass1_lm2,mass1_lm3,mass1_lm4,mass1_lm5,mass1_lm6,mass1_lm7) # remove candidate models

tab_model(mass_model,
          p.val = "kr",
          file = "model_summaries/Mass_Experiment_1.html") # output as html table

# ---- Agar Experiment ----

agar <- prop.test(x = c(0, 6), n = c(9, 9))
agar

# ---- End ----