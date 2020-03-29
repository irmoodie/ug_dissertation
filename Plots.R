# Thesis Script - Data Visualization
# Plots
# IR Moodie
# Last Modified = 29/03/2020

# ---- Start ----

rm(list=ls()) # clear environment

# ---- Packages ----

library(survival) # survival objects
library(survminer) # ggplot2 plotting of surv objects
library(cowplot) # plot themes
library(tidyverse) # the good stuff

# ---- Data Import ----

exp1 <- read_csv("exp1.csv", na = "na") # Experiment 1 dataset
exp2 <- read_csv("exp2.csv", na = "na") # Experiment 2 dataset
exp3 <- read_csv("exp2_effective_f.csv", na = "na") # Modified Exp 2 dataset so day is now equal to feeding day so a day where no food was give is not counted (gave up trying to do this in R so done manually)

# ---- Data tidying ----

consump1 <- gather(exp1, key = "day", value = "consumption", f0:f4) %>%
  separate(day,into=c("junk", "day"),-1) %>% 
  mutate(day = as.numeric(day)) %>%
  select(-mass_i, -mass_f, -junk, -survival, -censored) # transform consumption dataset from wide to long

consump1$day_f <- factor(consump1$day) %>%
  fct_recode("Day -1" = "0", "Day 0" = "1", "Day 1" = "2", "Day 2" = "3", "Day 3" = "4") # change day names

consump1$fungus_ordered <- factor(consump1$fungus, levels = c("dead", "alive")) # change order of fungal treatments

consump2 <- gather(exp2, key = "day", value = "consumption", f0:f4) %>%
  separate(day,into=c("junk", "day"),-1) %>% 
  mutate(day = as.numeric(day)) %>%
  select(-junk, -survival, -censored) %>%
  filter(feeding == "constant") # wide to long for exp 2

consump2$day_f <- factor(consump2$day) %>%
  fct_recode("Day -1" = "0", "Day 1" = "1", "Day 2" = "2", "Day 3" = "3", "Day 4" = "4")

consump2$method_f <- fct_recode(consump2$method, "Injection" = "injection", "Spray" = "spray") # capitalise for increased pretty-ness

consump2$fungus_ordered <- factor(consump2$fungus, levels = c("control", "alive")) # match order with previous graphs

consump3 <- gather(exp3, key = "day", value = "consumption", f0:f4) %>%
  separate(day,into=c("junk", "day"),-1) %>% 
  mutate(day = as.numeric(day)) %>%
  select(-junk, -survival, -censored) %>%
  filter(method == "spray", day != 3, day != 4)
  
consump3$day_f <- factor(consump3$day) %>%
  fct_recode("Feeding day -1" = "0", "Feeding day 0" = "1", "Feeding day 1" = "2") # set day to effective feeding day
consump3$feeding_f <- fct_recode(consump3$feeding, "Constant" = "constant", "Enforced anorexia" = "starved") # better titles
consump3$fungus_ordered <- factor(consump3$fungus, levels = c("control", "alive")) # sort order

survival1 <- exp1 %>%
  select(-mass_i, -mass_f, -f0, -f1, -f2, -f3, -f4) # survival data of interest in exp1

survival2 <- exp2 %>%
  filter(feeding == "constant") %>%
  select(-f0, -f1, -f2, -f3, -f4) # survival data of interest in exp2

survival2$method_f <- fct_recode(survival2$method, "Injection" = "injection", "Spray" = "spray")

survival3 <- exp2 %>%
  filter(method == "spray") %>%
  select(-f0, -f1, -f2, -f3, -f4) # survival data of interest in exp3

survival3$feeding_f <- fct_recode(survival3$feeding, "Constant" = "constant", "Enforced anorexia" = "starved")

rm(exp1, exp2, exp3) # remove original datasets

# ---- Consumption Plot Exp 1 ----

ggconsump1 <- ggplot(consump1, aes(y = consumption, x = fungus_ordered, colour = fungus_ordered)) + # setup objecy
  geom_jitter(height = 0, width = 0.1, alpha = 0.5) + # add points with jitter and alpha values
  scale_color_manual(values=c("#998ec3", "#f1a340")) + # set colours to match theme
  facet_grid(spore~day_f) + # facet by spore and day
  theme_half_open() + # cowplot pretty +1
  panel_border() + # more cowplot pretty
  theme(legend.position = "none") + # get rid of legend (added after)
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank()) # remove x axis text, ticks and title as well as Y title (titles added after)

ggconsump1 # show graph

ggsave("Consumption_Experiment_1.png",
       plot = ggconsump1,
       dpi = "retina",
       width = 10,
       height = 8,
       type = "cairo") # export graph

# ---- Consumption Plot Exp 2 ----

ggconsump2 <- ggplot(consump2, aes(y = consumption, x = fungus_ordered, colour = fungus_ordered)) + # setup object
  geom_jitter(height = 0, width = 0.1, alpha = 0.5) + # add points with jitter and alpha values
  scale_color_manual(values=c("#998ec3", "#f1a340")) + # add theme colours
  facet_grid(method_f~day_f) + # facet by method and day
  theme_half_open() + # pretty +1
  panel_border() + # pretty +2
  theme(legend.position = "none") + # get rid of legend (added after)
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank()) # remove x elements and y axis title (added after)

ggconsump2 # show graph (DAY 0 is added after)

ggsave("Consumption_Experiment_2.png",
       plot = ggconsump2,
       dpi = "retina", 
       width = 10,
       height = 4,
       type = "cairo") # export graph

# ---- Consumption Plot Exp 3 ----

ggconsump3 <- ggplot(consump3, aes(y = consumption, x = fungus_ordered, colour = fungus_ordered)) + # setup object
  geom_jitter(height = 0, width = 0.1, alpha = 0.5) + # add points with jitter and alpha
  scale_color_manual(values=c("#998ec3", "#f1a340")) + # add theme colours
  facet_grid(feeding_f~day_f) + # facet by feeding regime and day
  theme_half_open() + # cowplot pretty-ness
  panel_border() + # more pretty-ness
  theme(legend.position = "none") + # remove legend (added later)
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank()) # x ticks and titles remove and added later

ggconsump3 # show graph

ggsave("Consumption_Experiment_3.png",
       plot = ggconsump2b,
       dpi = "retina",
       width = 6,
       height = 4,
       type = "cairo") # export

# ---- Surivial Plot Exp 1 ----

spore.labs <- c("0","160","1447","13025") # create labels so that plots appear in order
names(spore.labs) <- c("a", "b", "c", "d") # add abcd so that label column can be used later

ggsurv1 <- ggsurvplot(
  fit = survfit(Surv(survival, censored) ~ fungus+label, data = survival1), # surv obj
  data = survival1,
  ylab = "Survival probability",
  xlab = "Days after injection",
  conf.int = FALSE,
  ylim = c(0, 1),
  palette = alpha(c("#f1a340", "#f1a340", "#f1a340", "#f1a340", "#998ec3", "#998ec3", "#998ec3", "#998ec3"),
                  c(0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7)) # colour scheme + alpha (v messy but works)
)

ggsurv1 <- ggsurv1$plot + # use ggplot to theme using cowplot
  theme_half_open() +
  background_grid() +
  panel_border() +
  scale_x_continuous(breaks = seq(from = 0, to = 10, by = 2)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
  facet_grid(cols = vars(label), labeller = labeller(label = spore.labs)) + # use labeller feature set up previously
  theme(legend.position = "none")

ggsurv1 # show plot

ggsave("Survival_Experiment_1.png", # export plot
       plot = ggsurv1,
       dpi = "retina",
       width = 10, 
       height = 3,
       type = "cairo")

# ---- Survival Plot Exp 2 ----

ggsurv2 <- ggsurvplot(
  fit = survfit(Surv(survival, censored) ~ method_f+fungus, data = survival2), # surv obj
  data = survival2,
  xlab = "Days after injection", 
  ylab = "Survival probability",
  conf.int = FALSE,
  ylim = c(0, 1),
  palette = alpha(c("#f1a340", "#998ec3", "#f1a340", "#998ec3"), # themeing
                  c(0.7,0.7,0.7,0.7))
)

ggsurv2 <- ggsurv2$plot + # use ggplot for themeing
  theme_half_open() +
  background_grid() +
  panel_border() +
  scale_x_continuous(breaks = seq(from = 0, to = 10, by = 2)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
  facet_grid(cols = vars(method_f)) +
  theme(legend.position = "none")

ggsurv2 # show plot

ggsave("Survival_Experiment_2.png", # export plot
       plot = ggsurv2,
       dpi = "retina",
       width = 5.33,
       height = 3,
       type = "cairo")

# ---- Survival Plot Exp 3 ----

ggsurv3 <- ggsurvplot(
  fit = survfit(Surv(survival, censored) ~ feeding_f+fungus, data = survival3),
  data = survival3,
  xlab = "Days after injection", 
  ylab = "Survival probability",
  conf.int = FALSE,
  ylim = c(0, 1),
  palette = alpha(c("#f1a340", "#998ec3", "#f1a340", "#998ec3"),
                  c(0.7,0.7,0.7,0.7))
)

ggsurv3 <- ggsurv3$plot +
  theme_half_open() +
  background_grid() +
  panel_border() +
  scale_x_continuous(breaks = seq(from = 0, to = 10, by = 2)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
  facet_grid(cols = vars(feeding_f)) +
  theme(legend.position = "none")

ggsurv3

ggsave("Survival_Experiment_3.png",
       plot = ggsurv3,
       dpi = "retina",
       width = 5.33,
       height = 3,
       type = "cairo")

# ---- End ----
  