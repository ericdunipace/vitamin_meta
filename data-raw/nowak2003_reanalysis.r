#### re-analysis of Zhao 2024 ####
# Nowak, G.; Siwek, M.; Dudek, D.; Zieba, A.; Pilc, A.
# Pol J Pharmacol 2003;55(6):1143-7
# 2003
# raw data printed in paper

#### Load packages ####
library(dplyr)
library(tidyr)
library(forcats)

#### Setup data ####
data <- read.csv(text = 'id, sex, age, hdrs.base, hdrs.final, beck.base, beck.final, group\n
1, f, 48, 20, 23, 42, 46,placebo\n
2, f, 41, 23, 9, 37, 27,placebo\n
3, f, 28, 28, 0, 37, 4,placebo\n
4, f, 55, 21, 5, 34, 26,placebo\n
5, m, 46, 28, 15, 31, 27,placebo\n
6, f, 47, 17, 13, 31, 28,placebo\n
7, f, 35, 20, 12, 23, 31,placebo\n
8, m, 47, 24, 13, 31, 28,placebo\n
9, m, 25, 22, 0, 29, 2,zinc\n
10, f, 53, 25, 12, 43, 35,zinc\n
11, m, 43, 32, 8, 31, 17,zinc\n
12, m, 26, 22, 4, 29, 25,zinc\n
13, m, 57, 20, 5, 30, 15,zinc\n
14, f, 49, 26, 2, 30, 3,zinc')

#### original estimates ####

data %>%
  group_by(group) %>% 
  summarize(hdrs = mean(hdrs.base),
            beck = mean(beck.base),
            n = n(),
            sd.hdrs = sd(hdrs.base),
            sd.beck = sd(beck.base))


data %>%
  group_by(group) %>% 
  summarize(hdrs = mean(hdrs.final),
            beck = mean(beck.final),
            n = n(),
            sd.hdrs = sd(hdrs.final),
            sd.beck = sd(beck.final))

#### get gain scores ####
data %>%
  group_by(group) %>% 
  summarize(hdrs = mean(hdrs.final - hdrs.base),
            beck = mean(beck.final - beck.base),
            n = n(),
            sd.hdrs = sd(hdrs.final - hdrs.base),
            sd.beck = sd(beck.final - beck.base))


#### linear model ####
fit_hdrs <- lm(hdrs.final ~ hdrs.base + group + age + sex, data = data) 
fit_beck <- lm(beck.final ~ beck.base + group + age + sex, data = data)

rbind(hdrs = summary(fit_hdrs)$coefficients["groupzinc",c("Estimate","Std. Error")],
      beck = summary(fit_beck)$coefficients["groupzinc",c("Estimate","Std. Error")])


length(fit_hdrs %>% coef)
length(fit_beck %>% coef)

sd(data$hdrs.final)
sd(data$beck.final)
