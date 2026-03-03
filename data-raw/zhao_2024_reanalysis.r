#### re-analysis of Zhao 2024 ####
# data for https://doi.org/10.1017/S0033291724000539
#osf: https://osf.io/dntb2/

#### Load packages ####
library(dplyr)
library(tidyr)
library(forcats)

#### Setup data ####
# load original data
zhao <- read.csv("data-raw/zhao2024_raw.csv") %>% 
  mutate(ID = gsub("-2","",ID),
         Group = forcats::fct_recode(as.factor(Group),
                                        Placebo = "NVD",
                                        "Vitamin D" = "VD",
                                        "Healthy" = "HC"),
         antdep_num = Antidepressants.1,
         Time = forcats::fct_recode(as.factor(Time),
                             Baseline = "Baseline2",
                             Baseline = "Baseline1",
                             "Follow-up" = "Follow-up1",
                             "Follow-up" = "Follow-up2")
         ) %>%
  filter(Group != "Healthy") 

data <- zhao %>%
  filter(Time == "Baseline") %>% 
  left_join(
    zhao %>%
      filter(Time == "Follow-up") %>%
      select(ID, HAMD, HAMA),
    by = "ID",
    suffix = c("_base", "_follow")
  ) 


#### Basic demo data
# make sure aligns with published paper
(checks <- data %>% group_by(Group) %>%
    summarise(n = n(),
              age = mean(age, na.rm = TRUE),
              gender = mean(gender, na.rm = TRUE),
              hamd = mean(HAMD_base, na.rm = TRUE),
              sd.hamd = sd(HAMD_base, na.rm = TRUE),
              hama = mean(HAMA_base, na.rm = TRUE),
              sd.hama = sd(HAMA_base)
    ))

# tests
stopifnot("sample sizes different" = all.equal(checks$n, c(26,20)))
stopifnot("ages different" = all.equal(round(checks$age,2), c(42.42,46.75)))
stopifnot("gender different" = all.equal(round(checks$gender*100,2), c(34.62,25)))
stopifnot("HAMA different" = all.equal(round(checks$hama,2), c(21.23,21.5))) # looks like error in writing of paper
stopifnot("HAMD different" = all.equal(round(checks$hamd,2), c(32.77,29.9))) # looks like error in paper
stopifnot("HAMD different" = all.equal(round(checks$sd.hamd,2), c(11.54,10.1))) # looks like error in paper
stopifnot("HAMA different" = all.equal(round(checks$sd.hama,2), c(7.24,8.24))) # looks like error in paper

# check missingness
sum(is.na(data$HAMA_base))
sum(is.na(data$HAMD_base))
sum(is.na(data$HAMA_follow))
sum(is.na(data$HAMD_follow))
nrow(data)

#### Endline means ####
data %>% group_by(Group) %>%
  summarise(n = n(),
            hamd = mean(HAMD_follow, na.rm = TRUE),
            sd.hamd = sd(HAMD_follow, na.rm = TRUE),
            hama = mean(HAMA_follow, na.rm = TRUE),
            sd.hama = sd(HAMA_follow, na.rm = TRUE)
  )

#### Gain Scores ####
data %>% group_by(Group) %>%
  summarise(n = sum(!is.na(HAMD_follow - HAMD_base)),
            hamd = mean(HAMD_follow - HAMD_base, na.rm = TRUE),
            sd.hamd = sd(HAMD_follow - HAMD_base, na.rm = TRUE),
            hama = mean(HAMA_follow - HAMA_base, na.rm = TRUE),
            sd.hama = sd(HAMA_follow - HAMA_base, na.rm = TRUE)
  )


#### HAMD model ####
fit <- data %>% lm(HAMD_follow ~ Group +
                     scale(Intervention.duration..months.) +
                     HAMD_base +
                     gender + age + scale(Illness.duration..months.)
                   + antdep_num
                   , data = .)
fit %>% summary
cbind(mean = c(vitD = coef(fit)["GroupVitamin D"]),
      se = c(sqrt(diag(vcov(fit)))[c("GroupVitamin D")]
      ))
sigma(fit)
length(coef(fit))


#### HAMA model ####
fit1 <- data %>% lm(HAMA_follow ~ Group +
                     scale(Intervention.duration..months.) +
                     HAMA_base +
                     gender + age + scale(Illness.duration..months.)
                   + antdep_num
                   , data = .)
fit1 %>% summary
cbind(mean = c(vitD = coef(fit1)["GroupVitamin D"]),
      se = c(sqrt(diag(vcov(fit1)))[c("GroupVitamin D")]
      ))
sigma(fit1)
length(coef(fit1))
