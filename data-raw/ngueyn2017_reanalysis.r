#### re-analysis of Nguyen 2017 ####
# data from https://doi.org/10.1186/s12905-017-0401-3
# adjusting for baseline values of CES-D 

#### Load packages ####
library(dplyr)
library(tidyr)
library(lme4)
library(forcats)

#### Setup data ####
# load original data
nguyen <- read.csv("data-raw/nguyen2017_raw.csv")

# pull out data included in trial
data <- nguyen %>%
  filter(!is.na(cesd_base)) %>% 
  mutate(id = 1:n())

# convert to long
long_EPDS <- data %>% 
  mutate(id = 1:n()) %>% 
  pivot_longer(cols = c("epds_prenatal_1", "epds_prenatal_2",
                        "epds_prenatal_3",
                        "epds_postpartum_3m"),
               names_to = "outcome",
               values_to = "epds") %>% 
  mutate(outcome = factor(outcome) %>% forcats::fct_relevel(
    c("epds_prenatal_1", "epds_prenatal_2",
      "epds_prenatal_3",
      "epds_postpartum_3m")
  ),
    time = rep(0:3, nrow(data))
  )
  

#### Basic demo data
# make sure aligns with published paper
(checks <- data %>% group_by(treatment_cat) %>%
  summarise(n = n(),
            mean_age = mean(age, na.rm = TRUE),
            mean_cesd = mean(cesd_baseline, na.rm = TRUE),
            
  ))

# tests
stopifnot("sample sizes different" = all.equal(checks$n, c(572, 526, 518)))
stopifnot("ages different" = all.equal(round(checks$mean_age,1), c(25.8, 25.9, 26.1)))
stopifnot("ces-d different" = all.equal(round(checks$mean_cesd,2), c(3.66, 3.33, 3.42)))
stopifnot("EPSD different" = all.equal(round(checks$mean_cesd,2), c(3.66, 3.33, 3.42)))

#### Baseline endline means ####
data %>% group_by(treatment_cat) %>%
  summarise(n = n(),
            mean_cesd = mean(Cesd_postpartum_12m, na.rm = TRUE),
            sd_cesd = sd(Cesd_postpartum_12m, na.rm = TRUE)
  )

data %>% group_by(treatment_cat) %>%
  summarise(n = n(),
            mean_cesd = mean(epds_postpartum_3m, na.rm = TRUE),
            sd_cesd = sd(epds_postpartum_3m, na.rm = TRUE)
  )

#### Check missing values ####
data %>% group_by(treatment_cat) %>% summarize(total = n(), n =sum(!is.na(epds_postpartum_3m)),  miss = sum(is.na(epds_postpartum_3m)), m = mean(epds_postpartum_3m, na.rm = TRUE))

data %>% group_by(treatment_cat) %>% summarize(total = n(), n =sum(!is.na(Cesd_postpartum_12m)),  miss = sum(is.na(Cesd_postpartum_12m)), m = mean(Cesd_postpartum_12m, na.rm = TRUE), gain = mean(Cesd_postpartum_12m - cesd_baseline, na.rm = TRUE), sd.gain = sd(Cesd_postpartum_12m - cesd_baseline, na.rm = TRUE))

#### Mixed model EPDS ####
mod <- long_EPDS %>% lmer(epds ~ (time) * treatment_cat +
                            scale(ses) + scale(cesd_baseline) +
                            scale(log(hb3)) + scale(logfer) + as.factor(mo_edu) +
                            anemia + scale(age)
                          + (1|id)
                   , data = .)
mod %>% summary
cbind(mean = c(multi = fixef(mod)["time:treatment_catMM"] * 3, 
        iron_folate = fixef(mod)["time:treatment_catIFA"] * 3),
      se = c(3*sqrt(diag(vcov(mod)))[c("time:treatment_catMM","time:treatment_catIFA")]
      ))
length(fixef(mod))
sigma(mod)
long_EPDS %>% group_by(id) %>% summarize(n = n()) %>% .$n %>% table()

#### model CES-D
fit <- data %>% lm(Cesd_postpartum_12m ~ treatment_cat +
                            scale(ses) + scale(cesd_baseline) +
                            scale(log(hb3)) + scale(logfer) + as.factor(mo_edu) +
                            anemia +  scale(age)
                          , data = .)
fit %>% summary
cbind(mean = c(multi = coef(fit)["treatment_catMM"], 
               iron_folate = coef(fit)["treatment_catIFA"]),
      se = c(sqrt(diag(vcov(fit)))[c("treatment_catMM","treatment_catIFA")]
      ))
sigma(fit)
length(coef(fit))
