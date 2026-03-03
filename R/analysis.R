# the purpose of this script is to run and save all necessary analyses
# then other code can be used to generate the figures and tables

#### Setup ####
suppressPackageStartupMessages({
  library(dplyr)
  library(multinma)
  library(brms)
  library(ggplot2)
  library(forcats)
  library(here)
  library(glue)
})


# reload box module
if(length(ls()) > 1 || length(ls()) ==0) {
  rm(list=ls())
  box::purge_cache()
  if(rlang::is_interactive()) {
    box::use(../R/vit)
  } else {
    box::use(R/vit)
  }
} else {
  box::reload(vit)
}

# set brms options
options(brms.backend = "cmdstanr")

# setup data holders
data <- NULL
outcome <- NULL
outcomes <- c("depression", "anxiety")

# null out all fits
fit_time <- fit_def <- fit_bias <- 
  fit_ssri_low <- fit_high <- fit_some <- 
  fit_low <- NULL

# setup brms arguments
refresh <- 0L
silent  <- 2L
threads = brms::threading(static = TRUE)

# seeds from random.org via random package
seeds   <- list(
  depression = c(
    "low" = 370549758,
    "some" = 236147528,
    "high" = 994407320,
    "some_only" = 661575489,
    "high_only" = 313509298,
    "ssri" = 283532347,
    "ssris"= 282362689,
    "bias" = 1293696306,
    "deficiency" = 1762515477,
    "deficiency_low" = 453698177
  ),
  anxiety = c(
    "low"  = 723753096,
    "some" = 841640689,
    "high" = 914570506,
    "some_only" = 966292209,
    "high_only" = 638547309,
    "ssri" = 23375810,
    "ssris"= 844769116,
    "bias" = 416662781,
    "deficiency" = 979046022,
    "deficiency_low" = 82389386
    )
)

# just run depression data in interactive mode
if (rlang::is_interactive()) {
  outcomes <- outcome <- c("depression")
  silent <- 1L
  refresh <- 100L
} 


#### Run Models and Save ####
for (outcome in outcomes) {
  outdir <- here::here("outputs", "saved_models", outcome)
  
  if (!file.exists(outdir)) {
    dir.create(outdir, 
               showWarnings = FALSE, recursive = TRUE)
  }
  
  if (!file.exists(here::here(outdir, "stan"))) {
    dir.create(here::here(outdir, "stan"), 
               showWarnings = FALSE, recursive = TRUE)
  }
  
  outcome_label <- switch(outcome,
                          depression = "depressed",
                          anxiety = "anxiety")
  
  data <-  vit$get_vitamin_data(outcome = outcome)
  
  cat("Running analyses for outcome:", outcome, "\n")
  #### Basic models ####
  cat("  Basic models...\n")
  fit_low <- vit$prep_brms_nma(data = data %>%
                               filter(bias == "low") %>% 
                               filter(control != "z.antidepressant")  # remove antidepressants b/c network not connected
                             , family = "student-t"
                             , confounders = confounders ~ 0  
                             , main = main ~ 0 + .trt
  ) %>% vit$fit_brms_nma(
    seed = seeds[[outcome]]["low"],
    refresh = refresh, silent = silent,
    threads = threads,
    file = here::here(outdir, "stan",
                          glue::glue("model_", 
                                     outcome, 
                                     "_low_bias_data")
    ),
    control = list(adapt_delta = 0.95, max_treedepth=15L))
  
  fit_some <- vit$prep_brms_nma(data = data %>% 
                      filter(bias != "high") %>% 
                      filter(control != "z.antidepressant")  # remove antidepressants b/c network not connected
                    , family = "student-t"
                    , confounders = confounders ~ 0  
                    , main = main ~ 0 + .trt
  ) %>% vit$fit_brms_nma(
    seed = seeds[[outcome]]["some"],
    warmup = 3000, iter = 4000,
    refresh = refresh, silent = silent,
    threads = threads,
    file = here::here(outdir, "stan",
                      glue::glue("model_", 
                                 outcome, 
                                 "_some_and_low_bias_data")
    ),
    control = list(
                   metric = "dense_e",
                   max_treedepth = 15
                   , adapt_delta = 0.99
                   , init_buffer = 100
                   , window = 50
                   , term_buffer = 100
    ))
  
  fit_high <- vit$prep_brms_nma(data = data
                    , family = "student-t"
                    , confounders = confounders ~ 0  
                    , main = main ~ 0 + .trt
  ) %>% vit$fit_brms_nma(
    seed = seeds[[outcome]]["high"],
    refresh = refresh, silent = silent,
    threads = threads,
    warmup = 3000, iter = 4000
    , control = list(
      metric = "dense_e",
      adapt_delta = 0.95, 
      max_treedepth=15L
      , init_buffer = 100
      , window = 50
      , term_buffer = 100
      ),
    file = here::here(outdir, "stan",
                         glue::glue("model_", 
                                    outcome, 
                                    "_all_data")
      )) 
  
  fit_some_only <- vit$prep_brms_nma(data = data %>% 
                                  filter(bias == "some concerns") %>% 
                                  filter(control != "z.antidepressant")  # remove antidepressants b/c network not connected
                                , family = "student-t"
                                , confounders = confounders ~ 0  
                                , main = main ~ 0 + .trt
  ) %>% vit$fit_brms_nma(
    seed = seeds[[outcome]]["some_only"],
    warmup = 3000, iter = 4000,
    refresh = refresh, silent = silent,
    threads = threads,
    file = here::here(outdir, "stan",
                      glue::glue("model_", 
                                 outcome, 
                                 "_some_data")
    ),
    control = list(
      metric = "dense_e",
      max_treedepth = 15
      , adapt_delta = 0.99
      , init_buffer = 100
      , window = 50
      , term_buffer = 100
    ))
  
  fit_high_only <- vit$prep_brms_nma(data = data %>% 
                                       filter(bias == "high")
                                , family = "student-t"
                                , confounders = confounders ~ 0  
                                , main = main ~ 0 + .trt
  ) %>% vit$fit_brms_nma(
    seed = seeds[[outcome]]["high_only"],
    refresh = refresh, silent = silent,
    threads = threads,
    warmup = 3000, iter = 4000
    , control = list(
      metric = "dense_e",
      adapt_delta = 0.95, 
      max_treedepth=15L
      , init_buffer = 100
      , window = 50
      , term_buffer = 100
    ),
    file = here::here(outdir, "stan",
                      glue::glue("model_", 
                                 outcome, 
                                 "_high_data")
    )) 
  
  if (outcome == "depression") { # no such studies for anxiety
    fit_ssri_low <- vit$prep_brms_nma(data = data %>%
                                        filter(bias == "low") %>% 
                                        filter(control == "z.antidepressant")
                                      , family = "student-t"
                                      , confounders = confounders ~ 0  
                                      , main = main ~ 0 + .trt
    ) %>% vit$fit_brms_nma(
      seed = seeds[[outcome]]["ssri"],
      refresh = refresh, silent = silent,
      threads = threads,
      warmup = 3000, iter = 4000,
      file = here::here(outdir, "stan",
                        glue::glue("model_", 
                                   outcome, 
                                   "_ssri_low_data")
      ),
      , control = list(metric = "dense_e", adapt_delta = 0.99, max_treedepth=15L
                       , init_buffer = 250
                       , window = 500
                       , term_buffer = 250
                       )) 
    
    saveRDS(fit_ssri_low, here::here(outdir, 
                                     glue::glue("model_", 
                                                outcome, 
                                                "_low_bias_ssri_data.rds")
    )
    )
    
    fit_ssri_some <- vit$prep_brms_nma(data = data %>%
                                        filter(bias != "high") %>% 
                                        filter(control == "z.antidepressant")
                                      , family = "student-t"
                                      , confounders = confounders ~ 0  
                                      , main = main ~ 0 + .trt
    ) %>% vit$fit_brms_nma(
      seed = seeds[[outcome]]["ssris"],
      refresh = refresh, silent = silent,
      threads = threads,
      warmup = 3000, iter = 4000,
      file = here::here(outdir, "stan",
                        glue::glue("model_", 
                                   outcome, 
                                   "_ssri_some_and_low_data")
      ),
      , control = list(metric = "dense_e", adapt_delta = 0.99, max_treedepth=15L
                       , init_buffer = 250
                       , window = 500
                       , term_buffer = 250
      )) 
    
    saveRDS(fit_ssri_some, here::here(outdir, 
                                     glue::glue("model_", 
                                                outcome, 
                                                "_low_and_some_bias_ssri_data.rds")
    )
    )
  }
  
  
  saveRDS(fit_high,     here::here(outdir, 
                                   glue::glue("model_", 
                                              outcome, 
                                        "_all_data.rds")
                                   )
          )
  saveRDS(fit_some,     here::here(outdir, 
                                   glue::glue("model_", 
                                              outcome, 
                                              "_some_and_low_bias_data.rds")
                                   )
          )
  saveRDS(fit_high_only,     here::here(outdir, 
                                   glue::glue("model_", 
                                              outcome, 
                                              "_high_bias_data.rds")
  )
  )
  saveRDS(fit_some_only,     here::here(outdir, 
                                   glue::glue("model_", 
                                              outcome, 
                                              "_some_bias_data.rds")
  )
  )
  saveRDS(fit_low,      here::here(outdir, 
                                   glue::glue("model_", 
                                              outcome, 
                                              "_low_bias_data.rds")
                                   )
          )
  
  #### Combined model for study biases ####
  cat("  Bias adjusted model...\n")
  fit_bias <- vit$prep_brms_nma(data = data 
                                    , confounders = confounders ~ 0  
                                    , interactions = interactions ~ 0 + (0 + .trt || bias)
                                    , main = main ~ 0 + .trt ,
                                family = "student-t"
  ) %>% 
    vit$fit_brms_nma( seed = seeds[[outcome]]["bias"],
                      refresh = refresh, silent = silent,
                      threads = threads,
                      file = here::here(outdir, "stan",
                                        glue::glue("model_", 
                                                   outcome, 
                                                   "_bias adjusted")
                      ),
                    control = list(adapt_delta = 0.95, max_treedepth=15L))
  
  saveRDS(fit_bias, here::here(outdir, 
                               glue::glue("model_", 
                                          outcome, 
                                          "_all_data_by_bias.rds")
                               )
          )
  
  #### vitamin deficiency ####
  cat("  Deficiency and depression adjusted model...\n")
  # interactions_formula <- glue::glue("interactions ~ 0 + (0 + .trt || bias/def/{outcome_label})") %>% as.formula()
  # main_formula <- glue::glue("main ~ 0  +.trt") %>% as.formula()
  # interactions_formula <- glue::glue("interactions ~ 0") %>% as.formula()
  # main_formula <- glue::glue("main ~ 0  + .trt + .trt:bias + .trt:bias:def:{outcome_label} + .trt:bias:def + .trt:bias:{outcome_label} + .trt:def + .trt:{outcome_label}") %>% as.formula()
  interactions_formula <- glue::glue("interactions ~ 0") %>% as.formula()
  main_formula <- glue::glue("main ~ 0  + .trt + .trt:def:{outcome_label} + .trt:def + .trt:{outcome_label} + (0 + .trt || bias)") %>% as.formula()
  fit_def <- vit$prep_brms_nma(data = data %>% 
                                          mutate(!!outcome_label :=  
                                                   factor(.data[[outcome_label]], ,
                                                                    levels = c(1:0),
                                                                    labels = c("yes","no")))
                                        , family = "student-t"
                                        , confounders = confounders ~ 0
                                        , interactions = interactions_formula
                                        , main = main_formula
  ) %>%  
    vit$fit_brms_nma(seed = seeds[[outcome]]["deficiency"],
                                   backend = "cmdstanr",
                     refresh = refresh, silent = silent,
                     threads = threads,
                     file = here::here(outdir, "stan",
                                       glue::glue("model_", 
                                                  outcome, 
                                                  "_def_and_dep_and_bias_adjusted")
                     ),
                     control = list(adapt_delta = 0.99, max_treedepth=15L))
  
  saveRDS(fit_def, here::here(outdir, 
            glue::glue("model_", 
                       outcome, 
                       "_all_data_by_deficiency_depression_and_bias.rds")
  )
  )  
  
  interactions_formula_low <- glue::glue("interactions ~ 0") %>% as.formula()
  main_formula_low <- glue::glue("main ~ 0  + .trt + .trt:def:{outcome_label} + .trt:def + .trt:{outcome_label}") %>% as.formula()
  fit_def_low <- vit$prep_brms_nma(data = data %>% 
                                 mutate(!!outcome_label :=  
                                          factor(.data[[outcome_label]], ,
                                                 levels = c(1:0),
                                                 labels = c("yes","no"))) %>% 
                                   filter(bias == "low") %>% 
                                   filter(control != "z.antidepressant")
                               , family = "student-t"
                               , confounders = confounders ~ 0
                               , interactions = interactions_formula_low
                               , main = main_formula_low
  ) %>%  
    vit$fit_brms_nma(seed = seeds[[outcome]]["deficiency_low"],
                     backend = "cmdstanr",
                     refresh = refresh, silent = silent,
                     threads = threads,
                     file = here::here(outdir, "stan",
                                       glue::glue("model_", 
                                                  outcome, 
                                                  "_def_and_dep_low_bias")
                     ),
                     control = list(adapt_delta = 0.99, max_treedepth=15L))
  
  saveRDS(fit_def_low, here::here(outdir, 
                              glue::glue("model_", 
                                         outcome, 
                                         "_low_bias_data_by_deficiency_depression.rds")
  )
  )  
  
}


