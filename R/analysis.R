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
  {
  # if (outcome == "depression") { # no such studies for anxiety
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

#### adverse events model ####
adv_mod_fun <- function(low_bias = TRUE, seed = NULL) {
  adverse_events_data <- vit$get_vitamin_data(outcome = NA_character_,
                                              simple_analysis = TRUE) %>% 
    distinct(study, intervention, .keep_all = TRUE) %>% 
    { if(low_bias == TRUE || low_bias == "main") {
      filter(.,bias == "low", control != "z.antidepressant") 
    } else if (low_bias == "ssri") {
      filter(.,bias == "low", control == "z.antidepressant")
    } else if (low_bias == FALSE || low_bias == "all") {
       .
      } 
      }%>%
    mutate(adverse.count = ifelse(is.na(adverse.count) & !is.na(adverse.censor.dir),
                                  adverse.censor.cut, adverse.count)) %>%
    mutate(cens = dplyr::case_when(
      is.na(adverse.censor.dir) ~ "none", # censored if no events but sample size reported
      adverse.censor.dir == ">=" ~ "right",
      adverse.censor.dir == "<=" ~ "left",
      TRUE ~ NA
    )) %>% 
    mutate(intervention = factor(intervention)) %>% 
    mutate(intervention = relevel(intervention, 
                                  if(low_bias != "ssri") "placebo" else "antidepressant"
                                  )) %>% 
    filter(!is.na(adverse.count))
    # filter(!is.na(adverse.count) | !is.na(adverse.censor.dir))
  
  adv_net <- multinma::set_agd_arm(
    data = adverse_events_data %>% filter(!is.na(adverse.count)),
    study = "study",
    trt = "intervention",
    r = "adverse.count",
    n = "final.N",
    trt_ref = if(low_bias == "ssri") "antidepressant" else "placebo"
  )
  
  adverse_prep <- vit$prep_brms_nma(data = adverse_events_data,
                                  family = "gaussian",
                                  confounders = confounders ~ 0,
                                  interaction = interactions ~ 0,
                                  main = main ~ 0 + .trt  )
  adverse_form <-  adverse_prep$formula$formula
  rlang::f_lhs(adverse_form) <- rlang::expr(adverse_count | trials(final.N) + cens(cens))
  adverse_prep$formula$formula <- adverse_form
  
  trt_mat <-  adverse_events_data %>% 
    rename(.trt = intervention) %>%
    model.matrix(~ .trt, data = .) %>% 
    .[,-1] %>% vit$special_clean_names()
  
  trt_form <- adverse_prep$formula$pforms$main %>% 
    rlang::f_rhs()
  
  rhs <- rlang::expr(
    0 + (!!trt_form) + .study + (1 | gr(.obs_re, cov = R))
  )
  
  full_formula <- rlang::new_formula(
    lhs = rlang::expr(adverse.count | trials(final.N) + cens(cens)),
    rhs = rhs
  )
  
  adverse_formula <- brms::bf(full_formula)
  
  # adverse_formula <- brms::bf(
  #   adverse.count | trials(final.N) + cens(cens) ~ 0 +
  #     .trtantidepressant + .trtantidepressant_b1_b2_b6 + .trtantidepressant_b12 +
  #     .trtantidepressant_b9 + .trtantidepressant_d + .trtantidepressant_magnesium +
  #     .trtantidepressant_zinc + .trtb1 + .trtb12 + .trtb12_b6_b9 + .trtb12_b9 +
  #     .trtb12_b9_d + .trtb9 + .trtd + .trtd_iron + .trtd_zinc + .trtmagnesium +
  #     .trtmagnesium_zinc + .trtselenium + .trtzinc +
  #     study + (1 | gr(.obs_re, cov = R))
  # )
  
  prior <- brms::prior(normal(0, 2.5), class = "b") + brms::prior(normal(0, 1), class = "sd", group = ".obs_re")
  
  R <- multinma::RE_cor(study = adverse_events_data$study,
                        trt = adverse_events_data$intervention, 
                        contrast = rep(FALSE, nrow(adverse_events_data)), 
                        type = "reftrt")
  which_R <- multinma::which_RE(study = adverse_events_data$study,
                              trt = adverse_events_data$intervention, 
                              contrast = rep(FALSE, nrow(adverse_events_data)), 
                              type = "reftrt")
  R_full <- diag(nrow(adverse_events_data))
  R_full[which_R, which_R] <- R
  rownames(R) <- colnames(R) <- paste0("re_",1:nrow(R))
  rownames(R_full) <- colnames(R_full) <- paste0("re_",1:nrow(R_full))
  dat <- cbind(adverse_events_data, trt_mat) %>%
    mutate(
      .obs_re = paste0("re_", 1:nrow(adverse_events_data))
    ) %>% 
    mutate(.obs_re = factor(.obs_re)) %>% 
    mutate(.study = study)
  
  temp_fit <- brms::brm(
    formula = adverse_formula,
    data = dat,
    data2 = list(R = R_full),
    prior = prior,
    family = binomial(),
    chains = 4, iter = 1000, warmup = 1000, empty = TRUE
    )
  
  sd <- brms::make_standata(
    formula = adverse_formula,
    data = dat,
    data2 = list(R = R_full),
    prior = prior,
    family = binomial()
  )
  
  stan_data <- c(
    sd,
    list(
      N_re = nrow(R),
      re_idx = as.array(which_R)  )
  )
  stan_data$Lcov_1 <- chol(R)
  stan_data$J_1 <- as.array(which_R) 
  
  {
    stan_code <- brms::make_stancode(
    formula = adverse_formula,
    data = dat,
    data2 = list(R = R_full),
    prior = prior,
    family = binomial())  %>% as.character()
  
  stan_code <- sub("array\\[M_1\\] vector\\[N_1\\] z_1;", "array[M_1] vector[N_re] z_1;", stan_code)
  stan_code <- sub("vector\\[N_1\\] r_1_1;", "vector[N_re] r_1_1;", stan_code)
  stan_code <- sub("mu[n] += r_1_1[J_1[n]] * Z_1_1[n];", "if(J_1[n] != 0) mu[n] += r_1_1[J_1[n]] * Z_1_1[n];" , stan_code, 
                   fixed = TRUE)
  stan_code <- sub("array\\[N\\] int<lower=1> J_1;","array[N] int<lower=0> J_1;",stan_code)
  stan_code <- sub("matrix[N_1, N_1] Lcov_1;", "matrix[N_re, N_re] Lcov_1;", stan_code,
                   fixed = TRUE)
  
  stan_code <-  stan_code %>% stringr::str_split("\n") %>% .[[1]]
  data_block_start <- sapply(stan_code, function(s) grepl("^data \\{", s)) %>% which()
  stan_code <- c(
    stan_code[1:(data_block_start)],
    "  int<lower=0> N_re;",
    stan_code[(data_block_start+1):length(stan_code)]
  )
  
  stan_code <- paste(stan_code, collapse = "\n")
  }
  
  stan_file <- cmdstanr::write_stan_file(stan_code)
  
  mod <- cmdstanr::cmdstan_model(stan_file)
  
  fit <- mod$sample(
    data = stan_data,
    seed = seed,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.99,
    max_treedepth = 15
  )
  
  adverse_fit <- temp_fit 
  adverse_fit$fit <-  brms::read_csv_as_stanfit(fit$output_files(), model = mod) 
  adverse_fit <- brms::rename_pars(adverse_fit)
  adverse_fit$prep <- list(
    formula = adverse_formula,
    data = dat,
    family = "binomial",
    priors = prior,
    X_nma = trt_mat,
    R_l = chol(R),
    stanvar = NULL,
    network = adv_net,
    orig_form = adverse_prep$orig_form,
    orig_family = "binomial",
    combos = adverse_prep$combos,
    reference_treatment = if(low_bias == "ssri") "antidepressant" else "placebo",
    class_ids = NULL
  ) 
  adverse_fit$data2 <- list(R = R)
  adverse_fit$model <- stan_code
  attr(adverse_fit$ranef, "levels") <- list(.obs_re = rownames(R))
  adverse_fit$basis$group_levels$.obs_re <- rownames(R)
  
  fit_adverse <- adverse_fit
  class(fit_adverse) <- c("vitfit", class(fit_adverse))
  
  return(fit_adverse) 
}
fit_adv_low <- adv_mod_fun(low_bias = TRUE, seed = 669030569)
fit_adv_low_ssri <- adv_mod_fun(low_bias = "ssri", seed = 787040255)
fit_adv <- adv_mod_fun(low_bias = FALSE, seed = 606213864 )

saveRDS(fit_adv_low, here::here("outputs", "saved_models", "adverse_events_low_bias.rds"))
saveRDS(fit_adv_low_ssri, here::here("outputs", "saved_models", "adverse_events_low_bias_ssri.rds"))
saveRDS(fit_adv, here::here("outputs", "saved_models", "adverse_events.rds"))


#### Secondary analysis of vitamin levels

# conceptual framework to figure out particulars of this analysis:
low_bias <- TRUE
level_data <- vit$get_vitamin_data(outcome = NA_character_,
                                   simple_analysis = TRUE) %>% 
  distinct(study, intervention, .keep_all = TRUE) %>% 
  group_by(study) %>% 
  mutate(end.D = end.D * 1e9, sd.end.D = sd.end.D * 1e9) %>%
  mutate(
    idx = which.min(priority),
    y.D = end.D - end.D[idx],
    sd.y.D = sqrt(sd.end.D^2/final.N + sd.end.D[idx]^2/final.N[idx]),
    y.D = ifelse(1:n() == idx, NA_real_, y.D),
    sd.y.D = ifelse(1:n() == idx, sd.y.D/sqrt(2), sd.y.D)
  ) %>%
  { if(low_bias == TRUE || low_bias == "main") {
    filter(.,bias == "low", control != "z.antidepressant") 
  } else if (low_bias == "ssri") {
    filter(.,bias == "low", control == "z.antidepressant")
  } else if (low_bias == FALSE || low_bias == "all") {
    .
  } 
  }

d.dat <- level_data %>% 
  filter(study %in% (level_data %>% filter(intervention == "D") %>% pull(study)) ) %>% 
  filter(!is.na(sd.y.D))

d.dat %>% 
  multinma::set_agd_contrast(
  study = "study",
  trt = "intervention",
  y = "y.D",
  se = "sd.y.D",
  sample_size = "final.N",
  trt_ref = "placebo"
) %>% plot() + ggtitle("Change in vitamin D levels: low bias studies only")

d.prep <- d.dat %>% 
  mutate(y = y.D, sd.y = sd.y.D) %>%
  vit$prep_brms_nma(
  family = "gaussian",
  confounders = confounders ~ 0,
  interactions = interactions ~ 0,
  main = main ~ 0 + .trt,
  trt_ref = "placebo"
) 
d.prep$priors <- d.prep$priors %>% 
  mutate(
  prior = if_else(nlpar == "main" & class == "b" & coef != "Intercept", "normal(0.0, 100)", as.character(prior))
) %>% 
  mutate(prior = if_else(class == "tau", "normal(0.0, 100)", as.character(prior)))

d.fit <- d.prep %>% 
  vit$fit_brms_nma(seed = 12345, refresh = 0L, silent = 2L, threads = threads)



# now cycle through all available endline data...
levels_run <- function(low_bias = TRUE) {
  level_data <- vit$get_vitamin_data(outcome = NA_character_,
                                     simple_analysis = TRUE) %>% 
    distinct(study, intervention, .keep_all = TRUE) %>% 
    { if(low_bias == TRUE || low_bias == "main") {
      filter(.,bias == "low", control != "z.antidepressant") 
    } else if (low_bias == "ssri_low") {
      filter(.,bias == "low", control == "z.antidepressant")
    } else if (low_bias == "ssri") {
      filter(., control == "z.antidepressant")
    } else if (low_bias == FALSE || low_bias == "all") {
      filter(., control != "z.antidepressant")
    } 
    }
  
  seed_levels <- switch(as.character(low_bias),
                        "TRUE" = 598999925,
                        "FALSE" = 359509057,
                        ssri_low = 338935590,
                        ssri = 283532347,
  )
  end.levels <- names(level_data)[startsWith(names(level_data), "end.") %>% which()] %>% setdiff("end.vitamin.sd")
  out <- vector("list", length(end.levels)) %>% setNames(end.levels)
  
  for(i in  end.levels) {
    micronut <- stringr::str_remove(i, "end.") %>% stringr::str_remove("\\.sd")
    sd.col   <- paste0("sd.", i)
    if( !(sd.col %in% names(level_data)) ) {
      next
    }
    
    dat <- level_data %>% 
      filter( study %in% 
                (level_data %>% filter(!is.na(.data[[i]]) ) %>% pull(study)) 
      ) %>% 
      filter( !is.na(.data[[i]]) ) %>% 
      filter( !is.na(.data[[sd.col]])) %>% 
      mutate( y = .data[[i]], sd.y = .data[[sd.col]] ) %>% 
      mutate(unit = vit$vector_convert_vitamin_from_molar(micronut,y)$unit ) %>% 
      mutate( y = vit$vector_convert_vitamin_from_molar(micronut,y)$value,
              sd.y = vit$vector_convert_vitamin_from_molar(micronut,sd.y)$value
      ) %>%
      group_by(study) %>% 
      mutate(
        idx = which.min(priority),
        y = y - y[idx],
        sd.y = sqrt(sd.y^2/final.N + sd.y[idx]^2/final.N[idx]),
        y = ifelse(1:n() == idx, NA_real_, y),
        sd.y = ifelse(1:n() == idx, sd.y/sqrt(2), sd.y)
      )  %>% 
      filter(n() > 1) %>%
      ungroup()
    
    if (nrow(dat) == 0) {
      next
    }
    trt_ref <- switch(i,
                      "end.B6" = "magnesium",
                      "end.C"  = "antidepressant",
                      "placebo")
    if(low_bias == "ssri" || low_bias == "ssri_low") trt_ref <- "antidepressant"
    prep <- dat %>% 
      vit$prep_brms_nma(
        family = "gaussian",
        confounders = confounders ~ 0,
        interactions = interactions ~ 0,
        main = main ~ 0 + .trt,
        trt_ref = trt_ref
      ) 
    
    prep$priors <- 
      if(isFALSE(low_bias) && i=="end.magnesium") {
        prep$priors %>% 
          mutate(
            prior = if_else(nlpar == "main" & class == "b" & coef != "Intercept", "normal(0.0, 1)", as.character(prior))
          ) %>% 
          mutate(prior = if_else(class == "tau", "normal(0.0, 1)", as.character(prior)))
      } else {
        prep$priors %>% 
      mutate(
        prior = if_else(nlpar == "main" & class == "b" & coef != "Intercept", "normal(0.0, 100)", as.character(prior))
      ) %>% 
      mutate(prior = if_else(class == "tau", "normal(0.0, 50)", as.character(prior)))
      }
    
    fit <- prep %>% 
      vit$fit_brms_nma(seed = seed_levels + which(end.levels == i), 
                       refresh = 0L, silent = 2L, 
                       threads = threads,
                       iter = 3000, warmup = 2000,
                       control = list(adapt_delta = 0.99, 
                                      max_treedepth=17L))
    fit$vitamin <- micronut
    fit$unit <- dat$unit[1]
    out[[i]] <- fit
    
  }
  return(out)
}

level_low      <- levels_run(TRUE )
level_ssri_low <- levels_run("ssri_low")
level_ssri     <- levels_run("ssri")
level_all      <- levels_run(FALSE)

saveRDS(level_low, here::here("outputs", "saved_models", "vitamin_levels_low_bias.rds"))
saveRDS(level_ssri_low, here::here("outputs", "saved_models", "vitamin_levels_low_bias_ssri.rds"))
saveRDS(level_ssri, here::here("outputs", "saved_models", "vitamin_levels_all_ssri.rds"))
saveRDS(level_all, here::here("outputs", "saved_models", "vitamin_levels_all.rds"))
