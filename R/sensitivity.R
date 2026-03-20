#### Setup ####
suppressPackageStartupMessages({
  library(dplyr)
  library(multinma)
  library(brms)
  library(ggplot2)
  library(forcats)
  library(glue)
  library(purrr)
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

# setup brms arguments
refresh <- 0L
silent  <- 2L
open_progress <- show_messages <- FALSE
threads = brms::threading(static = TRUE)

# seeds from random.org via random package
seeds   <- list(
  depression = c(multinma_fe = 69790750, 
                 multinma_re = 897755013, 
                 multinma_ume = 203868595, 
                 multinma_nodesplit = 812808054,  
                 multinma_ssri = 283532347,
                 multinma_high = 965529339, 
                 multinma_some = 899233348, 
                 multinma_nodehigh = 181523355, 
                 multinma_nodesome = 93084568,
                 ume_low = 518218714,
                 ume_some = 539544866,
                 ume_high = 791398097,
                 low = 1524405328, 
                 some = 2100830690, 
                 high = 994407320,
                 node_low = 853559197,
                 node_some = 240659605,
                 node_high = 184848027,
                 low_nogain = 437222330,
                 some_nogain = 645324198,
                 high_nogain = 271294043,
                 low_def_no_iran = 824002557,
                 def_no_iran = 372384114,
                 loo_low = 534263806, 
                 loo_some = 223221691, 
                 loo_high = 301721616,
                 loo_def_low = 948963238, 
                 loo_def_some = 183033670, 
                 loo_def_high = 944863566,
                 loo_low_gauss = 755392234, 
                 loo_some_gauss = 111584572, 
                 loo_high_gauss = 501002750,
                 loo_low_mna = 218098603,
                 loo_some_mna = 221657322,
                 loo_high_mna = 278781759,
                 loo_high_ume = 850838444,
                 loo_some_ume = 810467841,
                 loo_low_ume  = 154336884,
                 metareg = 236147528,
                 loo_metareg = 278892917,
                 "time" = 282669673,
                 "dose" = 440686404,
                 "dose2" = 673844248,
                 "time_cat" = 195682257,
                 "time_cat_low" = 218781043,
                 "pct" = 882723733),
  anxiety = c(multinma_fe = 44327420, 
              multinma_re = 551660263, 
              multinma_ume = 315681837, 
              multinma_nodesplit = 117567625,  
              multinma_ssri = 444869089,
              multinma_high = 159116459, 
              multinma_some = 195688763, 
              multinma_nodehigh = 277809323, 
              multinma_nodesome = 851630707,
              ume_low = 838782390,
              ume_some = 526198001,
              ume_high = 208987270,
              node_low = 361255549,
              node_some = 208851420,
              node_high = 559089640,
              low = 932005756, 
              some = 396567330, 
              high = 612224036,
              low_nogain = 121661413,
              some_nogain = 926110296,
              high_nogain = 52218309,
              low_def_no_iran = 336701971,
              def_no_iran = 137657926,
              loo_low = 506780338, 
              loo_some = 405652248, 
              loo_high = 266741761,
              loo_def_low = 294449283, 
              loo_def_some = 376626608, 
              loo_def_high = 45402035,
              loo_low_gauss = 859182597, 
              loo_some_gauss = 45028443, 
              loo_high_gauss = 332197140,
              loo_low_mna = 155603997,
              loo_some_mna = 692649922,
              loo_high_mna = 923649567,
              loo_high_ume = 884007497,
              loo_some_ume = 861011706,
              loo_low_ume  = 676296796,
              metareg = 706417974,
              loo_metareg = 566372216,
              "time" = 970531477,
              "dose" = 864543463,
              "dose2" = 806432633,
              "time_cat" = 68143789,
              "time_cat_low" = 448908120,
              pct = 271845657
)
)

# just run depression data in interactive mode
if (rlang::is_interactive()) {
  outcomes <- outcome <- c("depression")
  silent <- 1L
  refresh <- 100L
} 

# run sensitivity analyses for each outcome
cat("Running sensitivity analyses for outcomes:", paste(outcomes, collapse = ", "), "\n")
for(outcome in outcomes) {
  cat("Processing outcome:", outcome, "\n")
  
  # setup label for switching between outcomes in formulae
  outcome_label <- switch(outcome,
                          depression = "depressed",
                          anxiety = "anxiety")
  
  
  #### create output directory ####
  outdir <- here::here("outputs", "saved_models", outcome,
                       "sensitivity")
  
  if(!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  if(!dir.exists(here::here(outdir, "stan")) ) {
    dir.create(here::here(outdir, "stan"), recursive = TRUE)
  }
  
  #### Data ####
  data <- vit$get_vitamin_data(outcome = outcome)
  
  net <- vit$construct_nma_network(data %>%
                                       filter(bias == "low") %>%
                                       filter(control != "z.antidepressant"),
                                   trt_ref = "placebo")
  net_some <- vit$construct_nma_network(data %>%
                                            filter(bias != "high") %>%
                                            filter(control != "z.antidepressant"),
                                        trt_ref = "placebo")
  net_all <- vit$construct_nma_network(data,
                                       trt_ref = "placebo")

  #### run multinma models ####
  cat("  Running multinma models...\n")
  # run fixed vs random effects
  fit_FE <- nma(net,
                  trt_effects = "fixed",
                  prior_trt = normal(0, 2.5),
                  prior_het = normal(0, 0.5),
                  seed = seeds[[outcome]]["multinma_fe"],
                  cores = 4L,
                refresh = refresh,
                open_progress = open_progress,
                show_messages = show_messages
  )

  fit_RE <- nma(net,
                  trt_effects = "random",
                  prior_trt = normal(0, 2.5),
                  prior_het = normal(0, 0.5),
                  cores = 4L,
                  seed = seeds[[outcome]]["multinma_re"],
                  control = list(max_treedepth = 15L),
                refresh = refresh,
                open_progress = open_progress,
                show_messages = show_messages
  )

  # check unrelated means to see if there's inconsistency
  UME_RE <- nma(net,
                  trt_effects = "random",
                  consistency = "ume",
                  prior_trt = normal(0, 2.5),
                  prior_het = normal(0, 0.5),
                  cores = 4L,
                  seed = seeds[[outcome]]["multinma_ume"],
                  control = list(max_treedepth = 15L),
                refresh = refresh,
                open_progress = open_progress,
                show_messages = show_messages
  )
  if (outcome == "depression") {
    
    net_ssri <- vit$construct_nma_network(data %>%
                                              filter(bias == "low") %>%
                                              filter(control == "z.antidepressant"),
                                            trt_ref = "antidepressant"
    )
    
    fit_ssri <- nma(net_ssri,
                    trt_effects = "random",
                    prior_trt = normal(0, 2.5),
                    prior_het = normal(0, 0.5),
                    cores = 4L,
                    seed = seeds[[outcome]]["multinma_ssri"],
                    control = list(adapt_delta = 0.999, max_treedepth = 15L),
                    refresh = refresh,
                    open_progress = open_progress,
                    show_messages = show_messages
    )
    
    saveRDS(fit_ssri,
            file = here::here(outdir,
                              glue::glue("model_",
                                         outcome,
                                         "_mnma_ssri_low_bias.rds"))
    )
    
  }


  fit_RE_all <- nma(net_all,
                      trt_effects = "random",
                      prior_trt = normal(0, 2.5),
                      prior_het = normal(0, 0.5),
                      cores = 4L,
                      seed = seeds[[outcome]]["multinma_high"],
                      control = list(max_treedepth = 15L),
                    refresh = refresh,
                    open_progress = open_progress,
                    show_messages = show_messages
  )

  fit_RE_some <- nma(net_some,
                       trt_effects = "random",
                       prior_trt = normal(0, 2.5),
                       prior_het = normal(0, 0.5),
                       cores = 4L,
                       seed = seeds[[outcome]]["multinma_some"],
                       control = list(max_treedepth = 15L),
                     refresh = refresh,
                     open_progress = open_progress,
                     show_messages = show_messages
  )

  # cat("  Running multinma nodesplit models...\n")
  # nodesplit models
  # if (outcome == "depression") {
  #   node_split_low <- nma(net,
  #                         consistency = "nodesplit",
  #                         trt_effects = "random",
  #                         prior_intercept = normal(scale = 1),
  #                         prior_trt = normal(scale = 1),
  #                         prior_het = normal(0, 0.5),
  #                         cores = 4L,
  #                         seed = seeds[[outcome]]["multinma_nodesplit"],
  #                         control = list(max_treedepth = 15L),
  #                         refresh = refresh,
  #                         open_progress = open_progress,
  #                         show_messages = show_messages
  #   )
  #   
  #   node_split_some <- nma(net_some,
  #                          consistency = "nodesplit",
  #                          trt_effects = "random",
  #                          prior_intercept = normal(scale = 1),
  #                          prior_trt = normal(scale = 1),
  #                          prior_het = normal(0, 0.5),
  #                          cores = 4L,
  #                          seed = seeds[[outcome]]["multinma_nodesome"],
  #                          control = list(max_treedepth = 15L),
  #                          refresh = refresh,
  #                          open_progress = open_progress,
  #                          show_messages = show_messages
  #   )
  #   
  #   saveRDS(node_split_low,
  #           file = here::here(outdir,
  #                             glue::glue("model_",
  #                                        outcome,
  #                                        "_mnma_nodesplit_low_bias.rds"))
  #   )
  #   saveRDS(node_split_some,
  #           file = here::here(outdir,
  #                             glue::glue("model_",
  #                                        outcome,
  #                                        "_mnma_nodesplit_some_and_low_bias.rds"))
  #   )
  # }
  # 
  # # takes  a while
  # node_split_high <- nma(net_all,
  #                        consistency = "nodesplit",
  #                        trt_effects = "random",
  #                        prior_intercept = normal(scale = 1),
  #                        prior_trt = normal(scale = 1),
  #                        prior_het = normal(0, 0.5),
  #                        cores = 4L,
  #                        seed = seeds[[outcome]]["multinma_nodehigh"],
  #                        control = list(max_treedepth = 15L),
  #                        refresh = refresh,
  #                        open_progress = open_progress,
  #                        show_messages = show_messages
  # )

  saveRDS(fit_FE,
          file = here::here(outdir,
                            glue::glue("model_",
                                       outcome,
                                       "_mnma_fixed_effects.rds"))
  )

  saveRDS(fit_RE,
          file = here::here(outdir,
                            glue::glue("model_",
                                       outcome,
                                       "_mnma_random_effects.rds"))
  )

  saveRDS(UME_RE,
          file = here::here(outdir,
                            glue::glue("model_",
                                       outcome,
                                       "_mnma_ume_random_effects.rds"))
  )

  

  saveRDS(fit_RE_all,
          file = here::here(outdir,
                            glue::glue("model_",
                                       outcome,
                                       "_mnma_random_effects_all_data.rds"))
  )

  saveRDS(fit_RE_some,
          file = here::here(outdir,
                            glue::glue("model_",
                                       outcome,
                                       "_mnma_random_effects_some_and_low_data.rds"))
  )

  # saveRDS(node_split_high,
  #         file = here::here(outdir,
  #                           glue::glue("model_",
  #                                      outcome,
  #                                      "_mnma_nodesplit_all_data.rds"))
  # )

  #### BRMS Gaussian Models ####
  cat("  Running brms Gaussian models...\n")
  fit_brms_low_gaussian <- vit$prep_brms_nma(data = data %>%
                                               filter(bias == "low") %>%
                                               filter(control != "z.antidepressant")  # remove antidepressants b/c network not connected
                                             , main = main ~ 0 + .trt
                                             , family = "gaussian") %>%
    vit$fit_brms_nma(seed = seeds[[outcome]]["low"],
                     refresh = refresh, silent = silent,
                     threads = threads,
                     file = here::here(outdir, "stan",
                                            glue::glue("temp_brms_low_bias_data_gaussian_",
                                                       outcome)),
                     control = list(adapt_delta = 0.95, max_treedepth=15L))


  fit_brms_some_gaussian <- vit$prep_brms_nma(data = data %>%
                                                filter(bias != "high") %>%
                                                filter(control != "z.antidepressant")  # remove antidepressants b/c network not connected
                                              , main = main ~ 0+ .trt
                                              , family = "gaussian") %>%
    vit$fit_brms_nma(seed = seeds[[outcome]]["some"],
                     refresh = refresh, silent = silent,
                     threads = threads,
                     file = here::here(outdir, "stan",
                                       glue::glue("temp_brms_some_bias_data_gaussian_",
                                                  outcome)),
                     control = list(adapt_delta = 0.95, max_treedepth=15L))


  fit_brms_high_gaussian <- vit$prep_brms_nma(data = data
                                              , main = main ~ 0+ .trt
                                              , family = "gaussian") %>%
    vit$fit_brms_nma(seed = seeds[[outcome]]["high"],
                     refresh = refresh, silent = silent,
                     threads = threads,
                     file = here::here(outdir, "stan",
                                       glue::glue("temp_brms_all_data_gaussian_",
                                                  outcome)),
                     control = list(adapt_delta = 0.95, max_treedepth=15L))

  saveRDS(fit_brms_low_gaussian,
          file = here::here(outdir,
                            glue::glue("model_",
                                       outcome,
                                       "_low_bias_data_gaussian.rds"))
  )

  saveRDS(fit_brms_some_gaussian,
          file = here::here(outdir,
                            glue::glue("model_",
                                       outcome,
                                       "_some_and_low_bias_data_gaussian.rds"))
  )

  saveRDS(fit_brms_high_gaussian,
          file = here::here(outdir,
                            glue::glue("model_",
                                       outcome,
                                       "_all_data_gaussian.rds"))
  )
  
  ### UME model with student T ####
  fit_brms_low_ume <- vit$prep_brms_nma(data = data %>%
                                               filter(bias == "low") %>%
                                               filter(control != "z.antidepressant") %>% 
                                          mutate(ume = interaction(intervention,study)) %>% 
                                          droplevels()
                                             , main = main ~ 0 + ume
                                             , family = "student-t") %>%
    vit$fit_brms_nma(seed = seeds[[outcome]]["ume_low"],
                     refresh = refresh, silent = silent,
                     threads = threads,
                     file = here::here(outdir, "stan",
                                       glue::glue("temp_brms_low_bias_data_ume_",
                                                  outcome)),
                     control = list(adapt_delta = 0.95, max_treedepth=15L))
  
  
  fit_brms_some_ume <- vit$prep_brms_nma(data = data %>%
                                                filter(bias != "high") %>%
                                                filter(control != "z.antidepressant") %>% 
                                           mutate(ume = interaction(intervention,study))
                                         , main = main ~ 0 + ume
                                              , family = "student-t" ) %>%
    vit$fit_brms_nma(seed = seeds[[outcome]]["ume_some"],
                     refresh = refresh, silent = silent,
                     threads = threads,
                     file = here::here(outdir, "stan",
                                       glue::glue("temp_brms_some_bias_data_ume_",
                                                  outcome)),
                     control = list(adapt_delta = 0.95, max_treedepth=15L))
  
  
  fit_brms_high_ume <- vit$prep_brms_nma(data = data %>% 
                                           mutate(ume = interaction(intervention,study))
                                         , main = main ~ 0 + ume
                                              , family = "student-t") %>%
    vit$fit_brms_nma(seed = seeds[[outcome]]["ume_high"],
                     refresh = refresh, silent = silent,
                     threads = threads,
                     file = here::here(outdir, "stan",
                                       glue::glue("temp_brms_all_data_ume_",
                                                  outcome)),
                     control = list(adapt_delta = 0.95, max_treedepth=15L))
  
  #### Time course ####
  cat("  Time course model...\n")
  fit_time <- vit$prep_brms_nma(data = data
                                , confounders = confounders ~ 0  +
                                  scale(duration) + 
                                  I(scale(duration)^2),
                                interactions = interactions ~ 0 +
                                  .trt:scale(duration) +
                                  .trt:I(scale(duration)^2) + (.trt || bias)
                                , main = main ~ 0 +  .trt ,
                                family = "student-t"
  ) %>% 
    vit$fit_brms_nma(seed = seeds[[outcome]]["time"],
                     refresh = refresh, silent = silent,
                     threads = threads,
                     file = here::here(outdir, "stan",
                                       glue::glue("model_", 
                                                  outcome, 
                                                  "_time")
                     ),
                     control = list(adapt_delta = 0.95, max_treedepth = 15L))
  
  saveRDS(fit_time, here::here(outdir, 
                               glue::glue("model_", 
                                          outcome, 
                                          "_all_data_by_time.rds")
  )
  )
  
  #### models over time at rough quartiles ####
  # time_quartiles <- quantile(data$duration, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  # 8,12,26
  time_cat_mutate <- function(x) {
    x %>% 
      mutate(time_cat = cut(duration, breaks = c(-Inf, 8, 12, 26, Inf), labels = c("0-8 weeks", "8-12 weeks", "12-26 weeks", "26+ weeks"))) 
      # mutate(time_cat = cut(duration, breaks = c(-Inf, 4, 8, Inf), labels = c("0-4 weeks", "4-8 weeks", ">8 weeks")))
  }
  
  time_cat <- data %>% 
    time_cat_mutate () %>%
    vit$prep_brms_nma(confounders = confounders ~ 0,
                      interactions = interactions ~ 0 + (.trt || time_cat)
                      , main = main ~ 0 +  .trt ,
                      family = "student-t"
    ) %>% 
    vit$fit_brms_nma(seed = seeds[[outcome]]["time_cat"],
                     refresh = refresh, silent = silent,
                     threads = threads,
                     file = here::here(outdir, "stan",
                                       glue::glue("model_", 
                                                  outcome, 
                                                  "_time_cat")
                     ),
                     control = list(adapt_delta = 0.95, max_treedepth = 15L))
  
  
  time_cat_low <- data %>% 
    filter(bias == "low") %>% 
    filter(control != "z.antidepressant") %>%
    time_cat_mutate() %>% 
    vit$prep_brms_nma(confounders = confounders ~ 0,
                      interactions = interactions ~ 0 + (.trt || time_cat)
                      , main = main ~ 0 +  .trt ,
                      family = "student-t"
    ) %>% 
    vit$fit_brms_nma(seed = seeds[[outcome]]["time_cat_low"],
                     refresh = refresh, silent = silent,
                     threads = threads,
                     file = here::here(outdir, "stan",
                                       glue::glue("model_", 
                                                  outcome, 
                                                  "_time_cat_low")
                     ),
                     control = list(adapt_delta = 0.95, max_treedepth = 15L))
  
  saveRDS(time_cat, here::here(outdir, 
                               glue::glue("model_", 
                                          outcome, 
                                          "_all_data_by_time_cat.rds")
  )
  )
  
  saveRDS(time_cat_low, here::here(outdir, 
                                   glue::glue("model_", 
                                              outcome, 
                                              "_low_bias_data_by_time_cat.rds")
  )
  )
  
  #### dose  ####
  cat("  Dose response model...\n")
  interactions_formula <- glue::glue("interactions ~ 0") %>% as.formula()
  main_formula <- glue::glue("main ~ 0  + .trt + .dose(tdd.cat) + (0 + .trt || bias)") %>% as.formula()
  fit_dose <- purrr::quietly(vit$prep_brms_nma)(data = vit$get_vitamin_data(outcome, 
                                                                            simple_analysis = FALSE)
                                                , confounders = confounders ~ 0
                                                , interactions = interactions_formula
                                                , main = main_formula,
                                                family = "student-t"
  )$result %>% 
    vit$fit_brms_nma(seed = seeds[[outcome]]["dose"],
                     refresh = refresh, silent = silent,
                     threads = threads,
                     file = here::here(outdir, "stan",
                                       glue::glue("model_", 
                                                  outcome, 
                                                  "_dose")
                     ),
                     control = list(adapt_delta = 0.95, 
                                    max_treedepth = 15L))
  
  saveRDS(fit_dose, here::here(outdir, 
                               glue::glue("model_", 
                                          outcome,
                                          "_all_data_by_dose.rds")
  )
  )
  fit_dose_low <- purrr::quietly(vit$prep_brms_nma)(data = vit$get_vitamin_data(outcome, 
                                                                                simple_analysis = FALSE) %>% 
                                                      filter(bias == "low") %>% 
                                                      filter(control != "z.antidepressant" &
                                                               control != "z.antidepressant,z.B9")  # remove antidepressants b/c network not connected
                                                    , confounders = confounders ~ 0
                                                    , interactions = interactions_formula 
                                                    , main = glue::glue("main ~ 0  + .trt + .dose(tdd.cat)") %>% as.formula(),
                                                    family = "student-t"
  )$result %>% 
    vit$fit_brms_nma(seed = seeds[[outcome]]["dose2"],
                     refresh = refresh, silent = silent,
                     threads = threads,
                     file = here::here(outdir, "stan",
                                       glue::glue("model_", 
                                                  outcome, 
                                                  "_dose_low_bias")
                     ),
                     control = list(adapt_delta = 0.95, 
                                    max_treedepth = 15L))
  
  saveRDS(fit_dose_low, here::here(outdir, 
                                   glue::glue("model_", 
                                              outcome, 
                                              "_low_bias_data_by_dose.rds")
  )
  )
  
  #### check models 
  
  
  #### check data with pct dis and def ####
  cat("  Running brms models with pct dis and def...\n")
  # fit_pcts <- data %>% filter(bias == "low") %>% 
  #   filter(control != "z.antidepressant") %>% 
  # vit$prep_brms_nma(confounders = confounders ~ 0,
  #                   interactions = interactions ~ 0
  #                   , main = main ~ 0 +  .contrast(.trt : pct_dep + .trt : pct_deficient + .trt : pct_dep : pct_deficient) + .trt
  #                   , family = "student-t"
  # ) %>% 
  #   vit$fit_brms_nma(seed = seeds[[outcome]]["pct"],
  #                    refresh = refresh, silent = silent,
  #                    threads = threads,
  #                    control = list(adapt_delta = 0.95, max_treedepth = 15L))
  # grid <- tidyr::expand_grid(.trt = fit_pcts$prep$network$treatments, pct_dep = c(0,1), pct_deficient = c(0,1))
  # grid %>% model.matrix(~ .trt : pct_dep + .trt : pct_deficient + .trt : pct_dep : pct_deficient + .trt, data = .) %>%  as.data.frame() %>% vit$special_clean_names() %>% mutate(.study = fit_pcts$data$.study[1], .se = fit_pcts$data$.se[1], .cov = fit_pcts$data$.cov[1], pct_deficient = grid$pct_deficient, pct_dep = grid$pct_dep, .trt = grid$.trt) -> newdata
  # pct_sum <- vit$summary_brms_nma(fit_pcts, newdata = newdata, keep = c(".trt", 'pct_dep', "pct_deficient")) %>% 
  #   mutate(pct_dep = ifelse(pct_dep == 1, "yes", "no"),
  #          pct_deficient = ifelse(pct_deficient == 1, "yes", "no")) %>% 
  #   rename(depressed = pct_dep, deficiency = pct_deficient)
  # pct_sum %>% summary()
  # pct_sum %>% mutate(bias = "low") %>% vit$def_plot()
  pct_name <- switch(outcome,
                     "depression" = "pct_dep",
                     anxiety = "pct_anx")
  
  fit_pcts_approx <- data %>% filter(bias == "low") %>% 
    filter(control != "z.antidepressant") %>% 
    group_by(study) %>% 
    mutate(p_def = weighted.mean(pct_deficient, baseline.N)) %>% 
    mutate(p_dep = weighted.mean(.data[[pct_name]], baseline.N)) %>%
    # mutate(p_dep = ifelse(p_dep == 1, .999, p_dep)) %>% 
    ungroup() %>% 
    vit$prep_brms_nma(confounders = confounders ~ 0,
                      interactions = interactions ~ 0
                      , main = main ~ 0 +  .trt : mi(p_dep) + .trt : mi(p_def) + .trt : mi(p_dep) : mi(p_def) + .trt,
                      addl_form = list(
                        brms::bf(p_def | mi() + trunc(0,1) ~ scale(year) + def, family = "gaussian"),
                        brms::bf(p_dep | mi() + trunc(0,1) ~ scale(year) + depressed, family = "gaussian")
                      ), 
                      family = "student-t"
    ) %>% 
    vit$fit_brms_nma(seed = seeds[[outcome]]["pct"],
                     refresh = refresh, silent = silent,
                     threads = threads,
                     control = list(adapt_delta = 0.95, max_treedepth = 15L))
  
  saveRDS(fit_pcts_approx, here::here(outdir, 
                                   glue::glue("model_", 
                                              outcome, 
                                              "_low_bias_data_by_pct_def_and_dep.rds")
  )
  )
  
  #### load t-distribution models ####
  low_t <- readRDS(
    here::here("outputs", "saved_models", outcome, 
               glue::glue("model_", 
                          outcome, 
                          "_low_bias_data.rds"))
  ) #%>% brms::recompile_model()
  
  some_t <- readRDS(
    here::here("outputs", "saved_models", outcome,
               glue::glue("model_", 
                          outcome, 
                          "_some_and_low_bias_data.rds"))
  ) #%>% brms::recompile_model()
  
  high_t <- readRDS(
    here::here("outputs", "saved_models", outcome,
               glue::glue("model_", 
                          outcome, 
                          "_all_data.rds"))
  ) #%>% brms::recompile_model()
  
  fit_def_low <- readRDS(
    here::here("outputs", "saved_models", outcome,
               glue::glue("model_", 
                          outcome, 
                          "_low_bias_data_by_deficiency_depression.rds"))
  )
  
  # fit_def_some <- readRDS(
  #   here::here("outputs", "saved_models", outcome,
  #              glue::glue("model_", 
  #                         outcome, 
  #                         "_some_and_low_bias_data_by_deficiency_depression.rds"))
  # )
  
  # fit_def_high <- readRDS(
  #   here::here("outputs", "saved_models", outcome,
  #              glue::glue("model_", 
  #                         outcome, 
  #                         "_all_data_by_deficiency_depression.rds"))
  # )
  
  fit_def <- readRDS(
    here::here("outputs", "saved_models", outcome,
               glue::glue("model_", 
                          outcome, 
                          "_all_data_by_deficiency_depression_and_bias.rds"))
  )
  
  
  #### Check deficiency models without Iranian studies ####
  interactions_formula_low <- glue::glue("interactions ~ 0") %>% as.formula()
  main_formula_low <- glue::glue("main ~ 0  + .trt + .trt:def:{outcome_label} + .trt:def + .trt:{outcome_label}") %>% as.formula()
  
  fit_def_low_no_iran <- vit$prep_brms_nma(data = data %>% 
                                             mutate(!!outcome_label :=  
                                                      factor(.data[[outcome_label]], ,
                                                             levels = c(1:0),
                                                             labels = c("yes","no"))) %>% 
                                             filter(bias == "low") %>% 
                                             filter(control != "z.antidepressant") %>% 
                                             filter(country != "Iran")
                                           , family = "student-t"
                                           , confounders = confounders ~ 0
                                           , interactions = interactions_formula_low
                                           , main = main_formula_low
  ) %>%  
    vit$fit_brms_nma(seed = seeds[[outcome]]["low_def_no_iran"],
                     refresh = refresh, silent = silent,
                     threads = threads,
                     file = here::here(outdir, "stan",
                                       glue::glue("temp_def_no_iran_low_",
                                                  outcome)),
                     backend = "cmdstanr",
                     control = list(adapt_delta = 0.99, max_treedepth=15L))
  
  interactions_formula <- glue::glue("interactions ~ 0") %>% as.formula()
  main_formula <- glue::glue("main ~ 0  + .trt + .trt:def:{outcome_label} + .trt:def + .trt:{outcome_label} + (0 + .trt || bias)") %>% as.formula()
  fit_def_no_iran <- vit$prep_brms_nma(data = data %>% 
                                 mutate(!!outcome_label :=  
                                          factor(.data[[outcome_label]], ,
                                                 levels = c(1:0),
                                                 labels = c("yes","no"))) %>% 
                                 filter(country != "Iran")
                               , family = "student-t"
                               , confounders = confounders ~ 0
                               , interactions = interactions_formula
                               , main = main_formula
  ) %>%  
    vit$fit_brms_nma(seed = seeds[[outcome]]["def_no_iran"],
                     refresh = refresh, silent = silent,
                     threads = threads,
                     file = here::here(outdir, "stan",
                                       glue::glue("temp_def_no_iran_",
                                                  outcome)),
                     backend = "cmdstanr",
                     control = list(adapt_delta = 0.99, max_treedepth=15L))
  
  saveRDS(fit_def_low_no_iran,
          file = here::here(outdir,
                            glue::glue("model_",
                                       outcome,
                                       "_low_bias_data_by_deficiency_depression_no_iran.rds"))
  )
  
  saveRDS(fit_def_no_iran,
          file = here::here(outdir,
                            glue::glue("model_",
                                       outcome,
                                       "_all_data_by_deficiency_depression_and_bias_no_iran.rds"))
  )
  
  #### check data without gain scores added in ####
  cat("  Running brms models without gain scores...\n")
  no_gain <- vit$get_vitamin_data(outcome = outcome, no_gain_scores = TRUE)

  fit_brms_low_nogain <- vit$prep_brms_nma(data = no_gain %>%
                                               filter(bias == "low") %>%
                                               filter(control != "z.antidepressant")  # remove antidepressants b/c network not connected
                                             , main = main ~ 0 + .trt
                                             , family = "student-t") %>%
    vit$fit_brms_nma(seed = seeds[[outcome]]["low_nogain"],
                     refresh = refresh, silent = silent,
                     threads = threads,
                     file = here::here(outdir, "stan",
                                       glue::glue("temp_nogainscores_low_",
                                                  outcome)),
                     control = list(adapt_delta = 0.95, max_treedepth=15L))


  fit_brms_some_nogain <- vit$prep_brms_nma(data = no_gain %>%
                                                filter(bias != "high") %>%
                                                filter(control != "z.antidepressant")  # remove antidepressants b/c network not connected
                                              , main = main ~ 0+ .trt
                                              , family = "student-t") %>%
    vit$fit_brms_nma(seed = seeds[[outcome]]["some_nogain"],
                     refresh = refresh, silent = silent,
                     threads = threads,
                     warmup = 2000L, iter = 3000L,
                     file = here::here(outdir, "stan",
                                       glue::glue("temp_nogainscores_some_",
                                                  outcome)),
                     control = list(
                                    metric = "dense_e",
                                    adapt_delta = 0.99, max_treedepth=15L
                                    , init_buffer = 250
                                    , window = 500
                                    , term_buffer = 250))


  fit_brms_high_nogain <- vit$prep_brms_nma(data = no_gain
                                              , main = main ~ 0+ .trt
                                              , family = "student-t") %>%
    vit$fit_brms_nma(seed = seeds[[outcome]]["high_nogain"],
                     refresh = refresh, silent = silent,
                     threads = threads,
                     warmup = 2000L, iter = 3000L,
                     file = here::here(outdir, "stan",
                                       glue::glue("temp_nogainscores_high_",
                                                  outcome)),
                     control = list(
                       metric = "dense_e",
                       adapt_delta = 0.99, max_treedepth=15L
                       , init_buffer = 250
                       , window = 500
                       , term_buffer = 250))

  saveRDS(fit_brms_low_nogain,
          file = here::here(outdir,
                            glue::glue("model_",
                                       outcome,
                                       "_low_bias_data_nogain.rds"))
  )

  saveRDS(fit_brms_some_nogain,
          file = here::here(outdir,
                            glue::glue("model_",
                                       outcome,
                                       "_some_and_low_bias_data_nogain.rds"))
  )

  saveRDS(fit_brms_high_nogain,
          file = here::here(outdir,
                            glue::glue("model_",
                                       outcome,
                                       "_all_data_nogain.rds"))
  )
  
  #### Node splits on student-t model ####
  cat("  Running brms nodesplit models...\n")
  node_low <- vit$vit_nodesplit(low_t, seed = seeds[[outcome]]["node_low"])
  saveRDS(node_low, file = here::here(outdir, "nodesplit_low.rds"))
  
  node_some <- vit$vit_nodesplit(some_t, seed = seeds[[outcome]]["node_some"])
  saveRDS(node_some, file = here::here(outdir, "nodesplit_some.rds"))
  
  node_high <- vit$vit_nodesplit(high_t, seed = seeds[[outcome]]["node_high"])
  saveRDS(node_high, file = here::here(outdir, "nodesplit_high.rds"))
  
  #### LOO estimates t-dist ####
  cat("  Running LOO estimates for student-t models...\n")
  # warning, can take a long time to run, currently set to run on 8 cores
  # my computer has 8 high performance cores, hence the number
  loo_low <- vit$loo_estimates(
    fit = low_t,
    seed = seeds[[outcome]]["loo_low"],
    cores = 8L, reuse_metric = FALSE
  )
  saveRDS(loo_low, file = here::here(outdir, "loo_low.rds"))
  
  loo_def_low <- vit$loo_estimates(
    fit = fit_def_low,
    seed = seeds[[outcome]]["loo_def_low"],
    cores = 8L,
    keep = c(".trt","def", outcome_label),
    reuse_metric = TRUE
  )
  saveRDS(loo_def_low, file = here::here(outdir, "loo_def_low.rds"))
  
  # loo_def_some <- vit$loo_estimates(
  #   fit = fit_def_some,
  #   seed = seeds[[outcome]]["loo_def_some"],
  #   cores = 8L
  # )
  # saveRDS(loo_def_some, file = here::here(outdir, "loo_def_some.rds"))
  
  # warning, can take a long time to run
  loo_some <- vit$loo_estimates(
    fit = some_t,
    seed = seeds[[outcome]]["loo_some"],
    cores = 8L,
    reuse_metric = TRUE
  )
  saveRDS(loo_some, file = here::here(outdir, "loo_some.rds"))
  
  
  loo_def <- vit$loo_estimates(
    fit = fit_def,
    seed = seeds[[outcome]]["loo_def_high"],
    cores = 8L,
    keep = c(".trt","def","bias", outcome_label),
    reuse_metric = TRUE
  )
  saveRDS(loo_def, file = here::here(outdir, "loo_def.rds"))
  
  # warning, can take a long time to run
  loo_high <- vit$loo_estimates(
    fit = high_t,
    seed = seeds[[outcome]]["loo_high"],
    cores = 8L,
    reuse_metric = TRUE
  )
  saveRDS(loo_high, file = here::here(outdir, "loo_high.rds"))
  
  #### LOO estimates normal dist ####
  cat("  Running LOO estimates for gaussian models...\n")
  # warning, can take a long time to run, currently set to run on 8 cores
  # my computer has 8 high performance cores, hence the number
  loo_gauss_low <- vit$loo_estimates(
    fit = fit_brms_low_gaussian,
    seed = seeds[[outcome]]["loo_low_gauss"],
    cores = 8L,
    reuse_metric = FALSE
  )
  saveRDS(loo_gauss_low, file = here::here(outdir, "loo_gaussian_low.rds"))
  
  # warning, can take a long time to run
  loo_gauss_some <- vit$loo_estimates(
    fit = fit_brms_some_gaussian,
    seed = seeds[[outcome]]["loo_some_gauss"],
    cores = 8L,
    reuse_metric = FALSE
  )
  saveRDS(loo_gauss_some, , file = here::here(outdir, "loo_gaussian_some.rds"))
  
  # warning, can take a long time to run
  loo_gauss_high <- vit$loo_estimates(
    fit = fit_brms_high_gaussian,
    seed = seeds[[outcome]]["loo_high_gauss"],
    cores = 8L,
    verbose = FALSE,
    reuse_metric = TRUE
  )
  saveRDS(loo_gauss_high, file = here::here(outdir,
                                            "loo_gaussian_high.rds"))
  
  #### LOO estimates of UME ####
  cat("  Running LOO estimates for UME models...\n")
  # warning, can take a long time to run, currently set to run on 8
  # my computer has 8 high performance cores, hence the number
  loo_ume_low <- vit$loo_estimates(
    fit = fit_brms_low_ume,
    seed = seeds[[outcome]]["loo_low_ume"],
    cores = 8L,
    verbose = FALSE,
    get_estimates = FALSE,
    reuse_metric = FALSE
  )
  
  loo_ume_some <- vit$loo_estimates(
    fit = fit_brms_some_ume,
    seed = seeds[[outcome]]["loo_some_ume"],
    cores = 8L,
    verbose = FALSE,
    get_estimates = FALSE,
    reuse_metric = TRUE
  )
  
  loo_ume_high <- vit$loo_estimates(
    fit = fit_brms_high_ume,
    seed = seeds[[outcome]]["loo_high_ume"],
    cores = 8L,
    verbose = FALSE,
    get_estimates = FALSE,
    reuse_metric = TRUE
  )
  
  saveRDS(loo_ume_low, file = here::here(outdir,
                                            "loo_ume_low.rds"))
  
  saveRDS(loo_ume_some, file = here::here(outdir,
                                            "loo_ume_some.rds"))
  
  saveRDS(loo_ume_high, file = here::here(outdir,
                                            "loo_ume_high.rds"))
  
  # #### Loo estimates multnima models ####
  # l_low <- vit$multinma_loo(fit_RE, seed = seeds[[outcome]]["loo_low_mna"])
  # l_some <- vit$multinma_loo(fit_RE_some, seed = seeds[[outcome]]["loo_some_mna"])
  # l_all <- vit$multinma_loo(fit_RE_all, seed = seeds[[outcome]]["loo_high_mna"])
  # 
  # saveRDS(l_low, file = here::here(outdir, "loo_multinma_low.rds"))
  # saveRDS(l_some, file = here::here(outdir, "loo_multinma_some.rds"))
  # saveRDS(l_all, file = here::here(outdir, "loo_multinma_all.rds"))
  
}

cat("  Running meta-regression model...\n")
#### meta-regression analysis #####

{
  formula_bias_reg <- function(outcome_label) {
    def_inter <- glue::glue("+ .trt:def:{outcome_label} + .trt:def + .trt:{outcome_label}")
    if (outcome_label == "anxiety") {
      conf_sd<- "(s_se || country) + "
      inter_sd<- ""
      def_inter_main <- def_inter_conf <- ""
      def_inter_inter <- def_inter
    } else if (outcome_label == "depressed") {
      inter_sd <- "(s_se || country) + "
      conf_sd  <- ""
      def_inter_main <- def_inter
      def_inter_inter <- def_inter_conf <- ""
    } else {
      rlang::abort("Outcome not recognized")
    }
    main_formula <- glue::glue("main ~ 0 + .trt {def_inter_main}  + (0 + .trt || bias.overall)") %>% as.formula()
    interactions_formula <- glue::glue("interactions ~ {def_inter_inter} + {inter_sd} 0") %>% as.formula()
    confounders_formula <- glue::glue("confounders ~ 0 {def_inter_conf} +
    s_year + s_se +
                                      {conf_sd}
                          .contrast(mi(s_male) +  mi(s_age)) + 
                          s_duration + s_durationE2 + 
                          s_dropOut + s_assignedN + 
                          mi(s_baseline.y) +
                          mi(s_groupBaseline) + mi(s_groupBaseline) + 
                          mi(s_groupAge) + (1 + mi(s_groupBaseline) || scale) + 
                          mi(s_groupMale) +
                          bias.selective.outcome  +  
                          bias.other   + bias.overall +
                          bias.blinding.participants + bias.blinding.outcome + 
                          bias.incomplete.outcome +
                          bias.sequence + bias.concealment") %>% 
      as.formula()
    return(
      list(confounders = confounders_formula,
           interactions = interactions_formula,
           main = main_formula)
      )
  }
  bias_data_prep <- function(x) {
    x %>% 
    group_by(study) %>%
    arrange(desc(target)) %>% 
    distinct(study, intervention, .keep_all = TRUE) %>%
      # mutate( groupBaseline_anx = mean(baseline.outcome.mean[target == "anxiety"], na.rm = TRUE))  %>% 
      # filter(target == "depression") %>% 
      mutate(groupAge = mean(age, na.rm = TRUE),
             groupBaseline = mean(baseline.outcome.mean, na.rm = TRUE),
             groupMale = mean(male, na.rm = TRUE),
             dropOut = (1 - (final.N / baseline.N)),
      ) %>% 
      ungroup() %>% 
      mutate(s_duration = scale(duration),
             s_durationE2 = s_duration^2) %>% 
      mutate(
             s_year = scale(year),
             s_groupAge = scale(groupAge),
             s_groupBaseline = scale(groupBaseline),
             s_groupMale = scale(groupMale),
             s_baseline.y = scale(baseline.y),
             s_male = scale(male),
             s_age = scale(age),
             s_dropOut = scale(dropOut),
             s_assignedN = scale(assigned.N)
      ) #%>% 
      # mutate(#y_anx = ifelse(target == "anxiety", y, NA),
      #        y = ifelse(target == "anxiety", NA, y)
      # )
  }
  rewrite_bf_symbols <- function(bf, replacements) {
    rewrite_formula <- function(f, replacements) {
      txt <- deparse(f, width.cutoff = 500L)
      for (nm in names(replacements)) {
        txt <- gsub(
          pattern = nm,
          replacement = replacements[[nm]],
          x = txt,
          perl = TRUE
        )
      }
      ff <- as.formula(txt)
      environment(ff) <- environment(f)
      ff
    }
    
    bf_new <- bf
    for (nm in names(bf)) {
      if (inherits(bf[[nm]], "formula")) {
        bf_new[[nm]] <- rewrite_formula(bf[[nm]], replacements)
      }
    }
    bf_new
  }
  rewrite_brmsformula <- function(bf, replacements) {
    # 1) turn the whole brmsformula into a single string
    txt <- paste0(deparse(bf, width.cutoff = 500L), collapse = "\n")
    
    # 2) apply replacements in a controlled order (longer patterns first is safer)
    ord <- order(nchar(names(replacements)), decreasing = TRUE)
    pats <- names(replacements)[ord]
    
    for (pat in pats) {
      txt <- gsub(pat, replacements[[pat]], txt, perl = TRUE)
    }
    
    # 3) reconstruct brmsformula in the same environment as the original
    eval(parse(text = txt), envir = environment(bf))
  }
  rewrite_one_formula <- function(f, replacements) {
    txt <- paste(deparse(f, width.cutoff = 500L), collapse = "\n")
    
    # apply longer patterns first to reduce cascading surprises
    ord <- order(nchar(names(replacements)), decreasing = TRUE)
    for (pat in names(replacements)[ord]) {
      txt <- gsub(pat, replacements[[pat]], txt, perl = TRUE)
    }
    
    ff <- as.formula(txt)
    environment(ff) <- environment(f)
    ff
  }
  
  rewrite_brmsformula_parts <- function(bf, replacements, rename_pforms = TRUE) {
    bf2 <- bf
    
    # main response formula (includes y | vint() + vreal() ~ ...)
    bf2$formula <- rewrite_one_formula(bf2$formula, replacements)
    
    # sub-formulas (nonlinear parameters / distributional etc.)
    if (length(bf2$pforms)) {
      bf2$pforms <- lapply(bf2$pforms, rewrite_one_formula, replacements = replacements)
      
      if (rename_pforms) {
        # rename list entries too (e.g. main -> main_anx)
        nm <- names(bf2$pforms)
        for (pat in names(replacements)) {
          # only use replacements that are pure word matches for names
          # (i.e., patterns like \\bmain\\b)
          if (grepl("^\\\\b.*\\\\b$", pat)) {
            bare <- sub("^\\\\b", "", sub("\\\\b$", "", pat))
            nm[nm == bare] <- replacements[[pat]]
          }
        }
        names(bf2$pforms) <- nm
      }
    }
    
    # response name tracked separately
    if (!is.null(bf2$resp) && bf2$resp == "y") bf2$resp <- "y_anx"
    
    bf2
  }
  add_idx_to_mi <- function(formula, idx_name = ".obs") {
    
    add_idx_call <- function(expr) {
      if (!is.call(expr)) return(expr)
      
      # recurse first
      expr[] <- lapply(expr, add_idx_call)
      
      # if this is mi(...)
      if (identical(expr[[1]], quote(mi))) {
        args <- as.list(expr)[-1]
        has_idx <- any(names(args) == "idx")
        
        if (!has_idx) {
          args$idx <- as.name(idx_name)
          expr <- as.call(c(list(quote(mi)), args))
        }
      }
      
      expr
    }
    
    if (inherits(formula, "formula")) {
      formula[[3]] <- add_idx_call(formula[[3]])
    } else {
      stop("Input must be a formula")
    }
    
    formula
  }
  
  add_idx_to_bf <- function(bf_obj, idx_name = ".obs") {
    stopifnot(inherits(bf_obj, "brmsformula"))
    
    # 1) main formula
    main <- add_idx_to_mi(formula(bf_obj), idx_name = idx_name)
    
    # 2) distributional formulas, if present (sigma, zi, hu, etc.)
    dpars <- bf_obj$dpars
    if (length(dpars)) {
      for (nm in names(dpars)) {
        if (inherits(dpars[[nm]]$formula, "formula")) {
          dpars[[nm]]$formula <- add_idx_to_mi(dpars[[nm]]$formula, idx_name = idx_name)
        }
      }
    }
    
    # 3) nonlinear parameter formulas, if present
    nlpars <- bf_obj$nlpars
    if (length(nlpars)) {
      for (nm in names(nlpars)) {
        if (inherits(nlpars[[nm]]$formula, "formula")) {
          nlpars[[nm]]$formula <- add_idx_to_mi(nlpars[[nm]]$formula, idx_name = idx_name)
        }
      }
    }
    
    # rebuild a new brmsformula carrying over options
    out <- bf(
      main,
      family = bf_obj$family,
      subset = bf_obj$subset,
      nl = bf_obj$nl,
      center = bf_obj$center,
      sparse = bf_obj$sparse,
      decomp = bf_obj$decomp,
      autocor = bf_obj$autocor,
      cov_ranef = bf_obj$cov_ranef,
      loops = bf_obj$loops
    )
    
    # reattach updated dpars/nlpars
    out$dpars  <- dpars
    out$nlpars <- nlpars
    out
  }
  
  add_idx_to_bf2 <- function(bf_obj, idx_name = ".obs") {
    stopifnot(is.list(bf_obj), "formula" %in% names(bf_obj))
    
    # reuse your existing add_idx_to_mi()
    rewrite_if_formula <- function(x) {
      if (inherits(x, "formula")) add_idx_to_mi(x, idx_name = idx_name) else x
    }
    
    # 1) main formula (mean / location)
    bf_obj$formula <- rewrite_if_formula(bf_obj$formula)
    
    # 2) parameter formulas (sigma, zi, hu, etc.) stored in pforms
    if (!is.null(bf_obj$pforms) && length(bf_obj$pforms)) {
      bf_obj$pforms <- lapply(bf_obj$pforms, rewrite_if_formula)
    }
    
    bf_obj
  }
  add_index_term <- function(formula, idx_name = ".obs") {
    stopifnot(inherits(formula, "brmsformula"))
    
    lhs <- rlang::f_lhs(formula[[1]])
    if (!is.call(lhs)) return(formula)
    
    # We’re looking for: (|)(resp, <specials>)   i.e. resp | ...
    if (!identical(lhs[[1]], quote(`|`))) return(formula)
    
    # lhs[[2]] = response part, lhs[[3]] = specials part (often mi() or some call)
    specials <- lhs[[3]]
    
    # helper to see if specials already contains index(...)
    has_index <- function(expr) {
      if (!is.call(expr)) return(FALSE)
      if (identical(expr[[1]], quote(index))) return(TRUE)
      any(vapply(as.list(expr)[-1], has_index, logical(1)))
    }
    
    if (has_index(specials)) return(formula)
    
    # Add "+ index(.obs)" to the specials
    new_specials <- as.call(list(quote(`+`), specials, as.call(list(quote(index), as.name(idx_name)))))
    
    formula$formula[[2]] <- as.call(list(quote(`|`), lhs[[2]], new_specials))
    formula
  }
  
  add_subset_term <- function(formula, idx_name = ".idx_anx") {
    stopifnot(inherits(formula, "brmsformula"))
    
    lhs <- rlang::f_lhs(formula[[1]])
    if (!is.call(lhs)) return(formula)
    
    # We’re looking for: (|)(resp, <specials>)   i.e. resp | ...
    if (!identical(lhs[[1]], quote(`|`))) return(formula)
    
    # lhs[[2]] = response part, lhs[[3]] = specials part (often mi() or some call)
    specials <- lhs[[3]]
    
    # helper to see if specials already contains index(...)
    has_index <- function(expr) {
      if (!is.call(expr)) return(FALSE)
      if (identical(expr[[1]], quote(subset))) return(TRUE)
      any(vapply(as.list(expr)[-1], has_index, logical(1)))
    }
    
    if (has_index(specials)) return(formula)
    
    # Add "+ index(.idx)" to the specials
    new_specials <- as.call(list(quote(`+`), specials, as.call(list(quote(subset), as.name(idx_name)))))
    
    formula$formula[[2]] <- as.call(list(quote(`|`), lhs[[2]], new_specials))
    formula
  }
  
  
  
  replacements_anx <- c(
    "\\by\\b"                 = "y_anx",
    "\\.se\\b"                = ".se_anx",
    "\\.cov\\b"               = ".cov_anx",
    "\\.study\\b"             = ".study_anx",
    # "\\bzU\\b"                 = "zUAnx",
    # "s_groupBaseline"   = "s_groupBaseline_anx",
    # "s_baseline\\.y"    = "s_baseline.y_anx",
    # "s_se"             = "s_se_anx",
    "\\bs_groupBaseline\\b"   = "s_groupBaseline_anx",
    "\\bs_baseline\\.y\\b"    = "s_baseline.y_anx",
    "\\bs_se\\b"             = "s_se_anx",
    "\\bscale\\b"              = "scale_anx"
    # "\\bmain\\b"               = "mainAnx",
    # "\\binteractions\\b"       = "interactionsAnx",
    # "\\bconfounders\\b"        = "confoundersAnx"
  )
  
  formula_prep_out <- formula_bias_reg("depressed")
  
  depression_reg <- vit$prep_brms_nma(
    data = vit$get_vitamin_data(outcome = NULL, 
                                include_full_bias = TRUE) %>% 
      bias_data_prep() 
  , family = "student-t"
  , confounders = formula_prep_out$confounders
  , interactions = formula_prep_out$interactions
  , main = formula_prep_out$main 
  , addl_form = list(
    brms::bf(s_baseline.y | mi()  + subset(.idx_dep) + index(.obs) ~ 
               s_assignedN +
               .contrast(mi(s_age)) +
               .contrast(mi(s_male)) + 
               psych.dx + def + 
               s_se  ) + gaussian(),
    brms::bf(s_baseline.y_anx | mi() + subset(.idx_anx)  + index(.obs) ~
               s_assignedN +
               .contrast(mi(s_age)) +
               .contrast(mi(s_male)) +
               s_se_anx ) + gaussian(),
    brms::bf(.contrast(s_age) | mi() + index(.obs) ~  
               bias.overall) + gaussian(),
    brms::bf(.contrast(s_male) | mi() + index(.obs) ~    
               (1||country)    +
                bias.overall) + gaussian(),
    brms::bf(s_groupBaseline | mi() + subset(.idx_dep)  + index(.obs)~  
               mi(s_groupAge) + 
               mi(s_groupMale) +
               psych.dx + def + s_se ) + gaussian(),
    brms::bf(s_groupBaseline_anx | mi() + subset(.idx_anx)  + index(.obs) ~
               mi(s_groupAge) +
               mi(s_groupMale) +
               psych.dx + def + s_se_anx) + gaussian(),
    brms::bf( s_groupAge | mi() + index(.obs) ~   s_year + 
                (1||country) ) + gaussian(),
    brms::bf( s_groupMale | mi() + index(.obs) ~   s_year +
                (1||country) ) + gaussian()
  )
  )
  depression_reg$data <- depression_reg$data %>%
    mutate(
      y = ifelse(target == "anxiety", NA, y),
      .se = ifelse(target == "anxiety", NA, .se),
      sd.y = ifelse(target == "anxiety", NA, sd.y),
      sd_se = sd(.se, na.rm = TRUE),
      m_se  = mean(.se, na.rm = TRUE),
      s_se = (.se - m_se) / sd_se
    )
  dep_reg_main_form <- depression_reg$formula 
  dep_reg_main_form$forms$y$formula <- as.formula("y | vint(as.integer(.study)) + vreal(.se, .cov) + subset(.idx_dep) + index(.obs) ~ main + interactions + confounders + 0 ", env = environment(dep_reg_main_form$formula))
  attr(dep_reg_main_form$forms$y$formula, "nl") <- TRUE
  attr(dep_reg_main_form$forms$y$formula, "loop") <- TRUE
  # dep_reg_main_form$forms$sbaseliney <- dep_reg_main_form$forms$sbaseliney %>% add_subset_term(".idx_dep") %>% add_index_term()
  # dep_reg_main_form$forms$sbaseliney <- dep_reg_main_form$forms$sbaseliney %>% add_index_term()
  
  formula_prep_anx <- formula_bias_reg("anxiety")
  anx_reg <- vit$prep_brms_nma(
    data = vit$get_vitamin_data(outcome = "anxiety", 
                                include_full_bias = TRUE) %>% 
      bias_data_prep()
    , family = "student-t"
    , confounders = formula_prep_anx$confounders
    , interactions = formula_prep_anx$interactions
    , main = formula_prep_anx$main)
  
  anx_reg$data <- anx_reg$data %>%
    mutate(
      y = ifelse(target == "depression", NA, y),
      .se = ifelse(target == "depression", NA, .se),
      sd.y = ifelse(target == "depression", NA, sd.y),
      sd_se = sd(.se, na.rm = TRUE),
      m_se  = mean(.se, na.rm = TRUE),
      s_se = (.se - m_se) / sd_se
    )
  anx_reg_main_form <- anx_reg$formula %>% rewrite_brmsformula_parts(replacements_anx)
  anx_reg_main_form$formula <- as.formula("y_anx | vint(as.integer(.study_anx)) + vreal(.se_anx, .cov_anx) + subset(.idx_anx) + index(.obs) ~ main + interactions + confounders + 0 ", env = environment(anx_reg_main_form$formula))
  anx_reg_main_form$family$dpars <- c("mu","tau","zU")
  names(anx_reg_main_form$family$lb) <- 
    names(anx_reg_main_form$family$ub) <- c("mu","tau","zU")
  anx_reg_main_form$family$vars <- c("vint1", "vreal1", "vreal2", "Sigma_anx_L", "R_anx_L")
  
  attr(anx_reg_main_form$formula, "nl") <- TRUE
  attr(anx_reg_main_form$formula, "loop") <- TRUE
  anx_priors <- anx_reg$prior  %>% 
    mutate(resp = "yanx")
  anx_priors[anx_priors$nlpar == "main" & anx_priors$class == "sd", "prior"] <- "normal(0, 0.25)"
  # anx_priors <- anx_priors %>% mutate(
  #   class = ifelse(class == "tau", "tauAnx", class),
  #   dpar  = ifelse(dpar == "zU", "zUAnx", dpar),
  #   nlpar = ifelse(nlpar == "confounders", "confoundersAnx", nlpar),
  #   nlpar = ifelse(nlpar == "interactions", "interactionsAnx", nlpar),
  #   nlpar = ifelse(nlpar == "main", "mainAnx", nlpar)
  #   )
  anx_stanvar <- anx_reg$stanvar[[3]]
  anx_stanvar$scode <- gsub("Sigma","Sigma_anx", anx_stanvar$scode, fixed = TRUE)
  # anx_stanvar$scode <- gsub("Sigma_L","Sigma_L_anx", anx_stanvar$scode, fixed = TRUE)
  anx_stanvar$scode <- gsub("R","R_anx", anx_stanvar$scode, fixed = TRUE)
  # anx_stanvar$scode <- gsub("R_L","R_L_anx", anx_stanvar$scode, fixed = TRUE)
  anx_stanvar$scode <- gsub("N","N_yanx", anx_stanvar$scode, fixed = TRUE)
  anx_stanvar$scode <- gsub("vint1","vint1_yanx", anx_stanvar$scode, fixed = TRUE)
  anx_stanvar$scode <- gsub("vreal1","vreal1_yanx", anx_stanvar$scode, fixed = TRUE)
  anx_stanvar$scode <- gsub("vreal2","vreal2_yanx", anx_stanvar$scode, fixed = TRUE)
  # anx_stanvar$scode <- paste0(anx_stanvar$scode,  "\n", "print(\"first anx cov\");")
  anx_stanvar <- brms::stanvar(scode = anx_stanvar$scode, block = anx_stanvar$block,
                         position = anx_stanvar$position)
  
  dep_stanvar <- depression_reg$stanvar
  dep_stanvar[[3]]$scode <- gsub("N","N_y", dep_stanvar[[3]]$scode, fixed = TRUE)
  # dep_stanvar[[3]]$scode <- paste0(dep_stanvar[[3]]$scode,  "\n", "print(\"first cov\");")
  # dep_stanvar[[3]]$scode <- gsub("vint1","vint1_y", dep_stanvar[[3]]$scode, fixed = TRUE)
  # dep_stanvar[[3]]$scode <- gsub("vreal1","vreal1_y", dep_stanvar[[3]]$scode, fixed = TRUE)
  # dep_stanvar[[3]]$scode <- gsub("vreal2","vreal1_y", dep_stanvar[[3]]$scode, fixed = TRUE)
  
  bias_df <- depression_reg$data %>% 
    full_join(
      anx_reg$data #%>%
        # rename(
        #   y_anx = y,
        #   sd.y_anx = sd.y,
        #   s_baseline.y_anx = s_baseline.y
        # ) %>% 
        # select(y_anx, sd.y_anx, s_baseline.y_anx, study, intervention,
        #        s_groupAge_anx = s_groupAge,s_groupMale_anx = s_groupMale,s_groupBaseline_anx = s_groupBaseline,
        #        s_se_anx = s_se,
        #        .se_anx = .se,
        #        .cov_anx = .cov
        #        contrastSAge_anx = contrastSAge,
        #        
        #        )
        ,
      by = c(
        "study" = "study",
        "intervention" = "intervention"
      ),
      suffix = c("", "_anx")
      ) %>% 
    mutate(
      s_groupBaseline_anx = ifelse(is.nan(s_groupBaseline_anx), NA, s_groupBaseline_anx),
      s_baseline.y_anx = ifelse(is.nan(s_baseline.y_anx), NA, s_baseline.y_anx),
      s_groupAge = ifelse(is.na(s_groupAge), s_groupAge_anx, s_groupAge),
      s_groupMale = ifelse(is.na(s_groupMale), s_groupMale_anx, s_groupMale),
      contrastSAge = ifelse(is.na(contrastSAge), contrastSAge_anx, contrastSAge),
      contrastSMale = ifelse(is.na(contrastSMale), contrastSMale_anx, contrastSMale),
      s_assignedN = ifelse(is.na(s_assignedN), s_assignedN_anx, s_assignedN),
      s_age = ifelse(is.na(s_age), s_age_anx, s_age),
      s_male = ifelse(is.na(s_male), s_male_anx, s_male),
      # .se_anx = ifelse(is.na(.se_anx), 0, .se_anx),
      # .cov_anx = ifelse(is.na(.cov_anx), 0, .cov_anx),
      # .study_anx = as.integer(.study_anx),
      # .study_anx = ifelse(is.na(.study_anx), -1L, .study_anx),
      # # y_anx = ifelse(is.na(y_anx), 0, y_anx),
      # s_se_anx = ifelse(is.na(s_se_anx), 0, s_se_anx),
      # .idx = 1:n(),
      .idx_anx = !is.na(y_anx),
      .idx_dep = !is.na(y)
    )
  
  bias_reg_prep <- list(
    data = bias_df,
    family = "student-t",
    combos = depression_reg$combos %>% 
      full_join(anx_reg$combos,
                suffix = c("","_anx"),
                by =c(".trt", "def","bias.overall"), 
                relationship = "many-to-many") %>% distinct(),
    formula = dep_reg_main_form + anx_reg_main_form + set_rescor(FALSE),
    priors = rbind(depression_reg$prior %>% 
                     filter(class != "sd" &
                              !(resp %in% c("baseliney",
                                          "baselineyanx",
                                          "contrastSAge",
                                          "sgroupBaseline",
                                          "sgroupBaselineanx"))), 
                   anx_priors),
    class_ids = depression_reg$class_ids,
    reference_treatment = depression_reg$reference_treatment,
    stanvar = dep_stanvar + anx_stanvar #+
                # brms::stanvar(scode = "lprior += normal_lpdf(tau_y | 0,0.5); \n  lprior += normal_lpdf(tau_yanx | 0,0.5);", block = "tpar", position = "end")
  )
  names(bias_reg_prep$formula$responses)[length(bias_reg_prep$formula$responses)] <- "yanx"
  
  bias_reg_prep$formula$forms <- map(bias_reg_prep$formula$forms, add_idx_to_bf2)
  bias_reg_prep$formula$forms$contrastSAge <- bias_reg_prep$formula$forms$contrastSAge %>% add_index_term()
  bias_reg_prep$formula$forms$contrastSMale <- bias_reg_prep$formula$forms$contrastSMale %>% add_index_term()
  bias_reg_prep$formula$forms$sbaselineyanx <- bias_reg_prep$formula$forms$sbaselineyanx %>% add_subset_term() %>% add_index_term()
  bias_reg_prep$formula$forms$sgroupBaselineanx <- bias_reg_prep$formula$forms$sgroupBaselineanx %>% add_subset_term() %>% add_index_term()
  bias_reg_prep$formula$forms$sbaseliney <- bias_reg_prep$formula$forms$sbaseliney %>% add_subset_term(".idx_dep") %>%  add_index_term()
  bias_reg_prep$formula$forms$sgroupBaseline <- bias_reg_prep$formula$forms$sgroupBaseline %>% add_subset_term(".idx_dep") %>% add_index_term()
  bias_reg_prep$formula$
  bias_reg_prep$combos$.obs <- 1L:nrow(bias_reg_prep$combos)
  bias_reg_prep$combos$.idx_anx <- TRUE
  bias_reg_prep$combos$.idx_dep <- TRUE
  
  
  cd <- brms::make_stancode(bias_reg_prep$data, formula = bias_reg_prep$formula,
                      prior = bias_reg_prep$priors,
                      stanvar = bias_reg_prep$stanvar) %>% strsplit("\\n") %>% .[[1]]
  dat <- brms::make_standata(formula = bias_reg_prep$formula, bias_reg_prep$data)
  dat$N_y
  dat$N_yanx
  # s <- make_group_cov(dat$vint1_yanx, dat$vreal1_yanx, dat$vreal2_yanx)
  # r <- make_group_cor(dat$vint1_yanx)
  # eigen(r, only.values = TRUE)$values %>% min %>% `>`(0)
  # eigen(s, only.values = TRUE)$values %>% min %>% `>`(0)
  test <- vit$fit_brms_nma(prep = bias_reg_prep,
                           seed = seeds[[outcome]]["metareg"],
                           warmup = 2, 
                           iter = 3, cores = 1, chains = 1,
                           refresh = refresh, silent = silent,
                           threads = threads, backend = "cmdstanr",
                           control = list(max_treedepth = 15, 
                                          adapt_delta = 0.99)
  )
}

bias_reg <- vit$fit_brms_nma(prep = bias_reg_prep,
                 seed = seeds[["depression"]]["metareg"],
                 warmup = 1500L, 
                 iter = 2500L, cores = 4L, chains = 4L,
                 # init = 0.5,
                 refresh = refresh, silent = silent,
                 threads = threads, backend = "cmdstanr",
                 control = list(max_treedepth = 15, 
                                adapt_delta = 0.95)
)

saveRDS(bias_reg, file = here::here("outputs", "saved_models", "bias_regression_fit_overall.rds"))

# debugonce(vit$summary_brms_nma); 
# vit$summary_brms_nma(bias_reg,
#     keep = c(".trt",".idx_anx",".idx_dep","anxiety",
#              "depressed", "def","bias.overall"),
#     index = ".obs", resp = "y") %>% summary() %>%
#   filter(.trt == "D") %>% filter(bias.overall == "low") %>%
#   print(n = 100)
# 
# vit$summary_brms_nma(bias_reg,
#                      keep = c(".trt",".idx_anx","anxiety",
#                               "depressed", "def","bias.overall"),
#                      index = ".obs", resp = "yanx") %>% summary()  %>% filter(.trt == "D") -> anxD
# 
# anxD %>% filter(.trt == "D" & bias.overall == "low" & depressed == 1 ) %>%  arrange(.observed, anxiety, def) %>% print(n = 30) 

bl_test <- brms::loo(bias_reg, resp = "y")
saveRDS(bl_test, file = here::here("outputs", "saved_models", "bias_regression_loo_approx.rds"))

# cat("  Running LOO estimates for meta-regression model...\n")
# # b_l <- bias_reg %>% vit$reloo.vitfit(., loo::loo(.))
# b_l <- vit$loo_estimates(
#   fit = bias_reg,
#   seed = seeds[["depression"]]["loo_metareg"],
#   cores = 8L
# )
# saveRDS(b_l, file = here::here("outputs", "saved_models", "bias_regression_loo_overall.rds"))


# bias_reg <- 
#   vit$prep_brms_nma(data = vit$get_vitamin_data(outcome = outcome, 
#                                                 include_full_bias = TRUE) %>% 
#                       group_by(study) %>%
#                       mutate(groupAge = mean(age, na.rm = TRUE),
#                              groupBaseline = mean(baseline.outcome.mean, na.rm = TRUE),
#                              groupMale = mean(male, na.rm = TRUE),
#                              dropOut = (1 - (final.N / baseline.N))
#                       ) %>% 
#                       ungroup() %>% 
#                       mutate(s_duration = scale(duration),
#                              s_durationE2 = s_duration^2) %>% 
#                       mutate(sd_sd.y = sd(sd.y[!is.na(y)]),
#                              m_sd.y = mean(sd.y[!is.na(y)])) %>%
#                       mutate(s_se = (sd.y - m_sd.y) / sd_sd.y,
#                              s_year = scale(year),
#                              s_groupAge = scale(groupAge),
#                              s_groupBaseline = scale(groupBaseline),
#                              s_groupMale = scale(groupMale),
#                              s_baseline.y = scale(baseline.y),
#                              s_male = scale(male),
#                              s_age = scale(age),
#                              s_dropOut = scale(dropOut),
#                              s_assignedN = scale(assigned.N)
#                       ) 
#                     , family = "student-t"
#                     , confounders = formula_prep_anx$confounders
#                     , interactions = formula_prep_anx$interactions
#                     , main = formula_prep_anx$main 
#                     , addl_form = list(
#                       brms::bf(s_baseline.y | mi() ~ 
#                                  s_assignedN +
#                                  .contrast(mi(s_age)) +
#                                  .contrast(mi(s_male)) + 
#                                  psych.dx + def + 
#                                  (1 | scale) + 
#                                  s_duration + 
#                                  s_se +
#                                  mi(s_groupBaseline) +
#                                  (1 + s_se | country) +
#                                  bias.selective.outcome  +  
#                                  bias.other  + bias.overall +
#                                  bias.blinding.participants + 
#                                  bias.blinding.outcome + 
#                                  bias.incomplete.outcome +
#                                  bias.sequence + bias.concealment) + gaussian(),
#                       brms::bf(.contrast(s_age) | mi() ~  s_year +  
#                                  (1|country) + s_assignedN +
#                                  bias.selective.outcome  +  
#                                  bias.other  + bias.overall +
#                                  bias.blinding.participants + 
#                                  bias.blinding.outcome + 
#                                  bias.incomplete.outcome +
#                                  bias.sequence + bias.concealment) + gaussian(),
#                       brms::bf(.contrast(s_male) | mi() ~  s_year +  
#                                  (1|country)  + s_assignedN +
#                                  bias.selective.outcome  +  
#                                  bias.other  + bias.overall +
#                                  bias.blinding.participants + 
#                                  bias.blinding.outcome + 
#                                  bias.incomplete.outcome +
#                                  bias.sequence + bias.concealment) + gaussian(),
#                       brms::bf(s_groupBaseline | mi() ~  
#                                  mi(s_groupAge) + 
#                                  mi(s_groupMale) +
#                                  psych.dx + def + 
#                                  (1 | scale) + (1 + s_se|country) +
#                                  s_duration +
#                                  bias.selective.outcome  +  
#                                  bias.other  + bias.overall +
#                                  bias.blinding.participants + 
#                                  bias.blinding.outcome + 
#                                  bias.incomplete.outcome +
#                                  bias.sequence + bias.concealment) + gaussian(),
#                       brms::bf( s_groupAge | mi() ~  s_year + 
#                                   psych.dx + def + 
#                                   (1|country) ) + gaussian(),
#                       brms::bf( s_groupMale | mi() ~  s_year + 
#                                   psych.dx + def +
#                                   (1|country) ) + gaussian()
#                     )
#   ) %>%
#   vit$fit_brms_nma(
#     seed = seeds[[outcome]]["metareg"],
#     warmup = 2, 
#     iter = 3,
#     refresh = refresh, silent = silent,
#     threads = threads, backend = "cmdstanr",
#     control = list(max_treedepth = 15, 
#                    adapt_delta = 0.99)
#   )

# np <- nuts_params(bias_reg)
# draws <- as_draws_df(bias_reg)
# 
# # diagnostics wide: one row per chain/iter, columns divergent__/treedepth__/...
# np_wide <- np %>%
#   tidyr::pivot_wider(names_from = Parameter, values_from = Value) %>%
#   arrange(Chain, Iteration)
# 
# # sanity check: number of draws matches
# stopifnot(nrow(np_wide) == nrow(draws))
# 
# div_flag <- np_wide$divergent__ == 1
# td       <- np_wide$treedepth__
# max_td   <- max(td, na.rm = TRUE)
# td_flag  <- td >= max_td
# 
# param_cols <- grep("^(b_|sd_|cor_|r_|sigma$|nu$|Intercept)", names(draws), value = TRUE)
# 
# delta_div <- sapply(param_cols, function(p) {
#   x <- draws[[p]]
#   abs(mean(x[div_flag], na.rm = TRUE) - mean(x[!div_flag], na.rm = TRUE))
# })
# 
# head(sort(delta_div, decreasing = TRUE), 30)
# 
# delta_td <- sapply(param_cols, function(p) {
#   x <- draws[[p]]
#   abs(mean(x[td_flag], na.rm = TRUE) - mean(x[!td_flag], na.rm = TRUE))
# })
# 
# head(sort(delta_td, decreasing = TRUE), 30)
# 
# library(bayesplot)
# 
# top <- names(sort(delta_div, decreasing = TRUE))[1:12]
# 
# if(any(!is.na(top))) {
#   mcmc_pairs(
#     draws,
#     pars = colnames(draws)[colnames(draws) %in% top[1:5]],
#     np = np
#   )
# }
# 
# div <- np_wide$divergent__ == 1
# sus <- grep("^(sd_|cor_|tau)", names(draws), value = TRUE)
# 
# summ <- tibble(
#   par = sus,
#   mean_div = sapply(sus, \(p) mean(draws[[p]][div], na.rm=TRUE)),
#   mean_ok  = sapply(sus, \(p) mean(draws[[p]][!div], na.rm=TRUE)),
#   delta    = abs(mean_div - mean_ok)
# ) %>% arrange(desc(delta))
# 
# head(summ, 30)
# 
# aggregate(Value ~ Chain, data = subset(np, Parameter=="divergent__"), sum)
# aggregate(Value ~ Chain, data = subset(np, Parameter=="stepsize__"), median)
# aggregate(Value ~ Chain, data = subset(np, Parameter=="treedepth__"),
#           function(x) mean(x >= 15))
