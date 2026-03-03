# this code will generate all of the tables needed for the 
# main manuscript and supplementary materials

#### Load Packages ####
suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
  library(gt)
  library(stringr)
  library(tidyr)
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

outcomes <- c("depression","anxiety")

#### Load Data ####
# full data
data <- vit$get_vitamin_data(outcome = NULL, simple_analysis = FALSE,
                             include_full_bias = TRUE)

#### temp functions ####
gt_fmt <- function(x, ...) {
  mask_all_missing <- with(x, is.na(Estimate) & is.na(`q2.5%`) & is.na(`q97.5%`))
  
  x %>% 
    mutate(
      across(
        any_of("depressed"),
        ~ recode(factor(.x), "0" = "No", "1" = "Yes")
      )
    ) %>% 
    mutate(
      across(
        any_of("anxiety"),
        ~ recode(factor(.x), "0" = "No", "1" = "Yes")
      )
    ) %>% 
    mutate(
      across(
        any_of("depressed"),
        ~ recode(factor(.x), "no" = "No", "yes" = "Yes")
      )
    ) %>% 
    mutate(
      across(
        any_of("anxiety"),
        ~ recode(factor(.x), "no" = "No", "yes" = "Yes")
      )
    ) %>% 
    mutate(
      across(
        any_of("def"),
             ~ recode(factor(.x),
                    "NA" = "Unknown",
                    "no" = "No",
                    "yes" = "Yes"
             ))
      ) %>%
    rename(Intervention = any_of("intervention")) %>%
    rename(Intervention = any_of(".trt")) %>% 
    rename(Dose = any_of(".dose")) %>% 
    rename(Depressed = any_of("depressed")) %>% 
    rename(Anxious = any_of("anxiety")) %>% 
    rename(Bias = any_of("bias")) %>%
    rename(Bias = any_of("bias.overall")) %>%
    rename(Deficient = any_of("def")) %>% 
    mutate(
      across(
        any_of("Bias"),
        ~ forcats::fct_recode(.x, 
                              "low" = "low",
                              "unclear" = "some concerns",
                              "high" = "high")
      )
    ) %>%
    mutate(
      across(
        any_of("Bias"),
        ~ forcats::fct_relevel(.x, "low", "unclear", "high")
      )
    ) %>% 
    arrange(across(any_of(c("Intervention", "Bias","Dose", "Depressed","Anxious",  "Deficient")))) %>% 
    select(any_of(c("Intervention", "Bias","Dose", "Depressed","Anxious",  "Deficient")), 
           everything()) %>% 
    select(-c(S.E., `q50%`)) %>%
    gt(...) %>% 
    fmt_number(
      columns = c(Estimate, `q2.5%`, `q97.5%`),
      decimals = 3
    ) %>%
    sub_missing() %>% 
    cols_merge(
      columns = c(Estimate, `q2.5%`, `q97.5%`),
      pattern = "{1} ({2}, {3})"
    ) %>%
    text_transform(
      locations = cells_body(columns = Estimate, rows = mask_all_missing),
      fn = function(x) rep("—", length(x))
    ) %>%
    cols_label(
      Intervention = "Intervention",
      Estimate = "Effect Size (95% Credible Interval)"
    ) 
}


#### Baseline Demographic Data ####
# setup summary functions
{
  count_vars <- function(x, total) {
    r <- as.integer(round(x))
    pct_calc <- r/total
    pct_calc <- ifelse(r == 0, 0, pct_calc)
    pct <- sprintf('%.2f', round(pct_calc*100,digits = 2))
    glue::glue("{format(r,digits = 2,big.mark = ',')} ({format(pct, digits = 2)}%)")
  }
  cont_vars <- function(x) {
    glue::glue("{sprintf('%.2f', round(mean(x, na.rm = TRUE), digits = 2))} ({sprintf('%.2f', sd(x, na.rm = TRUE))})")
  }
  wt_cont_vars <- function(x, sds = NULL, weights) {
    m <- format(weighted.mean(x, w = weights, na.rm = TRUE), digits = 2,big.mark = ",")
    # mean_text <- "{sprintf('%.2f', round(m,digits = 2))}"
    mean_text <- "{m}"
    if( !is.null(sds) ) {
      w_var <- if ( length(x) > 1 ) {
        matrixStats::weightedVar(x = x, w = weights, na.rm = TRUE)
      } else {
        0
      }
      sd_w <- format(sqrt( w_var +
                             weighted.mean(sds^2, weights, na.rm = TRUE) 
      ), digits = 2,big.mark = ",")
      # sd_text <- " ({sprintf('%.2f', round(sd_w,digits = 2))})"
      sd_text <- " ({sd_w})"
    } else {
      sd_text <- ""
    }
    
    glue::glue(paste0(mean_text, sd_text))
  }
  clean_vit <- function(x) {
    sub("^(d\\.|min_d\\.|max_d\\.)", "", x)
  }
  clean_dose <- function(x) {
    sub("^(dose_d\\.|min_d\\.|max_d\\.)", "", x)
  }
  demo_create <- function(x) {
    x %>% 
    distinct(study, intervention, .keep_all = TRUE) %>%
      group_by(intervention) %>% 
      mutate(b.weights = baseline.N,
             f.weights = final.N,
             temp.sd = 0) %>% 
      # mutate(baseline.SMD = baseline.outcome.mean / sd(baseline.outcome.mean, na.rm = TRUE)) %>%
      summarize(
        n.study = n(),
        b.n = sum(baseline.N),
        n = format(sum(final.N), digits = 2, big.mark = ","), #count_vars(sum(final.N), ntotal), 
        age = wt_cont_vars(age, NULL, b.weights),
        male = count_vars(sum(male * b.weights, na.rm = TRUE), sum(b.weights)),
        depression = count_vars(sum(depressed), n()),
        anxiety = count_vars(sum(anxiety), n()),
        deficiency = count_vars(sum(def == "yes"), n()),
        # depression = paste0(wt_cont_vars(depressed*100, NULL, b.weights),"%"),
        # anxiety = paste0(wt_cont_vars(anxiety*100, NULL, b.weights),"%"),
        # depression = paste0(round(mean(depressed*100, na.rm = TRUE),digits = 2),"%"),
        # anxiety = paste0(round(mean(anxiety*100, na.rm=TRUE),digits = 2),"%"),
        # bias.low = glue::glue("{(mean(bias.overall == 'low', na.rm = TRUE) * 100) %>% round(digits = 2)}", "%"),
        # bias.unclear = glue::glue("{(mean(bias.overall == 'some concerns', na.rm = TRUE) * 100) %>% round(digits = 2)}", "%"),
        # bias.high = glue::glue("{(mean(bias.overall == 'high', na.rm = TRUE) * 100) %>% round(digits = 2)}", "%"),
        bias.low = count_vars(sum(bias.overall == 'low', na.rm = TRUE), n()),
        bias.unclear = count_vars(sum(bias.overall == 'some concerns', na.rm = TRUE),n()),
        bias.high = count_vars(sum(bias.overall == 'high', na.rm = TRUE), n()),
        # duration = wt_cont_vars(duration,temp.sd,f.weights),
        duration = wt_cont_vars(duration,NULL,rep(1,n())),
        
        across(starts_with("d."), 
               ~ 
                 if(n() > 1) {
                   if(any(.x != 0)) {
                     glue::glue(format(vit$convert_g_to_dose(clean_vit(cur_column()), weighted.mean(.x, b.weights)), digits = 2),
                                " [",
                                format(vit$convert_g_to_dose(clean_vit(cur_column()), min(.x, na.rm = TRUE)), digits = 2),
                                ", ",
                                format(vit$convert_g_to_dose(clean_vit(cur_column()), max(.x, na.rm = TRUE)), digits = 2),
                                "]"
                     )
                   } else {
                     NA_character_
                   }
                   
                 } else {
                   if(all(.x == 0)) {
                     NA_character_
                   } else {
                     glue::glue(format(vit$convert_g_to_dose(clean_vit(cur_column()), weighted.mean(.x, b.weights)), digits = 2)
                     )
                   }
                 }
               ,
               .names = "dose_{col}"),
        
        baseline.y = wt_cont_vars(baseline.y, baseline.sd.y, b.weights)
      ) %>% 
      select(
        -contains("aerobic_exercise"),
        -contains("antidepressant")
      ) %>% 
      ungroup() %>% 
      mutate(dose = 
               purrr::pmap_chr(
                 select(., starts_with("dose_d.") & !matches("_(min|max)_")), 
                 function(...) {
                   doses <- c(...)
                   dose_strs <- c()
                   # browser()
                   for(i in seq_along(doses)) {
                     dose_val <- doses[[i]]
                     if(!is.na(dose_val) && dose_val != "NA") {
                       col_name <- names(doses)[i]
                       vit_name <- clean_dose(col_name)
                       dose_strs <- c(dose_strs, paste0(vit_name, ": ", dose_val))
                     }
                   }
                   if(!is.null(dose_strs)) {
                     paste(dose_strs, collapse = ", ")
                   } else {
                     NA_character_
                   }
                 }
               )
      ) %>% 
      select(
        -starts_with("dose_d.")
      )
  }
  
  study_table <- function(x) {
    
    is_gt <- inherits(x, "gt_tbl")
    
    bias_pres <- switch(is_gt %>% as.character(),
                        "FALSE" = ("bias.low" %in% colnames(x)),
                        "TRUE" = ("bias.low" %in% colnames(x$`_data`))
                        )
    
    if ( !is_gt ) {
      x <- x %>% 
      mutate(intervention = stringr::str_to_title(as.character(intervention))) %>% 
      arrange(intervention) %>% 
      gt()
    }
    
    x <- x  %>% 
      # tab_header(
      #   title = "Study-Level Demographics (Low Bias Studies)"
      # ) %>% 
      gt::cols_label(
        intervention = "Intervention",
        n.study = "Number of Studies",
        depression = "Num. of Studies in Depressed Pop.",
        anxiety = "Num. of Studies in Anxious Pop.",
        deficiency = "Num. of Studies in Deficient Pop.",
        duration = "Duration in weeks",
        dose = "Dose*"
      ) %>% 
      tab_caption(
        "Study-Level Demographics. Values are n (%) or mean (SD) where appropriate. *Dose is reported as mean [min, max] across included studies."
      ) %>%
      sub_missing(
        columns = everything(),
        missing_text = ""
      )
    
    if(bias_pres) {
      x <- x %>% 
        tab_spanner(
          label = "Overall Study Bias (%)",
          columns = c(bias.low, bias.unclear, bias.high)
        ) %>% 
        gt::cols_label(
          bias.low = "Low",
          bias.unclear = "Unclear",
          bias.high = "High"
        ) 
    }
    if(is_gt) {
      x <- x %>% 
        overall_labeling()
    }
    return(x )
  }
  
  pop_table <- function(x) {
    is_gt <- inherits(x, "gt_tbl")
    
    if(!is_gt ){
      x <- x %>% 
        mutate(intervention = stringr::str_to_title(as.character(intervention))) %>% 
        arrange(intervention) %>% 
        gt::gt()
    }
    x <- x  %>% 
      gt::cols_label(
        intervention = "Intervention",
        n = "N",
        age = "Age",
        male = "Male"
        # ,baseline.y = "SMD at baseline"
      ) %>% 
      tab_caption(
        "Population-Level Demographics. Values are n (%) or mean (SD) where appropriate."
      ) %>%
      sub_missing(
        columns = everything(),
        missing_text = ""
      )
    
    if(is_gt) {
      x <- x %>% 
        overall_labeling()
    }
    
    return(x)
  }
  
  overall_labeling <- function(x) {
    bias_pres <- ("bias.low" %in% colnames(x$`_data`))
    study_cols <- if (bias_pres) {
      c("n.study", "depression", "anxiety", "deficiency", "duration", "dose",
        "bias.low", "bias.unclear", "bias.high")
    } else {
      c("n.study", "depression", "anxiety",  "deficiency", "duration", "dose")
    }
    
    bias_comment <- if(!bias_pres) {
      " in low bias studies"
    } else {
      ""
    }
    
    x %>% 
      tab_spanner(
        columns = c(
                    n,
                    age,
                    male
                    # , baseline.y
                    ),
        label = "Population-Level Characteristics"
      ) %>% 
      tab_spanner(label = "Study-Level Characteristics",
                    columns =  all_of(study_cols)) %>% 
      tab_caption(
        glue::glue("Baseline demographics{bias_comment}. Values are n (%) or mean (SD) where appropriate. *Dose is reported as mean [min, max] across included studies.")
      )
  }
}

# get counts
ntrt <- data$intervention %>% unique() %>% length()
nstudy <- data$study %>% unique() %>% length()
ntotal <- sum(data %>% distinct(study,
                                intervention, .keep_all = TRUE) %>% .$final.N)

# select baseline columns
baseline.cols <- rlang::quos(
  "final.N",
  "duration",
  "age",
  "male",
  "depressed",
  "anxiety",
  "def",
  matches("baseline\\.(?!outcome)", perl = TRUE),
  starts_with("d."),
  "bias.overall"
)

baseline.vars <- data %>% select(c("study","intervention", !!!baseline.cols))

# create demo tables
{
  # demographics tab overall
  demo_tab <- baseline.vars %>% 
    demo_create()
  
  # demographics tab low bias
  demo_tab_low <- baseline.vars %>% 
    filter(bias.overall == "low") %>% 
    demo_create()
  
  study_demo <- demo_tab %>% 
    select(
      intervention,
      n.study,
      depression,
      anxiety,
      deficiency,
      duration,
      dose,
      bias.low,
      bias.unclear,
      bias.high
    )
  
  study_demo_low <- demo_tab_low %>% 
    select(
      intervention,
      n.study,
      depression,
      anxiety,
      deficiency,
      duration,
      dose
    )
  
  
  pop_demo <- demo_tab %>% select(
    intervention,
    n,
    age,
    male
  )
    # , baseline.y
  # ) %>% 
  #   mutate(baseline.y = 
  #            ifelse(grepl("[0-9 ]*\\(NA\\)$", baseline.y, perl =TRUE) |
  #                     grepl("^NaN", baseline.y, perl = TRUE), 
  #                   "Ref.", baseline.y))
  
  pop_demo_low <- demo_tab_low %>% select(
    intervention,
    n,
    age,
    male
    # , baseline.y
  )
  # ) %>% 
  #   mutate(baseline.y = 
  #            ifelse(grepl("[0-9 ]*\\(NA\\)$", baseline.y, perl =TRUE) |
  #                     grepl("^NaN", baseline.y, perl = TRUE), 
  #                   "Ref.", baseline.y))
  
  dose_ranges <-   baseline.vars %>%
    distinct(study, intervention, .keep_all = TRUE) %>%
    mutate(b.weights = baseline.N,
           f.weights = final.N) %>% 
    mutate(across(starts_with("d."), 
                  ~ as.numeric(ifelse(.x == 0, NA_real_, .x))
    )
    ) %>%
    summarize(
      across(starts_with("d."), 
             ~ vit$convert_g_to_dose(clean_vit(cur_column()), min(.x, na.rm = TRUE)),
             .names = "min_{col}"
      ),
      across(starts_with("d."), 
             ~ vit$convert_g_to_dose(clean_vit(cur_column()), max(.x, na.rm = TRUE)),
             .names = "max_{col}"
      ),
      across(starts_with("d."), 
             ~ vit$convert_g_to_dose(clean_vit(cur_column()), weighted.mean(.x, b.weights, na.rm = TRUE)),
             .names = "mean_{col}"
      )
    )
  
  dose_ranges_low <-   baseline.vars %>%
    filter(bias.overall == "low") %>% 
    distinct(study, intervention, .keep_all = TRUE) %>%
    mutate(b.weights = baseline.N,
           f.weights = final.N) %>% 
    mutate(across(starts_with("d."), 
                  ~ as.numeric(ifelse(.x == 0, NA_real_, .x))
    )
    ) %>%
    summarize(
      across(starts_with("d."), 
             ~ vit$convert_g_to_dose(clean_vit(cur_column()), min(.x, na.rm = TRUE)),
             .names = "min_{col}"
      ),
      across(starts_with("d."), 
             ~ vit$convert_g_to_dose(clean_vit(cur_column()), max(.x, na.rm = TRUE)),
             .names = "max_{col}"
      ),
      across(starts_with("d."), 
             ~ vit$convert_g_to_dose(clean_vit(cur_column()), weighted.mean(.x, b.weights, na.rm = TRUE)),
             .names = "mean_{col}"
      )
    )
  
  if (!dir.exists(here::here("outputs","overall","tables","raw"))) {
    dir.create(here::here("outputs","overall","tables","raw"), recursive = TRUE)
  }
  saveRDS(dose_ranges,
          here::here("outputs","overall","tables","raw","dose_ranges.rds"))
  saveRDS(pop_demo,
          here::here("outputs","overall","tables","raw","pop_demo.rds"))
  saveRDS(study_demo,
          here::here("outputs","overall","tables","raw","study_demo.rds"))
  
  
  saveRDS(dose_ranges_low,
          here::here("outputs","overall","tables","raw","dose_ranges_low_bias.rds"))
  saveRDS(pop_demo_low,
          here::here("outputs","overall","tables","raw","pop_demo_low_bias.rds"))
  saveRDS(study_demo_low,
          here::here("outputs","overall","tables","raw","study_demo_low_bias.rds"))
  
  # create printable tables
  study_demo_low %>% 
    study_table() %>% 
    saveRDS(
      here::here("outputs","overall","tables",
                 "study_demo_low_bias_table.rds")
    )
    
  study_demo %>% 
    study_table() %>% 
    saveRDS(
      here::here("outputs","overall","tables",
                 "study_demo_table.rds")
    )
  
  pop_low_gt <- pop_demo_low %>% 
    pop_table()
  
  pop_gt <- pop_demo %>% 
    pop_table()
  
  pop_low_gt %>% 
    saveRDS(
      here::here("outputs","overall","tables",
                 "pop_demo_low_bias_table.rds")
    )
  
  pop_gt %>% 
    saveRDS(
      here::here("outputs","overall","tables",
                 "pop_demo_table.rds")
    )
  
  combine_tab <- dplyr::left_join(
    pop_demo,
    study_demo,
    by = "intervention") %>% 
    pop_table() %>% 
    study_table()
  
  combine_tab_low <- dplyr::left_join(
    pop_demo_low,
    study_demo_low,
    by = "intervention") %>% 
    pop_table() %>% 
    study_table()
  
  combine_tab %>% 
    saveRDS(
      here::here("outputs","overall","tables",
                 "combined_demo_table.rds")
    )
  
  combine_tab_low %>% 
    saveRDS(
      here::here("outputs","overall","tables",
                 "combined_demo_low_bias_table.rds")
    )
  
}

# create table by comparison to assess for overlap
{
  baseline.vars.comp <- data %>% 
    mutate(control = gsub("z\\.","",as.character(control)) %>% factor()) %>% 
    mutate(control = gsub("_"," ",as.character(control)) %>% factor()) %>% 
     mutate(control = as.character(control) %>% stringr::str_to_title()) %>% 
    mutate(intervention = as.character(intervention) %>% stringr::str_to_title() ) %>% 
    mutate(intervention = glue::glue("{intervention}:{control}")) %>% 
    mutate(intervention = gsub("z\\.","",as.character(intervention)) %>% factor()) %>% 
    mutate(intervention = gsub("_"," ",as.character(intervention)) %>% factor()) %>%
    mutate(intervention = gsub(" - ",":",as.character(intervention)) %>% factor()) %>%
    filter(!is.na(contrast)) %>% 
    select(c("study","intervention", !!!baseline.cols))
  
  # demographics tab overall
  demo_tab.comp <- baseline.vars.comp %>% 
    demo_create()
  
  # demographics tab low bias
  demo_tab_low.comp <- baseline.vars.comp %>% 
    filter(bias.overall == "low") %>% 
    demo_create()
  
  study_demo.comp <- demo_tab.comp %>% 
    select(
      intervention,
      n.study,
      depression,
      anxiety,
      deficiency,
      duration,
      dose,
      bias.low,
      bias.unclear,
      bias.high
    )
  
  study_demo_low.comp <- demo_tab_low.comp %>% 
    select(
      intervention,
      n.study,
      depression,
      anxiety,
      deficiency,
      duration,
      dose
    )
  
  
  pop_demo.comp <- demo_tab.comp %>% select(
    intervention,
    n,
    age,
    male
  )
  
  pop_demo_low.comp <- demo_tab_low.comp %>% 
  select(
    intervention,
    n,
    age,
    male
    # , baseline.y
  )
  
  dose_ranges.comp <-   baseline.vars.comp %>%
    distinct(study, intervention, .keep_all = TRUE) %>%
    mutate(b.weights = baseline.N,
           f.weights = final.N) %>% 
    mutate(across(starts_with("d."), 
                  ~ as.numeric(ifelse(.x == 0, NA_real_, .x))
    )
    ) %>%
    summarize(
      across(starts_with("d."), 
             ~ vit$convert_g_to_dose(clean_vit(cur_column()), min(.x, na.rm = TRUE)),
             .names = "min_{col}"
      ),
      across(starts_with("d."), 
             ~ vit$convert_g_to_dose(clean_vit(cur_column()), max(.x, na.rm = TRUE)),
             .names = "max_{col}"
      ),
      across(starts_with("d."), 
             ~ vit$convert_g_to_dose(clean_vit(cur_column()), weighted.mean(.x, b.weights, na.rm = TRUE)),
             .names = "mean_{col}"
      )
    )
  
  dose_ranges_low.comp <-   baseline.vars.comp %>%
    filter(bias.overall == "low") %>% 
    distinct(study, intervention, .keep_all = TRUE) %>%
    mutate(b.weights = baseline.N,
           f.weights = final.N) %>% 
    mutate(across(starts_with("d."), 
                  ~ as.numeric(ifelse(.x == 0, NA_real_, .x))
    )
    ) %>%
    summarize(
      across(starts_with("d."), 
             ~ vit$convert_g_to_dose(clean_vit(cur_column()), min(.x, na.rm = TRUE)),
             .names = "min_{col}"
      ),
      across(starts_with("d."), 
             ~ vit$convert_g_to_dose(clean_vit(cur_column()), max(.x, na.rm = TRUE)),
             .names = "max_{col}"
      ),
      across(starts_with("d."), 
             ~ vit$convert_g_to_dose(clean_vit(cur_column()), weighted.mean(.x, b.weights, na.rm = TRUE)),
             .names = "mean_{col}"
      )
    )
  
  if (!dir.exists(here::here("outputs","overall","tables","raw"))) {
    dir.create(here::here("outputs","overall","tables","raw"), recursive = TRUE)
  }
  saveRDS(dose_ranges.comp,
          here::here("outputs","overall","tables","raw","dose_ranges_comparison.rds"))
  saveRDS(pop_demo.comp,
          here::here("outputs","overall","tables","raw","pop_demo_comparison.rds"))
  saveRDS(study_demo.comp,
          here::here("outputs","overall","tables","raw","study_demo_comparison.rds"))
  
  
  saveRDS(dose_ranges_low.comp,
          here::here("outputs","overall","tables","raw","dose_ranges_low_bias_comparison.rds"))
  saveRDS(pop_demo_low.comp,
          here::here("outputs","overall","tables","raw","pop_demo_low_bias_comparison.rds"))
  saveRDS(study_demo_low.comp,
          here::here("outputs","overall","tables","raw","study_demo_low_bias_comparison.rds"))
  
  # create printable tables
  study_demo_low.comp %>% 
    study_table() %>% 
    gt::cols_label(intervention = "Comparison") %>% 
    saveRDS(
      here::here("outputs","overall","tables",
                 "study_demo_low_bias_table_comparison.rds")
    )
  
  study_demo.comp %>% 
    study_table() %>% 
    gt::cols_label(intervention = "Comparison") %>% 
    saveRDS(
      here::here("outputs","overall","tables",
                 "study_demo_table_comparison.rds")
    )
  
  pop_low_gt.comp <- pop_demo_low.comp %>% 
    pop_table() %>% 
    gt::cols_label(intervention = "Comparison") 
  
  pop_gt.comp <- pop_demo.comp %>% 
    pop_table() %>% 
    gt::cols_label(intervention = "Comparison") 
  
  pop_low_gt.comp %>% 
    saveRDS(
      here::here("outputs","overall","tables",
                 "pop_demo_low_bias_table_comparison.rds")
    )
  
  pop_gt.comp %>% 
    saveRDS(
      here::here("outputs","overall","tables",
                 "pop_demo_table_comparison.rds")
    )
  
  combine_tab.comp <- dplyr::left_join(
    pop_demo.comp,
    study_demo.comp,
    by = "intervention") %>% 
    pop_table() %>% 
    study_table() %>% 
    gt::cols_label(intervention = "Comparison") 
  
  combine_tab_low.comp <- dplyr::left_join(
    pop_demo_low.comp,
    study_demo_low.comp,
    by = "intervention") %>% 
    pop_table() %>% 
    study_table() %>% 
    gt::cols_label(intervention = "Comparison") 
  
  combine_tab.comp %>% 
    saveRDS(
      here::here("outputs","overall","tables",
                 "combined_demo_table_comparison.rds")
    )
  
  combine_tab_low.comp %>% 
    saveRDS(
      here::here("outputs","overall","tables",
                 "combined_demo_low_bias_table_comparison.rds")
    )
}

# create table by study for appendix
{
  data %>% 
    distinct(study, intervention,
             final.N, duration, depressed,anxiety,
             vit.def, bias.overall,
             .keep_all = TRUE) %>% 
    mutate(depressed = ifelse(depressed == 1, "Yes", "No"),
           anxiety = ifelse(anxiety == 1, "Yes", "No")) %>%
    select(Study = study, Intervention = intervention,
           N = final.N, `Duration in weeks` = duration, 
           Depressed = depressed,
           Anxious = anxiety,
           Deficiency = vit.def, Bias = bias.overall) %>% 
    arrange(tolower(Study)) %>% 
    gt::gt(groupname_col = "Study") %>% 
    gt::fmt_number(decimals = 0) %>% 
    gt::sub_missing() %>% 
    saveRDS(
      here::here("outputs","overall","tables",
                 "study_level_demo_table.rds")
    )
}

{
  vit$get_vitamin_data(outcome = NULL, simple_analysis = TRUE,
                       include_full_bias = TRUE) %>% 
    distinct(study, intervention,
             final.N, duration, depressed,anxiety,
             vit.def, bias.overall,
             .keep_all = TRUE) %>% 
    filter(bias.overall == "low") %>%
    mutate(depressed = ifelse(depressed == 1, "Yes", "No"),
           anxiety = ifelse(anxiety == 1, "Yes", "No")) %>%
    select(Study = study, Intervention = intervention,
           N = final.N, `Duration in weeks` = duration, 
           Depressed = depressed,
           Anxious = anxiety,
           Deficiency = vit.def, Bias = bias.overall) %>% 
    arrange(tolower(Study)) %>% 
    gt::gt(groupname_col = "Study") %>% 
    gt::fmt_number(decimals = 0) %>% 
    gt::sub_missing() %>% 
    saveRDS(
      here::here("outputs","overall","tables",
                 "study_level_demo_table_low.rds")
    )
}

#### Load models and tables for regressions ####
for (outcome in outcomes) {
  sensitivity_outdir <- here::here("outputs", "saved_models", outcome,
                                   "sensitivity")
  
  output_dir <- here::here("outputs",outcome,"table")
  outcome_label <- switch(outcome,
                          "depression" = "depressed",
                          "anxiety" = "anxiety")
  
  not_outcome_label <- switch(outcome,
                              "anxiety" = "depressed",
                              "depression" = "anxiety")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  # load models
  {
    low_t <- readRDS(
      here::here("outputs", "saved_models", outcome, 
                 glue::glue("model_", 
                            outcome, 
                            "_low_bias_data.rds"))
    ) 
    
    some_t <- readRDS(
      here::here("outputs", "saved_models", outcome,
                 glue::glue("model_",
                            outcome,
                            "_some_and_low_bias_data.rds"))
    )
    
    high_t <- readRDS(
      here::here("outputs", "saved_models", outcome,
                 glue::glue("model_",
                            outcome,
                            "_all_data.rds"))
    )
    
    # some_t <- readRDS(
    #   here::here("outputs", "saved_models", outcome,
    #              glue::glue("model_", 
    #                         outcome, 
    #                         "_some_bias_data.rds"))
    # ) 
    # 
    # high_t <- readRDS(
    #   here::here("outputs", "saved_models", outcome,
    #              glue::glue("model_", 
    #                         outcome, 
    #                         "_high_bias_data.rds"))
    # ) 
    
    
    fit_def_low <- readRDS(
      here::here("outputs", "saved_models", outcome,
                 glue::glue("model_", 
                            outcome, 
                            "_low_bias_data_by_deficiency_depression.rds"))
    )
    
    fit_def <- readRDS(
      here::here("outputs", "saved_models", outcome,
                 glue::glue("model_", 
                            outcome, 
                            "_all_data_by_deficiency_depression_and_bias.rds"))
    )
    
    fit_def_low_no_iran <- readRDS(
      here::here("outputs", "saved_models", outcome,
                 "sensitivity",
                 glue::glue("model_",
                                       outcome,
                                       "_low_bias_data_by_deficiency_depression_no_iran.rds"))
    )
    
    fit_def_no_iran <- readRDS(
      here::here("outputs", "saved_models", outcome,
                 "sensitivity",
                 glue::glue("model_",
                            outcome,
                            "_all_data_by_deficiency_depression_and_bias_no_iran.rds"))
    )
    
    fit_dose <- readRDS(
      here::here("outputs", "saved_models", outcome,
                 glue::glue("model_", 
                            outcome, 
                            "_low_bias_data_by_dose.rds"))
    )
    
    # fit_dose_full <- readRDS(
    #   here::here("outputs", "saved_models", outcome,
    #              glue::glue("model_", 
    #                         outcome, 
    #                         "_all_data_by_dose.rds"))
    # )
    
    fit_time <- readRDS(
      here::here("outputs", "saved_models", outcome, "sensitivity",
                 glue::glue("model_", 
                            outcome, 
                            "_low_bias_data_by_time_cat.rds"))
    )
    
    bias_reg <- readRDS(here::here("outputs", 
                                  "saved_models", 
                                  "bias_regression_fit_overall.rds"))
    
  }
  
  # summary tibbles
  {
    low_res <- vit$summary_brms_nma(low_t) %>% 
      summary() %>% 
      filter(.observed == TRUE) %>% 
      select(-.observed)
    
    some_res <- vit$summary_brms_nma(some_t) %>% 
      summary() %>% 
      filter(.observed == TRUE) %>% 
      select(-.observed)
    
    high_res <- vit$summary_brms_nma(high_t) %>% 
      summary() %>% 
      filter(.observed == TRUE) %>% 
      select(-.observed)
    
    dose_res <- vit$summary_brms_nma(fit_dose,
                                     keep = c(".trt",".dose")) %>% 
      summary() %>% 
      filter(.observed == TRUE)  %>% 
      select(-.observed) %>% 
      mutate(.dose = as.factor(.dose),
             .dose = forcats::fct_relevel(.dose,
                                         "B1 (> 100 mg)",
                                         "B12 (250 - 1000 µg), B6 (10 - 50 mg), B9 (0.5 - 5 mg)",
                                         "B12 (< 250 µg), B9 (< 0.5 mg)",
                                         "B12 (250 - 1000 µg), B9 (< 0.5 mg)",
                                         "B12 (250 - 1000 µg), B9 (< 0.5 mg), D (< 1000 IU)",
                                         "B6 (10 - 50 mg)",
                                         "B9 (< 0.5 mg)",
                                         "B9 (0.5 - 5 mg)",
                                         "C (200 - 1000 mg)",
                                         "D (< 1000 IU)",
                                         "D (1000 - 3000 IU)",
                                         "D (> 3000 IU)",
                                         "D (1000 - 3000 IU), iron (18 - 45 mg)",
                                         "magnesium (200 - 400 mg)",
                                         "magnesium (> 400 mg)",
                                         "magnesium (< 200 mg), zinc (> 30 mg)",
                                         "selenium (100 - 200 µg)",
                                         "zinc (> 30 mg)"
                                         ))
    
    
    time_res <- fit_time %>% 
      vit$summary_brms_nma(keep = c(".trt","time_cat")) %>% 
      summary() %>% 
      filter(.observed == TRUE)  %>% 
      select(-.observed) %>% 
      mutate(time_cat = factor(time_cat,
                                levels = c("0-8 weeks", "8-12 weeks", "12-26 weeks", "26+ weeks"),
                                labels = c("0-8 weeks", "8-12 weeks", "12-26 weeks", "> 26 weeks"), ordered = TRUE)) %>% 
      arrange(.trt, time_cat) %>% 
      rename(Duration = time_cat)
    
    def_res <- vit$summary_brms_nma(fit_def,
                                    keep = c(".trt",outcome_label, "def","bias")) %>% 
      summary() %>% 
      filter(.observed == TRUE)  %>% 
      select(-.observed) 
    
    def_total_res <- fit_def %>% vit$marginalize(block = c("bias"), margins = c(outcome_label, "def")) %>% 
      summary() 
    
    def_res_low <- vit$summary_brms_nma(fit_def_low,
                                    keep = c(".trt",outcome_label, "def")) %>% 
      summary() %>% 
      filter(.observed == TRUE)  %>% 
      select(-.observed) 
    
    def_total_res_low <- fit_def_low %>% vit$marginalize(margins = c(outcome_label, "def")) %>% 
      summary() 
    
    bias_reg_res <- vit$summary_brms_nma(bias_reg,
                                         keep = c(".trt",".idx_anx",".idx_dep",outcome_label, "def","bias.overall"), index = ".obs", resp = "y"
                                         ) %>% 
      # filter(.idx_anx & .idx_dep) %>% 
      # filter(.data[[not_outcome_label]] == 0) %>% 
      filter(.observed == TRUE) %>% 
      summary() %>% 
      select(-c(.observed,.idx_anx, .idx_dep,!!not_outcome_label))
  }
  
  # save basic results tables
  {
    low_res %>% 
      gt_fmt() %>% 
      saveRDS(
        # gt::gtsave(
        file = here::here(output_dir,
                              glue::glue("nma_{outcome}_low_bias_table.rds"))
      )
    
    some_res %>% 
      gt_fmt() %>% 
      saveRDS(
        # gt::gtsave(
        file = here::here(output_dir,
                              glue::glue("nma_{outcome}_some_and_low_bias_table.rds"))
      )
    
    high_res %>% 
      gt_fmt() %>%
      saveRDS(
      # gt::gtsave(
        file = here::here(output_dir,
                              glue::glue("nma_{outcome}_all_data_table.rds"))
      )
    
    dose_res %>% 
      mutate(
        .dose = factor(
          .dose,
          levels = levels(.dose),
          labels = stringr::str_wrap(levels(.dose), width = 20)
        )
      ) %>%
      gt_fmt(groupname_col = "Intervention") %>% 
      saveRDS(
        # gt::gtsave(
        file = here::here(output_dir,
                              glue::glue("nma_{outcome}_dose_response_table_low_bias.rds"))
      )
    
    def_res_low %>% 
      gt_fmt(groupname_col = "Intervention") %>% 
      saveRDS(
        # gt::gtsave(
        file = here::here(output_dir,
                              glue::glue("nma_{outcome}_deficiency_table_low_bias.rds"))
      )
    
    time_res %>% 
      gt_fmt(groupname_col = "Intervention") %>% 
      saveRDS(
        # gt::gtsave(
        file = here::here(output_dir,
                              glue::glue("nma_{outcome}_time_response_table_low_bias.rds"))
      )
    
    b_def_t <- bias_reg_res %>% 
      gt_fmt(groupname_col = "Intervention")  %>% 
      saveRDS(
        # gt::gtsave(
        file = here::here(output_dir,
                              glue::glue("nma_{outcome}_meta_regression_table.rds"))
      )
    
    def_res_low %>% 
      group_by(.trt) %>%
      tidyr::complete(
        !!outcome_label := c("yes", "no"),
        def = c("yes", "no", "NA")
      ) %>%
      gt_fmt(groupname_col = "Intervention") %>% 
      saveRDS(
        file = here::here(output_dir,
                          glue::glue("nma_{outcome}_def_low_bias_table.rds"))
      )
    
    def_res %>% 
      gt_fmt(groupname_col = "Intervention") %>%
      saveRDS(
        file = here::here(output_dir,
                          glue::glue("nma_{outcome}_def_all_data_table.rds"))
      )
    
    grid_apx <- tidyr::expand_grid(.trt = fit_pcts_approx$prep$network$treatments, p_dep = c(0,1), p_def = c(0,1))
    grid_apx %>% model.matrix(~ .trt : p_dep + .trt : p_def + .trt : p_dep : p_def + .trt, data = .) %>%  as.data.frame() %>% vit$special_clean_names() %>% mutate(.study = fit_pcts_approx$data$.study[1], .se = fit_pcts_approx$data$.se[1], .cov = fit_pcts_approx$data$.cov[1], p_def = grid_apx$p_def, p_dep = grid_apx$p_dep, .trt = grid_apx$.trt) -> newdata_apx
    pct_sum_apx <- vit$summary_brms_nma(fit_pcts_approx, newdata = newdata_apx, keep = c(".trt", 'p_dep', "p_def")) %>% 
      mutate(p_dep = ifelse(p_dep == 1, "yes", "no") %>% as.factor(),
             p_def = ifelse(p_def == 1, "yes", "no") %>% as.factor()) %>% 
      rename(!!outcome_label := p_dep, def = p_def)
    
    sel_trt <- fit_pcts_approx$prep$data %>% group_by(intervention) %>% 
      mutate(intervention = stringr::str_to_title(intervention)) %>% 
      summarize(n = n())
    pct_sum_apx %>% summary() %>% 
      filter(.trt %in% (sel_trt %>% filter(n > 1) %>% pull(intervention))) %>%
      gt_fmt(groupname_col = "Intervention") %>%
      saveRDS( 
        file = here::here(output_dir,
                          glue::glue("nma_{outcome}_pct_table_low_bias.rds"))
      )
  }
  
  #loo estimates 
  {
    loo_gauss_low <- readRDS(
      here::here(sensitivity_outdir,
                 "loo_gaussian_low.rds")
    )
    
    loo_gauss_some <- readRDS(
      here::here(sensitivity_outdir,
                 "loo_gaussian_some.rds")
    )
    
    loo_gauss_high <- readRDS(
      here::here(sensitivity_outdir,
                 "loo_gaussian_high.rds")
    )
    
    loo_low <- readRDS(
      here::here(sensitivity_outdir,
                 "loo_low.rds")
    )
    
    loo_some <- readRDS(
      here::here(sensitivity_outdir,
                 "loo_some.rds")
    )
    
    loo_high <- readRDS(
      here::here(sensitivity_outdir,
                 "loo_high.rds")
    )
    
    loo_meta_reg <- readRDS(
      here::here(sensitivity_outdir,
                 "bias_regression_loo.rds")
    )
    
    loo_def_low <- readRDS(
      here::here(sensitivity_outdir,
                 "loo_def_low.rds")
    )
    
    loo_def     <- readRDS(
      here::here(sensitivity_outdir,
                 "loo_def.rds")
    )
    
    loo_ume_low <- readRDS(
      here::here(sensitivity_outdir,
                 "loo_ume_low.rds"))
    
    loo_ume_some <- readRDS(
      here::here(sensitivity_outdir,
                 "loo_ume_some.rds"))
    
    loo_ume_high <- readRDS(
      here::here(sensitivity_outdir,
          "loo_ume_high.rds"))
    
    
  }
  
  # create loo lists
  {
    ord <- c("student-t","gaussian")
    compare_gauss <-
      list(low = loo::loo_compare(list("student-t" = vit$loo.vit_loo(loo_low), "gaussian" = vit$loo.vit_loo(loo_gauss_low)))[ord,],
           some = loo::loo_compare(list("student-t" = vit$loo.vit_loo(loo_some), "gaussian" = vit$loo.vit_loo(loo_gauss_some)))[ord,],
           high = loo::loo_compare(list("student-t" = vit$loo.vit_loo(loo_high), "gaussian" = vit$loo.vit_loo(loo_gauss_high)))[ord,]
           # ,high_w_meta_reg = loo::loo_compare(list("student-t" = vit$loo.vit_loo(loo_high), "gaussian" = vit$loo.vit_loo(loo_gauss_high), "meta regression" = vit$loo.vit_loo(loo_meta_reg)))[c("meta regression",ord),]
      )
    
    compare_mod <- list(
      low = loo::loo_compare(list("simple model" = vit$loo.vit_loo(loo_low),  
                                  "by disease/deficiency" = vit$loo.vit_loo(loo_def_low),                          
                                  "UME" = vit$loo.vit_loo(loo_ume_low))),
      some = loo::loo_compare(list("simple model" = vit$loo.vit_loo(loo_some), "UME" = vit$loo.vit_loo(loo_ume_some))),
      high = loo::loo_compare(list("simple model" = vit$loo.vit_loo(loo_high), "by disease/deficiency" = vit$loo.vit_loo(loo_def), 
                                   "UME" = vit$loo.vit_loo(loo_ume_high)))
    )
  }
  
  # create loo table
  {
    loo_tab <- dplyr::tibble(
      Likelihood = rownames(compare_gauss$high) %>% 
        as.factor() %>%
        recode("student-t" = "Student-t",
               "gaussian" = "Gaussian"
               # ,"meta regression" = "Student-t (meta-regression)"
        ),
      low = c(
        # NA_real_, 
        compare_gauss$low[, "elpd_diff"]
      ),
      lowse = c(
        # NA_real_, 
        compare_gauss$low[, "se_diff"]
      ),
      'unclear/low' = c(
        # NA_real_,
        compare_gauss$some[, "elpd_diff"]
      ),
      'unclear/lowse' = c(
        # NA_real_,
        compare_gauss$some[, "se_diff"]
      ),
      'high/unclear/low' = compare_gauss$high[, "elpd_diff"],
      'high/unclear/lowse' = compare_gauss$high[, 'se_diff']
      
    ) %>% gt() %>% 
      gt::cols_merge(
        columns = c(low, lowse),
        pattern = "{1} ({2})"
      ) %>% 
      gt::cols_merge(
        columns = c('unclear/low', 'unclear/lowse'),
        pattern = "{1} ({2})"
      ) %>%
      gt::cols_merge(
        columns = c('high/unclear/low', 'high/unclear/lowse'),
        pattern = "{1} ({2})"
      ) %>%
      # tab_spanner(
      #   label = "ELPD",
      #   columns = c(diff, se)
      # ) %>%
      fmt_number(decimals = 3) %>% 
      cols_label(
        Likelihood = "Likelihood",
        low = "Low Bias",
        'unclear/low' = "Unclear/Low Bias",
        'high/unclear/low' = "High/Unclear/Low Bias"
        # , diff = "Difference",
        # se   = "S.E."
      ) %>% 
      tab_header(
        title = "Comparison of Model Fit by Likelihood Function. Differences in expected log point-wise predictive densities (ELPD) estimated by leave-one-out cross validation. The second through third columns list the difference in ELPDs separated by the bias of the included studies. Standard error of difference in parentheses.") 
    
    # save gaussian comparisons
    saveRDS(loo_tab,
            here::here(output_dir,
                       "loo_compare_tab.rds"))
    }
  
  # other loo tabs
  {
    compare_mod$some <- rbind(compare_mod$some,
                              rep(NA_real_,
                                  ncol(compare_mod$some)) %>% setNames(colnames(compare_mod$some)))
    rownames(compare_mod$some)[nrow(compare_mod$some)] <- "by disease/deficiency"
    loo_def_tab <- compare_mod %>%
      purrr::imap(function(df, nm) {
        dplyr::as_tibble(df, rownames = "Model") %>%
          dplyr::mutate(bias = nm) %>%
          dplyr::mutate(elpd_diff = as.numeric(elpd_diff),
                        se_diff = as.numeric(se_diff)) %>%
          dplyr::mutate(Model = recode(Model,
                                     "simple model" = "Simple model",
                                     "by disease/deficiency" = "By Disease/Deficiency",
                                     "UME" = "Unrestricted Model of Effects")) %>%
          dplyr::select(Model, elpd_diff, se_diff, bias)
      }) %>%
      dplyr::bind_rows() %>% 
      tidyr::pivot_wider(names_from = bias, values_from = c(elpd_diff, se_diff)) %>%
      gt() %>% 
      gt::sub_missing() %>% 
      gt::cols_merge(
        columns = c(elpd_diff_low, se_diff_low),
        pattern = "{1} ({2})"
      ) %>% 
      gt::cols_merge(
        columns = c(elpd_diff_some, se_diff_some),
        pattern = "{1} ({2})"
      ) %>%
      gt::cols_merge(
        columns = c(elpd_diff_high, se_diff_high),
        pattern = "{1} ({2})"
      ) %>%
      fmt_number(decimals = 3) %>% 
      cols_label(
        Model = "Model",
        elpd_diff_low = "Low Bias",
        elpd_diff_some = "Unclear/Low Bias",
        elpd_diff_high = "High/Unclear/Low Bias"
      ) %>% 
      gt::fmt_number(decimals = 1) %>% 
      tab_header(
        title = "Comparison of Model Fit by Bias. Differences in expected log point-wise predictive densities (ELPD) estimated by leave-one-out cross validation. The second through third columns list the difference in ELPDs separated by the bias of the included studies. Standard error of difference in parentheses.") 
    
    # save model comparisons
    saveRDS(loo_def_tab,
            here::here(output_dir,
                       "loo_compare_model_tab.rds"))
    
  }
  
  # sucra and rank probs; leaving out given varying effects depending on population
  {
    rank_probs_low <- low_t %>% 
      vit$rank_probs_vit(sucra = TRUE,
                         cumulative = FALSE) %>% 
      rename(
        Treatment = .trt,
        SUCRA = sucra
      ) %>%
      rename_with(
        ~ stringr::str_replace(.x, "p_rank\\[(\\d+)\\]", "P(rank = \\1)"),
        starts_with("p_rank")
      ) %>% 
      select(Treatment, SUCRA, everything()) %>% 
      gt() %>% 
      fmt_number()
    
    rank_probs_all <- high_t %>% 
      vit$rank_probs_vit(sucra = TRUE,
                         cumulative = FALSE) %>% 
      rename(
        Treatment = .trt,
        SUCRA = sucra
      ) %>%
      rename_with(
        ~ stringr::str_replace(.x, "p_rank\\[(\\d+)\\]", "P(rank = \\1)"),
        starts_with("p_rank")
      ) %>% 
      select(Treatment, SUCRA, everything()) %>% 
      gt() %>% 
      fmt_number()
    
    rank_probs_def <-  fit_def_low %>%
      vit$rank_probs_vit(sucra = TRUE,
                         cumulative = FALSE,
                         filter = rlang::quo(!!rlang::sym(outcome_label) == "yes" & def == "yes"),
                         keep = c(".trt","def",outcome_label)) %>% 
      mutate(
        across(
          starts_with("p_rank"),
          ~tidyr::replace_na(.x, 0)
        )
      ) %>% 
      rename(
        Treatment = .trt,
        SUCRA = sucra
      ) %>%
      rename_with(
        ~ stringr::str_replace(.x, "p_rank\\[(\\d+)\\]", "P(rank = \\1)"),
        starts_with("p_rank")
      ) %>% 
      select(Treatment, SUCRA, everything()) %>% 
      gt() %>% 
      fmt_number()
    
    rank_probs_def_not_def <-  fit_def_low %>%
      vit$rank_probs_vit(sucra = TRUE,
                         cumulative = FALSE,
                         filter = rlang::quo(!!rlang::sym(outcome_label) == "yes" & def == "no"),
                         keep = c(".trt","def",outcome_label)) %>% 
      mutate(
        across(
          starts_with("p_rank"),
          ~tidyr::replace_na(.x, 0)
        )
      ) %>% 
      rename(
        Treatment = .trt,
        SUCRA = sucra
      ) %>%
      rename_with(
        ~ stringr::str_replace(.x, "p_rank\\[(\\d+)\\]", "P(rank = \\1)"),
        starts_with("p_rank")
      ) %>% 
      select(Treatment, SUCRA, everything()) %>% 
      gt() %>% 
      fmt_number()
    
    
    rank_probs_def_unknown_def <-  fit_def_low %>%
      vit$rank_probs_vit(sucra = TRUE,
                         cumulative = FALSE,
                         filter = rlang::quo(!!rlang::sym(outcome_label) == "yes" & def == "NA"),
                         keep = c(".trt","def",outcome_label)) %>% 
      mutate(
        across(
          starts_with("p_rank"),
          ~tidyr::replace_na(.x, 0)
        )
      ) %>% 
      rename(
        Treatment = .trt,
        SUCRA = sucra
      ) %>%
      rename_with(
        ~ stringr::str_replace(.x, "p_rank\\[(\\d+)\\]", "P(rank = \\1)"),
        starts_with("p_rank")
      ) %>% 
      select(Treatment, SUCRA, everything()) %>% 
      gt() %>% 
      fmt_number()
  }
  
  # contribution matrices
  {
    
    # yhat_lowt <- fitted.values(low_t); 
    # tx_lev <- low_t$prep$network$treatments %>% levels()
    # cinema <- low_t$prep$data %>% 
    #   mutate(id = as.integer(.study),  
    #          indirectness =1L, 
    #          rob = ifelse(bias == "low","1",ifelse(bias == "some concerns","2","3")) %>% as.numeric, 
    #          t2 = gsub(",z\\.",", ",as.character(control)), 
    #          t2 = gsub("z.","",t2)) %>% 
    #   filter(!is.na(y)) %>% 
    #   mutate(t2 = t2 %>% factor(levels = tx_lev) %>% 
    #            as.integer() %>% as.character()%>% paste0("t",.)
    #          ) %>% 
    #   mutate(t1 = intervention %>% as.character() %>% factor(levels = tx_lev) %>% as.integer() %>% as.character()%>% paste0("t",.)) %>% 
    #   select(id,t1 = t1, t2, rob,indirectness) %>% 
    #   mutate(effect = yhat_lowt[,1], se = yhat_lowt[,2]) %>% 
    #     bind_rows(
    #       low_t %>% vit$summary_brms_nma() %>% vit$contrasts("Magnesium","Zinc") %>% select(value) %>% summarize(id = 12, t1 = factor("magnesium",levels = tx_lev) %>% as.integer() %>% as.character() %>% paste0("t",.),t2 = factor("zinc",levels = tx_lev) %>% as.integer() %>% as.character() %>% paste0("t",.), rob = 1, indirectness = 1, effect = mean(value), se =sd(value))
    #       )
    #   write.csv(cinema, file = "cinema.csv")
    cm_low  <- vit$contribution_matrix(loo_low, fit = low_t)
    cm_high <- vit$contribution_matrix(loo_high, fit = high_t)
    
    cm_low_approx  <- vit$contribution_matrix(low_t) %>% 
      vit$alter_netcontrib(low_t)
    cm_high_approx <- vit$contribution_matrix(high_t) %>% 
      vit$alter_netcontrib(high_t)
    
    cm_low %>% 
      filter(trt2 == "Placebo") %>%  
      arrange(trt1) %>% 
      mutate(contribution = contribution * 100) %>% 
      tidyr::pivot_wider(id_cols = study,
                         names_from = comparison,
                         values_from = contribution) %>% 
      rename(Study = study) %>% 
      gt() %>%  
      fmt_number()  %>% 
      # gtsave(
      #   filename = here::here(output_dir,
      #                         glue::glue("contribution_matrix_low_bias.tex"))
      # )
       saveRDS(
         here::here(output_dir,
                    glue::glue("contribution_matrix_low_bias.rds"))
       )
       
       cm_low_approx %>% 
         filter(trt2 == "Placebo") %>%  
         arrange(trt1) %>% 
         mutate(contribution = contribution * 100) %>% 
         tidyr::pivot_wider(id_cols = study,
                            names_from = comparison,
                            values_from = contribution) %>% 
         rename(Study = study) %>% 
         gt() %>% 
         fmt_number() %>% 
         saveRDS(
           here::here(output_dir,
                      glue::glue("approx_contribution_matrix_low_bias.rds"))
         )
    
       cm_high %>% 
         filter(trt2 == "Placebo") %>%  
         arrange(trt1) %>% 
         mutate(contribution = contribution * 100) %>% 
         tidyr::pivot_wider(id_cols = study,
                            names_from = comparison,
                            values_from = contribution) %>% 
         left_join(
           high_t$prep$data %>% 
             select(study, bias) %>% 
             distinct(),
           by = "study"
         ) %>%
         select(study, bias, everything()) %>% 
         mutate(bias = factor(as.character(bias), 
                              labels = c("Low","Unclear","High"),
                              levels = c("low",
                                         "some concerns",
                                         "high"),
                              ordered = TRUE)) %>%
         rename(Bias = bias, Study = study) %>% 
         gt() %>% 
         fmt_number() %>% 
      # gtsave(
      #   filename = here::here(output_dir,
      #                         glue::glue("contribution_matrix_all_data.tex"))
      # )
      saveRDS(
        here::here(output_dir,
                   glue::glue("contribution_matrix_all_data.rds"))
      )
      
      cm_high_approx %>% 
        filter(trt2 == "Placebo") %>%  
        arrange(trt1) %>% 
        mutate(contribution = contribution * 100) %>% 
        tidyr::pivot_wider(id_cols = study,
                           names_from = comparison,
                           values_from = contribution) %>% 
        left_join(
          high_t$prep$data %>% 
            select(study, bias) %>% 
            distinct(),
          by = "study"
        ) %>%
        select(study, bias, everything()) %>% 
        mutate(bias = factor(as.character(bias), 
                             labels = c("Low","Unclear","High"),
                             levels = c("low",
                                        "some concerns",
                                        "high"),
                             ordered = TRUE)) %>%
        rename(Bias = bias, Study = study) %>% 
        gt() %>% 
        fmt_number() %>% 
      # gtsave(
      #   filename = here::here(output_dir,
      #                         glue::glue("contribution_matrix_all_data.tex"))
      # )
      saveRDS(
        here::here(output_dir,
                   glue::glue("approx_contribution_matrix_all_data.rds"))
      )
  }
  
  # cinema tables
  {
    cap_outcome_label <- stringr::str_to_title(outcome_label)
    
    cinema_filter <- function(x) {
      x %>% filter(!is.na(n_study) | trt2 == "Placebo" | overall_confidence %in% c("Low","Moderate","High","Very high"))
    }
    
    indirectness_low <- if(outcome == "depression") {
      data.frame(.trt =  c("B6", "B9, Iron", "C", "D, Iron", "Magnesium, Zinc", "B1", "Antidepressant, B12, B6, B9", "Antidepressant, D"),
                 rating = "Some Concerns")
    } else {
      data.frame(.trt = low_res$.trt %>% unique() %>% setdiff(c("Antidepressant","D","Magnesium","Placebo","Zinc")),
                 rating = "Some Concerns")
    }
    
    # low bias model
      cin_low  <- vit$cinema(low_t, 
                             contribution.matrix = cm_low,
                             indirectness = indirectness_low,
                             keep = c(".trt")) 
    
    # some bias model
      cin_some <- vit$cinema(some_t, 
                             contribution.matrix = loo_some,
                             keep = c(".trt")) 
      
    # high bias model
      cin_high <- vit$cinema(high_t, 
                             contribution.matrix = cm_high,
                             keep = c(".trt")) 
    
    # low bias deficiency and disease status
      cin_low_def <- vit$cinema(fit_def_low, 
                             contribution.matrix = loo_def_low,
                             indirectness = indirectness_low,
                             keep = c(".trt",outcome_label, "def")) %>% 
      mutate(def = forcats::fct_recode(
        as.factor(def), "Unknown" = "NA",
        "Yes" = "yes","No" = "no")) %>% 
      mutate(def = forcats::fct_relevel(def, "Unknown", "No","Yes")) %>%
      rename(Deficiency = def) %>%
      mutate(!!outcome_label := forcats::fct_recode(.data[[outcome_label]] %>% as.factor(),
                                             "Yes" = "yes",
                                             "No" = "no")) %>% 
      rename(!!cap_outcome_label := !!outcome_label) 
    
    cin_low_def_no_iran <- vit$cinema(fit_def_low_no_iran, 
                                      indirectness = indirectness_low,
                              keep = c(".trt",outcome_label, "def")) %>% 
      mutate(def = forcats::fct_recode(
        as.factor(def), "Unknown" = "NA",
        "Yes" = "yes","No" = "no")) %>% 
      mutate(def = forcats::fct_relevel(def, "Unknown", "No","Yes")) %>%
      rename(Deficiency = def) %>%
      mutate(!!outcome_label := forcats::fct_recode(.data[[outcome_label]] %>% as.factor(),
                                                    "Yes" = "yes",
                                                    "No" = "no")) %>% 
      rename(!!cap_outcome_label := !!outcome_label) 
    
    # all bias deficiency and disease status
    cin_def <- vit$cinema(fit_def, 
                          contribution.matrix = loo_def,
                          keep = c(".trt",outcome_label, "def","bias")) %>% 
      mutate(def = forcats::fct_recode(
        as.factor(def), "Unknown" = "NA",
        "Yes" = "yes","No" = "no")) %>% 
      mutate(def = forcats::fct_relevel(def, "Unknown", "No","Yes")) %>%
      rename(Deficiency = def) %>%
      mutate(!!outcome_label := forcats::fct_recode(.data[[outcome_label]] %>% as.factor(),
            "Yes" = "yes",
            "No" = "no")) %>% 
      rename(!!cap_outcome_label := !!outcome_label) %>% 
      mutate(bias = factor(as.character(bias), 
           labels = c("Low","Unclear","High"),
           levels = c("low",
                      "some concerns",
                      "high"),
           ordered = TRUE)) %>%
      rename(Bias = bias)
    
    
    cin_def_no_iran <- vit$cinema(fit_def_no_iran, 
                          keep = c(".trt",outcome_label, "def","bias")) %>% 
      mutate(def = forcats::fct_recode(
        as.factor(def), "Unknown" = "NA",
        "Yes" = "yes","No" = "no")) %>% 
      mutate(def = forcats::fct_relevel(def, "Unknown", "No","Yes")) %>%
      rename(Deficiency = def) %>%
      mutate(!!outcome_label := forcats::fct_recode(.data[[outcome_label]] %>% as.factor(),
                                                    "Yes" = "yes",
                                                    "No" = "no")) %>% 
      rename(!!cap_outcome_label := !!outcome_label) %>% 
      mutate(bias = factor(as.character(bias), 
                           labels = c("Low","Unclear","High"),
                           levels = c("low",
                                      "some concerns",
                                      "high"),
                           ordered = TRUE)) %>%
      rename(Bias = bias)
    
    # cin_def_no_bias <- vit$cinema(fit_def, 
    #             contribution.matrix = loo_def %>% 
    #               group_by(iter,.trt,
    #                        def, !!outcome_label,
    #                        left_out_study,
    #                        loo_iter
    #                        ) %>% 
    #               summarize(value = mean(value)), # everything but bias
    #             keep = c(".trt",outcome_label, "def"), re_formula = NA) %>% 
    #   mutate(def = forcats::fct_recode(
    #     as.factor(def), "Unknown" = "NA",
    #     "Yes" = "yes","No" = "no")) %>% 
    #   mutate(def = forcats::fct_relevel(def, "Unknown", "No","Yes")) %>%
    #   rename(Deficiency = def) %>%
    #   mutate(!!outcome_label := forcats::fct_recode(.data[[outcome_label]] %>% as.factor(),
    #                                                 "Yes" = "yes",
    #                                                 "No" = "no")) %>% 
    #   rename(!!cap_outcome_label := !!outcome_label)
    # 
    
    # save gt tables
    cin_low %>% cinema_filter() %>% 
      vit$cinema_table() %>% 
      saveRDS(
        here::here(output_dir,
                   glue::glue("cinema_table_low_bias.rds"))
      )
    
    cin_low %>% 
      filter(trt2 == "Placebo") %>% 
      vit$cinema_table() %>% 
      saveRDS(
        here::here(output_dir,
                   glue::glue("cinema_table_low_bias_placebo.rds"))
      )
    
    cin_some %>% cinema_filter() %>% 
      vit$cinema_table() %>% 
      saveRDS(
        here::here(output_dir,
                   glue::glue("cinema_table_some_and_low_bias.rds"))
      )
    
    cin_some %>% 
      filter(trt2 == "Placebo") %>% 
      vit$cinema_table() %>% 
      saveRDS(
        here::here(output_dir,
                   glue::glue("cinema_table_some_and_low_bias_placebo.rds"))
      )
    
    cin_high %>% cinema_filter() %>% 
      vit$cinema_table() %>% 
      saveRDS(
        here::here(output_dir,
                   glue::glue("cinema_table_all_data.rds"))
      )
    
    cin_high %>% 
      filter(trt2 == "Placebo") %>% 
      vit$cinema_table() %>% 
      saveRDS(
        here::here(output_dir,
                   glue::glue("cinema_table_all_data.rds"))
      )
    
    cin_low_def %>% 
      cinema_filter() %>% 
      vit$cinema_table() %>% 
      saveRDS(
        here::here(output_dir,
                   glue::glue("cinema_table_deficiency_low_bias.rds"))
      )
    
    cin_low_def %>% 
      filter(trt2 == "Placebo") %>% 
      vit$cinema_table() %>% 
      saveRDS(
        here::here(output_dir,
                   glue::glue("cinema_table_deficiency_low_bias_placebo.rds"))
      )
    
    cin_low_def_no_iran %>% 
      filter(trt2 == "Placebo") %>% 
      vit$cinema_table() %>% 
      saveRDS(
        here::here(output_dir,
                   glue::glue("cinema_table_deficiency_low_bias_no_iran_placebo.rds"))
      )
    
    cin_def %>% 
      cinema_filter() %>% 
      vit$cinema_table()  %>% 
      saveRDS(
        here::here(output_dir,
                   glue::glue("cinema_table_deficiency_all_data.rds"))
      )
    
    cin_def %>% 
      filter(trt2 == "Placebo") %>% 
      vit$cinema_table() %>% 
      saveRDS(
        here::here(output_dir,
                   glue::glue("cinema_table_deficiency_all_data_placebo.rds"))
      )
    
    # cin_def_no_bias %>% 
    #   filter(trt2 == "Placebo") %>% 
    #   vit$cinema_table()  %>% 
    #   saveRDS(
    #     here::here(output_dir,
    #                glue::glue("cinema_table_deficiency_all_data_no_bias_placebo.rds"))
    #   )
    
  }
  
}

