# this code will generate all of the tables needed for the 
# main manuscript and supplementary materials

#### Load Packages ####
suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
  library(gt)
  library(stringr)
  library(tidyr)
  library(flextable)
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

mem.maxVSize(1024*1024)

#### Load Data ####
# full data
data <- vit$get_vitamin_data(outcome = NULL, simple_analysis = FALSE,
                             include_full_bias = TRUE)

main_data <- vit$get_vitamin_data(outcome = NULL, simple_analysis = TRUE,
                                  include_full_bias = FALSE) %>% 
  filter(bias == "low")

main   <- main_data %>%  
  distinct(study) %>% pull(study)

#### setup directories ####
{
  if (!dir.exists(here::here("outputs","overall","tables","raw"))) {
    dir.create(here::here("outputs","overall","tables","raw"), recursive = TRUE)
  }
  
  if (!dir.exists(here::here("outputs","depression","tables","raw"))) {
    dir.create(here::here("outputs","depression","tables","raw"), recursive = TRUE)
  }
  
  if (!dir.exists(here::here("outputs","anxiety","tables","raw"))) {
    dir.create(here::here("outputs","anxiety","tables","raw"), recursive = TRUE)
  }
  
  if (!dir.exists(here::here("outputs","overall","tables","gt"))) {
    dir.create(here::here("outputs","overall","tables","gt"), recursive = TRUE)
  }
  
  if (!dir.exists(here::here("outputs","depression","tables","gt"))) {
    dir.create(here::here("outputs","depression","tables","gt"), recursive = TRUE)
  }
  
  if (!dir.exists(here::here("outputs","anxiety","tables","gt"))) {
    dir.create(here::here("outputs","anxiety","tables","gt"), recursive = TRUE)
  }
  
  if (!dir.exists(here::here("outputs","overall","tables","ft"))) {
    dir.create(here::here("outputs","overall","tables","ft"), recursive = TRUE)
  }
  
  if (!dir.exists(here::here("outputs","depression","tables","ft"))) {
    dir.create(here::here("outputs","depression","tables","ft"), recursive = TRUE)
  }
  
  if (!dir.exists(here::here("outputs","anxiety","tables","ft"))) {
    dir.create(here::here("outputs","anxiety","tables","ft"), recursive = TRUE)
  }
  
  if (!dir.exists(here::here("outputs","overall","tables","docx"))) {
    dir.create(here::here("outputs","overall","tables","docx"), recursive = TRUE)
  }
  
  if (!dir.exists(here::here("outputs","depression","tables","docx"))) {
    dir.create(here::here("outputs","depression","tables","docx"), recursive = TRUE)
  }
  
  if (!dir.exists(here::here("outputs","anxiety","tables","docx"))) {
    dir.create(here::here("outputs","anxiety","tables","docx"), recursive = TRUE)
  }
}


#### temp functions ####
{
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
                                "Low" = "low",
                                "Unclear" = "some concerns",
                                "High" = "high")
        )
      ) %>%
      mutate(
        across(
          any_of("Bias"),
          ~ forcats::fct_relevel(.x, "Low", "Unclear", "High")
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
  
  flextable_fmt <- function(x, digits = 2, ...) {
    
    fmt <- function(v, d) formatC(v, format = "f", digits = d)
    
    safe_merge_hierarchy <- function(ft,
                                     keys = c("Intervention", "Depressed", "Deficient"),
                                     targets = c("Deficient", "Depressed", "Intervention")) {
      # ft: a flextable made from `data`
      # data: the data.frame used to build ft (used only to detect columns robustly)
      # keys: ordered left->right hierarchy (most general to most specific)
      # targets: merge targets in order most-specific->most-general (right->left)
      
      if (!length(keys) || !length(targets)) return(ft)
      
      # helper: keys up to and including a target, e.g.
      # target="Deficient" -> keys=c("Intervention","Depressed","Deficient")
      key_upto <- function(target) {
        if (!target %in% keys) return(keys[1]) # fallback
        keys[seq_len(match(target, keys))]
      }
      
      # merge right-to-left (most specific first)
      for (tgt in targets) {
        kj <- key_upto(tgt)
        
        # If we're merging Depressed but Deficient exists, we *don't* include Deficient in the key
        # because we want Depressed merged across all Deficient levels within Intervention.
        # (That’s the usual journal-table look.)
        if (tgt == "Depressed" && "Depressed" %in% keys) {
          kj <- keys[seq_len(match("Depressed", keys))]
        }
        if (tgt == "Intervention" && "Intervention" %in% keys) {
          kj <- "Intervention"
        }
        
        ft <- flextable::merge_v(
          x = ft,
          j = kj,
          target = tgt,
          combine = TRUE
        )
      }
      
      ft
    }
    
    df <- x %>% 
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
                                "Low" = "low",
                                "Unclear" = "some concerns",
                                "High" = "high")
        )
      ) %>%
      mutate(
        across(
          any_of("Bias"),
          ~ forcats::fct_relevel(.x, "Low", "Unclear", "High")
        )
      ) %>% 
      arrange(across(any_of(c("Intervention", "Bias","Dose", "Depressed","Anxious",  "Deficient")))) %>% 
      select(any_of(c("Intervention", "Bias","Dose", "Depressed","Anxious",  "Deficient")), 
             everything()) %>% 
      select(-c(S.E., `q50%`)) %>%
      mutate(
        est_ci = dplyr::case_when(
          # only compute if these columns exist; otherwise leave NA and we won't drop them
          !all(c("Estimate", "q2.5%", "q97.5%") %in% names(.)) ~ NA_character_,
          is.na(Estimate) & is.na(`q2.5%`) & is.na(`q97.5%`) ~ "—",
          TRUE ~ glue::glue("{fmt(Estimate, digits)} ({fmt(`q2.5%`, digits)}, {fmt(`q97.5%`, digits)})")
        )
      ) %>%
      select(-c(Estimate, `q2.5%`, `q97.5%`)) 
    
    merge_cols <- intersect(c("Intervention", "Bias","Dose", "Depressed","Anxious",  "Deficient"), colnames(df))
    
    df %>% 
      # { if (length(merge_cols)>1) flextable::as_grouped_data(., groups = "Intervention")  else . } %>%
      flextable::flextable() %>% 
      # { if (length(merge_cols)>1) merge_v(., j = merge_cols) else . } %>%
      { if (length(merge_cols)>1) {
        safe_merge_hierarchy(
          ft = .,
          keys = merge_cols,
          targets = rev(merge_cols)  # right-to-left
        )
      }  else . } %>%
      flextable::set_header_labels(
        Intervention = "Intervention",
        est_ci = "Effect Size (95% Credible Interval)"
      ) %>%
      flextable::autofit() %>% 
      flextable::hline(border = officer::fp_border(color = "gray70", 
                                                   width = 0.5), 
                       part = "body")
    
  }
  
  save_gt_mult_formats <- function(x, file.name) {
    stopifnot(inherits(x, "gt_tbl"))
    # raw gt
    saveRDS(x, paste0(file.name, ".rds"))
    
    x %>% gt::gtsave(paste0(file.name, ".docx"))
  }
  
  save_mult_formats <- function(x, outdir, file.name, group = NULL) {
    
    # save raw data
    if(!dir.exists(here::here(outdir, "raw")) ) {
      dir.create(here::here(outdir, "raw"), recursive = TRUE)
    }
    saveRDS(x, here::here(outdir, "raw", paste0(file.name, ".rds")))
    
    # save gt
    ggt <- x %>% gt_fmt(groupname_col = group)
    stopifnot(inherits(ggt, "gt_tbl"))
    
    # raw gt
    if(!dir.exists(here::here(outdir, "gt")) ) {
      dir.create(here::here(outdir, "gt"), recursive = TRUE)
    }
    saveRDS(ggt, here::here(outdir, "gt", paste0(file.name, ".rds")))
    
    # gen flextable 
    ft <- x %>% flextable_fmt() 
    
    # save flextable
    if(!dir.exists(here::here(outdir, "ft")) ) {
      dir.create(here::here(outdir, "ft"), recursive = TRUE)
    }
    saveRDS(ft, here::here(outdir,  "ft", paste0(file.name, ".rds")))
    
    # save flextable word
    if(!dir.exists(here::here(outdir, "docx")) ) {
      dir.create(here::here(outdir, "docx"), recursive = TRUE)
    }
    ft %>% flextable::save_as_docx(path = here::here(outdir, "docx", paste0(file.name, ".docx")))
    
  }
  
  save_cm_mult_formats <- function(x, outdir, file.name, group = NULL) {
    
    # save raw data
    if(!dir.exists(here::here(outdir, "raw")) ) {
      dir.create(here::here(outdir, "raw"), recursive = TRUE)
    }
    saveRDS(x, here::here(outdir, "raw", paste0(file.name, ".rds")))
    
    # save gt
    ggt <- x %>% gt(groupname_col = group) %>% 
      gt::fmt_number(decimals = 2)
    
    # raw gt
    if(!dir.exists(here::here(outdir, "gt")) ) {
      dir.create(here::here(outdir, "gt"), recursive = TRUE)
    }
    saveRDS(ggt, here::here(outdir, "gt", paste0(file.name, ".rds")))
    
    # gen flextable 
    ft <- x %>% 
      flextable::flextable() 
    
    num_cols <- names(x)[sapply(x, is.numeric)]
    
    ft <- flextable::colformat_double(
      ft,
      j = num_cols,
      digits = 2,       # adjust as needed
      big.mark = ",",   # thousands separator
      decimal.mark = "."
    )
    
    # save flextable
    if(!dir.exists(here::here(outdir, "ft")) ) {
      dir.create(here::here(outdir, "ft"), recursive = TRUE)
    }
    saveRDS(ft, here::here(outdir,  "ft", paste0(file.name, ".rds")))
    
    # save flextable word
    if(!dir.exists(here::here(outdir, "docx")) ) {
      dir.create(here::here(outdir, "docx"), recursive = TRUE)
    }
    ft %>% flextable::save_as_docx(path = here::here(outdir, "docx", paste0(file.name, ".docx")))
    
  }
  
  save_cinema_mult_formats <- function(x, outdir, file.name) {
    x %>% 
      vit$cinema_table(format = "gt") %>% 
      saveRDS(here::here(outdir,
                         "gt",
                         glue::glue("{file.name}.rds")))
    
    ft <- x %>% 
      vit$cinema_table(format = "flextable") 
    
    ft %>% 
      saveRDS(here::here(outdir,
                         "ft",
                         glue::glue("{file.name}.rds")))
    
    ft %>% 
      flextable::save_as_docx(path = here::here(outdir,
                                         "docx",
                                         glue::glue("{file.name}.docx")))
  }
  
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
        n.study = "No. of Studies",
        depression = "No. of Studies in Depressed Populations",
        anxiety = "No. of Studies in Anxious Populations",
        deficiency = "No. of Studies in Deficient Populations",
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
  
  # ---- flextable demographic helpers ----
  .parse_num <- function(x) {
    suppressWarnings(as.numeric(gsub(",", "", trimws(as.character(x)))))
  }
  
  .fmt_count <- function(x) {
    x_chr <- as.character(x)
    out <- x_chr
    miss <- is.na(x_chr) | trimws(x_chr) == ""
    val <- .parse_num(x_chr)
    idx <- !miss & !is.na(val)
    out[idx] <- format(round(val[idx]), big.mark = ",", scientific = FALSE, trim = TRUE)
    out[miss] <- ""
    out
  }
  
  .fmt_pct <- function(x, digits = 1) {
    x_chr <- as.character(x)
    out <- x_chr
    miss <- is.na(x_chr) | trimws(x_chr) == ""
    val <- .parse_num(gsub("%", "", x_chr))
    idx <- !miss & !is.na(val)
    if (any(idx)) {
      adj <- ifelse(abs(val[idx]) <= 1, val[idx] * 100, val[idx])
      out[idx] <- paste0(formatC(adj, format = "f", digits = digits), "%")
    }
    out[miss] <- ""
    out
  }
  
  .fmt_count_pct <- function(x, pct_digits = 1) {
    x_chr <- as.character(x)
    out <- x_chr
    miss <- is.na(x_chr) | trimws(x_chr) == ""
    
    m <- regexec(
      "^\\s*([0-9,]+(?:\\.[0-9]+)?)\\s*\\(([-+]?[0-9]+(?:\\.[0-9]+)?)\\s*%\\)\\s*$",
      x_chr
    )
    reg <- regmatches(x_chr, m)
    idx <- lengths(reg) == 3
    if (any(idx)) {
      n_vals <- suppressWarnings(as.numeric(gsub(",", "", vapply(reg[idx], `[`, character(1), 2))))
      p_vals <- suppressWarnings(as.numeric(vapply(reg[idx], `[`, character(1), 3)))
      out[idx] <- paste0(
        format(round(n_vals), big.mark = ",", scientific = FALSE, trim = TRUE),
        " (",
        formatC(p_vals, format = "f", digits = pct_digits),
        "%)"
      )
    }
    
    # fallback: numeric or percent-only values
    rem <- !miss & !idx
    rem_num <- .parse_num(gsub("%", "", x_chr))
    rem_idx <- rem & !is.na(rem_num)
    if (any(rem_idx)) {
      adj <- ifelse(abs(rem_num[rem_idx]) <= 1, rem_num[rem_idx] * 100, rem_num[rem_idx])
      out[rem_idx] <- paste0(formatC(adj, format = "f", digits = pct_digits), "%")
    }
    
    out[miss] <- ""
    out
  }
  
  .fmt_mean_sd <- function(x, digits = 1) {
    x_chr <- as.character(x)
    out <- x_chr
    miss <- is.na(x_chr) | trimws(x_chr) == ""
    
    m <- regexec(
      "^\\s*([-+]?[0-9,]+(?:\\.[0-9]+)?)\\s*\\(([-+]?[0-9,]+(?:\\.[0-9]+)?)\\)\\s*$",
      x_chr
    )
    reg <- regmatches(x_chr, m)
    idx <- lengths(reg) == 3
    if (any(idx)) {
      mean_vals <- suppressWarnings(as.numeric(gsub(",", "", vapply(reg[idx], `[`, character(1), 2))))
      sd_vals <- suppressWarnings(as.numeric(gsub(",", "", vapply(reg[idx], `[`, character(1), 3))))
      out[idx] <- paste0(
        formatC(mean_vals, format = "f", digits = digits),
        " (",
        formatC(sd_vals, format = "f", digits = digits),
        ")"
      )
    }
    
    rem <- !miss & !idx
    rem_num <- .parse_num(x_chr)
    rem_idx <- rem & !is.na(rem_num)
    if (any(rem_idx)) {
      out[rem_idx] <- formatC(rem_num[rem_idx], format = "f", digits = digits)
    }
    
    out[miss] <- ""
    out
  }
  
  .format_demo_df_for_flex <- function(df) {
    count_cols <- intersect(c("n", "n.study"), names(df))
    pct_cols <- intersect(c("male", "depression", "anxiety", "deficiency",
                            "bias.low", "bias.unclear", "bias.high"), names(df))
    mean_sd_cols <- intersect(c("age", "duration"), names(df))
    
    if (length(count_cols) > 0) {
      df <- df %>% mutate(across(all_of(count_cols), ~ .fmt_count(.x)))
    }
    if (length(pct_cols) > 0) {
      df <- df %>% mutate(across(all_of(pct_cols), ~ .fmt_count_pct(.x, pct_digits = 1)))
    }
    if (length(mean_sd_cols) > 0) {
      df <- df %>% mutate(across(all_of(mean_sd_cols), ~ .fmt_mean_sd(.x, digits = 1)))
    }
    
    df %>% mutate(across(everything(), ~ ifelse(is.na(.x), "", as.character(.x))))
  }
  
  study_table_flex <- function(x) {
    is_flex <- inherits(x, "flextable")
    df <- if (is_flex) x$body$dataset else x
    
    bias_pres <- "bias.low" %in% colnames(df)
    
    df <- df %>%
      mutate(intervention = stringr::str_to_title(as.character(intervention))) %>%
      arrange(intervention) %>%
      .format_demo_df_for_flex()
    
    ft <- if (is_flex) {
      x$body$dataset <- df
      x
    } else {
      flextable::flextable(df)
    }
    
    ft <- ft %>%
      flextable::set_header_labels(
        intervention = "Intervention",
        n.study = "Number\nof\nStudies",
        depression = "No. of Studies\nin Depressed Populations",
        anxiety = "No. of Studies\nin Anxious Populations",
        deficiency = "No. of Studies\nin Deficient Populations",
        duration = "Duration (weeks)",
        dose = "Dose*",
        bias.low = "Low",
        bias.unclear = "Unclear",
        bias.high = "High"
      ) %>%
      flextable::set_caption(
        "Study-Level Demographics. Values are n (%) or mean (SD) where appropriate. *Dose is reported as mean [min, max] across included studies."
      ) %>%
      flextable::align(j = colnames(df), align = "center", part = "header") %>%
      flextable::align(j = "intervention", align = "left", part = "body") %>%
      flextable::align(j = setdiff(colnames(df), "intervention"), align = "right", part = "body") %>%
      flextable::empty_blanks() %>%
      flextable::autofit()
    
    if (bias_pres && all(c("bias.low", "bias.unclear", "bias.high") %in% colnames(df))) {
      cols_before_bias <- which(colnames(df) %in% "bias.low") - 1
      vals <- if (cols_before_bias > 0) c("", "Overall Study Bias (%)") else c("Overall Study Bias (%)")
      widths <- if (cols_before_bias > 0) c(cols_before_bias, 3) else c(3)
      
      ft <- ft %>%
        flextable::add_header_row(values = vals, colwidths = widths, top = TRUE) %>%
        flextable::align(align = "center", part = "header")
    }
    
    if (is_flex) {
      ft <- overall_labeling_flex(ft)
    }
    
    ft
  }
  
  pop_table_flex <- function(x) {
    is_flex <- inherits(x, "flextable")
    df <- if (is_flex) x$body$dataset else x
    
    df <- df %>%
      mutate(intervention = stringr::str_to_title(as.character(intervention))) %>%
      arrange(intervention) %>%
      .format_demo_df_for_flex()
    
    ft <- if (is_flex) {
      x$body$dataset <- df
      x
    } else {
      flextable::flextable(df)
    }
    
    ft <- ft %>%
      flextable::set_header_labels(
        intervention = "Intervention",
        n = "N",
        age = "Age",
        male = "Male"
      ) %>%
      flextable::set_caption(
        "Population-Level Demographics. Values are n (%) or mean (SD) where appropriate."
      ) %>%
      flextable::align(j = colnames(df), align = "center", part = "header") %>%
      flextable::align(j = "intervention", align = "left", part = "body") %>%
      flextable::align(j = setdiff(colnames(df), "intervention"), align = "right", part = "body") %>%
      flextable::empty_blanks() %>%
      flextable::autofit()
    
    if (is_flex) {
      ft <- overall_labeling_flex(ft)
    }
    
    ft
  }
  
  overall_labeling_flex <- function(x) {
    stopifnot(inherits(x, "flextable"))
    
    x_cols <- colnames(x$body$dataset)
    bias_pres <- "bias.low" %in% x_cols
    
    study_cols <- if (bias_pres) {
      c("n.study", "depression", "anxiety", "deficiency", "duration", "dose",
        "bias.low", "bias.unclear", "bias.high")
    } else {
      c("n.study", "depression", "anxiety", "deficiency", "duration", "dose")
    }
    
    pop_cols <- c("n", "age", "male")
    study_cols <- study_cols[study_cols %in% x_cols]
    pop_cols <- pop_cols[pop_cols %in% x_cols]
    
    hdr <- rep("", length(x_cols))
    names(hdr) <- x_cols
    if (length(pop_cols) > 0) hdr[pop_cols] <- "Population-Level Characteristics"
    if (length(study_cols) > 0) hdr[study_cols] <- "Study-Level Characteristics"
    
    bias_comment <- if (!bias_pres) {
      " in low bias studies"
    } else {
      ""
    }
    
    x %>%
      flextable::add_header_row(values = unname(hdr), colwidths = rep(1, length(hdr)), top = TRUE) %>%
      flextable::merge_h(part = "header") %>%
      flextable::align(align = "center", part = "header") %>%
      flextable::set_caption(
        glue::glue("Baseline demographics{bias_comment}. Values are n (%) or mean (SD) where appropriate. *Dose is reported as mean [min, max] across included studies.")
      )
  }
}

#### Baseline Demographic Data ####
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
    filter(study %in% main) %>% 
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
    n = b.n,
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
    n = b.n,
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
    filter(study %in% main) %>% 
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
      here::here("outputs","overall","tables", "gt",
                 "study_demo_low_bias_table.rds")
    )
  
  study_demo_low %>% 
    study_table_flex() %>% 
    saveRDS(
      here::here("outputs","overall","tables", "ft",
                 "study_demo_low_bias_table.rds")
    )
  
  study_demo %>% 
    study_table() %>% 
    saveRDS(
      here::here("outputs","overall","tables", "gt",
                 "study_demo_table.rds")
    )
  
  study_demo %>% 
    study_table_flex() %>% 
    saveRDS(
      here::here("outputs","overall","tables", "ft",
                 "study_demo_table.rds")
    )
  
  pop_low_gt <- pop_demo_low %>% 
    pop_table()
  
  pop_low_ft <- pop_demo_low %>% 
    pop_table_flex() 
  
  pop_gt <- pop_demo %>% 
    pop_table()
  
  pop_ft <- pop_demo %>% 
    pop_table_flex()
  
  pop_low_gt %>% 
    saveRDS(
      here::here("outputs","overall","tables", "gt",
                 "pop_demo_low_bias_table.rds")
    )
  
  pop_low_ft %>% 
    saveRDS(
      here::here("outputs","overall","tables", "ft",
                 "pop_demo_low_bias_table.rds")
    )
  
  pop_gt %>% 
    saveRDS(
      here::here("outputs","overall","tables", "gt",
                 "pop_demo_table.rds")
    )
  
  pop_ft %>% 
    saveRDS(
      here::here("outputs","overall","tables","ft",
                 "pop_demo_table.rds")
    )
  
  combine_tab <- dplyr::left_join(
    pop_demo,
    study_demo  %>% mutate(dose = stringr::str_replace_all(dose, "], ", "]\n"),
                           dose = stringr::str_replace_all(dose, "g, ", "g,\n")),
    by = "intervention") %>% 
    pop_table() %>% 
    study_table()
  
  combine_tab_ft <- dplyr::left_join(
    pop_demo,
    study_demo  %>% mutate(dose = stringr::str_replace_all(dose, "], ", "]\n"),
                           dose = stringr::str_replace_all(dose, "g, ", "g,\n")),
    by = "intervention") %>% 
    pop_table_flex() %>% 
    study_table_flex()
  
  combine_tab_low <- dplyr::left_join(
    pop_demo_low,
    study_demo_low,
    by = "intervention") %>% 
    pop_table() %>% 
    study_table()
  
  
  combine_tab_low_ft <- dplyr::left_join(
    pop_demo_low,
    study_demo_low %>% mutate(dose = stringr::str_replace_all(dose, "], ", "]\n"),
                              dose = stringr::str_replace_all(dose, "g, ", "g,\n")),
    by = "intervention") %>% 
    pop_table_flex() %>% 
    study_table_flex()
  
  combine_tab %>% 
    saveRDS(
      here::here("outputs","overall","tables", "gt",
                 "combined_demo_table.rds")
    )
  
  combine_tab_ft %>% 
    saveRDS(
      here::here("outputs","overall","tables", "ft",
                 "combined_demo_table.rds")
    )
  
  combine_tab_ft %>% 
    set_table_properties(layout = "fixed") %>%                # good to do before scaling
    fit_to_width(max_width = 10) %>%    # because margins will be smaller
    fontsize(size = 6, part = "all") %>% 
    valign(j = "dose", valign = "top") %>% 
    flextable::hline(border = officer::fp_border(color = "gray70", width = 0.5), part = "body") %>% 
    flextable::vline(
      j = c("intervention","male"),
      border = officer::fp_border(color = "gray40", width = 1)
    ) %>% 
    flextable::save_as_docx(
      path = here::here("outputs","overall","tables", "docx",
                        "combined_demo_table.docx"),
      pr_section = officer::prop_section(
        page_size = officer::page_size(width = 11, height = 8.5, orient = "landscape"),
        page_margins = officer::page_mar(left = 0.5, right = 0.5, top = 0.5, bottom = 0.5),
        type = "continuous"
      )
    )
  
  combine_tab_low %>% 
    saveRDS(
      here::here("outputs","overall","tables", "gt",
                 "combined_demo_low_bias_table.rds")
    )
  
  combine_tab_low_ft %>% 
    saveRDS(
      here::here("outputs","overall","tables", "ft",
                 "combined_demo_low_bias_table.rds")
    )
  
  
  combine_tab_low_ft %>% 
    set_table_properties(layout = "fixed") %>%                # good to do before scaling
    fit_to_width(max_width = 10) %>%    # because margins will be smaller
    fontsize(size = 7, part = "all") %>% 
    valign(j = "dose", valign = "top") %>% 
    flextable::hline(border = officer::fp_border(color = "gray70", width = 0.5), part = "body") %>% 
    flextable::vline(
      j = c("intervention","male"),
      border = officer::fp_border(color = "gray40", width = 1)
    ) %>% 
    flextable::save_as_docx(
      path = here::here("outputs","overall","tables","docx","combined_demo_low_bias_table.docx"),
      pr_section = officer::prop_section(
        page_size = officer::page_size(width = 11, height = 8.5, orient = "landscape"),
        page_margins = officer::page_mar(left = 0.5, right = 0.5, top = 0.5, bottom = 0.5),
        type = "continuous"
      )
    )
  
}

# create table by comparison to assess for overlap
{
  baseline.vars.comp <- data %>% 
    mutate(control = gsub("z\\.","",as.character(control)) %>% factor()) %>% 
    mutate(control = gsub("_"," ",as.character(control)) %>% factor()) %>% 
    mutate(control = as.character(control) %>% stringr::str_to_title()) %>% 
    mutate(control = gsub(",", ", ", as.character(control)) %>% factor()) %>% 
    mutate(intervention = as.character(intervention) %>% stringr::str_to_title() ) %>% 
    mutate(treatment = intervention) %>% 
    mutate(intervention = glue::glue("{intervention}:{control}")) %>% 
    mutate(intervention = gsub("z\\.","",as.character(intervention)) %>% factor()) %>% 
    mutate(intervention = gsub("_"," ",as.character(intervention)) %>% factor()) %>%
    mutate(intervention = gsub(" - ",":",as.character(intervention)) %>% factor()) %>%
    filter(!is.na(contrast)) %>% 
    select(c("study","intervention", !!!baseline.cols))
  
  # demographics tab overall
  demo_tab.comp <- baseline.vars.comp %>% 
    demo_create() %>% 
    tidyr::separate(
      intervention,
      into = c("intervention", "Comparator"),
      sep = ":"
    )
  
  # demographics tab low bias
  demo_tab_low.comp <- baseline.vars.comp %>% 
    filter(bias.overall == "low") %>% 
    filter(study %in% main) %>% 
    demo_create() %>% 
    tidyr::separate(
      intervention,
      into = c("intervention", "Comparator"),
      sep = ":"
    )
  
  study_demo.comp <- demo_tab.comp %>% 
    select(
      intervention,
      Comparator,
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
      Comparator,
      n.study,
      depression,
      anxiety,
      deficiency,
      duration,
      dose
    )
  
  
  pop_demo.comp <- demo_tab.comp %>% select(
    intervention,
    Comparator,
    n = b.n,
    age,
    male
  )
  
  pop_demo_low.comp <- demo_tab_low.comp %>% 
  select(
    intervention,
    Comparator,
    n = b.n,
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
    filter(study %in% main) %>% 
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
    saveRDS(
      here::here("outputs","overall","tables", "gt",
                 "study_demo_low_bias_table_comparison.rds")
    )
  
  study_demo_low.comp %>% 
    study_table_flex() %>% 
    saveRDS(
      here::here("outputs","overall","tables", "ft",
                 "study_demo_low_bias_table_comparison.rds")
    )
  
  study_demo.comp %>% 
    study_table() %>% 
    saveRDS(
      here::here("outputs","overall","tables", "gt",
                 "study_demo_table_comparison.rds")
    )
  
  study_demo.comp %>% 
    study_table_flex() %>% 
    saveRDS(
      here::here("outputs","overall","tables", "ft",
                 "study_demo_table_comparison.rds")
    )
  
  pop_demo_low.comp %>% 
    pop_table()  %>% 
    saveRDS(
      here::here("outputs","overall","tables", "gt",
                 "pop_demo_low_bias_table_comparison.rds")
    )
  
  pop_demo_low.comp %>% 
    pop_table_flex()  %>% 
    saveRDS(
      here::here("outputs","overall","tables", "ft",
                 "pop_demo_low_bias_table_comparison.rds")
    )
  
  pop_demo_low.comp %>% 
    pop_table_flex()  %>% 
    flextable::save_as_docx(
      path = here::here("outputs","overall","tables", "docx",
                 "pop_demo_low_bias_table_comparison.docx")
      
    )
  
  pop_demo.comp %>% 
    pop_table()  %>% 
    saveRDS(
      here::here("outputs","overall","tables", "gt",
                 "pop_demo_table_comparison.rds")
    )
  
  pop_demo.comp %>% 
    pop_table_flex()  %>% 
    saveRDS(
      here::here("outputs","overall","tables", "ft",
                 "pop_demo_table_comparison.rds")
    )
  
  pop_demo.comp %>% 
    pop_table_flex()  %>% 
    flextable::save_as_docx(
      path = here::here("outputs","overall","tables", "docx",
                 "pop_demo_table_comparison.docx")
    )
  
  combine_tab.comp <- dplyr::left_join(
    pop_demo.comp,
    study_demo.comp %>% mutate(dose = stringr::str_replace_all(dose, "], ", "]\n"),
                               dose = stringr::str_replace_all(dose, "g, ", "g,\n")),
    by = c("intervention", "Comparator")
    ) %>% 
    pop_table() %>% 
    study_table() 
  
  combine_tab_low.comp <- dplyr::left_join(
    pop_demo_low.comp,
    study_demo_low.comp %>% mutate(dose = stringr::str_replace_all(dose, "], ", "]\n"),
                                   dose = stringr::str_replace_all(dose, "g, ", "g,\n")),
    by = c("intervention", "Comparator") ) %>% 
    pop_table() %>% 
    study_table() 
  
  combine_tab.comp_flex <- dplyr::left_join(
    pop_demo.comp,
    study_demo.comp %>% mutate(dose = stringr::str_replace_all(dose, "], ", "]\n"),
                               dose = stringr::str_replace_all(dose, "g, ", "g,\n")),
    by = c("intervention", "Comparator") ) %>% 
    pop_table_flex() %>% 
    study_table_flex() 
  
  combine_tab_low.comp_flex <- dplyr::left_join(
    pop_demo_low.comp,
    study_demo_low.comp %>% mutate(dose = stringr::str_replace_all(dose, "], ", "]\n"),
                                   dose = stringr::str_replace_all(dose, "g, ", "g,\n")),
    by = c("intervention", "Comparator")) %>% 
    pop_table_flex() %>% 
    study_table_flex() 
  
  combine_tab.comp %>% 
    saveRDS(
      here::here("outputs","overall","tables", "gt",
                 "combined_demo_table_comparison.rds")
    )
  
  combine_tab_low.comp %>% 
    saveRDS(
      here::here("outputs","overall","tables", "gt",
                 "combined_demo_low_bias_table_comparison.rds")
    )
  
  combine_tab.comp_flex %>% 
    saveRDS(
      here::here("outputs","overall","tables", "ft",
                 "combined_demo_table_comparison.rds")
    )
  
  combine_tab_low.comp_flex %>% 
    saveRDS(
      here::here("outputs","overall","tables", "ft",
                 "combined_demo_low_bias_table_comparison.rds")
    )
  
  combine_tab.comp_flex %>% 
    set_table_properties(layout = "fixed") %>%                # good to do before scaling
    fit_to_width(max_width = 10) %>%    # because margins will be smaller
    fontsize(size = 6, part = "all") %>% 
    valign(j = "dose", valign = "top") %>% 
    flextable::hline(border = officer::fp_border(color = "gray70", width = 0.5), part = "body") %>% 
    flextable::vline(
      j = c("Comparator","male"),
      border = officer::fp_border(color = "gray40", width = 1)
    ) %>% 
    flextable::save_as_docx(
      path = here::here("outputs","overall","tables", "docx",
                 "combined_demo_table_comparison.docx"),
      pr_section = officer::prop_section(
        page_size = officer::page_size(width = 11, height = 8.5, orient = "landscape"),
        page_margins = officer::page_mar(left = 0.5, right = 0.5, top = 0.5, bottom = 0.5),
        type = "continuous"
      )
    )
  
  combine_tab_low.comp_flex %>% 
    set_table_properties(layout = "fixed") %>%                # good to do before scaling
    fit_to_width(max_width = 10) %>%    # because margins will be smaller
    fontsize(size = 6, part = "all") %>% 
    valign(j = "dose", valign = "top") %>% 
    flextable::hline(border = officer::fp_border(color = "gray70", width = 0.5), part = "body") %>% 
    flextable::vline(
      j = c("Comparator","male"),
      border = officer::fp_border(color = "gray40", width = 1)
    ) %>% 
    flextable::save_as_docx(
      path = here::here("outputs","overall","tables", "docx",
                 "combined_demo_low_bias_table_comparison.docx"),
        pr_section = officer::prop_section(
          page_size = officer::page_size(width = 11, height = 8.5, orient = "landscape"),
          page_margins = officer::page_mar(left = 0.5, right = 0.5, top = 0.5, bottom = 0.5),
          type = "continuous"
        )
    )
}

# create table by study for appendix
{
  study_demo_all <- data %>% 
    distinct(study, intervention,
             final.N, duration, depressed,anxiety,
             vit.def, bias.overall,
             .keep_all = TRUE) %>% 
    mutate(depressed = ifelse(depressed == 1, "Yes", "No"),
           anxiety = ifelse(anxiety == 1, "Yes", "No")) %>%
    mutate(
      across(
        any_of("bias.overall"),
        ~ forcats::fct_recode(.x, 
                              "Low" = "low",
                              "Unclear" = "some concerns",
                              "High" = "high")
      )
    ) %>% 
    select(Study = study, Intervention = intervention,
           N = final.N, `Duration in weeks` = duration, 
           Depressed = depressed,
           Anxious = anxiety,
           Deficiency = vit.def, Bias = bias.overall) %>% 
    arrange(tolower(Study)) 
  
  study_demo_all %>%
    saveRDS(
      here::here("outputs","overall","tables", "raw",
                 "study_level_demo_table.rds")
    )
  
  study_demo_all %>%
    gt::gt(groupname_col = "Study") %>% 
    gt::fmt_number(decimals = 1) %>% 
    gt::sub_missing() %>% 
    saveRDS(
      here::here("outputs","overall","tables", "gt",
                 "study_level_demo_table.rds")
    )
  
  study_demo_all %>%
    flextable() %>% 
    flextable::merge_v(j = c("Study", "Intervention")) %>%
    flextable::colformat_double(digits = 1) %>% 
    flextable::colformat_char(na_str = "—") %>%
    flextable::autofit() %>% 
    flextable::fit_to_width(max_width = 6.5) %>%    # because margins will be smaller
    flextable::fontsize(size = 6, part = "all") %>% 
    flextable::hline(border = officer::fp_border(color = "gray70", width = 0.5), part = "body") %>%
    saveRDS(
      here::here("outputs","overall","tables", "ft",
                 "study_level_demo_table.rds")
    )
  
  study_demo_all %>%
    flextable::flextable() %>% 
    flextable::merge_v(j = c("Study", "Intervention")) %>%
    flextable::colformat_double(digits = 1) %>% 
    flextable::colformat_char(na_str = "—") %>%
    flextable::autofit() %>% 
    flextable::fit_to_width(max_width = 6.5) %>%    # because margins will be smaller
    flextable::fontsize(size = 6, part = "all") %>% 
    flextable::hline(border = officer::fp_border(color = "gray70", width = 0.5), part = "body") %>% 
    flextable::save_as_docx(
      path = here::here("outputs","overall","tables", "docx",
                 "study_level_demo_table.docx"),
      pr_section = officer::prop_section(
        page_size = officer::page_size(width = 8.5, height = 11, orient = "portrait"),
        page_margins = officer::page_mar(left = 0.5, right = 0.5, top = 0.5, bottom = 0.5),
        type = "continuous"
      )
    )
}

{
  study_demo <- vit$get_vitamin_data(outcome = NULL, simple_analysis = TRUE,
                       include_full_bias = TRUE) %>% 
    distinct(study, intervention,
             final.N, duration, depressed,anxiety,
             vit.def, bias.overall,
             .keep_all = TRUE) %>% 
    filter(bias.overall == "low") %>%
    filter(study %in% main) %>% 
    mutate(depressed = ifelse(depressed == 1, "Yes", "No"),
           anxiety = ifelse(anxiety == 1, "Yes", "No")) %>%
    select(Study = study, Intervention = intervention,
           N = final.N, `Duration in weeks` = duration, 
           Depressed = depressed,
           Anxious = anxiety,
           Deficiency = vit.def, Bias = bias.overall) %>% 
    arrange(tolower(Study))
  
  study_demo %>% 
    saveRDS(
      here::here("outputs","overall","tables", "raw",
                 "study_demo_low_bias.rds")
    )
  
  study_demo %>% 
    gt::gt(groupname_col = "Study") %>% 
    gt::fmt_number(decimals = 0) %>% 
    gt::sub_missing() %>% 
    saveRDS(
      here::here("outputs","overall","tables", "gt",
                 "study_level_demo_table_low.rds")
    )
    
  study_demo %>% 
    # flextable::as_grouped_data(groups = "Study") %>%
    flextable::flextable() %>% 
    flextable::merge_v(j = c("Study", "Intervention")) %>%
    flextable::colformat_double(digits = 0) %>% 
    flextable::colformat_char(na_str = "—") %>%
    flextable::autofit() %>% 
    flextable::fit_to_width(max_width = 7.5) %>%
    saveRDS(
      here::here("outputs","overall","tables", "ft",
                 "study_level_demo_table_low.rds")
    )
}

#### adverse events table ####
fit_adv <- here::here("outputs", "saved_models", "adverse_events.rds") %>% readRDS()
newdata <- vit$nma_newdata_for_summary(fit_adv)
newdata$final.N <- 1
newdata$.study <- fit_adv$data$.study[1]
newdata$.obs_re <- as.factor(newdata$.obs_re)
adverse_summ <- vit$summary_brms_nma(fit_adv, newdata = newdata)

adverse_summ %>% 
  mutate(value = exp(value)) %>% 
  summary() %>% 
  filter(.observed) %>% 
  select(-.observed) %>% 
  save_mult_formats(outdir = here::here("outputs", "overall", "tables"),
                    file.name = glue::glue("adverse_outcomes_all_data_table"))

fit_adv_low <- here::here("outputs", "saved_models", "adverse_events_low_bias.rds") %>% readRDS()
newdata_low <- vit$nma_newdata_for_summary(fit_adv_low)
newdata_low$final.N <- 1
newdata_low$.study <- fit_adv_low$data$.study[1]
newdata_low$.obs_re <- as.factor(newdata_low$.obs_re)
adverse_summ_low <- vit$summary_brms_nma(fit_adv_low, newdata = newdata_low)

adverse_summ_low %>% 
  mutate(value = exp(value)) %>% 
  summary() %>% 
  filter(.observed) %>% 
  select(-.observed) %>% 
  save_mult_formats(outdir = here::here("outputs", "overall", "tables"),
                    file.name = glue::glue("adverse_outcomes_low_bias_table"))


fit_adv_low_ssri <- here::here("outputs", "saved_models", "adverse_events_low_bias_ssri.rds") %>% readRDS()
newdata_ssri <- vit$nma_newdata_for_summary(fit_adv_low_ssri)
newdata_ssri$final.N <- 1
newdata_ssri$.study <- fit_adv_low_ssri$data$.study[1]
newdata_ssri$.obs_re <- as.factor(newdata_ssri$.obs_re)
adverse_summ_ssri <- vit$summary_brms_nma(fit_adv_low_ssri, newdata = newdata_ssri)

adverse_summ_ssri %>% 
  mutate(value = exp(value)) %>% 
  summary() %>% 
  filter(.observed) %>% 
  select(-.observed) %>% 
  save_mult_formats(outdir = here::here("outputs", "overall", "tables"),
                    file.name = glue::glue("adverse_outcomes_low_bias_ssri_table"))


#### create vitamin level tables ####
{
  level_low <- readRDS(here::here("outputs", "saved_models", "vitamin_levels_low_bias.rds"))
  level_ssri_low <- readRDS(here::here("outputs", "saved_models", "vitamin_levels_low_bias_ssri.rds"))
  level_ssri <- readRDS(here::here("outputs", "saved_models", "vitamin_levels_all_ssri.rds"))
  level_all <- readRDS(here::here("outputs", "saved_models", "vitamin_levels_all.rds"))
  
  level_tab <- function(lev, nm) {
    llt <- vector("list", length(lev))
    
    for(i in seq_along(lev)) {
      if(is.null(lev[[i]])) next
      
      llt[[i]] <- lev[[i]] %>% 
        vit$summary_brms_nma(keep = c(".trt")) %>% summary()
      llt[[i]]$micro <- names(lev)[i] %>% stringr::str_remove("end\\.") %>% stringr::str_to_title()
      llt[[i]]$unit <- lev[[i]]$unit
      
    }
    fmt <- function(v, d) formatC(v, format = "f", digits = 2)
    
    llt_out <- llt %>% bind_rows() %>% 
      select(-.observed) %>% 
      mutate(micro = ifelse(micro == "Rbc.b9", "R.B.C. B9", micro)) %>% 
      select(Intervention = .trt, Micronutrient = micro, Unit = unit, everything()) %>% 
      select(-c(S.E., `q50%`)) %>%
      mutate(M.D. = fmt(Estimate, digits)) %>% 
      mutate(ci = glue::glue("({fmt(`q2.5%`, digits)}, {fmt(`q97.5%`, digits)})")) %>% 
      select(-c(Estimate, `q2.5%`, `q97.5%`)) %>% 
      rename("95% Credible Interval" = ci) %>% 
      mutate(Intervention = as.character(Intervention))
    
    data.table::setorder(llt_out, Intervention, Micronutrient)
    
    llt_out %>% saveRDS(
      here::here("outputs","overall","tables", "raw",
                 glue::glue("vitamin_levels_{nm}_tab.rds"))
    )
    
    llt_out %>% gt() %>%
      fmt_number() %>% 
      saveRDS(
        here::here("outputs","overall","tables", "gt",
                   glue::glue("vitamin_levels_{nm}_tab.rds"))
      )
    
    llt_out %>% flextable::flextable() %>% 
      flextable::autofit() %>% 
      saveRDS(
        here::here("outputs","overall","tables", "gt",
                   glue::glue("vitamin_levels_{nm}_tab.rds"))
      )
    
    llt_out %>% flextable::flextable() %>% 
      flextable::autofit() %>% 
      flextable::save_as_docx(
        path = here::here("outputs","overall","tables", "docx",
                          glue::glue("vitamin_levels_{nm}_tab.docx"))
      )
    
  }
  
  level_tab(level_low, "low_bias")
  level_tab(level_ssri_low, "low_bias_ssri")
  level_tab(level_ssri, "all_bias_ssri")
  level_tab(level_all, "all_data")
  
}

#### Regression tables by outcome ####
outcomes <- c("depression","anxiety")

if (rlang::is_interactive()) {
  outcomes <- outcome <- c("depression")
} 

for (outcome in outcomes) {
  sensitivity_outdir <- here::here("outputs", "saved_models", outcome,
                                   "sensitivity")
  
  output_dir <- here::here("outputs",outcome,"tables")
  outcome_label <- switch(outcome,
                          "depression" = "depressed",
                          "anxiety" = "anxiety")
  
  not_outcome_label <- switch(outcome,
                              "anxiety" = "depressed",
                              "depression" = "anxiety")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  if (! dir.exists(tmp.dir <- here::here(output_dir, "raw"))) {
    dir.create(tmp.dir, recursive = TRUE)
  }
  if (! dir.exists(tmp.dir <- here::here(output_dir, "gt"))) {
    dir.create(tmp.dir, recursive = TRUE)
  }
  if (! dir.exists(tmp.dir <- here::here(output_dir, "ft"))) {
    dir.create(tmp.dir, recursive = TRUE)
  }
  if (! dir.exists(tmp.dir <- here::here(output_dir, "docx"))) {
    dir.create(tmp.dir, recursive = TRUE)
  }
  # load models
  cat("  Loading models for", outcome, "...\n")
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
      here::here("outputs", "saved_models", outcome, "sensitivity",
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
    
    fit_pcts_approx <- readRDS(
      here::here("outputs", "saved_models", outcome, "sensitivity",
                 glue::glue("model_", 
                            outcome, 
                            "_low_bias_data_by_pct_def_and_dep.rds")
    ))
    
    bias_reg <- readRDS(here::here("outputs", 
                                  "saved_models", 
                                  "bias_regression_fit_overall.rds"))
    
  }
  
  cat("  Creating summary tables for", outcome, "...\n")
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
  
  cat("  Saving basic tables for", outcome, "...\n")
  # save basic results tables
  {
    
    low_res %>% 
      save_mult_formats(outdir = output_dir,
                        file.name = glue::glue("nma_{outcome}_low_bias_table"))
      
    some_res %>% 
      save_mult_formats(outdir = output_dir,
                        file.name = glue::glue("nma_{outcome}_some_and_low_bias_table"))
    
    high_res %>% 
      save_mult_formats(outdir = output_dir,
                        file.name = glue::glue("nma_{outcome}_all_data_table"))
    
    dose_res %>% 
      mutate(
        .dose = factor(
          .dose,
          levels = levels(.dose),
          labels = stringr::str_wrap(levels(.dose), width = 20)
        )
      ) %>%
      save_mult_formats(outdir = output_dir,
                        file.name = 
                              glue::glue("nma_{outcome}_dose_response_table_low_bias"),
                              group = "Intervention")
    
    def_res_low %>%
      save_mult_formats(outdir = output_dir,
                        file.name =
                              glue::glue("nma_{outcome}_deficiency_table_low_bias"),
                        group = "Intervention")
    
    time_res %>% 
      save_mult_formats(outdir = output_dir,
                        file.name =
                          glue::glue("nma_{outcome}_time_response_table_low_bias"),
                        group = "Intervention")
    
    b_def_t <- bias_reg_res %>% 
      save_mult_formats(outdir = output_dir,
                       file.name = glue::glue("nma_{outcome}_meta_regression_table"),
                       group = "Intervention")
    
    def_res_low %>% 
      group_by(.trt) %>%
      tidyr::complete(
        !!outcome_label := c("yes", "no"),
        def = c("yes", "no", "NA")
      ) %>%
      # gt_fmt(groupname_col = "Intervention") %>% 
      save_mult_formats(
        outdir = output_dir,
        file.name =glue::glue("nma_{outcome}_deficiency_table_low_bias_complete"),
        group = "Intervention")
    
    def_res %>% 
      save_mult_formats(outdir = output_dir,
                        file.name =
                          glue::glue("nma_{outcome}_deficiency_table_all_data"),
                        group = "Intervention")
    
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
      save_mult_formats(outdir = output_dir,
                        file.name =
                          glue::glue("nma_{outcome}_pct_table_low_bias"),
                        group = "Intervention")
  }
  
  cat("  Loading loo models for", outcome, "...\n")
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
    
    # loo_meta_reg <- readRDS(
    #   here::here(sensitivity_outdir,
    #              "bias_regression_loo.rds")
    # )
    
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
  
  cat("  Creating loo tables for", outcome, "...\n")
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
  
  cat("  Saving loo tables for", outcome, "...\n")
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
      
    ) 
    
    
    loo_tab %>% gt() %>% 
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
        title = "Comparison of Model Fit by Likelihood Function. Differences in expected log point-wise predictive densities (ELPD) estimated by leave-one-out cross validation. The second through third columns list the difference in ELPDs separated by the bias of the included studies. Standard error of difference in parentheses.") %>% 
      saveRDS(
        here::here(output_dir, "gt",
                   "loo_compare_gaussian_tab.rds")
      )
    
    loo_tab_ft <- loo_tab %>% 
      mutate(
        low = ifelse(is.na(low), NA_character_, sprintf("%.3f (%.3f)", low, lowse)),
        `unclear/low` = ifelse(is.na(`unclear/low`), NA_character_, sprintf("%.3f (%.3f)", `unclear/low`, `unclear/lowse`)),
        `high/unclear/low` = ifelse(is.na(`high/unclear/low`), NA_character_, sprintf("%.3f (%.3f)", `high/unclear/low`, `high/unclear/lowse`))
      ) %>%
      select(-any_of(c("se_diff_low",
                       "se_diff_some",
                       "se_diff_high",
                       "lowse",
                       "unclear/lowse",
                       "high/unclear/lowse"))) %>%
      flextable::flextable() %>% 
      flextable::set_header_labels(
        Likelihood = "Likelihood",
        low = "Low Bias",
        'unclear/low' = "Unclear/Low Bias",
        'high/unclear/low' = "High/Unclear/Low Bias"
        # , diff = "Difference",
        # se   = "S.E."
      ) 
    
    loo_tab_ft %>%
      saveRDS(
        here::here(output_dir, "ft",
                   "loo_compare_gaussian_tab.rds")
      )
    
    loo_tab_ft %>%
      flextable::save_as_docx(
        path = here::here(output_dir, "docx",
                   "loo_compare_gaussian_tab.docx")
      )
    
    # other loo tabs
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
      tidyr::pivot_wider(names_from = bias, values_from = c(elpd_diff, se_diff))
    
    loo_def_tab %>% saveRDS(
      here::here(output_dir, "raw",
                 "loo_compare_model_tab.rds")
    )
    
    loo_def_tab %>%
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
        title = "Comparison of Model Fit by Bias. Differences in expected log point-wise predictive densities (ELPD) estimated by leave-one-out cross validation. The second through third columns list the difference in ELPDs separated by the bias of the included studies. Standard error of difference in parentheses.") %>% 
      saveRDS(
          here::here(output_dir, "gt",
                     "loo_compare_model_tab.rds")
        )
    
    loo_def_tab_ft <- loo_def_tab %>%
      mutate(
        elpd_diff_low = ifelse(is.na(elpd_diff_low), NA_character_, sprintf("%.3f (%.3f)", elpd_diff_low, se_diff_low)),
        elpd_diff_some = ifelse(is.na(elpd_diff_some), NA_character_, sprintf("%.3f (%.3f)", elpd_diff_some, se_diff_some)),
        elpd_diff_high = ifelse(is.na(elpd_diff_high), NA_character_, sprintf("%.3f (%.3f)", elpd_diff_high, se_diff_high))
      ) %>%
      select(-any_of(c("se_diff_low",
                       "se_diff_some",
                       "se_diff_high",
                       "lowse",
                       "unclear/lowse",
                       "high/unclear/lowse"))) %>%
      flextable() %>% 
      flextable::set_header_labels(
        Model = "Model",
        elpd_diff_low = "Low Bias",
        elpd_diff_some = "Unclear/Low Bias",
        elpd_diff_high = "High/Unclear/Low Bias"
      ) %>% 
      flextable::colformat_char(na_str = "—")
    
    
    loo_def_tab_ft %>%
      saveRDS(
        here::here(output_dir, "ft",
                   "loo_compare_model_tab.rds")
      )
    
    loo_def_tab_ft %>%
      flextable::save_as_docx(
        path = here::here(output_dir, "docx",
                   "loo_compare_model_tab.docx")
      )
    
  }
  
  # sucra and rank probs; leaving out given varying effects depending on population
  # {
  #   rank_probs_low <- low_t %>% 
  #     vit$rank_probs_vit(sucra = TRUE,
  #                        cumulative = FALSE) %>% 
  #     rename(
  #       Treatment = .trt,
  #       SUCRA = sucra
  #     ) %>%
  #     rename_with(
  #       ~ stringr::str_replace(.x, "p_rank\\[(\\d+)\\]", "P(rank = \\1)"),
  #       starts_with("p_rank")
  #     ) %>% 
  #     select(Treatment, SUCRA, everything()) %>% 
  #     gt() %>% 
  #     fmt_number()
  #   
  #   rank_probs_all <- high_t %>% 
  #     vit$rank_probs_vit(sucra = TRUE,
  #                        cumulative = FALSE) %>% 
  #     rename(
  #       Treatment = .trt,
  #       SUCRA = sucra
  #     ) %>%
  #     rename_with(
  #       ~ stringr::str_replace(.x, "p_rank\\[(\\d+)\\]", "P(rank = \\1)"),
  #       starts_with("p_rank")
  #     ) %>% 
  #     select(Treatment, SUCRA, everything()) %>% 
  #     gt() %>% 
  #     fmt_number()
  #   
  #   rank_probs_def <-  fit_def_low %>%
  #     vit$rank_probs_vit(sucra = TRUE,
  #                        cumulative = FALSE,
  #                        filter = rlang::quo(!!rlang::sym(outcome_label) == "yes" & def == "yes"),
  #                        keep = c(".trt","def",outcome_label)) %>% 
  #     mutate(
  #       across(
  #         starts_with("p_rank"),
  #         ~tidyr::replace_na(.x, 0)
  #       )
  #     ) %>% 
  #     rename(
  #       Treatment = .trt,
  #       SUCRA = sucra
  #     ) %>%
  #     rename_with(
  #       ~ stringr::str_replace(.x, "p_rank\\[(\\d+)\\]", "P(rank = \\1)"),
  #       starts_with("p_rank")
  #     ) %>% 
  #     select(Treatment, SUCRA, everything()) %>% 
  #     gt() %>% 
  #     fmt_number()
  #   
  #   rank_probs_def_not_def <-  fit_def_low %>%
  #     vit$rank_probs_vit(sucra = TRUE,
  #                        cumulative = FALSE,
  #                        filter = rlang::quo(!!rlang::sym(outcome_label) == "yes" & def == "no"),
  #                        keep = c(".trt","def",outcome_label)) %>% 
  #     mutate(
  #       across(
  #         starts_with("p_rank"),
  #         ~tidyr::replace_na(.x, 0)
  #       )
  #     ) %>% 
  #     rename(
  #       Treatment = .trt,
  #       SUCRA = sucra
  #     ) %>%
  #     rename_with(
  #       ~ stringr::str_replace(.x, "p_rank\\[(\\d+)\\]", "P(rank = \\1)"),
  #       starts_with("p_rank")
  #     ) %>% 
  #     select(Treatment, SUCRA, everything()) %>% 
  #     gt() %>% 
  #     fmt_number()
  #   
  #   
  #   rank_probs_def_unknown_def <-  fit_def_low %>%
  #     vit$rank_probs_vit(sucra = TRUE,
  #                        cumulative = FALSE,
  #                        filter = rlang::quo(!!rlang::sym(outcome_label) == "yes" & def == "NA"),
  #                        keep = c(".trt","def",outcome_label)) %>% 
  #     mutate(
  #       across(
  #         starts_with("p_rank"),
  #         ~tidyr::replace_na(.x, 0)
  #       )
  #     ) %>% 
  #     rename(
  #       Treatment = .trt,
  #       SUCRA = sucra
  #     ) %>%
  #     rename_with(
  #       ~ stringr::str_replace(.x, "p_rank\\[(\\d+)\\]", "P(rank = \\1)"),
  #       starts_with("p_rank")
  #     ) %>% 
  #     select(Treatment, SUCRA, everything()) %>% 
  #     gt() %>% 
  #     fmt_number()
  # }
  
  cat("  Creating contribution matrices for", outcome, "...\n")
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
      save_cm_mult_formats(outdir = output_dir,
                           file.name = glue::glue("contribution_matrix_low_bias"))
       
    cm_low_approx %>% 
         filter(trt2 == "Placebo") %>%  
         arrange(trt1) %>% 
         mutate(contribution = contribution * 100) %>% 
         tidyr::pivot_wider(id_cols = study,
                            names_from = comparison,
                            values_from = contribution) %>% 
         rename(Study = study) %>% 
      save_cm_mult_formats(outdir = output_dir,
                           file.name = glue::glue("approx_contribution_matrix_low_bias"))
    
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
         save_cm_mult_formats(outdir = output_dir,
                        file.name = glue::glue("contribution_matrix_all_data"))
      
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
        save_cm_mult_formats(outdir = output_dir,
                          file.name = glue::glue("approx_contribution_matrix_all_data"))
      
  }
  
  cat("  Creating cinema tables for", outcome, "...\n")
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
                          # contribution.matrix = loo_def,
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
      save_cinema_mult_formats(outdir = output_dir, file.name ="cinema_table_low_bias")

    cin_low %>% 
      filter(trt2 == "Placebo") %>% 
      save_cinema_mult_formats(outdir = output_dir, file.name ="cinema_table_low_bias_placebo")
    
    cin_some %>% cinema_filter() %>% 
      save_cinema_mult_formats(outdir = output_dir, file.name ="cinema_table_some_and_low_bias")
    
    cin_some %>% 
      filter(trt2 == "Placebo") %>% 
      save_cinema_mult_formats(outdir = output_dir, file.name ="cinema_table_some_and_low_bias_placebo")
    
    cin_high %>% cinema_filter() %>% 
      save_cinema_mult_formats(outdir = output_dir, file.name ="cinema_table_all_data")
    
    cin_high %>% 
      filter(trt2 == "Placebo") %>% 
      save_cinema_mult_formats(outdir = output_dir, file.name ="cinema_table_all_data")
    
    cin_low_def %>% 
      cinema_filter() %>% 
      save_cinema_mult_formats(outdir = output_dir, 
                               file.name ="cinema_table_deficiency_low_bias")
    
    cin_low_def %>% 
      filter(trt2 == "Placebo") %>% 
      save_cinema_mult_formats(outdir = output_dir, file.name ="cinema_table_deficiency_low_bias_placebo")
    
    cin_low_def_no_iran %>% 
      filter(trt2 == "Placebo") %>% 
      save_cinema_mult_formats(outdir = output_dir, file.name ="cinema_table_deficiency_low_bias_no_iran_placebo")
    
    cin_def %>% 
      cinema_filter() %>% 
      save_cinema_mult_formats(outdir = output_dir, file.name = "cinema_table_deficiency_all_data")
    
    cin_def %>% 
      filter(trt2 == "Placebo") %>% 
      save_cinema_mult_formats(outdir = output_dir, file.name ="cinema_table_deficiency_all_data_placebo")
    
    
  }
  
}
