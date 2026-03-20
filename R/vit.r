# R/vit.r (a box module)

box::use(
  magrittr[`%>%`],
  dplyr[filter, mutate, group_by, ungroup,summarize,select,all_of,arrange,distinct,rename,left_join, case_when],
  stringr[str_replace_all,str_extract,str_match,str_split],
  janitor[clean_names],
  rstantools[log_lik]
)
box::use(brms)
box::use(stats)
box::use(utils)
box::use(multinma)
box::use(tidyr)
box::use(rlang)
box::use(Matrix)
box::use(forcats)
box::use(ggplot2)
box::use(ggh4x)
box::use(mvtnorm)
box::use(loo)
box::use(snakecase)
box::use(foreach[`%dopar%`])
box::use(here)
box::use(doRNG[`%dorng%`])
box::use(doFuture)
box::use(future)
box::use(glue)
box::use(rstan)
box::use(purrr)
box::use(matrixStats)
box::use(cli)

mom_log_normal <- function(mu, sigma) {
  m = log(mu) - 0.5 * log(1 + (sigma^2 / mu^2))
  list(mean = m,
       sd = sqrt(log(1 + (sigma^2 / mu^2)))
       )
}

log_normal_prob_dis <- function(mu, sigma, scale, lower = FALSE) {
  cutoff <- cutoff_lookup[scale]
  moments <- mom_log_normal(mu, sigma)
  moments$mean[scale == "Log BAI"] <- mu[scale == "Log BAI"]
  moments$sd[scale == "Log BAI"] <- sigma[scale == "Log BAI"]
  cutoff[scale == "Log BAI"] <- exp(cutoff[scale == "Log BAI"])
  
  complete <- !is.na(cutoff) & !is.na(moments$mean) & !is.na(moments$sd)
  out <- ifelse(complete, stats::pnorm(q=log(cutoff), mean = moments$mean, sd = moments$sd, lower.tail = lower), NA)
  return(out)
}

normal_prob_dis <- function(mu, sigma, scale, lower = FALSE) {
  cutoff <- cutoff_lookup[scale]
  mu[scale == "Log BAI"] <- exp(mu[scale == "Log BAI"] + 0.5 * sigma[scale == "Log BAI"]^2)
  sigma[scale == "Log BAI"] <- sqrt((exp(sigma[scale == "Log BAI"]^2) - 1) * exp(2 * mu[scale == "Log BAI"] + sigma[scale == "Log BAI"]^2))
  cutoff[scale == "Log BAI"] <- exp(cutoff[scale == "Log BAI"])
  
  complete <- !is.na(cutoff) & !is.na(mu) & !is.na(sigma)
  out <- ifelse(complete, stats::pnorm(q=cutoff, mean = mu, sd = sigma, lower.tail = lower), NA)
  return(out)
}


log_normal_prob_from_mom <- function(mu, sigma, cutoff, lower = FALSE) {
  moments <- mom_log_normal(mu, sigma)
  
  complete <- !is.na(cutoff) & !is.na(moments$mean) & !is.na(moments$sd)
  out <- ifelse(complete, stats::pnorm(q=log(cutoff), mean = moments$mean, sd = moments$sd, lower.tail = lower), NA)
  return(out)
}

normal_prob_from_mom <- function(mu, sigma, cutoff, lower = FALSE) {
  complete <- !is.na(cutoff) & !is.na(mu) & !is.na(sigma)
  out <- ifelse(complete, stats::pnorm(q=cutoff, mean = mu, sd = sigma, lower.tail = lower), NA)
  return(out)
}

#' Gets J function correction for hedge's g
J_fun <- function(m) {
  exp(lgamma(m/2)- 0.5 * log(m/2)  - lgamma((m-1)/2))
}

psych_scale_cutoffs <- dplyr::tribble(
  ~scale,                                  ~cutoff,      ~dx,          ~transform,
  "Beck Depression Inventory",              14,           "depression", NA,
  "von Zersen mood scale (ZMS)",            27,           "depression", NA,
  "PHQ-9",                                  10,           "depression", NA,
  "MADRS",                                  7,            "depression", NA,
  "BDI-II",                                 14,           "depression", NA,
  "HAM-D",                                  8,            "depression", NA,
  "HADS-D",                                 8,            "depression", NA,
  "DASS-21, depression",                    5,            "depression", NA,
  "DASS-21, anxiety",                       4,            "anxiety",    NA,
  "CES-D",                                  16,           "depression", NA,
  "HAM-A",                                  8,            "anxiety",    NA,
  "SIGH-SAD",                               20,           "depression", NA,
  "BAI",                                    8,            "anxiety",    NA,
  "Log BAI",                                log(8),            "anxiety",    "log",
  "EPDS",                                   10,           "depression", NA,
  "HADS-A",                                 8,            "anxiety",    NA,
  "GDS",                                    10,           "depression", NA,
  "STAI-S",                                 40,           "anxiety",    NA,
  "depression per HADS-D",                  8,            "depression", NA,
  "SCAARED",                                25,           "anxiety",    NA,
  "MFQ",                                    27,           "depression", NA,
  "GAD-7",                                  10,           "anxiety",    NA,
  "AUC BDI-II",                             14,           "depression", NA,
  "AUC MADRS",                              7,            "depression", NA,
  "HRS",                                    18,           "anxiety",    NA,
  "Yasavage and Brink Depression Scale",    10,           "depression", NA,
  "BDI",                                    10,           "depression", NA,
  "GDS-15",                                 5,            "depression", NA,
  "MDI",                                    20,           "depression", NA,
  "STAI-T",                                 40,           "anxiety",    NA,
  "SF-36 MH",                              -42,           "depression", NA,
  "SF-12 MCS",                             -45.6,         "depression", NA,
  "DASS-42, depression",                    10,           "depression", NA,
  "DASS-42, anxiety",                       8,            "anxiety",    NA,
  "ΗDRS",                                   8,            "depression", NA,
  "HSCL-25",                                1.75,         "depression", NA,
  "HADS-D >= 8",                            8,            "depression", NA,
  "HADS-A >=8",                             8,            "anxiety",    NA
)

cutoff_lookup <- psych_scale_cutoffs %>%
  select(scale, cutoff) %>%
  tibble::deframe()

psych_dx_scale <- function(scale,value) {
  case_when(
    scale == "Beck Depression Inventory" & value >= 14 ~ "depression",
    scale == "von Zersen mood scale (ZMS)" & value >= 27 ~ "depression",
    scale == "PHQ-9" & value >= 10 ~ "depression",
    scale == "MADRS" & value >= 7 ~ "depression",
    scale == "BDI-II" & value >= 14 ~ "depression",
    scale == "HAM-D" & value >= 8 ~ "depression",
    scale == "HADS-D" & value >= 8 ~ "depression",
    scale == "DASS-21, depression" & value >= 5 ~ "depression",
    scale == "DASS-21, anxiety" & value >= 4 ~ "anxiety",
    scale == "CES-D" & value >= 16 ~ "depression",
    scale == "HAM-A" & value >= 8 ~ "anxiety",
    scale == "SIGH-SAD" & value >= 20 ~ "depression",
    scale == "BAI" & value >= 8 ~ "anxiety",
    scale == "Log BAI" & value >= log(8) ~ "anxiety",
    scale == "EPDS" & value >= 10 ~ "depression",
    scale == "HADS-A" & value >= 8 ~ "anxiety",
    scale == "GDS" & value >= 10 ~ "depression",
    scale == "STAI-S" & value >= 40 ~ "anxiety",
    scale == "depression per HADS-D" & value >= 8 ~ "depression",
    scale == "SCAARED" & value >= 25 ~ "anxiety",
    scale == "MFQ" & value >= 27 ~ "depression",
    scale == "GAD-7" & value >= 10 ~ "anxiety",
    scale == "AUC BDI-II" & value >= 14 ~ "depression",
    scale == "AUC MADRS" & value >= 7 ~ "depression",
    scale == "HRS" & value >= 18 ~ "anxiety",   # assumes HRS = HAM-A
    scale == "Yasavage and Brink Depression Scale" & value >= 10 ~ "depression",
    scale == "BDI" & value >= 10 ~ "depression",
    scale == "GDS-15" & value >= 5 ~ "depression",
    scale == "MDI" & value >= 20 ~ "depression",
    scale == "STAI-T" & value >= 40 ~ "anxiety",
    scale == "SF-36 MH" & value >= -42 ~ "depression",
    scale == "SF-12 MCS" & value >= -45.6 ~ "depression",
    scale == "DASS-42, depression" & value >= 10 ~ "depression",
    scale == "DASS-42, anxiety" & value >= 8 ~ "anxiety",
    scale == "ΗDRS" & value >= 8 ~ "depression",
    scale == "HSCL-25" & value >= 1.75 ~ "depression",
    scale == "HADS-D >= 8" & value >= 8 ~ "depression",
    scale == "HADS-A >=8" & value >= 8 ~ "anxiety",
    is.na(value) ~ NA_character_,
    TRUE ~ "none"
  )
}

vitamin_cutoff_tbl <- data.frame(matrix(c("zinc", "70", "ug/dL",
                                       "vitamin B9", "3", "ng/mL", # per WHO https://iris.who.int/server/api/core/bitstreams/781a6c61-8263-438b-a982-f0c57cb879ef/content
                                       "homocysteine", "15", "umol/L",
                                       "vitamin B1", "70", "nmol/L", #https://pmc.ncbi.nlm.nih.gov/articles/PMC6392124/
                                       "vitamin B2", "26.5", "nmol/L", #https://www.cambridge.org/core/journals/proceedings-of-the-nutrition-society/article/plasma-riboflavin-concentration-as-novel-indicator-for-vitaminb2-status-assessment-suggested-cutoffs-and-its-association-with-vitaminb6-status-in-women/B67BBD2C54B48221676A23565136C574   # Plasma <26.5 nmol/L: Strongly correlates with an EGRAC of ≥1.25, which is the IOM's threshold for marginal deficiency.
                                       "vitamin B6", "20", "nmol/L", #https://www.ncbi.nlm.nih.gov/books/NBK470579/#:~:text=The%20plasma%20PLP%20concentration%20is,short%2Dterm%20vitamin%20B6%20status.
                                       "vitamin D", "20", "ng/mL", #https://ods.od.nih.gov/factsheets/VitaminD-HealthProfessional/#:~:text=Assessing%20vitamin%20D%20status,analyses%20%5B5%2C6%5D.
                                       "vitamin B12", "200", "pg/mL", # https://ods.od.nih.gov/factsheets/VitaminB12-HealthProfessional/#:~:text=Another%20marker%20is%20total%20plasma,nutrient%20intakes%20of%20healthy%20people.
                                       "magnesium", "2.0", "mg/dL", #1.7 per IOM or 2.0 per https://pmc.ncbi.nlm.nih.gov/articles/PMC9186275/#:~:text=Conclusions,;%201.7%20mEq/L).
                                       "ferritin", "50", "ng/mL", # per https://pmc.ncbi.nlm.nih.gov/articles/PMC11817370/   # who recommends 15; had previously used 30 in last updates
                                       "vitamin C", "23", "umol/L" , #per WHO 11.4, other study of 23 https://pmc.ncbi.nlm.nih.gov/articles/PMC10450057/#:~:text=Hypovitaminosis%20C%20definition.,asthenia%2C%20muscular%20fatigue%20and%20depression.
                                       "selenoprotein P", "2.6", "mg/L" , #closes citation I can find https://www.sciencedirect.com/science/article/abs/pii/S0946672X07001708
                                       "selenoproteinP", "2.6", "mg/L" ,
                                       "selenium", "70", "ug/L" , # https://www.ncbi.nlm.nih.gov/sites/books/NBK482260/#:~:text=Tests%20commonly%20used%20to%20assess,increase%20the%20risk%20of%20fractures.
                                       "cobalamin", "200", "pg/mL" , # see B12 above
                                       "methylmalonate", "0.271", "umol/L" , #https://ods.od.nih.gov/factsheets/VitaminB12-HealthProfessional/#:~:text=Assessing%20Vitamin%20B12%20Status,6,9,10%5D.
                                       "Red Blood Cell B9", "140", "ng/mL", # per WHO 100 but IOM 140 ng/mL
                                       "rbc.B9", "140", "ng/mL"
), 
ncol=3, byrow=TRUE))
colnames(vitamin_cutoff_tbl) <- c("vitamin", "cutoff_value", "unit")

#' convert back to wacky concentration units
convert_from_molar <- function(molar, mw, to_unit) {
factor <- dplyr::case_when(
  to_unit == "ug/dL"  ~ mw * 1e6 / 10,  # mol/L → µg/dL
  to_unit == "ug/L"   ~ mw * 1e6,       # mol/L → µg/L
  to_unit == "µg/L"   ~ mw * 1e6,       # mol/L → µg/L
  to_unit == "µg/dL"  ~ mw * 1e6 / 10,  # mol/L → µg/dL
  to_unit == "ng/mL"  ~ mw * 1e9 / 1e3,   # mol/L → ng/mL
  to_unit == "pg/mL"  ~ mw * 1e12/ 1e3,      # mol/L → pg/mL
  to_unit == "nmol/L" ~ 1e9,            # mol/L → nmol/L
  to_unit == "umol/L" ~ 1e6,            # mol/L → µmol/L
  to_unit == "mg/dL"  ~ mw * 1e3 / 10,  # mol/L → mg/dL
  to_unit == "mg/L"   ~ mw * 1e3,       # mol/L → mg/L
  TRUE ~ NA_real_
)
molar * factor
}

convert_IU_to_ug <- function(dose) {
  
  # out <- if (  treatment == "vitamin D") {
  #   if (grepl("g", dose)) {
  #     40 * dose # convert g to IU for vitamin D
  #   } else if (grepl("IU", dose) | grepl("iu", dose)) {
  #     dose # already in IU
  #   } else {
  #     NA_real_
  #   }
  # }
  
  dose/40.0
  
  # return (out)
  
         
}

#' converting based upon specific vitamin
convert_vitamin_from_molar <- function(vitamin, value) {
  
  if(is.na(value)) return(sprintf("%s: %s", vitamin, "NA"))
  mw <- molecular_mass[vitamin]
  
  new_val <- convert_from_molar(value, mw, 
                      to_unit = vitamin_cutoff_tbl$unit[vitamin_cutoff_tbl$vitamin == vitamin])
  
  return(sprintf("%s: %.2f %s", vitamin, new_val, vitamin_cutoff_tbl$unit[vitamin_cutoff_tbl$vitamin == vitamin]))
}

#' convert based on the vitamin name, and also convert to more human-friendly units if the value is large
numeric_convert_vitamin_from_molar <- function(vitamin, value) {
  search_string <- glue::glue("{vitamin}$")
  mw <- molecular_mass[grepl(search_string, names(molecular_mass))]
  if(vitamin == "rbc.B9") {
    search_string <- "Red Blood Cell B9$"
    vitamin <- "Red Blood Cell B9"
  }
  new_val <- convert_from_molar(value, mw, 
                                to_unit = vitamin_cutoff_tbl$unit[grepl(search_string, vitamin_cutoff_tbl$vitamin)])
  
  unit <- vitamin_cutoff_tbl$unit[vitamin_cutoff_tbl$vitamin == vitamin]
  
  if (any(as.numeric(new_val)/1000 >= 1)) {
    new_val <- new_val/1000
    unit <- dplyr::case_when(
      unit == "ug/dL" ~ "mg/dL",
      unit == "ug/L" ~ "mg/L",
      unit == "µg/L" ~ "mg/L",
      unit == "µg/dL" ~ "mg/dL",
      unit == "ng/mL" ~ "µg/mL",
      unit == "pg/mL" ~ "ng/mL",
      unit == "nmol/L" ~ "µmol/L",
      unit == "umol/L" ~ "mmol/L",
      unit == "mg/dL" ~ "g/dL",
      unit == "mg/L" ~ "g/L",
      TRUE ~ unit
    )
  }
  
  
  return(sprintf("%.2f %s", new_val, unit))
}


vector_convert_vitamin_from_molar <- function(vitamin, value) {
  search_string <- glue::glue("{vitamin}$")
  mw <- molecular_mass[grepl(search_string, names(molecular_mass))]
  if(vitamin == "rbc.B9") {
    search_string <- "Red Blood Cell B9$"
    vitamin <- "Red Blood Cell B9"
  }
  if( grepl("B9", vitamin) &  !grepl("rbc", vitamin, ignore.case = TRUE)) {
    search_string <- "vitamin B9$"
  }
  unit <- vitamin_cutoff_tbl$unit[grepl(search_string, vitamin_cutoff_tbl$vitamin)]
  
  new_val <- convert_from_molar(value, mw, 
                                to_unit = unit)
  
  if (any(as.numeric(new_val)/1000 >= 1)) {
    new_val <- new_val/1000
    unit <- dplyr::case_when(
      unit == "ug/dL" ~ "mg/dL",
      unit == "ug/L" ~ "mg/L",
      unit == "µg/L" ~ "mg/L",
      unit == "µg/dL" ~ "mg/dL",
      unit == "ng/mL" ~ "µg/mL",
      unit == "pg/mL" ~ "ng/mL",
      unit == "nmol/L" ~ "µmol/L",
      unit == "umol/L" ~ "mmol/L",
      unit == "mg/dL" ~ "g/dL",
      unit == "mg/L" ~ "g/L",
      TRUE ~ unit
    )
  }
  
  
  return(list(value = new_val, unit = unit))
}

#' convert from mass per volume to molar
convert_conc <- function(mass_per_volume, molar_mass, unit,
                         direction = c("toMolarity", "toMass")) {
  direction <- match.arg(direction)
  
  # unit conversion factors
  mass_factors <- c(g = 1, mg = 1e-3, mcg = 1e-6, ug = 1e-6, µg = 1e-6, ng = 1e-9, pg = 1e-12)
  volume_factors <- c(L = 1, l = 1, dL = 1e-1, dl = 1e-1, ml = 1e-3, mL = 1e-3, ul = 1e-6, uL = 1e-6)
  
  # split the unit into mass and volume parts
  unit_parts <- str_split(unit, "/")
  
  mass_unit <- sapply(unit_parts,`[`,1) 
  volume_unit <- sapply(unit_parts,`[`,2)
  
  mass_included <- (mass_unit %in% names(mass_factors)) | is.na(mass_unit)
  volume_included <- (volume_unit %in% names(volume_factors)) | is.na(volume_unit)
  
  if(!all(mass_included)) {
    warning("Some mass unit names not found. Converting to NA")
    mass_unit[!mass_included] <- NA
  }
  if(!all(volume_included)) {
    warning("Some volume unit names not found. Converting to NA")
    volume_unit[!volume_included] <- NA
  }
  
  # normalize to g/L
  conc_gL <- mass_per_volume * mass_factors[mass_unit] / volume_factors[volume_unit]
  
  if (direction == "toMolarity") {
    # molarity (mol/L)
    return(conc_gL / molar_mass)
  } else {
    # go from molarity to g/L, then back to chosen mass/vol units
    mass_per_L <- mass_per_volume * molar_mass  # mol/L * g/mol = g/L
    return(mass_per_L / mass_factors[mass_unit] * volume_factors[volume_unit])
  }
}

#' convert from various mass units to grams
convert_mass_to_g <- function(mass, unit) {
  mass_factors <- c(g = 1, mg = 1e-3, mcg = 1e-6, ug = 1e-6, µg = 1e-6, ng = 1e-9, pg = 1e-12,
                    IU = 1e-6/40, iu = 1e-6/40)
  
  mass*mass_factors[unit]
}

convert_g_to_dose <- function(vitamin, value) {
  vitamin <- stringr::str_replace_all(vitamin, "vitamin ", "")
  dose <- dplyr::case_when(
  vitamin == "A" ~ sprintf("%2.f %s", value/1e-6/0.3, "IU"),
  vitamin == "D" ~ sprintf("%2.f %s", value/1e-6/0.025, "IU"),
  vitamin == "B1" ~ sprintf("%2.f %s", value/1e-3, "mg"),
  vitamin == "B2" ~ sprintf("%2.f %s", value/1e-3, "mg"),
  vitamin == "B3" ~ sprintf("%2.f %s", value/1e-3, "mg"),
  vitamin == "B5" ~ sprintf("%2.f %s", value/1e-3, "mg"),
  vitamin == "B6" ~ sprintf("%2.f %s", value/1e-3, "mg"),
  vitamin == "B7" ~ sprintf("%2.f %s", value/1e-6, "ug"),
  vitamin == "B9" ~ sprintf("%2.f %s", value/1e-3, "mg"),
  vitamin == "B12" ~ sprintf("%2.f %s", value/1e-6, "ug"),
  vitamin == "C" ~ sprintf("%2.f %s", value/1e-3, "mg"),
  vitamin == "E" ~ sprintf("%2.f %s", value/1e-3/0.67, "IU"),
  vitamin == "iron" ~ sprintf("%2.f %s", value/1e-3, "mg"),
  vitamin == "magnesium" ~ sprintf("%2.f %s", value/1e-3, "mg"),
  vitamin == "selenium" ~ sprintf("%2.f %s", value/1e-6, "ug"),
  vitamin == "zinc" ~ sprintf("%2.f %s", value/1e-3, "mg"),
  TRUE ~ NA_character_
  )
  
  dose
}

#' convert from different molarity units to SI units
convert_to_molar <- function(conc, unit = "mol") {
  # factors relative to mol/L
  unit_factors <- c(
    mol  = 1,
    mmol = 1e-3,
    umol = 1e-6,
    µmol = 1e-6,
    nmol = 1e-9,
    pmol = 1e-12
  )
  
  volume_factors <- c(L = 1, l = 1, dL = 1e-1, mL = 1e-3, uL = 1e-6)
  
  
  units <- str_split(unit, "/")
  
  mol_units <- sapply(units, `[`, 1) %>% 
    stringr::str_remove_all("\\s")
  vol_units <- sapply(units, `[`, 2) %>% 
    stringr::str_remove_all("\\s")
  
  mol_included <- mol_units %in% c(names(unit_factors),NA)
  volume_included <- vol_units %in% c(names(volume_factors),NA)
  
  if(!all(mol_included)) {
    warning("Some molar unit names not found. Converting to NA")
    mol_units[!mol_included] <- NA
  }
  if(!all(volume_included)) {
    warning("Some volume unit names not found. Converting to NA")
    vol_units[!volume_included] <- NA
  }
  
  return(conc * unit_factors[mol_units] / volume_factors[vol_units])
}

#' clean vitamin data to SI units
vitamin_data_to_SI <- function(data, colname) {
  data.frame(
    id = data$study,
    intervention = data$intervention,
    dose = data$tdd,
    name = data[[colname]]) %>% 
  distinct(id, intervention, dose, .keep_all = TRUE) %>% 
    tidyr::separate_rows(name, sep = ",\\s*") %>%
    tidyr::separate_wider_delim(name, delim = ": ", 
                         names = c("name", "value"),
                         too_few = "align_start") %>% 
    mutate(name = ifelse(stringr::str_starts(name, "normal|insufficient|deficient|deficiency"),
                         "",name)) %>% 
    mutate(value = ifelse(stringr::str_starts(value, "deficiency"),
                          "",value)) %>% 
    mutate(name = forcats::fct_recode(as.character(name), "vitamin D" = "25(OH)D",
                             "ferritin" = "Ferritin",
                             "vitamin B9" ="folate",
                             "vitamin D" = "vitamin d",
                             "vitamin B12" = "vitamin b12",
                             "vitamin B12" = "vitaminb12",
                             "vitamin B2" = "Vitamin B2",
                             "selenoproteinP" = "selenoprotein P",
                             "rbc.B9" = "RBC vitamin B9"),
           name = forcats::fct_na_level_to_value(name, extra_levels = "")) %>% 
    mutate(value = ifelse(name %>% is.na(), NA, value)) %>% 
    # tidyr::separate_wider_regex(
    #   value,
    #   patterns = c(
    #     "\\s*",
    #     value = "\\S+",   # first non-space chunk
    #     "\\s+",           # one or more spaces (dropped; unnamed)
    #     units = ".*"      # the rest
    #   ),
    #   too_few = "align_start"
    # ) %>% 
    tidyr::separate_wider_regex(
      value,
      patterns = c(
        "\\s*",
        value = "[-+]?(?:\\d*\\.\\d+|\\d+\\.?\\d*)(?:[eE][-+]?\\d+)?",
        "\\s*",
        units = ".*"
      ),
      too_few = "align_start"
    ) %>% 
    mutate(value = as.numeric(value)) %>% 
    mutate(std_value = 
             ifelse(
               # units %in% c("ng/mL", "md/dL","ug/L","pg/mL","mg/dL", "ng/mL", "µg/L","µg/l"), 
               grepl("g",units),
               convert_conc(value, molecular_mass[as.character(name)], units),
               convert_to_molar(value, units))) %>% 
    distinct() %>% 
    mutate(name = stringr::str_remove(as.character(name), "^vitamin\\s+") %>% as.factor())
}

#' get hedge's g from cohen's d
gfromd <- function(d, n1, n2, p = 1) {
  df <- n1 + n2 - p
  
  J_fun(df) * d
}

# function for OLS to hedge's g
ols_to_g <- function(b, seb, sdy, n1, n2, p) { 
  p <- as.numeric(p)
  n1 <- as.numeric(n1)
  n2 <- as.numeric(n2)
  seb <- as.numeric(seb)
  sdy <- as.numeric(sdy)
  b <- as.numeric(b)
  
  totaln = n1 + n2 ;
  df = totaln  - p ;
  sdpooled = sqrt((sdy*sdy*(totaln - 1)-(b*b*(n1*n2)/(n1+n2)))/(totaln - 2)) ;
  d  = b/sdpooled ;
  
  J_val <- J_fun(df)
  g  = d * J_val
  
  vd = seb^2/sdpooled^2  + d^2/(2*df);
  vg = vd * J_val^2
  se = sqrt(vg) 
  
  list(g = g, sd = se)
}

#' LMM to hedge's g
lmer_to_g <- function(b, seb, sdy, n1, n2, extra) { 
  extra <- strsplit(extra,",")
  
  p     <- sapply(extra,  `[`, 1) %>% as.numeric # correlation
  if (any(is.na(p))) p[is.na(p)] <- 0.2
  
  nc     <- sapply(extra,  `[`, 2) %>% as.numeric  # number of clusters
  
  n1 <- as.numeric(n1) 
  n2 <- as.numeric(n2)
  n1_ <- n1 * nc
  n2_ <- n2 * nc
  seb <- as.numeric(seb)
  sdy <- as.numeric(sdy)
  b <- as.numeric(b)
  
  n = n1 + n2 # num obs, also should be equal to totaln/nc
  totaln <- n1_ + n2_ #rows of data
  
  lambda = 1 - (2*(n-1)*p)/(totaln - 1)
  
  df <- (((totaln-2)-2*(n-1)*p)^2)/((totaln-2)*((1-p)^2)+ n*(totaln-2*n)*((p)^2)+2*(totaln-2*n)*p*(1-p))
  
  sdpooled = sqrt((sdy*sdy*(totaln - 1)-(b*b*(n1_*n2_)/(totaln)))/(totaln - 2))
  
  J_val <- J_fun(df)
  
  d  = b/sdpooled * sqrt(lambda);
  g  = d *  J_val;
  
  vd = lambda * (seb/sdpooled)^2  + d^2/(2*df);
  vg = vd * J_val^2 ;
  se = sqrt(vg);
  
  list(g = g, sd = se)
  
  # var b	  = parseFloat(frm.b.value) ;
  # var seb	  = parseFloat(frm.seb.value) ;
  # var sdy	  = parseFloat(frm.sdy.value) ;
  # var grp1n	  = parseFloat(frm.grp1n.value) ;
  # var grp2n	  = parseFloat(frm.grp2n.value) ;
  # var p         = parseFloat(frm.p.value) ;
  # var nc        = parseFloat(frm.nc.value) ;
  # totaln = grp1n + grp2n ;
  # n = totaln/nc ;
  # h = Math.pow((totaln-2)-2*(n-1)*p,2)/((totaln-2)*Math.pow(1-p,2)+ n*(totaln-2*n)*Math.pow(p,2)+2*(totaln-2*n)*p*(1-p)) ;
  # c = 1 - 3/(4*h-1) ;
  # lambda = 1 - (2*(n-1)*p)/(totaln - 1) ;
  # sdpooled = Math.sqrt((sdy*sdy*(totaln - 1)-(b*b*(grp1n*grp2n)/(grp1n+grp2n)))/(totaln - 2)) ;
  # d = (b/sdpooled)*Math.sqrt(lambda) ;
  # g = (c*b/sdpooled)*Math.sqrt(lambda) ;
  # v   = Math.pow(seb/sdpooled,2)*lambda + Math.pow(d,2)/(2*h) ;
  # vg  = Math.pow(c,2)*Math.pow(seb/sdpooled,2)*lambda + Math.pow(g,2)/(2*h) ;
  # se  = Math.sqrt(v) ;
  # seg = Math.sqrt(vg) ;
  # doutput(d,v,se,g,vg,seg) ; 
}

#' logistic estimates to hedge's g
logistic_to_g <- function(b, seb, n1, n2, p) {
  
  p <- as.numeric(p)
  n1 <- as.numeric(n1)
  n2 <- as.numeric(n2)
  seb <- as.numeric(seb)
  b <- as.numeric(b)
  
  loggedor  = b;
  vloggedor = seb*seb ;
  
  dl  =  sqrt(3)/pi * loggedor ;
  gl  = gfromd(dl,n1,n2,p);
  vl  =  3/(pi^2) * vloggedor ;
  vgl = vl * J_fun(n1 + n2 - p) ;
  segl = sqrt(vgl) ;
  
  list(g = gl, sd = segl)
}

#' molecular mass used in mass conversions
molecular_mass <- c("zinc" = 65.38,
                    "vitamin B9" = 441.4,
                    "homocysteine" = 135.18,
                    "vitamin B1" = 265.355,
                    "vitamin B2" = 376.36,
                    "vitamin B6" = 169.18,
                    "vitamin D" = 400.64,
                    "vitamin B12" = 1355.38,
                    "magnesium" = 24.305,
                    "ferritin" = 474000,
                    "vitamin C" = 176.124,
                    "selenoproteinP" = 98000,
                    selenium = 78.96,
                    cobalamin = 1355.38,
                    methylmalonate = 118.091,
                    rbc.B9 = 441.4
)

quick_convert_vit <- function(x, conc = "conc") {
  if (conc == "conc") {
    convert_conc(vitamin_cutoff_tbl %>% filter(vitamin == x) %>% dplyr::pull(cutoff_value) %>% as.numeric(), molecular_mass[x], vitamin_cutoff_tbl %>% filter(vitamin == x) %>% dplyr::pull(unit)) %>% as.numeric()
  } else {
    convert_to_molar(vitamin_cutoff_tbl %>% filter(vitamin == x) %>% dplyr::pull(cutoff_value) %>% as.numeric(), vitamin_cutoff_tbl %>% filter(vitamin == x) %>% dplyr::pull(unit)) %>% as.numeric()
  }
}

#' vitamin cutoffs for deficient
cutoff_vitamins <- c("zinc" = quick_convert_vit("zinc"),
                     "vitamin B9" = quick_convert_vit("vitamin B9"),
                     "homocysteine" = quick_convert_vit("homocysteine","molar"),
                     "vitamin B1" = quick_convert_vit("vitamin B1","molar"),
                     "vitamin B2" = quick_convert_vit("vitamin B2","molar"),
                     "vitamin B6" = quick_convert_vit("vitamin B6","molar"),
                     "vitamin D" = quick_convert_vit("vitamin D"),
                     "vitamin B12" = quick_convert_vit("vitamin B12"),
                     "magnesium" = quick_convert_vit("magnesium"),
                     "ferritin" = quick_convert_vit("ferritin"),
                     "vitamin C" = quick_convert_vit("vitamin C","molar"),
                     "selenoproteinP" = quick_convert_vit("selenoproteinP"),
                     selenium = quick_convert_vit("selenium"),
                     cobalamin = quick_convert_vit("cobalamin"),
                     methylmalonate = quick_convert_vit("methylmalonate","molar"),
                     rbc.B9 = quick_convert_vit("rbc.B9")
)

here::here("data","vitamins.rds") -> vitamin_data_path
here::here("data","quality.rds") -> bias_path

#' Load and clean vitamin data
get_vitamin_data <- function(outcome = c(NA_character_, "depression","anxiety"),
                             simple_analysis = TRUE,
                             include_full_bias = FALSE,
                             no_gain_scores = FALSE,
                             additive_tx = FALSE) {
  outcome_to <- match.arg(outcome)
  
  
  data <- readRDS(vitamin_data_path) %>% # load cleaned data
    {if(simple_analysis)
      filter(., !(study %in% c(
      "Gugger/2019/Switzerland",
      "Narula/2017/Canada",
      "Venkatasubramanian/2013/India",
      "Mozaffari-Khosravi/2013/Iran",
      "Penckofer/2022/USA"
    ))) else . } %>%
    {if(isTRUE(additive_tx)) {
      .
    } else {
      filter(., !(study %in% "Nguyen/2009/Guatemala"))
    }
    } %>%
    # no NA in data but just in case
    filter(!is.na(bias)) %>% 
    # clean name so easier to read
    mutate(intervention = as.character(intervention) %>% 
             gsub("z\\.","",.) %>% 
             gsub(",",", ",.) %>% 
             as.factor()) %>% 
    # remove underscore to make easier to read
    mutate(intervention = as.character(intervention) %>% 
             gsub("aerobic_exercise","aerobic exercise",.) ,
           intervention = intervention %>% 
             gsub("SSRI_TCA","SSRI/TCA",.) %>% 
             as.factor()) %>%
    mutate(treatment = as.character(treatment) %>% 
             gsub("aerobic_exercise","aerobic exercise",.) ,
           treatment = treatment %>% 
             gsub("SSRI_TCA","SSRI/TCA",.) %>% 
             as.factor()) %>%
    mutate(control = as.character(control) %>% 
             gsub("aerobic_exercise","aerobic exercise",.) ,
           control = control %>% 
             gsub("SSRI_TCA","SSRI/TCA",.) %>% 
             as.factor()) %>%
    # prefer gain scores and model based estimates for increased precision
    {if(!no_gain_scores) {
      mutate(., y = ifelse(!is.na(y.gain) & !is.na(sd.y.gain), 
                        y.gain, 
                        y),
             # only for not missing in SD since baseline category has no contrast (is NA)
             sd.y = ifelse(!is.na(sd.y.gain), 
                           sd.y.gain, 
                           sd.y))
    } else {
      .
    } } %>% 
    # remove rows without variance estimates
    filter(!is.na(sd.y)) %>%
    # prefer scales with higher preference (more common in data)
    group_by(study, intervention, target) %>%
    filter(scale.preference == min(scale.preference)) %>% 
    ungroup() %>%
    # drop unused levels for easier reading
    droplevels() %>% 
    mutate(depressed = ifelse(grepl("depression", psych.dx), 1, 0)) %>%
    mutate(anxiety = ifelse(grepl("anxiety", psych.dx), 1, 0)) %>%
    mutate(def = forcats::fct_na_value_to_level(ifelse(!is.na(vit.def) & !grepl("none",vit.def), 1, 0 * as.numeric(vit.def)) %>% factor(levels = c(0,1),labels = c("no","yes")) ,"NA")) %>% 
    mutate(bias = forcats::fct_relevel(bias, "low", "some concerns", "high")) %>%
    group_by(study,target) %>%
    filter(dplyr::n() > 1) %>%
    ungroup()
  
  if( isTRUE(include_full_bias)) {
    bias_ratings <- bias_path %>% 
      readRDS() %>%
      filter(study %in% data$study)
    data %>% dplyr::left_join(bias_ratings, by = "study") -> data
    
  }
  
  if (is.na(outcome_to[1])) {
    return(data)
  } else if (outcome_to == "depression") {
    return(data %>% filter(target == "depression"))
  } else if (outcome_to == "anxiety") {
    return(data %>% filter(target == "anxiety"))
  } else {
    stop("target must be one of NA, 'depression', or 'anxiety'")
  }
  
}

construct_nma_network <- function(data, trt_ref = NULL) {
  network <- multinma::set_agd_contrast(data
                                 , study = "study",
                                 trt = "intervention",
                                 y = "y",
                                 se = "sd.y",
                                 sample_size = "final.N",
                                 trt_ref = trt_ref
  )
  if(multinma::is_network_connected(network) == FALSE) {
    warning("The treatment network is not connected. Please check your data.")
  }
  return(network)
}

network_disconnect_check <- function(network, outcome = NULL, error.on.true = TRUE) {
  
  connected <- multinma::is_network_connected(network)
  
  if(connected) {
    msg  <- glue::glue("full ",outcome, " target network is connected")
  } else if (!connected) {
    msg  <- glue::glue("full ",outcome, " target network is *not* connected")
    if(isTRUE(error.on.true)) rlang::abort(msg)
  } else {
    stop("connected status unknown")
  }
  names(connected) <- msg
  
  
  return(invisible(connected))
}

mvn_stan_fun <- "
real mvn_nma_lpdf(vector y, vector mu,
                      real tau,
                      vector z_u, // standard normal random effects
                      array[] int study,
                      array[] real se,
                      array[] real cov,
                      matrix S_L, // Cholesky factor of covariance of data
                      matrix R_L // cholesky factor of correlation of random effects
                      ) {
    vector[rows(y)] u = tau * (R_L * z_u); // study level random effects
    real lp = multi_normal_cholesky_lpdf(y | mu + u, S_L); // Sigma = S_L*S_L'
    return lp;
}
  
real mvn_nma_2_lpdf(vector y, vector mu,
                      real tau,
                      array[] int study,
                      array[] real se,
                      array[] real cov,
                      matrix S_L, // Cholesky factor of covariance of data
                      matrix R_L // cholesky factor of correlation of random effects
                      ) {
    matrix[rows(y),rows(y)] Sigma = S_L*S_L' + tau*tau * R_L*R_L'; // study level random effects
    real lp = multi_normal_lpdf(y | mu, Sigma); // Sigma = S_L*S_L'
    return lp;
  }
  
  vector mvn_nma_rng(vector mu, 
                      real tau,
                      vector z_u, // standard normal random effects
                      array[] int study,
                      array[] real se,
                      array[] real cov,
                      matrix S_L, // Cholesky factor of covariance of data
                      matrix R_L // cholesky factor of correlation of random effects
                      ) {
    vector[rows(mu)] u = tau * (R_L * z_u); // study level random effects
    vector[rows(mu)] y_sim;
    y_sim = multi_normal_cholesky_rng(mu + u, S_L); // Sigma = S_L*S_L'
    return y_sim;
  }
"

mvn_nma <- brms::custom_family(
  "mvn_nma",
  dpars = c("mu", "tau", "zU"),  
  lb = c(NA,  0, NA),
  ub = c(NA, NA, NA),
  vars = c(
    "vint1","vreal1","vreal2",
    "Sigma_L","R_L"),
  links = c("identity","identity","identity"),
  type  = "real",
  loop = FALSE
)

log_lik_mvn_nma <- function(i, prep) {
  study <- prep$data$vint1[i]
  study_idx <- which(prep$data$vint1 == study)
  n_j   <- length(study_idx)
  
  # likelihood
  out <- if (i == min(study_idx)) {
    se    <- prep$data$vreal1[study_idx]
    cov   <- prep$data$vreal2[study_idx]
    
    mu  <- brms::get_dpar(prep, "mu",  i = study_idx)
    tau <- brms::get_dpar(prep, "tau", i = study_idx)
    z_u <- brms::get_dpar(prep, "zU",  i = study_idx)
    
    S <- matrix(cov, n_j, n_j)
    diag(S) <- se^2
    R <- matrix(0.5, n_j, n_j)
    diag(R) <- 1
    
    # response value for obs i
    y_i <- c(prep$data$Y[study_idx])
    if(n_j > 1) {
      R_L <- eigen_sqrt(R)
      sapply(1:prep$ndraws, function(j) {
        u   <- (R_L %*% z_u[j,,drop = TRUE]) * tau[j]
        mvtnorm::dmvnorm(y_i, mean = mu[j,] + u, sigma = S,
                      log = TRUE,
                      checkSymmetry = FALSE)
      })
      
    } else {
      u <- sqrt(c(R)) * tau * z_u
      sigma <- sqrt(c(S))
      stats::dnorm(y_i, mean = (mu + u), sd = sigma, log = TRUE)
    }
  } else {
    NA_real_
  }
  
  return(out)
}

posterior_epred_mvn_nma <- function(prep, ...) {
  
  mu <- brms::get_dpar(prep, "mu")
  
  return(mu)
}

posterior_predict_mvt_nma <- function(i, prep, ...) {
  mu    <- brms::get_dpar(prep, "mu", i = i)
  tau   <- brms::get_dpar(prep, "tau", i = i)
  z_u   <- brms::get_dpar(prep, "zU", i = i)
  
  # study <- prep$data$vint1
  # study_idx_list <- split(1:nrow(prep$data), study)
  
  y_sim <- matrix(NA, nrow = prep$ndraws, ncol = 1)
  
  # for(study_i in names(study_idx_list)) {
  study_idx <- i # study_idx <- study_idx_list[[study_i]]
  n_j   <- 1 #length(study_idx)
  
  se    <- prep$data$vreal1[study_idx]
  cov   <- prep$data$vreal2[study_idx]
  
  S <- matrix(cov, n_j, n_j)
  diag(S) <- se^2
  R <- matrix(0.5, n_j, n_j)
  diag(R) <- 1
  
  
  
  for (j in 1:prep$ndraws) {
    y_sim[j,] <- if (n_j > 1) {
      R_L <- eigen_sqrt(R)
      u   <- (R_L %*% z_u[j,drop = TRUE]) * tau[j]
      mvtnorm::rmvnorm(1, mean = mu[j] + u,
                       sigma = S,
                       checkSymmetry = FALSE)
    } else {
      u   <- z_u[j,drop = TRUE] * tau[j]
      stats::rnorm(1, mean = mu[j] + c(u),
                   sd = sqrt(c(S)))
    }
    
  }
  
  return(y_sim)
}


{covariance_functions_stan <- "
  // Build an N x N covariance matrix where:
  //  - diag elements are individual variances var[i]
  //  - off-diagonals for observations in the same group g are cov_group[g]
  //  - off-diagonals for different groups are 0
  //
  // Arguments:
  //   group: integer group indicator for each obs (1..G)
  //   var:   length-N vector of individual variances
  //   cov:   length-N vector of within-group covariances
  matrix make_group_cov(array[] int group,
                        array[] real se,
                        array[] real cov) {
    int N = size(group);
    matrix[N, N] Sigma;
    vector[N] V;

    // start with diagonal matrix of variances
    for(i in 1:N) {
      V[i] = se[i] * se[i];
    }
    Sigma = diag_matrix(V);

    // fill in off-diagonals
    if (N > 1) {
      for (i in 1:(N - 1)) {
        for (j in (i + 1):N) {
          if (group[i] == group[j]) {
            // same group: use that group's covariance
            real c = cov[i];
            Sigma[i, j] = c;
            Sigma[j, i] = c;
          } else {
            // different groups: covariance is zero
            Sigma[i, j] = 0;
            Sigma[j, i] = 0;
          }
        }
      }
    }

    return Sigma;
  }

matrix make_group_cor(array[] int group) {
    int N = size(group);
    matrix[N, N] R;
    real c = 0.5;

    // start with diagonal matrix of variances
    R = diag_matrix(ones_vector(N));

    // fill in off-diagonals
    if (N > 1) {
      for (i in 1:(N - 1)) {
        for (j in (i + 1):N) {
          if (group[i] == group[j]) {
            // same group: use that group's covariance
            R[i, j] = c;
            R[j, i] = c;
          } else {
            // different groups: covariance is zero
            R[i, j] = 0;
            R[j, i] = 0;
          }
        }
      }
    }

    return R;
  }

"}

{mvt_stan_fun <- "
real mvt_nma_lpdf(vector y, vector mu,
                      real tau,
                      vector z_u, // standard normal random effects
                      array[] int study,
                      array[] real se,
                      array[] real cov,
                      matrix S_L, // Cholesky factor of covariance of data
                      matrix R_L // cholesky factor of correlation of random effects
                      ) {
    vector[rows(y)] u = tau * (R_L * z_u); // study level random effects
    real lp = multi_student_t_cholesky_lpdf(y | 30.0,  mu + u, S_L); // Sigma = S_L*S_L'
    return lp;
  }

  vector mvt_nma_rng(vector mu, 
                      real tau,
                      vector z_u, // standard normal random effects
                      array[] int study,
                      array[] real se,
                      array[] real cov,
                      matrix S_L, // Cholesky factor of covariance of data
                      matrix R_L // cholesky factor of correlation of random effects
                      ) {
    vector[rows(mu)] u = tau * (R_L * z_u); // study level random effects
    vector[rows(mu)] y_sim;
    y_sim = multi_student_t_cholesky_rng(30.0, mu + u, S_L); // Sigma = S_L*S_L'
    return y_sim;
  }
"}

posterior_epred_mvt_nma <- function(prep) {
  
  mu <- brms::get_dpar(prep, "mu")
  
  return(mu)
}

posterior_predict_mvt_nma <- function(i, prep, ...) {
  mu    <- brms::get_dpar(prep, "mu", i = i)
  tau   <- brms::get_dpar(prep, "tau", i = i)
  z_u   <- brms::get_dpar(prep, "zU", i = i)
  
  # study <- prep$data$vint1
  # study_idx_list <- split(1:nrow(prep$data), study)
  
  y_sim <- matrix(NA, nrow = prep$ndraws, ncol = 1)
  
  # for(study_i in names(study_idx_list)) {
    study_idx <- i # study_idx <- study_idx_list[[study_i]]
    n_j   <- 1 #length(study_idx)
    
    se    <- prep$data$vreal1[study_idx]
    cov   <- prep$data$vreal2[study_idx]
    
    S <- matrix(cov, n_j, n_j)
    diag(S) <- se^2
    R <- matrix(0.5, n_j, n_j)
    diag(R) <- 1
    
    R_L <- eigen_sqrt(R)
    
    for (j in 1:prep$ndraws) {
      u   <- (R_L %*% z_u[j,drop = TRUE]) * tau[j]
      y_sim[j,] <- mvtnorm::rmvt(1, delta = mu[j] + u,
                                         sigma = S,
                                         df = 30.0,
                                         checkSymmetry = FALSE)
    }
  # }
  
  return(y_sim)
}

eigen_sqrt <- function(mat) {
  e_list <- eigen(mat)
  root_v <- e_list$values %>% sqrt()
  root_mat <- e_list$vectors %*% diag(root_v) %*% t(e_list$vectors)
  return(root_mat)
}

log_lik_mvt_nma <- function(i, prep) {
  study <- prep$data$vint1[i]
  study_idx <- which(prep$data$vint1 == study)
  n_j   <- length(study_idx)
  
  # likelihood
  out <- if (i == min(study_idx)) {
    se    <- prep$data$vreal1[study_idx]
    cov   <- prep$data$vreal2[study_idx]
    
    mu  <- brms::get_dpar(prep, "mu",  i = study_idx)
    tau <- brms::get_dpar(prep, "tau", i = study_idx)
    z_u <- brms::get_dpar(prep, "zU",  i = study_idx)
    
    S <- matrix(cov, n_j, n_j)
    diag(S) <- se^2
    R <- matrix(0.5, n_j, n_j)
    diag(R) <- 1
    
    # response value for obs i
    y_i <- c(prep$data$Y[study_idx])
    if(n_j > 1) {
      R_L <- eigen_sqrt(R)
      sapply(1:prep$ndraws, function(j) {
        u   <- (R_L %*% z_u[j,,drop = TRUE]) * tau[j]
        mvtnorm::dmvt(y_i, delta = mu[j,] + u, sigma = S,
                      df = 30.0, 
                      log = TRUE,
                      checkSymmetry = FALSE)
      })
      
    } else {
      u <- sqrt(c(R)) * tau * z_u
      sigma <- sqrt(c(S))
      stats::dt((y_i - (mu + u))/sigma, df = 30.0, log = TRUE) - log(sigma)
    }
  } else {
    NA_real_
  }
  # u <- brms::as_draws_df(get("fit_brms_low", envir = .GlobalEnv))[paste0("u_draw[",study_idx,"]")] %>% as.matrix()
  # out <- if (i == min(study_idx)) {
  #   if(n_j > 1) {
  #     sapply(1:prep$ndraws, function(j) {
  #       mvtnorm::dmvnorm(y_i, mean = mu[j,] + u[j,], sigma = S, log = TRUE,
  #                        checkSymmetry = FALSE)
  #     })
  #     
  #   } else {
  #     stats::dnorm(y_i, mean = mu + c(u), sd = sqrt(S), log = TRUE)
  #   }
  # } else {
  #   NA_real_
  # }
  return(out)
}

mvt_nma <- brms::custom_family(
  "mvt_nma",
  dpars = c("mu", "tau","zU"),  
  lb = c(NA, 0, NA),
  ub = c(NA, NA, NA),
  vars = c(
    "vint1", "vreal1", "vreal2",
           "Sigma_L","R_L"),
  links = c("identity","identity","identity"),
  type  = "real",
  loop = FALSE,
  posterior_epred = posterior_epred_mvt_nma,
  posterior_predict = posterior_predict_mvt_nma,
  log_lik = log_lik_mvt_nma
)

{binomial_stan_fun <- "
int binomial_logit_nma_lpdf(int y, vector mu,
                      real tau,
                      vector z_u, // standard normal random effects
                      array[] int study,
                      array[] int trials,
                      matrix R_L // cholesky factor of correlation of random effects
                      ) {
    vector[rows(y)] u = tau * (R_L * z_u); // study level random effects
    real lp = binomial_logit_lpmf(y | trials,  mu + u); 
    return lp;
  }

  vector binomial_logit_nma_rng(vector mu,
                      real tau,
                      vector z_u, // standard normal random effects
                      array[] int study,
                      array[] int trials,
                      matrix R_L // cholesky factor of correlation of random effects
                      ) {
    vector[rows(mu)] u = tau * (R_L * z_u); // study level random effects
    vector[rows(mu)] prob = inv_logit(mu + u);
    vector[rows(mu)] y_sim;
    y_sim = binomial_rng(trials, prob); 
    return y_sim;
  }
"}

binomial_nma <- brms::custom_family(
  "binomial_logit_nma",
  dpars = c("mu", "tau","zU"),  
  lb = c(NA, 0, NA),
  ub = c(NA, NA, NA),
  vars = c(
    "vint1", "vint2",
    "R_L"),
  links = c("identity","identity","identity"),
  type  = "int",
  loop = FALSE,
  posterior_epred = posterior_epred_mvt_nma,
  posterior_predict = NULL,
  log_lik = NULL
)


log_lik.vitfit <- function(object, newdata = NULL, re_formula = NULL,
                            resp = NULL, ndraws = NULL, draw_ids = NULL,
                            pointwise = FALSE, combine = TRUE,
                            add_point_estimate = FALSE,
                            cores = NULL, ...) {
  pointwise <- brms:::as_one_logical(pointwise)
  combine <- brms:::as_one_logical(combine)
  add_point_estimate <- brms:::as_one_logical(add_point_estimate)
  brms:::contains_draws(object)
  object <- brms:::restructure(object)
  prep <- brms:::prepare_predictions(
    object, newdata = newdata, re_formula = re_formula, resp = resp,
    ndraws = ndraws, draw_ids = draw_ids, check_response = TRUE, ...
  )
  if (add_point_estimate) {
    # required for the loo_subsample method
    # Computing a point estimate based on the full prep object is too
    # difficult due to its highly nested structure. As an alternative, a second
    # prep object is created from the point estimates of the draws directly.
    attr(prep, "point_estimate") <- brms:::prepare_predictions(
      object, newdata = newdata, re_formula = re_formula, resp = resp,
      ndraws = ndraws, draw_ids = draw_ids, check_response = TRUE,
      point_estimate = "median", ...
    )
  }
  if (pointwise) {
    stopifnot(combine)
    log_lik <- log_lik_pointwise
    # names need to be 'data' and 'draws' as per ?loo::loo.function
    attr(log_lik, "data") <- data.frame(i = seq_len(brms:::choose_N(prep)))
    attr(log_lik, "draws") <- prep
  } else {
    log_lik <- log_lik(prep, combine = combine, cores = cores)
    study <- if ( !is.null(prep$resps$y$data$vint1) ) {
      prep$resps$y$data$vint1
    } else {
      prep$data$vint1
    }
    study_unique <- !(duplicated(study))
    log_lik <- log_lik[,study_unique, drop = FALSE]
    if (anyNA(log_lik)) {
      brms:::warning2(
        "NAs were found in the log-likelihood. Possibly this is because ",
        "some of your responses contain NAs. If you use 'mi' terms, try ",
        "setting 'resp' to those response variables without missing values. ",
        "Alternatively, use 'newdata' to predict only complete cases."
      )
    }
  }
  log_lik
}

reloo.vitfit <- function(x, loo = NULL, k_threshold = 0.7, newdata = NULL,
         resp = NULL, check = TRUE, recompile = NULL,
         future_args = list(), ...) {
  stopifnot(brms:::is.brmsfit(x), is.list(future_args))
  if (brms:::is.brmsfit_multiple(x)) {
    brms:::warn_brmsfit_multiple(x)
    class(x) <- c("vitfit", "brmsfit")
  }
  loo <- loo %||% x$criteria[["loo"]]
  if (is.null(loo)) {
    brms:::stop2("No 'loo' object was provided and none is stored within the model.")
  } else if (!loo::is.loo(loo)) {
    brms:::stop2("Inputs to the 'loo' argument must be of class 'loo'.")
  }
  if (is.null(newdata)) {
    mf <- stats::model.frame(x)
  } else {
    mf <- as.data.frame(newdata)
  }
  mf <- brms:::rm_attr(mf, c("terms", "brmsframe"))
  study_name <- mf[[".study"]]
  study <- study_name %>% as.integer()
  n_study   <- length(unique(study))
  unique_study <- unique(study)
  unique_sn    <- unique(study_name)
  
  if (n_study != NROW(loo$pointwise)) {
    brms:::stop2("Number of observations in 'loo' and 'x' do not match.")
  }
  check <- brms:::as_one_logical(check)
  if (check) {
    yhash_loo <- attr(loo, "yhash")
    yhash_fit <- brms:::hash_response(x, newdata = newdata)
    if (!brms:::is_equal(yhash_loo, yhash_fit)) {
      brms:::stop2(
        "Response values used in 'loo' and 'x' do not match. ",
        "If this is a false positive, please set 'check' to FALSE."
      )
    }
  }
  if (is.null(loo$diagnostics$pareto_k)) {
    brms:::stop2("No Pareto k estimates found in the 'loo' object.")
  }
  obs <- loo::pareto_k_ids(loo, k_threshold)
  J <- length(obs)
  if (J == 0L) {
    message(
      "No problematic observations found. ",
      "Returning the original 'loo' object."
    )
    return(loo)
  }
  
  # ensure that the model can be run in the current R session
  x <- brms:::recompile_model(x, recompile = recompile)
  
  # split dots for use in log_lik and update
  dots <- list(...)
  ll_arg_names <- brms:::arg_names("log_lik")
  ll_arg_names <- intersect(names(dots), ll_arg_names)
  ll_args <- dots[ll_arg_names]
  ll_args$allow_new_levels <- TRUE
  ll_args$sample_new_levels <-
    brms:::first_not_null(ll_args$sample_new_levels, "gaussian")
  ll_args$resp <- resp
  ll_args$combine <- TRUE
  # cores is used in both log_lik and update
  up_arg_names <- setdiff(names(dots), setdiff(ll_arg_names, "cores"))
  up_args <- dots[up_arg_names]
  up_args$object <- x
  up_args$refresh <- 0
  
  .reloo <- function(j) {
    message(
      "\nFitting model ", j, " out of ", J,
      " (leaving out study ", obs[j], ": ",unique_sn[j],")"
    )
    drop_study <- unique_study[obs[j]]
    omitted <- which(study == drop_study)
    mf_omitted <- mf[-omitted, , drop = FALSE]
    up_args$newdata <- mf_omitted
    up_args$data2 <- brms:::subset_data2(x$data2, -omitted)
    
    # cov_mats <- list(R = x$prep$R[-omitted, -omitted],
    #                  S = x$prep$Sigma[-omitted, -omitted],
    #                  R_l = chol(x$prep$R[-omitted, -omitted]),
    #                  Sigma_l = chol(x$prep$Sigma[-omitted, -omitted]))
    # up_args$stanvars <- make_stanvar(cov_mats, x$prep$stan_fun)
    fit_j <- brms:::SW(brms:::do_call(stats::update, up_args))
    class(fit_j) <- c("vitfit", "brmsfit")
    ll_args$object <- fit_j
    ll_args$newdata <- mf[omitted, , drop = FALSE]
    ll_args$newdata2 <- brms:::subset_data2(x$data2, omitted)
    ll_args$allow_new_levels <- TRUE
    ll_args$skip_validation <- TRUE
    return(brms::do_call(log_lik, ll_args))
  }
  
  message(
    J, " problematic observation(s) found.",
    "\nThe model will be refit ", J, " times."
  )
  # TODO: separate parallel and non-parallel code to enable better printing?
  # future_args$X <- seq_len(J)
  # future_args$FUN <- .reloo
  # future_args$future.seed <- TRUE
  # browser(); #.reloo(2)
  # lls <- brms:::do_call("future_lapply", future_args, pkg = "future.apply")
  # lls <- do.call("lapply",future_args)
  lls <- lapply(seq_len(J), .reloo)
  
  # most of the following code is taken from rstanarm:::reloo
  # compute elpd_{loo,j} for each of the held out observations
  elpd_loo <- brms:::ulapply(lls, brms:::log_mean_exp)
  # compute \hat{lpd}_j for each of the held out observations (using log-lik
  # matrix from full posterior, not the leave-one-out posteriors)
  drop_studies <- unique_study[obs]
  omitted <- which(study %in% drop_studies)
  mf_obs <- mf[omitted, , drop = FALSE]
  data2_obs <- brms:::subset_data2(x$data2, omitted)
  ll_x <- log_lik(x, newdata = mf_obs, newdata2 = data2_obs)
  hat_lpd <- apply(ll_x, 2, brms:::log_mean_exp)
  # compute effective number of parameters
  p_loo <- hat_lpd - elpd_loo
  # replace parts of the loo object with these computed quantities
  sel <- c("elpd_loo", "p_loo", "looic")
  loo$pointwise[obs, sel] <- cbind(elpd_loo, p_loo, -2 * elpd_loo)
  new_pw <- loo$pointwise[, sel, drop = FALSE]
  loo$estimates[, 1] <- colSums(new_pw)
  loo$estimates[, 2] <- sqrt(nrow(loo$pointwise) * apply(new_pw, 2, stats::var))
  # what should we do about pareto-k? for now setting them to 0
  loo$diagnostics$pareto_k[obs] <- 0
  loo
}

special_clean_names <- function(data) {
  
  clean_chain <- function(x) {
    x %>% 
      str_replace_all(",","") %>% 
      str_replace_all("-","_") %>% 
      str_replace_all(" ","_") %>% 
      str_replace_all("/","_") %>%
      str_replace_all(":","..") %>%
      str_replace_all("\\(",".") %>%
      str_replace_all("\\)",".") %>%
      str_replace_all(">",".GT.") %>%
      str_replace_all("<",".LT.") %>%
      str_replace_all("__",".") %>%
      tolower()
  }
  
  if(is.data.frame(data) || is.matrix(data) ) {
    colnames(data) <- colnames(data) %>% clean_chain()
  } else if (!is.null(names(data))) {
    names(data) <- names(data) %>% clean_chain()
  } else if (is.character(data)) {
    data <- data %>% clean_chain()
  }
  
  data
}

restore_dose_names <- function(dosecol) {
  varname <- paste0(".dose", dosecol %>%
    special_clean_names())
  cbind(label = dosecol, varname = varname) %>% 
    as.data.frame() %>% 
    distinct() %>% 
    arrange(label)
}

construct_nma_formula <- function(main, interactions, confounders) {
  # this function adds any .contrast() terms to the nma formula
  # and will set up the contrast model matrix if ".trt" is used
  # in either a or b forms
  find_all_contrast_args <- function(expr) {
    out <- list()
    if (is.call(expr)) {
      if (identical(expr[[1]], as.name(".contrast"))) {
        out <- append(out, list(expr[[2]]))
      }
      for (e in as.list(expr)[-1]) {
        out <- append(out, Recall(e))
      }
    }
    out
  }
  
  
  # start with just study intercepts
  nma_form <- stats::as.formula(~-1 + .study)
  
  if(missing(main)) main <- main ~ 0 + .trt
  if(missing(interactions)) interactions <- interactions ~ 0
  if(missing(confounders)) confounders <- confounders ~ 0
  main <- stats::as.formula(main)
  interactions <- stats::as.formula(interactions)
  confounders <- stats::as.formula(confounders)
  
  # check if .contrast() is in either form
  # if so, extract any terms inside the .contrast and add to nma_form
  # eg .contrast(def + sd.y) would extract def + sd.y into nma_form
  if(any(grepl(".contrast",deparse(main)))) {
    # Pull RHS
    arhs <- rlang::f_rhs(main)
    
    # Find .contrast() call and extract argument
    find_all_contrast_args(arhs) -> a_contrast_terms
    a_contrast_terms <- gsub("mi\\(([^)]*)\\)", "\\1", a_contrast_terms)
    
    nma_form <- stats::update(nma_form, paste0("~. +", a_contrast_terms))
  }
  if(any(grepl(".contrast",deparse(interactions)))) {
    # Pull RHS
    brhs <- rlang::f_rhs(interactions)
    
    # Find .contrast() call and extract argument
    find_all_contrast_args(brhs) -> b_contrast_terms
    b_contrast_terms <- gsub("mi\\(([^)]*)\\)", "\\1", b_contrast_terms)
    
    
    nma_form <- stats::update(nma_form, paste0("~. +", b_contrast_terms))
  }
  if(any(grepl(".contrast",deparse(confounders)))) {
    # Pull RHS
    crhs <- rlang::f_rhs(confounders)
    
    # Find .contrast() call and extract argument
    find_all_contrast_args(crhs) -> c_contrast_terms
    c_contrast_terms <- gsub("mi\\(([^)]*)\\)", "\\1", c_contrast_terms)
    
    
    nma_form <- stats::update(nma_form, paste0("~. +", c_contrast_terms))
  }
  
  
  # check if .trt is in either form
  # if so, extract any terms with .trt and add to nma_form
  # this allows for interaction terms with .trt to be included
  # e.g., .trt:depressed or .trt:def:depressed
  # but also allows for other terms to be included in a or b
  # that do not involve .trt, e.g., age, baseline severity, etc
  if(any(grepl(".trt",deparse(main)))) {
    # stringr::str_split(deparse(main), "\\s\\+\\s")[[1]] %>%
    #   stringr::str_trim() %>%
    #   grep(".trt", ., value = TRUE)  %>% 
    #   paste0(., collapse = " + ")  -> a_trt_terms
    a_trt_terms <- ".trt"
    
    nma_form <- stats::update(nma_form, paste0("~. +", a_trt_terms))
  }
  
  if(any(grepl(".trt",deparse(interactions)))) {
    # stringr::str_split(deparse(interactions), "\\s\\+\\s")[[1]] %>%
    #   stringr::str_trim() %>%
    #   grep(".trt", ., value = TRUE)  %>% 
    #   paste0(., collapse = " + ") -> b_trt_terms
    b_trt_terms <- ".trt"
    
    nma_form <- stats::update(nma_form, paste0("~. +", b_trt_terms))
  }
  
  if(any(grepl(".trt",deparse(confounders)))) {
    # stringr::str_split(deparse(interactions), "\\s\\+\\s")[[1]] %>%
    #   stringr::str_trim() %>%
    #   grep(".trt", ., value = TRUE)  %>% 
    #   paste0(., collapse = " + ") -> b_trt_terms
    c_trt_terms <- ".trt"
    
    nma_form <- stats::update(nma_form, paste0("~. +", c_trt_terms))
  }
  
  attr(nma_form,"dose") <- FALSE
  
  if(any(grepl(".dose",deparse(main)))) {
    # stringr::str_split(deparse(main), "\\s\\+\\s")[[1]] %>%
    #   stringr::str_trim() %>%
    #   grep(".trt", ., value = TRUE)  %>% 
    #   paste0(., collapse = " + ")  -> a_trt_terms
    # a_dose_terms <- ".dose"
    # 
    # nma_form <- stats::update(nma_form, paste0("~. +", a_dose_terms))
    
    attr(nma_form,"dose") <- stringr::str_extract(paste0(deparse(main), collapse = ""), "\\.dose\\([^()]+\\)")
  }
  
  if(any(grepl(".dose",deparse(interactions)))) {
    # stringr::str_split(deparse(interactions), "\\s\\+\\s")[[1]] %>%
    #   stringr::str_trim() %>%
    #   grep(".trt", ., value = TRUE)  %>% 
    #   paste0(., collapse = " + ") -> b_trt_terms
    # b_dose_terms <- ".dose"
    # 
    # nma_form <- stats::update(nma_form, paste0("~. +", b_dose_terms))
    if(isFALSE(attr(nma_form,"dose"))) {
      attr(nma_form,"dose") <- stringr::str_extract(paste0(deparse(interactions), collapse = ""), "\\.dose\\([^()]+\\)")
    } else {
      attr(nma_form,"dose") <- c(attr(nma_form,"dose"),
                                 stringr::str_extract(paste0(deparse(interactions), collapse = ""), "\\.dose\\([^()]+\\)"))
    }
  }
  
  if(any(grepl(".dose",deparse(confounders)))) {
    # stringr::str_split(deparse(interactions), "\\s\\+\\s")[[1]] %>%
    #   stringr::str_trim() %>%
    #   grep(".trt", ., value = TRUE)  %>% 
    #   paste0(., collapse = " + ") -> b_trt_terms
    # c_dose_terms <- ".dose"
    # 
    # nma_form <- stats::update(nma_form, paste0("~. +", c_dose_terms))
    
    if(isFALSE(attr(nma_form,"dose"))) {
      attr(nma_form,"dose") <-  stringr::str_extract(paste0(deparse(confounders), collapse = ""), "\\.dose\\([^()]+\\)")
    } else {
      attr(nma_form,"dose")  <- c(attr(nma_form,"dose"), stringr::str_extract(paste0(deparse(confounders), collapse = ""), "\\.dose\\([^()]+\\)"))
    }
  }
  
  
  return(nma_form)
}

check_collinear_x <- function(X) {
  
  qr.x <- qr(X)
  stopifnot(length(qr.x$pivot) == ncol(X))
  
  return(X[,qr.x$pivot[c(1:qr.x$rank)], drop = FALSE])
  
}

expand_nma_terms <- function(form, var = c(".trt", ".z")) {
  # this function takes a formula and expands any terms with .trt or .z
  # into the full set of terms needed for the nma model matrix
  # e.g., .trt:depressed -> trt1:depressed + trt2:depressed + ...
  # .z -> z.age + z.baseline_severity + ...
  var <- match.arg(var)
  terms <- stats::terms(form)
  
  mapping <- attr(terms, "factors")
  if(length(mapping) == 0) return(form) # no terms to expand
  
  var_terms <-   names(mapping[var,][mapping[var,] == 0])
  trt_terms <-   names(mapping[var,][mapping[var,] == 1])
  if(length(trt_terms) == 0) return(form)
  
  if(attr(terms, "intercept") == 1) {
    var_terms <- c("1", var_terms)
  } else if(attr(terms, "intercept") ==  0) {
    var_terms <- c("0", var_terms)
  }
  
  return(
    rlang::new_formula(rlang::f_lhs(form), 
                            str2lang(
                              paste(c(var_terms,".trt"), collapse = " + ")
                              ))
         )
}

vit_make_nma_model_matrix <- function(nma_formula,
                                  dat_agd_contrast = tibble::tibble(),
                                  agd_contrast_bl = logical(),
                                  xbar = NULL,
                                  consistency = c("consistency", "nodesplit", "ume"),
                                  nodesplit = NULL,
                                  classes = FALSE,
                                  newdata = FALSE) {
  #copied from multinma but adapted to needs here
  # Checks
  if (!rlang::is_formula(nma_formula)) abort("`nma_formula` is not a formula")
  stopifnot(
            is.data.frame(dat_agd_contrast))
  consistency <- rlang::arg_match(consistency)
  if (!rlang::is_bool(classes)) abort("`classes` should be TRUE or FALSE")
  if (nrow(dat_agd_contrast) && !rlang::is_logical(agd_contrast_bl, n = nrow(dat_agd_contrast)))
    rlang::abort("`agd_contrast_bl` should be a logical vector of length nrow(agd_contrast)")
  if (!is.null(xbar) && (
    !(rlang::is_double(xbar) || rlang::is_integer(xbar)) || !rlang::is_named(xbar)))
    rlang::abort("`xbar` should be a named numeric vector")
  
  if (!consistency %in% c("consistency", "ume", "nodesplit")) {
    rlang::abort(glue::glue("Inconsistency '{consistency}' model not yet supported."))
  }
  
  if (consistency == "nodesplit") {
    if (!rlang::is_vector(nodesplit, 2))
      rlang::abort("`nodesplit` should be a vector of length 2.")
  }
  
  .has_agd_contrast <- if (nrow(dat_agd_contrast)) TRUE else FALSE
  
  # Sanitise factors
  if (.has_agd_contrast) {
    dat_agd_contrast <- dplyr::mutate_at(dat_agd_contrast,
                                         .vars = if (classes) c(".trt", ".study", ".trtclass") else c(".trt", ".study"),
                                         .funs = multinma:::fct_sanitise)
  }
  
  # Define contrasts for UME model
  if (consistency == "ume") {
    # For IPD and AgD (arm-based), take the first-ordered arm as baseline
    # So the contrast sign will always be positive
      contrs_arm <- tibble::tibble()
    
    # For AgD (contrast-based), take the specified baseline arm (with .y = NA)
    # Need to make sure to construct contrast the correct way around (d_12 instead of d_21)
    # and then make a note to change the sign if necessary
    if (.has_agd_contrast) {
      contrs_contr <- dat_agd_contrast %>%
        dplyr::distinct(.data$.study, .data$.trt) %>%  # In case integration data passed
        dplyr::left_join(dplyr::distinct(dat_agd_contrast[agd_contrast_bl, ], .data$.study, .data$.trt) %>%
                           dplyr::transmute(.data$.study, .trt_b = .data$.trt), by = ".study") %>%
        dplyr::distinct(.data$.study, .data$.trt, .data$.trt_b) %>%
        dplyr::mutate(.contr_sign = dplyr::if_else(as.numeric(.data$.trt) < as.numeric(.data$.trt_b), -1, 1),
                      .contr = dplyr::if_else(.data$.trt == .data$.trt_b,
                                              "..ref..",
                                              dplyr::if_else(.data$.contr_sign == 1,
                                                             paste0(.data$.trt, " vs. ", .data$.trt_b),
                                                             paste0(.data$.trt_b, " vs. ", .data$.trt))))
    } else {
      contrs_contr <- tibble::tibble()
    }
    
    # Make contrast info
    contrs_all <- dplyr::bind_rows(contrs_arm, contrs_contr)
    
    # Vector of all K(K-1)/2 possible contrast levels
    nt <- nlevels(contrs_all$.trt)
    ctr <- which(lower.tri(diag(nt)), arr.ind = TRUE)
    trt_lev <- levels(contrs_all$.trt)
    c_lev <- paste(trt_lev[ctr[, "row"]], trt_lev[ctr[, "col"]], sep = " vs. ")
    
    contrs_all <- dplyr::transmute(contrs_all,
                                   .data$.study, .data$.trt,
                                   .contr = forcats::fct_drop(factor(.data$.contr, levels = c("..ref..", c_lev))),
                                   .data$.contr_sign)
    
    # Join contrast info on to study data
    if (.has_agd_contrast) {
      dat_agd_contrast <- dplyr::left_join(dat_agd_contrast, contrs_all, by = c(".study", ".trt"))
    }
  }
  
  # Derive .omega indicator for node-splitting model
  if (consistency == "nodesplit") {
    if (.has_agd_contrast) {
      dat_agd_contrast <- dplyr::group_by(dat_agd_contrast, .data$.study) %>%
        dplyr::mutate(.omega = all(nodesplit %in% .data$.trt) & .data$.trt == nodesplit[2]) %>%
        dplyr::ungroup()
    }
  }
  
  # Construct design matrix all together then split out, so that same dummy
  # coding is used everywhere
  dat_all <- dat_agd_contrast
  
  # Check that required variables are present in each data set, and non-missing
  # check_regression_data(nma_formula,
  #                       dat_ipd = dat_ipd,
  #                       dat_agd_arm = dat_agd_arm,
  #                       dat_agd_contrast = dat_agd_contrast,
  #                       newdata = newdata)
  
  # Center
  if (!is.null(xbar)) {
    dat_all[, names(xbar)] <-
      purrr::map2(dat_all[, names(xbar), drop = FALSE], xbar, ~.x - .y)
  }
  
  # Explicitly set contrasts attribute for key variables
  fvars <- all.vars(nma_formula)
  if (".trt" %in% fvars) stats::contrasts(dat_all$.trt) <- "contr.treatment"
  if (".trtclass" %in% fvars) stats::contrasts(dat_all$.trtclass) <- "contr.treatment"
  if (".contr" %in% fvars) stats::contrasts(dat_all$.contr) <- "contr.treatment"
  if (".omega" %in% fvars) stats::contrasts(dat_all$.omega) <- "contr.treatment"
  # .study handled separately next (not always a factor)
  
  # Drop study to factor to 1L if only one study (avoid contrasts need 2 or
  # more levels error)
  if (".study" %in% fvars && dplyr::n_distinct(dat_all$.study) == 1) {
    
    # Save study label to restore
    single_study_label <- unique(dat_all$.study)
    dat_all$.study_temp <- dat_all$.study
    dat_all$.study <- 1L
    
    # Fix up model formula with an intercept (if .study is main effect, which is not usually the case for aux_regression)
    if (".study" %in% colnames(attr(stats::terms(nma_formula), "factors"))) nma_formula <- stats::update.formula(nma_formula, ~. + 1)
  } else {
    single_study_label <- NULL
    if (".study" %in% fvars) stats::contrasts(dat_all$.study) <- "contr.treatment"
  }
  
  # Apply NMA formula to get design matrix
  mf <- stats::model.frame(nma_formula, data = dat_all, na.action = stats::na.pass)
  X_all <-  mf%>% 
    stats::model.matrix(nma_formula, .)
  offsets <- stats::model.offset(mf)
  has_offset <- !is.null(offsets)
  
  disc_names <- setdiff(names(attr(X_all, "contrasts")), c(".study", ".trt", ".trtclass", ".contr", ".omega"))
  
  if (!is.null(single_study_label)) {
    # Restore single study label and .study column
    colnames(X_all) <- stringr::str_replace(colnames(X_all),
                                            "^\\.study$",
                                            paste0(".study", single_study_label))
    dat_all <- dat_all %>%
      dplyr::mutate(.study = .data$.study_temp) %>%
      dplyr::select(-".study_temp")
    
    # Drop intercept column from design matrix
    if ("(Intercept)" %in% colnames(X_all)) X_all <- X_all[, -1, drop = FALSE]
  }
  
  # Remove columns for reference level of .trtclass
  if (classes) {
    ref_class <- levels(dat_all$.trtclass)[1]
    col_trtclass_ref <- grepl(paste0(".trtclass\\Q", ref_class, "\\E"),
                              colnames(X_all), perl = TRUE)
    X_all <- X_all[, !col_trtclass_ref, drop = FALSE]
  }
  
  # Remove columns for reference levels of discrete covariates
  if (length(disc_names)) for (xvar in disc_names) {
    # Check contrast type
    ctype <- attr(dat_all[[xvar]], "contrasts")
    if (is.null(ctype)) ctype <- getOption("contrasts")[if (is.ordered(dat_all[[xvar]])) "ordered" else "unordered"]
    
    # Get reference level for treatment/SAS contrasts
    if (ctype == "contr.treatment") {
      x_ref <-
        if (is.factor(dat_all[[xvar]])) levels(dat_all[[xvar]])[1]
      else if (is.logical(dat_all[[xvar]])) FALSE
      else levels(as.factor(dat_all[[xvar]]))[1]
    } else if (ctype == "contr.SAS") {
      x_ref <-
        if (is.factor(dat_all[[xvar]])) rev(levels(dat_all[[xvar]]))[1]
      else if (is.logical(dat_all[[xvar]])) FALSE
      else rev(levels(as.factor(dat_all[[xvar]])))[1]
    } else {
      x_ref <- NULL
    }
    
    # Remove reference level columns if present
    if (!is.null(x_ref)) {
      col_x_ref <- grepl(paste0("(`?)\\Q", xvar, "\\E(`?)\\Q", x_ref, "\\E"), colnames(X_all), perl = TRUE)
      X_all <- X_all[, !col_x_ref, drop = FALSE]
    }
  }
  
  # Remove columns for interactions with reference level of .trt or .trtclass
  ref_trt <- levels(dat_all$.trt)[1]
  regex_int_ref <- paste0("\\:\\.trt\\Q", ref_trt, "\\E$|^\\.trt\\Q", ref_trt, "\\E\\:")
  if (classes)
    regex_int_ref <- paste0(regex_int_ref, "|",
                            "\\:\\.trtclass\\Q", ref_class, "\\E$|^\\.trtclass\\Q", ref_class, "\\E\\:")
  col_int_ref <- grepl(regex_int_ref, colnames(X_all), perl = TRUE)
  X_all <- X_all[, !col_int_ref, drop = FALSE]
  
  # Remove global intercept column (present for aux_regression)
  X_all <- X_all[, colnames(X_all) != "(Intercept)", drop = FALSE]
  
  if (consistency == "ume") {
    # Set relevant entries to +/- 1 for direction of contrast, using .contr_sign
    contr_cols <- grepl("^\\.contr", colnames(X_all))
    X_all[, contr_cols] <- sweep(X_all[, contr_cols, drop = FALSE], MARGIN = 1,
                                 STATS = dat_all$.contr_sign, FUN = "*")
  }
  
  if (.has_agd_contrast) {
    X_agd_contrast_all <- X_all[1:nrow(dat_agd_contrast), , drop = FALSE]
    offset_agd_contrast_all <-
      if (has_offset) {
        offsets[nrow(dat_ipd) + nrow(dat_agd_arm) + 1:nrow(dat_agd_contrast)]
      } else {
        NULL
      }
    
    # Difference out the baseline arms
    X_bl <- X_agd_contrast_all[agd_contrast_bl, , drop = FALSE]
    if (has_offset) offset_bl <- offset_agd_contrast_all[agd_contrast_bl]
    
    X_agd_contrast <- X_agd_contrast_all[!agd_contrast_bl, , drop = FALSE]
    offset_agd_contrast <-
      if (has_offset) {
        offset_agd_contrast_all[!agd_contrast_bl]
      } else {
        NULL
      }
    
    # Match non-baseline rows with baseline rows by study
    for (s in unique(dat_agd_contrast$.study)) {
      nonbl_id <- which(dat_agd_contrast$.study[!agd_contrast_bl] == s)
      bl_id <- which(dat_agd_contrast$.study[agd_contrast_bl] == s)
      
      bl_id <- rep_len(bl_id, length(nonbl_id))
      
      X_agd_contrast[nonbl_id, ] <- X_agd_contrast[nonbl_id, , drop = FALSE] - X_bl[bl_id, , drop = FALSE]
      if (has_offset) offset_agd_contrast[nonbl_id] <- offset_agd_contrast[nonbl_id] - offset_bl[bl_id]
    }
    
    # Remove columns for study baselines corresponding to contrast-based studies - not used
    s_contr <- unique(dat_agd_contrast$.study)
    bl_s_reg <- paste0("^\\.study(\\Q", paste0(s_contr, collapse = "\\E|\\Q"), "\\E)$")
    bl_cols <- grepl(bl_s_reg, colnames(X_agd_contrast), perl = TRUE)
    
    X_agd_contrast <- X_agd_contrast[, !bl_cols, drop = FALSE]
  } else {
    X_agd_contrast <- offset_agd_contrast <- NULL
  }
  
  return(list(
              X_agd_contrast = X_agd_contrast,
              offset_agd_contrast = offset_agd_contrast))
}

get_se_df <- function(x) {
  s_ij <- x[is.na(x$.y), ".se", drop = TRUE]^2  # Covariances
  s_ii <- x[!is.na(x$.y), ".se", drop = TRUE]^2  # Variances
  n_contr <- length(s_ii)
  se_mat <- data.frame(.cov = rep(s_ij, n_contr))
  return(se_mat)
}

get_se_values <- function(x) {
  return(unclass(
    by(x,
       forcats::fct_inorder(forcats::fct_drop(x$.study)),
       FUN = get_se_df,
       simplify = FALSE)) %>% unlist())
}


construct_nma_data <- function(network, nma_form) {
  
  if(missing(nma_form)) {
    nma_form <- stats::as.formula("~ -1 + .study + .trt")
  }
  
  if(!isFALSE(attr(nma_form,"dose"))) {
    nma_dose_form <- TRUE
  } else{
    nma_dose_form <- FALSE
  }
  
  contrast_dat <- network$agd_contrast
  
  # make model matrix with appropriate contrasts
  X_nma <- vit_make_nma_model_matrix(nma_form, 
                        dat_agd_contrast = contrast_dat,
                        agd_contrast_bl = contrast_dat$.y %>% is.na(),
                        consistency = "consistency")$X_agd_contrast %>% 
    # check_collinear_x() %>%
    # janitor::clean_names() %>% 
    special_clean_names() 
  
  # clean up contrast names
  cn_x <- colnames(X_nma)
  contrast_names <- cn_x %>% grep("^(?!\\.trt)", ., perl = TRUE, value = TRUE)
  colnames(X_nma)[cn_x %in% contrast_names] <- remove_contrast_terms(contrast_names) %>% 
    gsub("^","contrast ",.) %>% snakecase::to_lower_camel_case()
    
  # get data frame
  cols_to_drop_z <- contrast_dat %>% dplyr::filter(!is.na(y)) %>%
    dplyr::select(dplyr::starts_with("z.")) %>%
    dplyr::select(dplyr::where(~ (is.numeric(.x) || is.logical(.x)) &&
                   all(tidyr::replace_na(.x, 0) == 0))) %>%
    names()
  db <- contrast_dat %>% dplyr::filter(!is.na(y)) %>% 
    select(-all_of(cols_to_drop_z))
  # R <- multinma::RE_cor(db$.study, db$.trt, rep(TRUE, nrow(db)), "reftrt")
  # Sigma <- multinma:::make_Sigma(network$agd_contrast) %>% Matrix::bdiag()
  # 
  # R_l  <- t(chol(R))
  # S_l  <- t(chol(Sigma))
  
  db$.cov <- get_se_values(network$agd_contrast)
  db$.obs <- as.factor(1:nrow(db))
  
  if (isTRUE(nma_dose_form )) {
    dose_call <- attr(nma_form,"dose")
    dose_val <- stringr::str_extract(dose_call, "(?<=\\().*(?=\\))")
    
    new_network <- multinma::set_agd_contrast(
      network$agd_contrast,
      study = ".study",
      trt = !!dose_val,
      y = ".y",
      se = ".se",
      trt_class = ".trt",
      sample_size = "final.N",
      trt_ref = NULL
    )
    nma_form2 <- nma_form
    attr(nma_form2,"dose") <- FALSE
    
    dose_X <- construct_nma_data(new_network,
                       nma_form = nma_form2)$X_nma
    dose_X <- dose_X[, grep(".trt", colnames(dose_X)), drop = FALSE]
    colnames(dose_X) <- gsub("\\.trt",".dose", colnames(dose_X))
    X_nma <- cbind(X_nma, dose_X)
    
    dose_prior_ids <- which_CE(new_network$classes)
    dose_prior_ids$varname <- dose_val
  } else {
    dose_prior_ids <- NULL
  }
  
  return(list(X_nma = X_nma,
              dat = db %>% cbind(X_nma),
              # R_l = R_l,
              # Sigma_l = S_l,
              # R = R,
              # S = Sigma,
              dose_prior_ids = dose_prior_ids))
}

remove_contrast_terms <- function(term) {
  if( !any(grepl(".contrast",term)) ) {
    return(term)
  }
  term %>% 
    as.character() %>% 
    gsub("^.contrast","contrast ",.) %>%
    gsub("[(]"," ",.) %>% 
    gsub("[)]"," ",.) %>% 
    snakecase::to_lower_camel_case()
}

split_plus <- function(expr) {
  if (rlang::is_call(expr, "+")) {
    # Get all args of this `+` call
    args <- rlang::call_args(expr)
    # Recursively split each arg and flatten the result
    unlist(lapply(args, split_plus), recursive = FALSE)
  } else {
    list(expr)
  }
}

filter_formula_for_contrast <- function(form) {
  
  if(!any(grepl(".contrast",deparse(form)))) {
    return(form)
  }
  
  if(isbrmsform <- brms::is.brmsformula(form)) {
    orig_form <- form
    form <- form$formula
  }
  
  lhs <- rlang::f_lhs(form) %>% deparse()
  rhs <- rlang::f_rhs(form) %>% split_plus() %>% sapply(deparse) 
  con <- mark_each_contrast(form)
  
  mapply(function(terms, contrast) {
    if (contrast) {
      gsub("mi\\(([^)]*)\\)", "\\1", terms)
    } else {
      terms
    }
  }, terms = rhs, contrast = con, SIMPLIFY = TRUE) -> rhs
    
  
  add_back_mi_flag <- FALSE
  if (any(grepl("| mi()",lhs))) {
    add_back_mi_flag <- TRUE
    lhs <- gsub(" \\| mi\\(.*\\)","", lhs)
  }

  lhs_clean <- remove_contrast_terms(lhs) 
  
  if (add_back_mi_flag) {
    lhs_clean <- paste0(lhs_clean, " | mi()")
  }
  lhs_clean <- lhs_clean %>% str2lang()
  rhs_clean <- sapply(rhs, remove_contrast_terms) %>% paste0(collapse = " + ") %>% str2lang()
  
  if (isbrmsform) {
     out <-  orig_form
     out$formula <- rlang::new_formula(lhs_clean, 
                       rhs_clean)
  } else {
    out <- rlang::new_formula(lhs_clean, 
                       rhs_clean)
  }
  return(
    out
  )
}

formula_terms <- function(form) {
  split_plus(rlang::f_rhs(form))
}

is_mi <- function(term_expr) {
  rlang::is_call(term_expr, "mi")
}

strip_contrast_calls <- function(expr) {
  # If it's .contrast(x) or contrast(x), return x (recursively cleaned)
  if (rlang::is_call(expr, ".contrast") || rlang::is_call(expr, "contrast")) {
    strip_contrast_calls(expr[[2]])
    
    # If it's a "+" call, clean its arguments and rebuild the call
  } else if (rlang::is_call(expr, "+")) {
    args <- rlang::call_args(expr)
    # for "+" this should be length 2, but this is a bit more general
    cleaned <- lapply(args, strip_contrast_calls)
    # rebuild with +
    Reduce(function(x, y) rlang::call2("+", x, y), cleaned)
    
  } else {
    # anything else: leave as-is
    expr
  }
}

strip_contrast_in_formula <- function(form) {
  rhs <- rlang::f_rhs(form)
  new_rhs <- strip_contrast_calls(rhs)
  rlang::new_formula(rlang::f_lhs(form), new_rhs, env = rlang::f_env(form))
}

mark_mi_formula <- function(form) {
  strip_contrast <- function(x) {
    sub("^\\.contrast\\((.*)\\)$", "\\1", x)
  }
  
  if(brms::is.brmsformula(form)) {
    form <- form$formula
  }
  
  terms <- formula_terms(strip_contrast_in_formula(form))
  vapply(terms, is_mi, logical(1))
}

add_back_mi <- function(form, mi_flags) {
  
  if(isbrms <- brms::is.brmsformula(form)) {
    orig_form <- form
    form <- form$formula
  }
  
  terms <- formula_terms(strip_contrast_in_formula(form))
  
  new_terms <- mapply(function(term, is_mi) {
    if (is_mi) {
      rlang::call2("mi", term)
    } else {
      term
    }
  }, terms, mi_flags, SIMPLIFY = FALSE)
  
  new_rhs <- Reduce(function(x, y) rlang::call2("+", x, y), new_terms)
  new_form <- rlang::new_formula(rlang::f_lhs(form), new_rhs, env = rlang::f_env(form))
  
  if (isbrms) {
    out <- orig_form
    out$formula <- new_form
    return(out)
  } else {
    return(new_form)
  }
}

mark_contrast_formula <- function(form) {
  
  if(brms::is.brmsformula(form)) {
    form <- form$formula
  }
  
  rhs <- rlang::f_rhs(form)
  
  terms <- rhs %>% strip_contrast_calls() %>% split_plus() %>% as.character() %>% unlist()
  cf  <- paste0(deparse(rhs), collapse = "")
  # contrast_contents <- stringr::str_match(cf,
  #   stringr::regex("\\.contrast\\((.*)\\)", dotall = TRUE))[, 2] %>% 
  #   strsplit("\\s*\\+\\s*") %>% unlist()
  
  m <- regexpr("\\.contrast\\(", cf)
  if (m[1] == -1L) return(rep(FALSE, length(terms)))
  
  start <- as.integer(m)
  # position just after the opening "("
  pos <- start + attr(m, "match.length")
  
  depth <- 1L
  i <- pos
  n <- nchar(cf)
  
  while (i <= n && depth > 0L) {
    ch <- substr(cf, i, i)
    if (ch == "(") {
      depth <- depth + 1L
    } else if (ch == ")") {
      depth <- depth - 1L
    }
    i <- i + 1L
  }
  
  # `i` is now 1 past the closing ")", so the interior is [pos, i-2]
  inner <- substr(cf, pos, i - 2L)
  
  # Split on "+" into individual terms
  contrast_contents <- strsplit(inner, "\\s*\\+\\s*")[[1]]
  
  terms %in% contrast_contents
}

collapse_contrast <- function(form) {
  if(is_brms <- brms::is.brmsformula(form)) {
    orig_form <- form
    form <- form$formula
  }
  
  
  rhs <- rlang::f_rhs(form)
  new_rhs <- gsub("contrast\\(.*\\)","contrast", paste0(deparse(rhs), collapse = "")) %>%
    str2lang()
  new_form <- rlang::new_formula(rlang::f_lhs(form), new_rhs, env = rlang::f_env(form))
  
  if (is_brms) {
    out <- orig_form
    out$formula <- new_form
    return(out)
  } else {
    return(new_form)
  }
}

mark_each_contrast <- function(form) {
  
  if(brms::is.brmsformula(form)) {
    form <- form$formula
  }
  rhs <- rlang::f_rhs(form)
  x <- paste0(deparse(rhs), collapse = "")
  
  terms <- unlist(strsplit(x, "\\s*\\+\\s*"))
  is_contrast <- grepl("^\\s*\\.contrast\\s*\\(", terms)
  is_contrast
  
}

replace_dotfun_call <- function(expr, dot_fun = ".contrast") {
  if (is.call(expr)) {
    # match either .contrast(...) or contrast(...)
    fn <- expr[[1]]
    
    if (identical(fn, as.name(dot_fun)) || identical(fn, as.name(dot_fun) ) ) {
      # replace the whole call with the bare symbol
      return(as.name(dot_fun))  # keep the dot, per your desired output
    }
    
    # recurse into call arguments (skip [[1]] which is the function name)
    for (i in 2:length(expr)) {
      expr[[i]] <- replace_dotfun_call(expr[[i]], dot_fun = dot_fun)
    }
  }
  expr
}

construct_brms_formula <- function(data_list, main, interactions, confounders, addl_form = NULL, family) {
  X <- data_list$X_nma
  dat <- data_list$dat
  # trt terms
  trt_form <- if(any(grepl("^.trt", colnames(X)))) {
    paste0("`", c(
      grep("^.trt", colnames(X),value = TRUE)
      # , paste0(grep("^trt", colnames(X_nma),value = TRUE),":depressed")
      # , paste0(grep("^trt", colnames(X_nma),value = TRUE),":def:depressed")
    ), "`",
    collapse = " + ") %>% str2lang()
  } else {
    NULL
  }
  
  dose_form <- if(any(grepl("^.dose", colnames(X)))) {
    paste0("`", c(
      grep("^.dose", colnames(X),value = TRUE)
      # , paste0(grep("^trt", colnames(X_nma),value = TRUE),":depressed")
      # , paste0(grep("^trt", colnames(X_nma),value = TRUE),":def:depressed")
    ), "`",
    collapse = " + ") %>% str2lang()
  } else {
    NULL
  }
  
  add_trt_form <- paste0("`",  c(
    grep("^z.", colnames(dat),value = TRUE)
  ),  "`",
  collapse = " + ") %>% str2lang()
  
  #contrast terms
  contrast_form <- if(any(grepl("^contrast", colnames(X)))) {
    paste0("`", c(
      grep("^^contrast", colnames(X),value = TRUE)
      # , paste0(grep("^trt", colnames(X_nma),value = TRUE),":depressed")
      # , paste0(grep("^trt", colnames(X_nma),value = TRUE),":def:depressed")
    ), "`",
    collapse = " + ") %>% str2lang()
  } else {
    NULL
  }
  
  # formula
  bf_mu <- brms::bf(y | vint(as.integer(.study)) + vreal(.se, .cov) ~ main + interactions + confounders + 0, 
                    zU ~ 0 + (1 | .obs),
                    nl = TRUE) + family
  
  # main terms
  main <- stats::as.formula(main) 
  main_mi   <- mark_mi_formula(main)
  main_cont <- mark_contrast_formula(main)
  main <- main %>% replace_dotfun_call(".contrast") %>% 
    replace_dotfun_call(".dose")
  
  main <- main %>% 
    # expand_nma_terms(".trt") %>% 
    list(., list(
      .contrast = contrast_form)) %>%
    do.call(substitute, .) %>% 
    add_back_mi(main_mi & main_cont) %>% 
    list(., list(.trt = trt_form,
                 .z = add_trt_form,
                 .dose = dose_form)) %>% 
    do.call(substitute, .)
  
  bf_mu <- bf_mu + brms::lf(main, cmc = TRUE, center = FALSE)
  
  # interaction terms
  interactions <- stats::as.formula(interactions) 
  interactions_mi   <- mark_mi_formula(interactions) 
  interactions_cont <- mark_contrast_formula(interactions)
  interactions <- interactions %>% replace_dotfun_call(".contrast") %>% 
    replace_dotfun_call(".dose")
  
  interactions <- interactions %>% 
    list(., list(
      .contrast = contrast_form)) %>%
    do.call(substitute, .) %>% 
    add_back_mi(interactions_mi & interactions_cont) %>% 
    list(., list(.trt = trt_form,
                 .z = add_trt_form,
                 .dose = dose_form)) %>% 
    do.call(substitute, .)
  
  bf_mu <- bf_mu +  brms::lf(interactions, cmc = TRUE, center = FALSE)
  
  
  # confounders terms
  confounders <- stats::as.formula(confounders) 
  confounders_mi   <- mark_mi_formula(confounders) 
  confounders_cont <- mark_contrast_formula(confounders)
  confounders <- confounders %>% replace_dotfun_call(".contrast") %>% 
    replace_dotfun_call(".dose")
  # confounders <- gsub("contrast\\([^)]*\\)","contrast", paste0(deparse(confounders), collapse = "")) %>%
  #   stats::as.formula()
  # confounders <- gsub("\\.dose\\([^)]*\\)",".dose", paste0(deparse(confounders), collapse = "")) %>%
  #   stats::as.formula()
  
  confounders <- confounders %>% 
    # expand_nma_terms(".trt") %>% 
    list(., list(
                 .contrast = contrast_form)) %>%
    do.call(substitute, .) %>% 
    add_back_mi(confounders_mi & confounders_cont) %>% 
    list(., list(.trt = trt_form,
                 .z = add_trt_form,
                 .dose = dose_form)) %>% 
    do.call(substitute, .)
  
  bf_mu <- bf_mu +  brms::lf(confounders, cmc = TRUE, center = FALSE)
  
  
  if(!is.null(addl_form)) {
    for (i in seq_along(addl_form)) {
      addl_mi <- mark_mi_formula(addl_form[[i]])
      addl_cont<- mark_each_contrast(addl_form[[i]])
      addl_form[[i]] <- filter_formula_for_contrast(addl_form[[i]]) %>% 
        add_back_mi(addl_mi & addl_cont)
      if(!brms::is.brmsformula(addl_form[[i]])) {
        addl_form[[i]] <- brms::bf(addl_form[[i]])
      }
      bf_mu <- bf_mu + addl_form[[i]]
      
    }
  }
  
  if(!is.null(dose_form)) {
    bf_mu$dose_prior_ids <- data_list$dose_prior_ids
    bf_mu$dose_prior_ids$vars <- grep("^.dose", colnames(X),value = TRUE)
  }
  
  return(bf_mu)
}

get_form_part <- function(bf_mu, part = c("main","interactions","confounders")) {
  part <- match.arg(part)
  form_check <- if (is.null(bf_mu$pforms[[part]])) {
    bf_mu$forms$y$pforms[[part]]
  } else {
    bf_mu$pforms[[part]]
  }
  
  form_check
}



construct_brms_prior <- function(bf_mu) {
  
  check_intercept <- function(bf_mu, part = c("main","interactions","confounders")) {
    get_form_part(bf_mu, part)  %>% stats::terms() %>%  attributes() %>% .$intercept %>%  `==`( 0)
  }
  
  # if(!brms::is.brmsformula(bf_mu)) {
  #   stop("bf_mu must be a brmsformula object")
  # }
  
  if(!is.null(bf_mu$forms)) {
    resp <- "y"
  } else {
    resp <- ""
  }
  
  # heterogeneity prior
  bf_prior <- brms::set_prior("normal(0,0.5)",  class = "tau", resp = resp) # heterogeneity prior; always need this
  
  # re priors
  bf_prior <- bf_prior + brms::set_prior(
    "constant(1)",
    class = "sd",
    group = ".obs",
    dpar = "zU",
    resp = resp
  )
  
  # confounders terms
  confounders_part <- get_form_part(bf_mu, "confounders")
  bf_prior <- bf_prior + # add horsehoe if a terms exist
    if(!is.null(confounders_part) && length(brms:::brmsterms(confounders_part)$dpars$mu$fe %>% all.vars())>0) {
      brms::set_prior("horseshoe(df = 3, scale_global = 0.1, scale_slab = 1)",class = "b", nlpar = "confounders", resp = resp)
      # brms::set_prior("normal(0,0.25)",class = "b", nlpar = "confounders", resp = resp)
    } else {NULL}
  bf_prior <- bf_prior + # add sd priors if re terms exist
    if(!is.null(confounders_part) && length(brms:::brmsterms(confounders_part)$dpars$mu$re)>0) {
      brms::set_prior("normal(0.0,0.25)",  class = "sd", nlpar = "confounders", resp = resp)
    } else {NULL}
  
  # interactions terms
  interactions_part <- get_form_part(bf_mu, "interactions")
  bf_prior <- bf_prior + # add horsehoe if a terms exist
    if(!is.null(interactions_part) && length(brms:::brmsterms(interactions_part)$dpars$mu$fe %>% all.vars())>0) {
      # brms::set_prior("horseshoe(df = 1, scale_global = .05, scale_slab = 1)",class = "b", nlpar = "interactions", resp = resp)
      brms::set_prior("normal(0,0.25)",class = "b", nlpar = "interactions", resp = resp)
    } else {NULL}
  bf_prior <- bf_prior + # add sd priors if re terms exist
    if(!is.null(interactions_part) && length(brms:::brmsterms(interactions_part)$dpars$mu$re)>0) {
      brms::set_prior("normal(0.0,0.5)",  class = "sd", nlpar = "interactions", resp = resp)
    } else {NULL}
  if (check_intercept(bf_mu, "interactions")) {
    bf_prior <- bf_prior + brms::set_prior("constant(0)", coef = "Intercept", nlpar = "interactions", resp = resp)
    # bf_prior <- bf_prior + brms::set_prior(
    #   "normal(0, 1e-6)", coef = "Intercept", nlpar = "b", resp = resp
    # )
    if (!is.null(bf_mu$pforms$interactions)) {
      bf_mu$pforms$interactions <- stats::update.formula(bf_mu$pforms$interactions, . ~ . + 1)
    } else {
      bf_mu$forms$y$pforms$interactions<- stats::update.formula(bf_mu$forms$y$pforms$interactions, . ~ . + 1)
    }
  } 
  
  # main terms
  main_part <- get_form_part(bf_mu, "main")
  bf_prior <- bf_prior + # add normal prior if b fe terms exist
    if (!is.null(main_part) && 
       length(brms:::brmsterms(main_part)$dpars$mu$fe %>% all.vars())>0) {
      brms::set_prior("normal(0.0,2.5)", class = "b",nlpar = "main", resp = resp)
    } else {NULL}
  bf_prior <- bf_prior + # add sd priors if re terms exist in b
    if (!is.null(main_part) && 
       length(brms:::brmsterms(main_part)$dpars$mu$re)>0) {
      brms::set_prior("normal(0.0,0.5)",  class = "sd", nlpar = "main", resp = resp)
    } else {NULL}
  # bf_prior <- bf_prior + brms::set_prior("constant(0)", coef = "Intercept", nlpar = "b", resp = resp)
  if (check_intercept(bf_mu, "main")) {
    bf_prior <- bf_prior + brms::set_prior("constant(0)", coef = "Intercept", nlpar = "main", resp = resp)
    # bf_prior <- bf_prior + brms::set_prior(
    #   "normal(0, 1e-6)", coef = "Intercept", nlpar = "b", resp = resp
    # )
    if (!is.null(bf_mu$pforms$main)) {
      bf_mu$pforms$main <- stats::update.formula(bf_mu$pforms$main, . ~ . + 1)
    } else {
      bf_mu$forms$y$pforms$main<- stats::update.formula(bf_mu$forms$y$pforms$main, . ~ . + 1)
    }
  } 
  # else {
  #   bf_prior <- bf_prior + brms::set_prior(
  #     "normal(0, 1e-6)", coef = "Intercept", nlpar = "b", resp = resp
  #   )
  # }
  
  
  # if(a_part %>% stats::terms() %>%  attributes() %>% .$intercept %>%  `==`( 0)) {
  #   bf_prior <- bf_prior + brms::set_prior("constant(0)", coef = "Intercept", nlpar = "a")
  #   bf_mu$pforms$a <- stats::update.formula(bf_mu$pforms$a, a ~ . + 1)
  # }
  
  # addl terms
  if(length(bf_mu$forms) > 0) {
    nms <- names(bf_mu$forms)
    nms <- nms[ nms != "y"]
    
    for (nm in nms) {
      bf_prior <- bf_prior + 
        brms::set_prior("normal(0,0.25)", class = "b",resp = nm) +
        brms::set_prior("normal(0,0.25)", class = "Intercept", resp = nm)
      
      if ( !is.null(brms::brmsterms(bf_mu$forms[[nm]])$dpars$mu$re) ) {
        bf_prior <- bf_prior + brms::set_prior("normal(0,0.25)",  class = "sd", resp = nm)
      }
      
      if ( brms::brmsterms(bf_mu$forms[[nm]])$family$family == "gaussian") {
        bf_prior <- bf_prior +
          brms::set_prior("exponential(1)", class = "sigma", resp = nm)
      }
    }
    
  }
  
  if (!is.null(bf_mu$dose_prior_ids) & !is.null(main_part)) {
    dose_vars <- bf_mu$dose_prior_ids$vars
    ids       <- bf_mu$dose_prior_ids$id
    
    # will go along class names and add in specific priors for dose effects
    for(i in seq_along(ids)) {
      class_id <- ids[i]
      bf_prior <- bf_prior + 
        if(class_id == 0) {
        brms::set_prior("constant(0)", coef = dose_vars[i], nlpar = "main", resp = resp)
      } else {
        brms::set_prior(
          paste0("normal(0,sd_class[", class_id, "])"),
          coef = dose_vars[i],
          nlpar = "main",
          resp = resp)
      }
    }
    
  }
  
  
  return(list(bf_prior = bf_prior, bf_mu = bf_mu))
}

construct_stanvar <- function(data_list, bf_mu, family_list) {
  if(!is.null(bf_mu$forms)) {
      resp <- "_y"
    } else {
      resp <- ""
    }
    
  stan_fun <- family_list$stan_fun
  stanvar <- brms::stanvar(scode = stan_fun, block = "functions") +
    # brms::stanvar(scode = "vector[N] z_u;  // latent standard-normal random effects (stacked)",
    #               block = "parameters", position = "end")
    brms::stanvar(scode = covariance_functions_stan, block = "functions")
    
  stanvar <- stanvar + if(!is.null(data_list$dose_prior_ids)) {
      brms::stanvar(
        scode = "
        vector<lower = 0>[N_class] sd_class_raw;
        ",
        block = "parameters", position = "end") +
      brms::stanvar(
        scode = "lprior += std_normal_lpdf(sd_class_raw);", block = "tparameters", position = "end") +
      brms::stanvar(
        scode = "vector[N_class] sd_class = sd_class_raw * 0.25;", block = "tparameters", position = "start")
    } else {
      NULL
    }
    
  stanvar <- stanvar + 
    # brms::stanvar(
    #   scode = "lprior += std_normal_lpdf(z_u);", block = "tparameters", position = "end") +
    # brms::stanvar(scode = "
    #   vector[N] u_draw = tau * (R_L * z_u);",
    #               block = "genquant", position = "start"
    # ) +
    brms::stanvar(
      scode = glue::glue("matrix[N,N] Sigma = make_group_cov(vint1{resp}, vreal1{resp}, vreal2{resp});
      matrix[N,N] R = make_group_cor(vint1{resp});
      matrix[N,N] Sigma_L = cholesky_decompose(Sigma);
      matrix[N,N] R_L = cholesky_decompose(R);"),
      block = "tdata", position = "end"
    ) +
    # brms::stanvar(
    #   x =  data_list$Sigma_l,
    #   name = "Sigma_L", block = "data") +
    # brms::stanvar(
    #   x =  data_list$R_l,
    #   name = "R_L", block = "data") +
#     brms::stanvar(
#       x =  as.matrix(data_list$S),
#       name = "Sigma2", block = "data") +
#     brms::stanvar(
#       x =  as.matrix(data_list$R),
#       name = "R2", block = "data") +
#     brms::stanvar(
#       scode = '
#       // suppose A and B are two matrices you want to compare
# {
#   real max_diff = max(abs(Sigma2 - Sigma));
# 
#    print("max diff = ", max_diff);
# 
#   if (max_diff > 1e-12) {
#     print("Matrices differ! max diff = ", max_diff);
#   } else {
#     print("Matrices are equal within tolerance");
#   }
# }
# 
# {
#   real max_diff = max(abs(R2 - R));
# 
#   print("max diff = ", max_diff);
# 
#   if (max_diff > 1e-12) {
#     print("Matrices differ! max diff = ", max_diff);
#   } else {
#     print("Matrices are equal within tolerance");
#   }
# }
#       ', block = "tdata", position = "end"
#     ) +
    if(!is.null(data_list$dose_prior_ids)) {
      brms::stanvar(
        x =  as.integer(max(data_list$dose_prior_ids$id)),
        name = "N_class", block = "data")
    } else {
      NULL
    }
  
  # ll_fun <- switch(family_list$family_name,
  #                  "gaussian" = glue::glue(
  #                    "mvn_nma_2_lpdf(Y{resp} | mu{resp}, tau{resp}, vint1{resp}, Sigma_L, R_L);"),
  #                  "student-t" = glue::glue("mvt_nma_2_lpdf(Y{resp} | mu{resp}, tau{resp}, vint1{resp}, Sigma_L, R_L);")
  # )
  # 
  # stanvar <- stanvar + brms::stanvar(
  #   scode = " real log_likelihood = log_lik_full;",
  #                 block = "genquant", position = "start"
  #   ) +
  #   brms::stanvar(scode = glue::glue("log_lik_full = {ll_fun}"),
  #                                    block = "likelihood", position = "end"
  #   ) +
  #   brms::stanvar(scode = "real log_lik_full;",
  #                         block = "tparameters", position = "end"
  #   )
  
  return(stanvar)
}

construct_combo_ledger <- function(data_list, network, main, interactions, confounders, addl_form) {
  
  get_interaction_vars <- function(trt, formula) {
    terms         <- stats::terms(formula)
    term.factors  <- attr(terms, "factors")
    if(length(term.factors) == 0) return(NULL)
    trt_found     <- which(term.factors[trt,] > 0)
    interaction_l <- apply(term.factors[,trt_found,drop=FALSE],2,function(x) x>0) %>% 
                        rowSums() %>% 
                        `>`(0)
    setdiff(all.vars(terms)[interaction_l],trt)
  }
  
  get_vars <- function(formula) {
    brms::brmsterms(formula)$allvars %>% all.vars()
  }
  
  dat  <- data_list$dat #%>% filter(!is.na(y))
  original_data <- network$agd_contrast #%>% filter(!is.na(y))
  
  # if(any(grepl(".trt",colnames(dat)))) {
  #   trt  <- grep("^\\.trt[^.]*$", colnames(dat), value = TRUE)
  #   trt_var <- ".trt"
  # } else {
  #   trt  <- grep("^z\\.", colnames(dat), value = TRUE)
  #   trt_var <- ".z"
  # }
  # interactions <- c(
  #   get_interaction_vars(trt_var,main),
  #   get_interaction_vars(trt_var,interactions)
  # ) %>% unique()
  
  check_trt <- c(
    get_vars(main),
    get_vars(interactions),
    get_vars(confounders),
    unlist(lapply(addl_form, get_vars))
  ) %>% unique()
  
  if(any(grepl(".trt",check_trt))) {
    trt  <- grep("^\\.trt[^.]*$", colnames(dat), value = TRUE)
    trt_var <- ".trt"
  } else {
    trt  <- grep("^z\\.", colnames(dat), value = TRUE)
    trt_var <- ".z"
    check_trt <- check_trt[check_trt != ".z"]
    check_trt <- c(check_trt, grep("^z\\.", colnames(dat), value = TRUE))
  }
  
  keep <- intersect(check_trt, names(dat))
  out <- original_data %>% 
    dplyr::mutate(
      dplyr::across(dplyr::any_of(keep), \(x)
                    if (is.factor(x)) as.character(x) else x),
      .keep = "none"
    ) %>%
    # cbind(original_data[trt_var]) %>% 
    dplyr::distinct()  %>% 
    { if (".trt" %in% names(.) && trt_var == ".trt") filter(., .trt != "placebo") else . } 
  
  if(trt_var == ".z") {
    z_cols <- grep("^z\\.", colnames(dat), value = TRUE)
    out <- out %>%
      mutate(
        .orig_row = dplyr::row_number(),
        .n_on = rowSums(dplyr::across(dplyr::all_of(z_cols)), na.rm = TRUE)
      ) %>%
      {                               # <- prevent implicit LHS into bind_rows()
        dplyr::bind_rows(
          # keep rows with 0 or 1 "on"
          dplyr::filter(., .n_on <= 1) %>% dplyr::select(-.n_on),
          
          # expand rows with >1 "on"
          dplyr::filter(., .n_on > 1) %>%
            tidyr::pivot_longer(dplyr::all_of(z_cols), names_to = "var", values_to = "value") %>%
            dplyr::filter(value == 1) %>%
            dplyr::select(-value, -.n_on) %>%
            dplyr::mutate(..id = 1:dplyr::n(), 
                          value = 1L,
                   var = factor(var, levels = z_cols)) %>%   # ensure all z_cols exist
            tidyr::pivot_wider(
              # id_cols = c(..id, .orig_row),
              names_from = var, values_from = value,
              values_fill = 0, names_expand = TRUE
            )
        )
      } %>%
      dplyr::arrange(.orig_row) %>%
      dplyr::select(-c(..id, .orig_row)) %>% 
      dplyr::distinct()
    
  }

  factors <- sapply(
    original_data %>% dplyr::select(all_of(colnames(out))), is.factor
  )
  if (any(factors)) {
    out <- out %>%
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(names(factors)[factors]),
          \(x) factor(x,
                      levels = levels(original_data[[dplyr::cur_column()]]))
        )
      )
  }
  return(out)
}

construct_family_nma <- function(family =c("gaussian","student-t")) {
  family_name <- match.arg(family, choices = c("gaussian","student-t"))
  if(family_name == "gaussian") {
    family <- mvn_nma
    stan_fun <- mvn_stan_fun
  } else if(family_name == "student-t") {
    family <- mvt_nma
    stan_fun <- mvt_stan_fun
  }
  return(list(family = family, stan_fun = stan_fun, family_name = family_name))
}

which_CE <- function(classes)   {
  # from multinma
  # Class vector, without network reference treatment
  x <- classes[-1]
  
  # Identify sole occupancy classes
  solo_classes <- levels(x)[table(x) == 1]
  
  # Set sole occupancy classes to NA (no class effects) and drop unused levels
  x <- droplevels(x, exclude = solo_classes)
  
  # Create numeric ID vector (0 = no class effect)
  id <- as.numeric(x)
  id[is.na(id)] <- 0
  
  # Create class labels
  label <- levels(x)
  
  return(list(id = id, label = label))
}

prep_brms_nma <- function(data, 
                          family = c("gaussian","student-t"),
                         main = main ~ 0 + .trt,
                         interactions = interactions ~ 0,
                         confounders = confounders ~ 0,
                         addl_form = NULL,
                         trt_ref = NULL) {
  
  this_call    <- match.call()
  if(missing(data)) stop("argument `data` is required")
  orig_form    <- list(main = main, 
                       interactions = interactions, 
                       confounders = confounders,
                       addl_form = addl_form)
  
  # set up data
  network      <- construct_nma_network(data, trt_ref = trt_ref)
  
  # make model matrix from multinma
  nma_form     <- construct_nma_formula(main, 
                                        interactions, 
                                        confounders)
  
  # get data list
  data_list    <- construct_nma_data(network, nma_form)
  
  # set family
  family_list  <- construct_family_nma(family)
  
  # make brms formula
  bf_mu        <- construct_brms_formula(data_list, 
                                         main, 
                                         interactions, 
                                         confounders,
                                         addl_form, 
                                         family_list$family)
  
  # priors  
  prior_list   <- construct_brms_prior(bf_mu)
  bf_prior     <- prior_list$bf_prior
  bf_mu        <- prior_list$bf_mu
  
  # stanvars
  stanvar      <- construct_stanvar(data_list, bf_mu, family_list)
    
  # check combos of covariates
  combo_ledger <- construct_combo_ledger(data_list,
                                    network,
                                    main, interactions,
                                    confounders,
                                    addl_form)
  
  return(
    list(
      data = data_list$dat,
      family = family_list$family,
      formula = bf_mu,
      priors = bf_prior,
      X_nma = data_list$X_nma,
      # R_l = data_list$R_l,
      # Sigma_l = data_list$Sigma_l,
      # Sigma = data_list$Sigma,
      # R = data_list$R,
      stanvar = stanvar,
      network = network,
      orig_form = orig_form,
      orig_family = family,
      combos = combo_ledger,
      stan_fun = family_list$stan_fun,
      reference_treatment = levels(network$treatments)[1],
      class_ids = data_list$dose_prior_ids,
      call = this_call
    )
  )
}

dose_hack <- function(...) {
  dots <- list(...)
  
  dots$empty <- TRUE
  dots$backend <- "rstan"
  
  fit <- do.call(brms::brm, dots)
  
  if(!is.null(dots$silent) && dots$silent != 2) message("Fitting required rstan backend with non-centered parameterization for sd_class priors for dose response.")
  
  # verbosity
  if (!is.null(dots$silent) && dots$silent == 2) {
    open_progress <- FALSE
    show_messages <- FALSE
  } else if (!is.null(dots$silent) && dots$silent == 1) {
    open_progress <- FALSE
    show_messages <- TRUE
  } else {
    open_progress <- TRUE
    show_messages <- TRUE
  }
  
  # rewrite stan code to non-centered parameterization for sd_class priors
  stanc <- fit %>% brms::make_stancode()
  stand <- fit %>% brms::make_standata()
  mod_sc <- rewrite_stan_noncenter_sd_class(stanc)
  sm <- rstan::stan_model(model_code = mod_sc)
  stan_fit <- rstan::sampling(sm, data = stand, chains = dots$chains, 
                              iter = dots$iter,
                              warmup = dots$warmup,
                              cores = dots$cores,
                              seed = dots$seed,
                              control = fit$stan_args$control,
                              open_progress = open_progress,
                              show_messages = show_messages,
                              refresh = dots$refresh)
  fit$fit <- stan_fit
  fit$model <- mod_sc
  
  return(brms::rename_pars(fit))
}

#' Fit Bayesian multivariate meta-analysis model
fit_brms_nma <- function(prep, 
                         iter = 2000, warmup = 1000, chains = 4, cores = 4,
                         control = list(adapt_delta = 0.95, max_treedepth = 15),
                         seed = NA,
                         ...) {
  
  
  # fit model
  fit <- if (is.null(prep$class_ids ) ) {
      brms::brm(
        formula = prep$formula,
        data = prep$data,
        prior = prep$priors,
        stanvars = prep$stanvar,
        iter = iter,
        warmup = warmup,
        chains = chains,
        cores = cores,
        control = control,
        seed = seed,
        ...
      )
  } else {
    dose_hack(
        formula = prep$formula,
        data = prep$data,
        prior = prep$priors,
        stanvars = prep$stanvar,
        iter = iter,
        warmup = warmup,
        chains = chains,
        cores = cores,
        control = control,
        seed = seed,
        ...
      )
  }
  
  fit$prep <- prep
  
  class(fit) <- c("vitfit", class(fit))
  return(fit)
}

rewrite_stan_noncenter_sd_class <- function(code) {
  # If code is one big string, split into lines
  if (length(code) == 1L && grepl("\n", code, fixed = TRUE)) {
    code <- strsplit(code, "\n", fixed = TRUE)[[1]]
  }
  
  # Helper to escape regex metacharacters
  escape_regex <- function(x) {
    gsub("([][(){}+*.^$|?\\\\])", "\\\\\\1", x)
  }
  
  # Regex for priors of the form:
  # lprior += normal_lpdf(<coef> | 0, sd_class[k]);
  prior_rx <- "^\\s*lprior\\s*\\+=\\s*normal_lpdf\\(([^|]+)\\|\\s*0(?:\\.0)?\\s*,\\s*(sd_class\\[\\d+\\])\\s*\\);"
  
  m <- regexec(prior_rx, code)
  matches <- regmatches(code, m)
  has_match <- lengths(matches) > 0L
  if (!any(has_match)) {
    message("No sd_class-based normal priors found to rewrite.")
    return(code)
  }
  
  match_lines <- which(has_match)
  
  for (line_idx in match_lines) {
    parts <- matches[[line_idx]]
    # parts: full_match, coef_expr, scale_expr
    coef_expr_raw <- parts[2]            # e.g. "b_main[60] "
    scale_expr    <- parts[3]            # e.g. "sd_class[6]"
    
    coef_expr <- trimws(coef_expr_raw)
    
    # Find assignment line:
    # <coef_expr> = <param_name>;
    assign_rx <- sprintf(
      "^\\s*%s\\s*=\\s*([A-Za-z][A-Za-z0-9_\\.]*)\\s*;",
      escape_regex(coef_expr)
    )
    assign_matches <- regexec(assign_rx, code)
    assign_m <- regmatches(code, assign_matches)
    
    assign_lines <- which(lengths(assign_m) > 0L)
    # We only care about lines where coef_expr is on the LHS
    if (length(assign_lines) == 0L) {
      warning(sprintf("No assignment line found for '%s'. Skipping this prior.", coef_expr))
      next
    }
    # Take the first matching assignment
    assign_line <- assign_lines[1L]
    assign_parts <- assign_m[[assign_line]]
    param_name <- assign_parts[2]   # captured RHS symbol
    
    # ---- 1) Rewrite prior line: std_normal on param_name ----
    code[line_idx] <- sprintf("  lprior += std_normal_lpdf(%s);", param_name)
    
    # ---- 2) Rewrite assignment: coef_expr = param_name * sd_class[k]; ----
    code[assign_line] <- sprintf("  %s = %s * %s;", coef_expr, param_name, scale_expr)
  }
  
  new_code <- paste(code, collapse = "\n")
  new_code
}

rewrite_stan_efficient_data <- function(code) {
  # If code is one big string, split into lines
  if (length(code) == 1L && grepl("\n", code, fixed = TRUE)) {
    code <- strsplit(code, "\n", fixed = TRUE)[[1]]
  }
  
  # find Sigma_L and R_L data declarations
  sigma_l_rx <- "^\\s*matrix.*\\s+Sigma_L\\s*;"
  r_l_rx     <- "^\\s*matrix.*\\s+R_L\\s*;"
  sigma_l_line <- which(grepl(sigma_l_rx, code))
  r_l_line     <- which(grepl(r_l_rx, code))
  
  # find Sigma and R data declarations
  sigma_rx <- "^\\s*matrix.*\\s+Sigma\\s*;"
  r_rx     <- "^\\s*matrix.*\\s+R\\s*;"
  sigma_line <- which(grepl(sigma_rx, code))
  r_line     <- which(grepl(r_rx, code))
  
  # change numbers to variable N
  if (length(sigma_l_line) > 0) {
    code[sigma_l_line] <- gsub("matrix\\s*\\[\\s*\\d+\\s*,\\s*\\d+\\s*\\]", "matrix[N, N]", code[sigma_l_line])
  }
  if (length(r_l_line) > 0) {
    code[r_l_line] <- gsub("matrix\\s*\\[\\s*\\d+\\s*,\\s*\\d+\\s*\\]", "matrix[N, N]", code[r_l_line])
  }
  if (length(sigma_line) > 0) {
    code[sigma_line] <- gsub("matrix\\s*\\[\\s*\\d+\\s*,\\s*\\d+\\s*\\]", "matrix[N, N]", code[sigma_line])
  }
  if (length(r_line) > 0) {
    code[r_line] <- gsub("matrix\\s*\\[\\s*\\d+\\s*,\\s*\\d+\\s*\\]", "matrix[N, N]", code[r_line])
  } 
  new_code <- paste(code, collapse = "\n")
  
  new_code
}

drop_desired_study <- function(fit, study_i) {
  drop.idx <- which(fit$prep$data$.study == study_i)
  
  data <- brms::standata(fit, 
                 newdata = fit$data %>% filter(.study != study_i),
                 allow_new_levels = TRUE)
  
  # data$Sigma_L <- data$Sigma_L[-drop.idx, -drop.idx, drop = FALSE]
  # data$R_L     <- data$R_L[-drop.idx, -drop.idx, drop = FALSE]
  # data$Sigma   <- data$Sigma[-drop.idx, -drop.idx, drop = FALSE]
  # data$R       <- data$R[-drop.idx, -drop.idx, drop = FALSE]
  
  return(data)
  
  # if(any(drop.idx > data$N) ){
  #   rlang::abort(sprintf("Study %s index %d exceeds number of studies %d", 
  #                   study_i, drop.idx, data$N))
  # }
  # if(length(drop.idx) == 0) {
  #   warning(sprintf("Study %s not found in data. Returning original data.", study_i))
  #   return(data)
  # }
  # data$N <- data$N - length(drop.idx)
  # data$Y <- data$Y[-drop.idx, drop = FALSE]
  # 
  # data$vint1 <- data$vint1[-drop.idx, drop = FALSE]
  # 
  # data$X_main <- data$X_main[-drop.idx, , drop = FALSE]
  # if (!is.null(data$X_interactions)) {
  #   data$X_interactions <- data$X_interactions[-drop.idx, , drop = FALSE]
  # }
  # if (!is.null(data$X_confounders)) {
  #   data$X_confounders <- data$X_confounders[-drop.idx, , drop = FALSE]
  # }
  # data$Sigma_L <- data$Sigma_L[-drop.idx, -drop.idx, drop = FALSE]
  # data$R_L     <- data$R_L[-drop.idx, -drop.idx, drop = FALSE]
  # data$Sigma   <- data$Sigma[-drop.idx, -drop.idx, drop = FALSE]
  # data$R       <- data$R[-drop.idx, -drop.idx, drop = FALSE]
  # 
  # return(data)
}

get_desired_study_raw_data <- function(fit, studies, study_i) {
  drop.idx <- which(studies == study_i)
  newdata     <- fit$prep$data[drop.idx,,drop = FALSE]
  
  return(newdata)
}

# recover_interaction_bases <- function(interaction_cols, formula, data) {
#   # 1) tokens after '..' from your column names
#   tokens <- unique(unlist(lapply(interaction_cols, function(s) {
#     strsplit(s, "\\.\\.")[[1]]
#   })))
#   
#   if (length(tokens) == 0) return(character(0))
#   
#   terms <- sapply(formula, function(x) stats::terms(x)  %>% 
#                     attr(.,"variables") %>% 
#                     as.character() %>% 
#                     .[-1:-2]
#                     ) %>% unlist() %>% unique()
#   
#   inter_terms <- terms[sapply(terms, function(x) any(grepl(x, tokens)))]
#   
#   if (length(inter_terms) == 0) {
#     list(terms = character(0),
#          levels = list())
#   }
#   
#   if( all(inter_terms %in% colnames(data))) {
#    levels <- sapply(inter_terms, function(x) {
#       if(is.factor(data[[x]])) {
#         levels(data[[x]])
#       } else if(is.character(data[[x]])) {
#         unique(data[[x]])
#       } else if(is.numeric(data[[x]])) {
#         NULL
#       } else {
#         stop(sprintf("Variable %s is of unsupported type: %s", x, class(data[[x]])))
#       }
#     }, simplify = FALSE)
#     
#   } else {
#     levels <- sapply(inter_terms, function(x) {
#       cur_token <- tokens[grepl(x, tokens)]
#       if(length(cur_token) == 0) return(NULL)
#       if(length(cur_token) == 1) return(c(".base_val",".other_val"))
#       str_split(cur_token, x) %>% 
#         lapply(function(y) y[-1]) %>% 
#         unlist() %>% 
#         unique() %>% 
#         c(".base_val",.) %>%
#         return()
#     })
#   }
#     
#     if(!is.null(levels$.trt )) {
#       levels$.trt <- levels$.trt[levels$.trt != ".base_val"]
#     }
#     if (!is.null(levels$`z.` )) {
#       levels$`z.`  <- levels$`z.` [levels$`z.`  != ".base_val"]
#   }
#   
#   
#   return(list(terms = inter_terms,
#               levels = levels)
#    )
#   
# }
# 
# all_possible_contrasts <- function(txvar, main_effects, actual_inter, interaction_levels_and_terms) {
#   
#   all_interactions <- function(vars, sep = "..") {
#     out <- character()
#     n <- length(vars)
#     for (k in 1:n) {
#       combos <- utils::combn(vars, k, simplify = FALSE)
#       out <- c(out, vapply(combos, function(x) paste(x, collapse = sep), character(1)))
#     }
#     unique(out)
#   }
#   
#   # trt_levels is a vector of treatment levels
#   # interaction_levels is a named list of interaction levels
#   # e.g., list(def = c("no","yes"), depressed = c("no","yes"))
#   interaction_levels <- interaction_levels_and_terms$levels
#   interaction_terms <- interaction_levels_and_terms$terms
#   if(length(interaction_levels) == 0) {
#     return(trt_levels)
#   }
#   
#   inter_names <- names(interaction_levels)
#   inter_combos <- expand.grid(interaction_levels, stringsAsFactors = FALSE)
#   
#   contrasts <- matrix(0, nrow =  nrow(inter_combos), ncol = length(main_effects) + length(actual_inter))  
#   colnames(contrasts) <- c(main_effects, actual_inter)
#   
#   for ( i in 1:nrow(inter_combos) ) {
#     cur_combo <- inter_combos[i,]
#     cur_combo <- cur_combo[,cur_combo != ".base_val",drop = FALSE]
#     vars      <- colnames(cur_combo)
#     levels    <- sub(".other_val","",cur_combo)
#     cur_name  <- paste0(vars, levels)
#     cur_var   <- all_interactions(cur_name)
#     if ( any(cur_var %in% colnames(contrasts)) && any(cur_var %in% actual_inter) ) {
#       cur_names <- unique(c(cur_name, cur_var))
#       cur_names <- cur_names[cur_names %in% colnames(contrasts)]
#       if(length(cur_names) > 0) contrasts[i, cur_names] <- 1
#     }
#   }
#   keep    <- rowSums(contrasts) > 0
#   
#   base_vals <- inter_combos %>% 
#     select(-starts_with(txvar)) %>% 
#     filter(if_all(everything(), ~ .x == ".base_val")) %>% 
#     .[1,]
#   rownames(base_vals) <- NULL
#   
#   factors <- rbind(inter_combos[keep, , drop = FALSE],
#                    cbind(data.frame(.trt = gsub(txvar, "", main_effects)), base_vals))
#   return(list(factors = factors,
#               matrix = contrasts[keep, , drop = FALSE]))
# }
# 
# verify_with_data <- function(interaction, data) {
#   # this function verifies the combinations that exist in the data
#   # are also present in the interaction matrix
#   # and calibrates the interaction matrix to only include those levels
#   # it also will check that the estimated effect matches the 
#   # combo in the data, ie, the main effect for one treatment
#   # may actually correspond to an interaction but there was 
#   # only one row in the data for that combo
#   # so we need to make sure the listed factors actually reflect
#   # the data
#   if(nrow(interaction$factors) == 0) return(interaction)
#   if(nrow(interaction$matrix) == 0) return(interaction)
#   
#   # make all factors lower case for matching
#   # this assumes data is not case sensitive
#   # if this is not the case, then need to modify
#   # this function to be case sensitive
#   # and ensure interaction$factors is also case sensitive
#   data_factors <- data %>%
#     select(all_of(colnames(interaction$factors))) %>%
#     distinct() %>% 
#     mutate(across(where(is.factor), function(x) x %>% as.character() %>% tolower())) %>% 
#     mutate(.trt = .trt %>%
#              str_replace_all(",","") %>% 
#              str_replace_all("-","_") %>% 
#              str_replace_all(" ","_") %>% 
#              str_replace_all(":",".."))
#   
#   interaction$factors <- interaction$factors %>%
#     mutate(across(where(is.factor), function(x) x %>% as.character()))
#   
#   # replace .base_val and .other_val with actual values from data
#   for (col in colnames(interaction$factors)) {
#     if (col %in% colnames(data_factors)) {
#       base_val <- data_factors %>%
#         select(!!col) %>% 
#         distinct() %>%
#         dplyr::pull() %>% 
#         setdiff(interaction$factors %>%
#                   select(!!col) %>% 
#                   distinct() %>%
#                   dplyr::pull()) %>% 
#         sort(decreasing = TRUE)
#       other_vals <- interaction$factors %>%
#         select(!!col) %>% 
#         distinct() %>%
#         # filter(.[[col]] == ".other_val") %>%
#         dplyr::pull() %>% 
#         setdiff(data_factors %>%
#                   select(!!col) %>% 
#                   distinct() %>%
#                   dplyr::pull())
#       if (length(base_val) == 1) {
#         interaction$factors[[col]][interaction$factors[[col]] == ".base_val"] <- base_val
#       }
#       if (length(other_vals) > 0) {
#         interaction$factors[[col]][interaction$factors[[col]] == ".other_val"] <- other_vals[1] # just take the first other_val
#       }
#     }
#   }
#   
#   # find rows in interaction$factors that match data_factors
#   match_rows <- apply(interaction$factors_lower, 1, function(row) {
#     any(apply(data_factors, 1, function(drow) all(row == drow)))
#   })
#   
#   return(interaction)
# }
# 
# list_main_and_inter_eff <- function(trt, formula) {
#   # 1) tokens after '..' from your column names
#   included_vars <- formula  %>% brms::brmsterms() %>% .$allvars %>% stats::terms() %>% attr(.,"term.labels")
#   if("y" %in% included_vars) included_vars <- setdiff(included_vars, "y")
#   
#   main_eff_lookup <- switch(trt,
#                             ".trt" = paste0("^\\.trt[^.]*$"),
#                             "z." = paste0("^z\\.[^.]*$"))
#   main_tx <- included_vars[grepl(main_eff_lookup, included_vars)]
#   if(length(main_tx) == 0) stop("No treatment terms found in model")
#   
#   inter_eff_lookup <- switch(trt,
#                              ".trt" = "^\\.trt.*\\.\\.",
#                              "z." = "^z\\..*\\.\\.")
#   interactions <- included_vars[grepl(inter_eff_lookup, included_vars)]
#   other_vars <- setdiff(included_vars, c(main_tx, interactions))
#   
#   return(list(main_tx = main_tx,
#               interactions = interactions,
#               other_vars = other_vars))
# }
# 
# recover_interaction_variables <- function(trt,fit) {
#   formula <- fit$orig_form
#   data <- fit$orig_data
#   
#   effects <- list_main_and_inter_eff(trt,fit$formula)
#   interactions <- effects$interactions
#   
#   tokens <- unique(unlist(lapply(interactions, function(s) {
#     strsplit(s, "\\.\\.")[[1]]
#   })))
#   
#   if (length(tokens) == 0) return(character(0))
#   
#   terms <- sapply(formula, function(x) stats::terms(x)  %>% 
#                     attr(.,"variables") %>% 
#                     as.character() %>% 
#                     .[-1:-2]
#   ) %>% unlist() %>% unique()
#   
#   inter_terms <- terms[sapply(terms, function(x) any(grepl(x, tokens)))]
#   
#   if (length(inter_terms) == 0) {
#     list(terms = character(0),
#          levels = list())
#   }
#   
#   if( all(inter_terms %in% colnames(data))) {
#     levels <- sapply(inter_terms, function(x) {
#       if(is.factor(data[[x]])) {
#         levels(data[[x]]) %>% tolower() %>% special_clean_names()
#       } else if(is.character(data[[x]])) {
#         unique(data[[x]]) %>% tolower()  %>% special_clean_names()
#       } else if(length(unique(data[[x]])) < 10) {
#         unique(data[[x]])
#       } else if(is.numeric(data[[x]])) {
#         NULL
#       } else {
#         stop(sprintf("Variable %s is of unsupported type: %s", x, class(data[[x]])))
#       }
#     }, simplify = FALSE)
#     
#   } else {
#     levels <- sapply(inter_terms, function(x) {
#       cur_token <- tokens[grepl(x, tokens)]
#       if(length(cur_token) == 0) return(NULL)
#       if(length(cur_token) == 1) return(c(".base_val",".other_val"))
#       str_split(cur_token, x) %>% 
#         lapply(function(y) y[-1]) %>% 
#         unlist() %>% 
#         unique() %>% 
#         c(".base_val",.) %>%
#         return()
#     })
#   }
#   
#   if(!is.null(levels$.trt )) {
#     levels$.trt <- levels$.trt[levels$.trt != "placebo"]
#   }
#   if (!is.null(levels$`z.` )) {
#     levels$`z.`  <- levels$`z.` [levels$`z.`  != ".base_val"]
#   }
#   
#   
#   return(list(terms = inter_terms,
#               levels = levels)
#   )
#   
# }
# 
# get_interaction_contrasts <- function(trt, fit) {
#   form <- fit$prep$orig_form
#   data <- fit$prep$data
#   
#   
#   interaction_terms_and_levels <- recover_interaction_variables(trt,fit)
#   
#   interaction_levels <- interaction_terms_and_levels$levels
#   interaction_terms <- interaction_terms_and_levels$terms
#   if(length(interaction_levels) == 0) {
#     return(trt_levels)
#   }
#   
#   inter_names <- names(interaction_levels)
#   inter_combos <- expand.grid(interaction_levels, stringsAsFactors = FALSE)
#   
#   
#   # skeleton
#   # [x] get all contrasts
#   # [ ] verify with data which combinations exist
#   # [ ] create design matrix to get predictions
#   # [ ] ensure proper labels align with each row (may be automatic given combos?)
#   
#   # following functions don't exist yet
#   actual_combos <- check_combos(inter_combos, data, trt)
#   # in real data
#     # filter for main effects == 1
#     # filter for each combination of interactions and for column name to exist
#   # see b12_b6_b9 at right: depressed == 0 and def unknown but interaction effect was dropped due to colinearity. Could redo main effect I suppose but hard to know how to
#   # automate the dropping to preserve weird combos like this
#   #   depressed def        .trt .trtb1 .trtb12_b6_b9 .trtb12_b9 .trtb12_b9_d .trtb6
#   # 24         0  NA B12, B6, B9      0             1          0            0      0
#   # .trtb9 .trtb9_iron .trtc .trtd .trtd_iron .trtmagnesium .trtmagnesium_zinc
#   # 24      0           0     0     0          0             0                  0
#   # .trtselenium .trtzinc .trtb12_b9..depressed .trtd..depressed
#   # 24            0        0                     0                0
#   # .trtmagnesium..depressed .trtselenium..depressed .trtd..defyes .trtd..defna
#   # 24                        0                       0             0            0
#   # .trtzinc..defna .trtd..depressed..defyes .trtd..depressed..defna
#   # 24               0                        0                       0
#   
#   design <- create_design(actual_combos, trt, data)
#   # given combos, create design matrix
#   
#   return(list(factors = actual_combos,
#               matrix = design))
# }

set_trt <- function(fit) {
  .trt_found <- any(grepl("^.trt",fit$data %>% colnames()))
  .z_found   <- any(grepl("^z.",fit$data %>% colnames()))
  
  if(.trt_found && !.z_found) {
    trt <- ".trt"
  } else if (.z_found && !.trt_found) {
    trt <- "z."
  } else {
    stop("Could not determine treatment variable in data. Both .trt and .z found")
  }
  return(trt)
}

nma_newdata_for_summary <- function(fit) {
  
  trt <- set_trt(fit)
  
  combos <- fit$prep$combos
  
  if(trt == ".trt") {
    main_tx <- stats::model.matrix( ~ 0 + .trt, data = fit$prep$combos) %>% 
      special_clean_names()
    combos <- combos %>%
      select(-.trt) %>% 
      cbind(main_tx)
  }  
  
  if (any(grepl(".dose", colnames(fit$data)))) {
    classes    <- fit$prep$class_ids
    dose_form  <- stats::as.formula(
      paste0("~ 0 + ", paste(classes$varname, collapse = " + "))
    )
    dose_tx    <- stats::model.matrix(dose_form, data = fit$prep$combos) %>% 
      special_clean_names()
    colnames(dose_tx) <- gsub(glue::glue("^{classes$varname}"),"\\.dose",colnames(dose_tx))
    
    combos <- cbind(combos, dose_tx)
  }
  
  combos$.study <- fit$data$.study[1]
  
  # gets other continuous covariates
  if (ncol(fit$data %>% select(-c(y, .study, starts_with(trt), starts_with(".dose")))) > 0) {
    other_covs <- fit$data %>%
      select(-c(y, .study, starts_with(trt), starts_with(".dose"), any_of(colnames(combos))))
    combos <- combos %>%
      dplyr::cross_join(
        other_covs
      )
    
  }
  
  return(combos %>% 
           distinct(across(-c(.se,.cov,.obs)), .keep_all = TRUE))
  
}

nma_newdata_for_summary <- function(fit, keep = NULL) {
  
  trt <- set_trt(fit)
  if(is.null(keep)){
    keep_q <- rlang::quo(get_trt_terms(trt, fit))
  } else if(length(keep) == 1 && keep == ".") {
    keep_q <- rlang::quo(dplyr::everything())
  } else {
    keep_q <- rlang::quo(get_terms(keep, fit))
  }
  
  grid <- make_cell_grid(fit, trt_var = trt, cell_vars = keep,
                 drop_placebo_from_trts = TRUE)
  extra_cov <- fit$data %>%
    select(-c(starts_with(trt), starts_with(".dose"), any_of(colnames(grid)))) %>%
    distinct() %>%
    mutate( dplyr::across(-dplyr::any_of(!!keep_q), ~ {
      if(is.character(.) || is.factor(.) || is.logical(.)) {
        Mode(.)
      } else if (is.numeric(.)) {
        mean(., na.rm = TRUE)
      }
    })) %>%
    distinct(dplyr::across(dplyr::any_of(!!keep_q)),
             .keep_all = TRUE)
  
  dat <- cbind(grid,extra_cov)
  dat$.study <- factor(dat$.study, levels(fit$data$.study))
  return(dat)
  
  combos <- fit$prep$combos
  
  if(trt == ".trt") {
    main_tx <- stats::model.matrix( ~ 0 + .trt, data = fit$prep$combos) %>% 
      special_clean_names()
    combos <- combos %>%
      select(-.trt) %>% 
      cbind(main_tx)
  }  
  
  if (any(grepl(".dose", colnames(fit$data)))) {
    classes    <- fit$prep$class_ids
    dose_form  <- stats::as.formula(
      paste0("~ 0 + ", paste(classes$varname, collapse = " + "))
    )
    dose_tx    <- stats::model.matrix(dose_form, data = fit$prep$combos) %>% 
      special_clean_names()
    colnames(dose_tx) <- gsub(glue::glue("^{classes$varname}"),"\\.dose",colnames(dose_tx))
    
    combos <- cbind(combos, dose_tx)
  }
  
  combos$.study <- factor(fit$data$.study[1], levels(fit$data$.study))
  
  # gets other continuous covariates
  if (ncol(fit$data %>% select(-c(y, .study, starts_with(trt), starts_with(".dose")))) > 0) {
    other_covs <- fit$data %>%
      select(-c(y, .study, starts_with(trt), starts_with(".dose"), any_of(colnames(combos))))
    combos <- combos %>%
      dplyr::cross_join(
        other_covs
      )
    
  }
  
  return(combos %>% 
           distinct(across(-c(.se,.cov,.obs)), .keep_all = TRUE))
  
}
make_cell_grid <- function(fit,
                           trt_var = ".trt",
                           cell_vars = NULL,
                           include_levels = list(),
                           drop_placebo_from_trts = TRUE) {

  dat <- fit$prep$combos
  placebo_level <- get_reference_treatment(fit)
  if (is.null(cell_vars)) cell_vars <- colnames(dat)
  # treatments (keep factor levels if they exist)
  trt <- dat[[trt_var]]
  trt_levels <- if (is.factor(trt)) levels(trt) else sort(unique(trt))
  
  if (drop_placebo_from_trts) trt_levels <- setdiff(trt_levels, placebo_level)
  
  dose_present <- FALSE
  if(any(cell_vars == ".dose") || all(cell_vars == ".")) {
    cell_vars <- setdiff(cell_vars, ".dose")
    dose_present <- TRUE
  }
  
  if(all(cell_vars == ".") ) {
    cell_vars <- colnames(fit$prep$data)
  }
  
  # for each cell var: either user-specified levels, or all factor levels / unique values
  cell_levels <- lapply(cell_vars, function(v) {
    if (!is.null(include_levels[[v]])) return(include_levels[[v]])
    x <- fit$prep$data[[v]]
    if (is.factor(x)) levels(x) else sort(unique(x))
  })
  names(cell_levels) <- cell_vars
  
  # make sure levels exist in fit data
  for (v in cell_vars) {
    if (!all(cell_levels[[v]] %in% fit$prep$data[[v]])) {
      cell_levels[[v]] <- cell_levels[[v]][cell_levels[[v]] %in% fit$prep$data[[v]]]
    }
  }
  
  # expand grid
  if (trt_var %in% cell_vars) {
    cell_vars <- setdiff(cell_vars, trt_var)
    cell_levels[[trt_var]] <- NULL
  }
  grid <- tidyr::expand_grid(
    !!trt_var := trt_levels,
    !!!cell_levels
  )
  
  # observed flag: does this combo exist in data (ignoring treatment or including it)
  # Here: including treatment in "observed"
  key_dat <- dat %>%
    distinct(across(any_of(c(trt_var, cell_vars)))) %>%
    mutate(.observed = TRUE)
  
  grid <- grid %>%
    left_join(key_dat, by = c(trt_var, cell_vars),
              relationship = "one-to-many") %>%
    mutate(.observed = dplyr::if_else(is.na(.observed), FALSE, .observed))
  
  if(trt_var == ".trt") {
    
    if (any(grepl(".dose", colnames(fit$data))) && dose_present) {
      classes    <- fit$prep$class_ids
      dose_form  <- stats::as.formula(
        paste0("~ 0 + ", paste(classes$varname, collapse = " + "))
      )
      dist_comb <- fit$prep$combos %>% 
        select(all_of(c(trt_var, classes$varname))) %>%
        distinct()
      dose_tx    <- stats::model.matrix(dose_form, data = dist_comb) %>% 
        special_clean_names()
      colnames(dose_tx) <- gsub(glue::glue("^{classes$varname}"),"\\.dose",colnames(dose_tx))
      dat2 <- cbind(dist_comb, dose_tx)
      
      grid <- grid %>% dplyr::full_join(
        dat2,
        by = c(trt_var),
        relationship = "many-to-many"
      ) %>% 
        mutate(across(all_of(colnames(dose_tx)), ~ ifelse(is.na(.), 0, .))) %>% 
        select(-any_of(classes$varname))
      
    }
    lvls <- fit$prep$network$treatments %>% levels()
    if(is.null(lvls)) lvls <- c(get_reference_treatment(fit), sort(unique(grid[[trt_var]])))
    grid[[trt_var]] <- factor(grid[[trt_var]], levels = lvls)
    main_tx <- stats::model.matrix( ~ 0 + .trt, data = grid)[,-1,drop = FALSE] %>% 
      special_clean_names()
    
    grid <- grid %>%
      select(-.trt) %>% 
      cbind(main_tx)
  }  
  
  grid
}

Mode <- function(InVec, mult = FALSE) {
  if (!is.factor(InVec)) InVec <- factor(InVec)
  A <- tabulate(InVec)
  if (isTRUE(mult)) {
    levels(InVec)[A == max(A)]
  } 
  else levels(InVec)[which.max(A)]
}

get_terms <- function(var, fit) {
  terms_fit <- stats::terms(stats::model.frame(fit) )
  lbls <- attr(terms_fit, "term.labels")
  fctr <- attr(terms_fit, "factors")
  
  all.vars <- rownames(fctr)
  fctr[unlist(lapply(var, function(v) grep(v,all.vars))), , drop = FALSE] %>% colSums() %>% `>`(0) -> cols_involving_var
  
  # pick only terms that involve .trt somewhere
  trt_terms <- fctr[,cols_involving_var, drop = FALSE] %>%
    rowSums() %>%
    `>`(0) %>%
    which() %>% names()
  
  trt_terms
}

get_trt_terms <- function(trt, fit) {
  unique(c(get_terms(trt, fit), get_terms(".dose", fit)))
}

get_trt_cols <- function(trt,fit) {
  grep(paste0(trt,"|.dose"), colnames(fit$data),value = TRUE)
}

# get treatment labels used in column names for pretty printing
get_trt_labels <- function(trt, fit) {
  trt_cols <- c(trt, fit$prep$class_ids$varname)
  trt_labels <- lapply(trt_cols, function(col) {
    if (col %in% colnames(fit$prep$combos)) {
      levels <- 
        if(is.factor(fit$prep$combos[[col]])) {
          fit$prep$combos %>%
            select(all_of(col)) %>%
            dplyr::pull() %>% 
            levels() %>% sort()
        } else {
          fit$prep$combos %>%
            select(all_of(col)) %>%
            dplyr::pull() %>% unique() %>% sort()
        }
    }   else {
      character(0)
    }
    levels
  })
  stats::setNames(trt_labels, trt_cols) -> trt_labels
  if(!is.null(trt_labels[[trt]])) {
    trt_labels[[trt]] <- trt_labels[[trt]][trt_labels[[trt]] != get_reference_treatment(fit)]
    trt_labels[[trt]] <-  trt_labels[[trt]] %>% stringr::str_to_title()
  }
  trt_labels
}

get_map <- function(lbl,fit) {
  dose_name <- fit$prep$class_ids$varname
  
  names(lbl)[names(lbl) == dose_name] <- ".dose"
  out <- vector("list", length(lbl))
  names(out) <- names(lbl)
  for (i in seq_along(lbl)) {
    vn <- names(lbl)[i]
    if(vn == ".trt") {
      lvls <- fit$prep$network$treatments %>% levels()
    } else{
      lvls <- NULL
    }
    if(is.null(lvls)) lvls <- c(get_reference_treatment(fit), sort(unique(lbl[[vn]])))
    df <- data.frame(x = factor(lbl[[vn]], levels =  lvls %>% stringr::str_to_title())) %>% 
      rename(!!vn := x)
    out[[i]] <- stats::model.matrix(~ . + 0, data = df)[,-1,drop = FALSE] %>%
      colnames() %>%
      special_clean_names()
  }
  
  out <- unlist(out)
  names(out) <- unlist(lbl)
  out
}

summary.vitfit <- function(object, probs = c(0.025, 0.5, 0.975),
                                     drop = NULL, ...) {
  check_vitfit(object)
  
  object %>% summary_brms_nma(...) %>% summary(probs = probs, drop = drop, ...)
}

#' summarize brms nma fit
summary_brms_nma <- function(fit, keep = NULL, 
                             newdata = NULL,
                             placebo = NULL,
                             re_formula = NULL, 
                             resp = NULL,
                             index = NULL,
                             ...) {
  call <- match.call()
  
  check_vitfit(fit)
  
  trt <- set_trt(fit)
  all_trt_cols   <- get_trt_cols(trt, fit) 
  all_trt_labels <- get_trt_labels(trt, fit)
  trt_label_map  <- get_map(all_trt_labels,fit) 
  
  stopifnot(all((trt_label_map %>%  stringi::stri_trans_nfkc()) %in% (all_trt_cols %>% stringi::stri_trans_nfkc())))
  
  if(is.null(keep)){
    keep_q <- rlang::quo(get_trt_terms(trt, fit))
  } else if(length(keep) == 1 && keep == ".") {
    keep_q <- rlang::quo(dplyr::everything())
  } else {
    cn <- c(get_terms(keep, fit), keep %>% setdiff(c(trt,".dose"))) %>% unique()
    keep_q <- rlang::quo(cn)
  }
  
  # get newdata for predictions
  if (is.null(newdata)) {
    newdata <- nma_newdata_for_summary(fit, keep = keep)
  }
  if (!is.null(index)) {
    newdata[[index]] <- 1:nrow(newdata)
  }
  if (is.null(placebo)) {
    placebo <- newdata %>%
      mutate(dplyr::across(dplyr::any_of(all_trt_cols), ~  0))
  }
  if (inherits(fit$formula, "mvbrmsformula") && is.null(resp)) {
    resp <- "y"
  }
  if (!is.null(resp) && resp == "yanx" && any(grepl(".idx_anx", colnames(newdata)))) {
    subset_idx <- newdata$.idx_anx == 1
  } else if (!is.null(resp) && resp == "y" && any(grepl(".idx_dep", colnames(newdata)))) {
    subset_idx <- newdata$.idx_dep == 1
  } else {
    subset_idx <- rep(TRUE, nrow(newdata))
  }
  
  # get predictions
  preds <- brms::posterior_linpred(
    fit, newdata = newdata,
    allow_new_levels = TRUE, re_formula = re_formula, resp = resp,
    ...
  ) -
    brms::posterior_linpred(
      fit, newdata = placebo,
      allow_new_levels = TRUE, re_formula = re_formula, resp = resp,
      ...
    )
  
  if(is.matrix(preds)) {
    preds <- t(preds)
  } else if (is.array(preds) && length(dim(preds)) == 3) {
    preds <- preds[,,"y"] %>% t()
  } else {
    stop("Unexpected format of posterior_linpred output")
  }
  
  colnames(preds) <- paste0("pred_", 1:ncol(preds))
  
  pred.df <- cbind(newdata %>% 
                     select(dplyr::all_of(!!keep_q), 
                            any_of(".observed")) %>% 
                     dplyr::filter(subset_idx), 
                   preds ) %>% as.data.frame()
  
  d_comb <- pred.df %>%
    tidyr::pivot_longer(cols = dplyr::starts_with(trt), 
                        names_to = ".trt", 
                        values_to = ".level") %>%
    filter(.level == 1) %>% 
    # drop helper column
    select(-c(.level)) %>% 
    { if(any(grepl(".dose", colnames(.)))) {
      tidyr::pivot_longer(., cols = dplyr::starts_with(".dose"), 
                        names_to = ".dose", 
                        values_to = ".level") %>% 
        
        filter(.level == 1) %>%  
        # drop helper column
        select(-c(.level)) 
      } else {
        .
      } 
      } %>%
    tidyr::pivot_longer(cols = dplyr::starts_with("pred_"), 
                        names_to = "iter", 
                        values_to = "value") %>% 
    # remove prefix for readability
    # mutate(.trt = stringr::str_remove(.trt, "^\\.trt")) %>%
    # mutate(.trt = stringr::str_remove(.trt, "^z\\.")) %>%
    mutate(.trt = factor(.trt, labels = names(trt_label_map), 
                         level = trt_label_map)) %>% 
    mutate(iter = stringr::str_remove(iter, "^pred_") %>% as.integer())
    
  if(".dose" %in% colnames(d_comb)) {
    d_comb <- d_comb %>%
      mutate(.dose = factor(.dose %>%  stringi::stri_trans_nfkc(), labels = names(trt_label_map), 
                           levels = trt_label_map %>%  stringi::stri_trans_nfkc()))
  }
  
  # if ( !is.null(margins) ) {
  #   if ( margins == "." ) margins <- rlang::quo(dplyr::everything())
  #     
  #   d_comb <- d_comb %>%
  #     group_by(dplyr::across(c(all_of(margins), iter))) %>%
  #     summarize(value = mean(value), .groups = "drop") %>%
  #     ungroup() 
  # }
  d_comb <- data.table::as.data.table(d_comb)
  
  
  if (inherits(d_comb, "data.table")) {
    data.table::setattr(d_comb, "class", c("summary_brms_nma", class(d_comb)))
    data.table::setattr(d_comb, "reference_treatment",  get_reference_treatment(fit))
    data.table::setattr(d_comb, "ndraws", brms::ndraws(fit))
    data.table::setattr(d_comb, "ntreatments", ntreatments(fit))
    data.table::setattr(d_comb, "treatments", fit$prep$network$treatments %>% levels())
  } else if(inherits(d_comb, "tbl_df")) {
    class(d_comb) <- c("summary_brms_nma", class(d_comb))
    attr(d_comb, "reference_treatment") <- get_reference_treatment(fit)
    attr(d_comb, "ndraws") <- brms::ndraws(fit)
    attr(d_comb, "ntreatments") <- ntreatments(fit)
    attr(d_comb, treatments) <- fit$prep$network$treatments %>% levels()
  } else {
    rlang::warn("Could not set attributes on summary object because it is not a data.table or tibble")
  }
  
  return(d_comb)
}

summary.summary_brms_nma <- function(object, probs = c(0.025, 0.5, 0.975),
                                     drop = NULL, 
                                     na.rm = FALSE,
                                     ...) {
  # object %>% group_by(across(-c(iter,value, dplyr::any_of(drop)))) %>%
  #   summarize(Estimate = mean(value), 
  #             `S.E.` = stats::sd(value), 
  #             `l-95% CI`= stats::quantile(value, .025), 
  #             `median`= stats::quantile(value, 0.5),
  #             `u-95% CI`=stats::quantile(value,.975),
  #             .groups = "drop")
  
  # 1. Create formatted names for the quantiles (e.g., "2.5%", "50%", "97.5%")
  q_names <- paste0("q",probs * 100, "%")
  
  # 2. Define the quantile functions using the provided probs
  # Use purrr::map or a loop to create a list of quantile calls
  q_funs <- purrr::map(probs, ~ ~stats::quantile(.x, probs = .y, na.rm = na.rm), .y = probs)
  # Alternatively, a simple list approach:
  qs <- stats::setNames(lapply(probs, function(p) function(x) stats::quantile(x, p, na.rm = na.rm)), q_names)
  
  # object %>%
  #   group_by(across(-c(iter, value, dplyr::any_of(drop)))) %>%
  #   summarize(
  #     Estimate = mean(value, na.rm = na.rm),
  #     `S.E.` = stats::sd(value, na.rm = na.rm),
  #     # 3. Use across() to apply the named list of quantile functions to 'value'
  #     across(value, qs, .names = "{fn}"),
  #     .groups = "drop"
  #   )
  DT <- data.table::as.data.table(object)
  
  by_cols <- setdiff(names(DT), c("iter", "value", drop))
  
  res <- DT[, c(
    list(
      Estimate = mean(value, na.rm = na.rm),
      `S.E.`   = stats::sd(value, na.rm = na.rm)
    ),
    lapply(qs, function(f) f(value))
  ), by = by_cols]
  
  # Give quantile columns the names you expect (the names of qs)
  data.table::setnames(res, old = names(qs), new = names(qs))
  
  res
}
  
print.summary_brms_nma <- function(x, n = 10, drop = NULL, summary = FALSE, ...) {
    if(!inherits(x, "summary_brms_nma")) stop("x must be of class 'summary_brms_nma'")
  n_tx <- length(unique(x$.trt))
  niter<- attr(x, "ndraws")
  
  other_grps <- setdiff(colnames(x), c(".trt","iter","value"))
  
  out <- if(isTRUE(summary)) {
      summary(x, probs = c(0.025, 0.5, 0.975), drop = drop) %>% 
        rename(`q2.5%` = `l-95% CI`,
               `q50%`  = "median",
               `q97.5%`= `u-95% CI`)
  } else {
    dplyr::as_tibble(x)[seq_len(min(n, nrow(x))), , drop = FALSE]
  }
  n_max<- min(n, nrow(out))
  
  info <- glue::glue("These are the ",
                    cli::style_underline("raw"),
                    " data. To see a summary of the results, use ",
                     cli::col_blue("summary = TRUE"), 
                    " in the ",
                     cli::col_blue("`print`"),
                     " call or use the ",
                     cli::col_blue("`summary()`"),
                     " function: ",
                     cli::col_blue("summary(vit$summary_brms_nma(x))"))
  
  cli::cli_h1("BRMS NMA fit")
  
  use_s_tx <- if(n_tx != 1){'s'}else{''}
  use_s_samp <- if(niter != 1){'s'}else{''}
  cli::cli_text(
    "Fit with {cli::col_yellow(n_tx)} treatment{use_s_tx} and {cli::col_yellow(niter)} posterior sample{use_s_samp}."
  )
  
  if (length(other_grps) > 0) {
    cli::cli_alert_info(
      "Results are stratified by: {.emph {paste(other_grps, collapse = ', ')}}"
    )
  }
  
  id <- cli::cli_div(
    class = "indent",
    theme = list(
      ".indent" = list(`margin-left` = 2)
    )
  )
  cli::cli_verbatim(
    utils::capture.output(
      print(out, n = n)
    )
  )
  
  cli::cli_end(id)
  
  if (n_max < nrow(x)) {
    cli::cli_alert_info(
      glue::glue("Showing first {n_max} of {nrow(x)} rows. To see more, increase ",
      cli::col_blue("n"), " in the ",
      cli::col_blue("`print`"),
      "call.")
    )
  }
  
  if(isFALSE(summary)) {
    cli::cli_alert_info("{info}")
  }
  
  invisible(out)
    
}

as.array.vitfit <- function(x, ...) {
  check_vitfit(x)
  
  trt <- set_trt(x)
  
  newdata <- nma_newdata_for_summary(x)
  placebo <- newdata %>%
    mutate(dplyr::across(dplyr::all_of(get_trt_cols(trt, x)), ~  0))
  
  preds <- brms::posterior_linpred(
    x, newdata = newdata,
    allow_new_levels = TRUE
  ) -
    brms::posterior_linpred(
      x, newdata = placebo,
      allow_new_levels = TRUE
    )
  
  return(preds)
}

as.matrix.vitfit <- function(x,...) {
  as.array.vitfit(x,...) %>% 
    {
      if(is.matrix(.)) {
        .
      } else if (is.array(.) && length(dim(.)) == 3) {
        .[,, "y"]
      } else {
        stop("Unexpected format of posterior_linpred output")
      }
    } %>%
    t()
}

all_contrasts <- function(summary_nma, trt_var = ".trt") {
  treatment <- c(get_reference_treatment(summary_nma) %>% 
                   stringr::str_to_title()
                   , levels(summary_nma$.trt))
  contr <- lapply(treatment[-1], function(a) {
    rows <- lapply(treatment, function(b) {
      return(contrasts(summary_nma, .x = a, .y = b) %>% 
               mutate(a = a, b = b))
    }) 
    dplyr::bind_rows(rows)
  }) %>% dplyr::bind_rows()
  return(contr)
}

upper_tri_contrasts <- function(summary_nma, trt_var = ".trt") {
  
  if(!inherits(summary_nma, "data.table")) {
    summary_nma <- data.table::as.data.table(summary_nma)
  }
  return(upper_tri_contrasts_dt(summary_nma, 
                                trt_var = trt_var))
  
  treatment <- c(get_reference_treatment(summary_nma) %>% 
                   stringr::str_to_title()
                 , levels(summary_nma$.trt))
  ntrt      <- length(treatment)
  contr <- lapply(1:(ntrt-1), function(i) {
    rows <- lapply((i+1):ntrt, function(j) {
      a <- treatment[i]
      b <- treatment[j]
      return(contrasts(summary_nma, .x = a, .y = b) %>% 
               mutate(a = a, b = b))
    }) 
    dplyr::bind_rows(rows)
  }) %>% dplyr::bind_rows()
  return(contr)
}

upper_tri_contrasts_dt <- function(summary_dt, trt_var = ".trt") {
  trts <- unique(as.character(summary_dt[[trt_var]]))
  ref_trt <- get_reference_treatment(summary_dt) %>% stringr::str_to_title()
  
  treatment <- unique(c(ref_trt, trts))
  if (length(treatment) < 2L) return(summary_dt)
  
  cmb <- utils::combn(treatment, 2L)
  # swap placebo if in first column to ensure it's always the reference in contrasts
  placebo_idx <- which(cmb == ref_trt, arr.ind = TRUE)[,2L]
  if (length(placebo_idx) > 0L) {
    for (i in placebo_idx) {
      cmb[,i] <- rev(cmb[,i])
    }
  }
  pairs <- data.table::data.table(a = cmb[1L,], b = cmb[2L,])
  
  res <- vector("list", nrow(pairs))
  for (k in seq_len(nrow(pairs))) {
    res[[k]] <- contrasts_dt(summary_dt, pairs$a[k], pairs$b[k], ref_trt)
    res[[k]][, `:=`(a = pairs$a[k], b = pairs$b[k])]
  }
  out <- data.table::rbindlist(res, use.names = TRUE, fill = TRUE)
  
  data.table::setattr(out,"class", c("summary_brms_nma", class(out)))
  data.table::setattr(out, "reference_treatment", ref_trt)
  data.table::setattr(out, "ndraws", attr(summary_dt, "ndraws"))
  data.table::setattr(out, "ntreatments", attr(summary_dt, "ntreatments"))
  data.table::setattr(out, "treatments", attr(summary_dt, "treatments"))
  return(out)
}


league_table <- function(fit, summary = TRUE) {
  
  check_vitfit(fit)
  
  newdata <- nma_newdata_for_summary(fit) %>% 
    distinct(across(-c(.se,.cov,.obs)),.keep_all = TRUE)
  
  sum <- summary_brms_nma(fit, newdata = newdata)
  
  contr <- all_contrasts(sum, ".trt")
  
  if (!summary) {
    return(contr)
  }
  
  league_summary_table <- contr %>%
    summary()
  
  league <- league_summary_table %>%
    mutate(across(where(is.numeric), as.numeric),
           est_ci = sprintf("%.2f (%.2f, %.2f)", Estimate, `q2.5%`, `q97.5%`)) %>%
    select(a, b, est_ci) %>%
    mutate(est_ci = ifelse(a == b, NA, est_ci)) %>%
    tidyr::pivot_wider(names_from = a, values_from = est_ci) %>%
    tibble::column_to_rownames("b")
  
  return(league)
}

theme_vit <- function(...) {
  ggplot2::theme_light(...) +
    ggplot2::theme(
      # panel.border = ggplot2::element_rect(colour = "grey70", fill = NA),
      # panel.grid.major = ggplot2::element_line(colour = "grey95"),
      # panel.grid.minor = ggplot2::element_line(colour = "grey95"),
      # strip.background = ggplot2::element_rect(colour = "grey70", fill = "grey90"),
      strip.text = ggplot2::element_text(colour = "black"),
      strip.text.y.left = ggplot2::element_text(angle = 0, 
                                   hjust = 1, vjust = 0.5, 
                                   size = 8,
                                   lineheight = 0.9,
                                   margin = ggplot2::margin(r = 0, l = 0)),
     strip.text.y.right = ggplot2::element_text(angle = 0, 
                                  hjust = 1, vjust = 0.5, 
                                  size = 8,
                                  lineheight = 0.9,
                                  margin = ggplot2::margin(r = 0, l = 0)),
     legend.position = "bottom",
     strip.background = ggplot2::element_rect(fill = "white", color = NA),
     strip.background.y = ggplot2::element_rect(fill = "white", color = NA),
     strip.background.x = ggplot2::element_rect(fill = "white", color = NA),
     strip.placement = "outside",
     panel.spacing.y = ggplot2::unit(0.02, "lines"),  # small within-treatment gap
     panel.border   = ggplot2::element_rect(fill = NA, color = "grey88", linewidth = 0.25),
     axis.line.x = ggplot2::element_line(color = "grey40", linewidth = 0.4),
     axis.line.y = ggplot2::element_line(color = "grey40", linewidth = 0.4),   
     panel.grid.major.y = ggplot2::element_blank(),
     panel.grid.minor.y = ggplot2::element_blank(),
     panel.grid.major.x = ggplot2::element_line(color = "grey90", linewidth = 0.2),
     strip.switch.pad.grid = ggplot2::unit(0.2, "lines"),
     legend.title = ggplot2::element_text(hjust = 0.5),
     legend.justification = "center",
     legend.box.margin = ggplot2::margin(t = -10, b = 5),
     plot.margin = ggplot2::margin(l = 1, r = 5, t = 5, b = 0.2),
     axis.text.x = ggplot2::element_text(
       margin = ggplot2::margin(t = 0.5)   # smaller = closer to axis
     )
    )
}
  
plot.summary_brms_nma <- function(x, color = NULL, 
                                  grouping = NULL,
                                  base_size = 11, 
                                  position = "identity",
                                  .width = c(0.66, 0.95),
                                  add_left_wall = TRUE,
                                  ...) {
    if(!inherits(x, "summary_brms_nma")) stop("x must be of class 'summary_brms_nma'")
  # .x <- x %>%
  #   mutate(group = 
  #            purrr::pmap_chr(dplyr::across(-c(.trt, iter, value)), 
  #                     ~ paste(names(list(...)), unlist(list(...)), sep = ": ", collapse = ", "))
  #   )
  add_left_wall <- rlang::is_true(add_left_wall)
  .x <- 
  if( ncol(x %>% select(-c(.trt, iter, value, !!color,
                           !!grouping)))> 0) {
      x %>%
      mutate(
        group = purrr::pmap_chr(
          dplyr::across(-c(.trt, iter, value, !!color, !!grouping)),
          ~ {
            vals <- lapply(list(...), as.character)   # preserve labels
            paste(names(vals), unlist(vals), sep = ": ", collapse = ", ")
          }
        )
      )
  } else {
    x
  }
  
  # position <- if (!is.null(color)) ggplot2::position_dodge() else "identity"
  
  y_var <- if (!is.null(.x[["group"]])) "group" else if (!is.null(color)) color else if (!is.null(grouping)) grouping else ".trt"
  
  mapping <- list(
    y = rlang::expr(.data[[!!y_var]]),
    x = rlang::expr(value)
  )
  
  if (!is.null(color)) {
    mapping$colour <- rlang::expr(.data[[!!color]])
  }
  if (!is.null(grouping)) {
    mapping$group <- rlang::expr(.data[[!!grouping]])
  }
  
  p <- .x %>% 
    group_by(.trt) %>%
    mutate(
      # order_key = ifelse(bias == "low", value, NA_real_),
      order_key = mean(value,na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(.trt = .trt %>% stringr::str_wrap(width = width)) %>% 
    mutate(.trt = forcats::fct_reorder(.trt, order_key)) %>% 
    ggplot2::ggplot(ggplot2::aes(!!!mapping)) + 
  ggplot2::geom_vline(xintercept = 0, lty = 2, col = "grey50") +
  ggdist::stat_pointinterval(position = position,
                             .width = .width) +
  ggplot2::xlab("SMD") + 
  theme_vit(base_size = base_size)
  
  if(is.null(.x[["group"]])) {
    p <- p + 
      ggplot2::facet_wrap(~.trt, ncol = 1, 
                          space = "free_y",
                          scales = "free_y", strip.position = "left") + 
      ggplot2::ylab(NULL) + ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      strip.text.y.left = ggplot2::element_text(size = base_size * 0.9)
    )
  } else {
    p <- p + 
      ggplot2::scale_y_discrete(position = "right") +
      ggplot2::facet_wrap(~.trt, ncol = 1, 
                          space = "free_y",
                          scales = "free_y", strip.position = "left") +
      ggplot2::theme(strip.text.y.right = ggplot2::element_text(size = base_size * 0.9))
    
    if (add_left_wall) {
      p <- p + 
        ggplot2::geom_segment(
          data = data.frame(x = -Inf, xend = -Inf, y = -Inf, yend = Inf),
          ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
          inherit.aes = FALSE,
          linewidth = 0.7,
          colour = "black"
        ) 
    }
                     
  }
  
  # print(p)
  return(p)
}

re_plot <- function(x, variable, varnames = NULL) {
  
  re <- brms::ranef(x,probs = c(0.025, 0.17,0.83,.975))[[variable]]
  
  if(!is.null(varnames)) {
    if(!(length(varnames) == length(dimnames(re)[[3]]))) rlang::abort("Variable names not the same length as number of coefficients")
    dimnames(re)[[3]] <- varnames
  }
  
  re_int <- re[,,1] %>% as.data.frame() %>% 
    mutate(!!variable := rownames(re),
           coef = dimnames(re)[[3]][1])
  re_coef<- re[,,2] %>% as.data.frame() %>% 
    mutate(!!variable := rownames(re),
           coef = dimnames(re)[[3]][2])
  
  dplyr::bind_rows(re_int,
            re_coef) %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x    = Estimate,
        y    = !!rlang::sym(variable)
      )
    ) +
    ggdist::geom_pointinterval(ggplot2::aes(xmin = Q2.5,
                                   xmax = Q97.5,),
                               linewidth = 0.6) +
    ggdist::geom_pointinterval(ggplot2::aes(xmin = Q17,
                                   xmax = Q83),
                               size = 4) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",color = "gray70") +
    ggplot2::facet_grid(cols = ggplot2::vars(coef), rows = ggplot2::vars(!!rlang::sym(variable)), scales = "free_y",
                        switch = "y")  +
    ggplot2::ylab(variable) +
    theme_vit() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
    )
}

def_plot <- function(x, filter_by = NULL, position = "identity", ...) {
  filter_quo <- rlang::enquo(filter_by)
  if(!("depressed" %in% colnames(x)) ) {
    outcome <- "anxiety"
  } else {
    outcome <- "depressed"
  }
  
  p <- x %>% 
    # group_by(.trt) %>%
    # mutate(
    #   order_key = ifelse(bias == "low", value, NA_real_),
    #   order_key = mean(value,na.rm = TRUE)) %>%
    # ungroup() %>%
    # mutate(.trt = forcats::fct_reorder(.trt, order_key)) %>% 
    mutate(.trt = .trt %>% stringr::str_wrap(width = width)) %>% 
    {
      expr <- rlang::quo_get_expr(filter_quo)
      
      if (is.null(expr) || rlang::quo_is_missing(filter_quo)) {
        .
      } else {
        dplyr::filter(., !!filter_quo)
      }
    } %>% 
    ggplot2::ggplot(
      ggplot2::aes( 
        y = deficiency, 
        x = value,
        color = deficiency,
        ...
        # !!!(
        #   if (rlang::is_symbol(alpha)) list(alpha = alpha) else NULL
        # )
      )) + 
    ggplot2::geom_vline(xintercept = 0, lty = 2, col = "grey50") +
    ggdist::stat_pointinterval(
      position = position
    ) +
    ggplot2::scale_color_manual(values = def_colors) +
    theme_vit() + 
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   legend.position = "bottom",
                   axis.text.y = ggplot2::element_blank()) +
    ggplot2::xlab("SMD") + 
    # ggplot2::facet_grid(rows = ggplot2::vars(.trt, bias),
    #                     cols =  ggplot2::vars(!!rlang::sym(outcome))) +
    ggplot2::scale_y_discrete(drop = FALSE, expand = ggplot2::expansion(mult = 0)) +
    ggplot2::theme(
      panel.spacing.y = ggplot2::unit(0.05, "lines"),  # small within-treatment gap
      strip.switch.pad.grid = ggplot2::unit(0.2, "lines"),
      strip.text.y.left = ggplot2::element_text(size = 8, margin = ggplot2::margin(0, 4, 0, 4)),
      strip.background.y = ggplot2::element_rect(fill = "grey95", color = NA),
      strip.placement = "outside"
    ) +
    # Add a visible line to separate treatments
    ggplot2::theme(ggh4x.facet.nestline = ggplot2::element_line(color = "grey60", linewidth = 0.5)) 
  
  p <-  if (length(unique(x$bias)) > 1) {
    if(x$.trt %>% dplyr::n_distinct() == 1) {
      p + 
        ggh4x::facet_nested(
          cols = ggplot2::vars(bias, !!rlang::sym(outcome)),
          # scales = "free_y",
          space  = "free_y",
          labeller = ggplot2::labeller(
            bias      = ggplot2::label_both,      # prints "bias = low"
            depressed = ggplot2::label_both,      # prints "depressed = yes"
            .multi_line = TRUE          # keep on one line
          ),
          strip =  ggh4x::strip_nested(           # optional: nicer strip backgrounds
            background_y = ggh4x::elem_list_rect(fill = NA, color = NA),  # no box
            background_x = ggh4x::elem_list_rect(fill = NA, color = NA),   # (optional) columns
            by_layer_y   = FALSE,
            bleed = FALSE
          )
        ) 
    } else {
    p + 
      ggh4x::facet_nested(
        rows = ggplot2::vars(.trt),
        cols = ggplot2::vars(bias, !!rlang::sym(outcome)),
        # scales = "free_y",
        space  = "free_y",
        labeller = ggplot2::labeller(
          bias      = ggplot2::label_both,      # prints "bias = low"
          depressed = ggplot2::label_both,      # prints "depressed = yes"
          .multi_line = TRUE          # keep on one line
        ),
        strip =  ggh4x::strip_nested(           # optional: nicer strip backgrounds
          background_y = ggh4x::elem_list_rect(fill = NA, color = NA),  # no box
          background_x = ggh4x::elem_list_rect(fill = NA, color = NA),   # (optional) columns
          by_layer_y   = FALSE,
          bleed = FALSE
        )
      ) 
    }
  } else {
    if(x$.trt %>% dplyr::n_distinct() == 1) {
      p + 
        ggh4x::facet_nested(
          cols = ggplot2::vars(!!rlang::sym(outcome)),
          # scales = "free_y",
          space  = "free_y",
          labeller = ggplot2::labeller(
            bias      = ggplot2::label_both,      # prints "bias = low"
            depressed = ggplot2::label_both,      # prints "depressed = yes"
            .multi_line = TRUE          # keep on one line
          ),
          strip =  ggh4x::strip_nested(           # optional: nicer strip backgrounds
            background_y = ggh4x::elem_list_rect(fill = NA, color = NA),  # no box
            background_x = ggh4x::elem_list_rect(fill = NA, color = NA),   # (optional) columns
            by_layer_y   = FALSE,
            bleed = FALSE
          )
        ) 
    } else {
     p + 
      ggh4x::facet_nested(
        rows = ggplot2::vars(.trt),
        cols = ggplot2::vars(!!rlang::sym(outcome)),
        # scales = "free_y",
        space  = "free_y",
        labeller = ggplot2::labeller(
          bias      = ggplot2::label_both,      # prints "bias = low"
          depressed = ggplot2::label_both,      # prints "depressed = yes"
          .multi_line = TRUE          # keep on one line
        ),
        strip =  ggh4x::strip_nested(           # optional: nicer strip backgrounds
          background_y = ggh4x::elem_list_rect(fill = NA, color = NA),  # no box
          background_x = ggh4x::elem_list_rect(fill = NA, color = NA),   # (optional) columns
          by_layer_y   = FALSE,
          bleed = FALSE
        )
      ) 
    }
  }
  
  
  p <- p +
    # make row strips span both columns and sit on the left
    ggplot2::theme(
      strip.text.y = ggplot2::element_text(angle = 0, hjust = 0, vjust = 0.5, size = 8,
                                  lineheight = 0.9),
      strip.background = ggplot2::element_rect(fill = "white", color = NA),
      strip.background.y = ggplot2::element_rect(fill = "white", color = NA),
      strip.background.x = ggplot2::element_rect(fill = "white", color = NA),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold"),
      ggh4x.facet.nestline = ggplot2::element_line(color = "grey60"),  # draws the “combining” line
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(l = 1, r = 10, t = 5, b = 0.2)
    ) +
    ggplot2::coord_cartesian(clip = "off")
  p
}

loo_plot <- function(x, loo_est) {
  if(missing(x)) stop("x is required")
  if(missing(loo_est)) stop("loo_est is required")
  check_vitfit(x)
  if(!inherits(loo_est, "summary_brms_nma") & inherits(loo_est, "vitfit")) loo_est <- loo_est %>% summary_brms_nma()
  if(!inherits(x, "summary_brms_nma")) x <- x %>% summary_brms_nma()
  loo_edit <- loo_est %>% 
    group_by(.trt, left_out_study, interventions_affected) %>% 
    summarize(E = stats::median(value)) %>% 
    mutate(.trt = .trt %>% stringr::str_wrap(width = width))
  
  plt <- x %>% 
    select(-.observed) %>% 
    plot() 
  
  plt +
    # ggdist::stat_pointinterval(color = "gray70")  +
    ggplot2::geom_point(data = loo_edit,
                 ggplot2::aes(x = E, y = .trt, color = interventions_affected,
                   alpha = interventions_affected,
                   stroke = interventions_affected,
                   shape = interventions_affected),
               size = 6) +
    ggplot2::scale_color_manual(name = "Dropped study targets intervention", 
                       values = c("FALSE" = "#6a51a3", "TRUE" = "orange"),
                       labels = c("FALSE" = "No", "TRUE" = "Yes")) +
    ggplot2::scale_alpha_manual(values = c("FALSE" = 0.2, "TRUE" = 1)) +
    ggplot2::scale_shape_manual(
      name = "Dropped study targets intervention",
      values = c("FALSE" = 1, "TRUE" = 4),
      labels = c("FALSE" = "No", "TRUE" = "Yes")) +
    # ggplot2::scale_discrete_manual(
    #   name = "Dropped study targets intervention",
    #   aesthetics = "stroke",
    #   values = c("FALSE" = 0.1, "TRUE" = 2)) +
    ggdist::stat_pointinterval(color = "gray50")  +
    ggplot2::facet_wrap(~.trt, ncol = 1, scales = "free_y", strip.position = "left") +
    ggplot2::guides(alpha = "none", stoke = "none") 
}

contrasts <- function(summary_brms, .x, .y) {
  ref_trt <- get_reference_treatment(summary_brms) %>% stringr::str_to_title()
  
  if(!inherits(summary_brms, "data.table")) {
    atts <- attributes(summary_brms)
    summary_brms <- data.table::as.data.table(summary_brms)
    
    for(i in names(atts)) {
      data.table::setattr(summary_brms, i, atts[[i]])
    }
  }
  
  return(contrasts_dt(summary_brms, 
              x_trt = .x, y_trt = .y, ref_trt = ref_trt))
  
  # get contrasts between .x and .y levels of .var
  d_filt_x <- if(.x != ref_trt) {
    summary_brms %>%
      filter(.trt %in% c(.x))
  } else {
    create_baseline_rows(summary_brms)
  }
  
  d_filt_y <- if(.y != ref_trt) {
    summary_brms %>%
      filter(.trt %in% c(.y))
  } else {
    create_baseline_rows(summary_brms)
  }
  
  grp <- setdiff(colnames(d_filt_x), c(".trt","iter","value"))
  names(grp) <- grp
  
  d_contr <- d_filt_x %>%
    dplyr::inner_join(
      d_filt_y %>%
        select(-.trt),
      by = c("iter" = "iter", grp)
    ) %>%
    mutate(value = value.x - value.y) %>%
    mutate(.trt = paste0("[",.trt, "] vs. [", .y,"]")) %>% 
    select(-c(value.x, value.y))
  
  return(d_contr)
}

contrasts_dt <- function(summary_dt, x_trt, y_trt, ref_trt) {
  # summary_dt: data.table with columns .trt, iter, value, plus grouping cols
  # x_trt/y_trt/ref_trt are character (Title-case already)
  
  grp <- setdiff(names(summary_dt), c(".trt","iter","value",".observed"))
  
  if(x_trt == ref_trt || y_trt == ref_trt) {
    # baseline rows: for ref treatment we want value = 0 for each (iter, grp)
    baseline_dt <- unique(summary_dt[, c(grp, "iter"), with = FALSE])
    baseline_dt[, `:=`(.trt = ref_trt, value = 0)]
  }
  
  dx <- if (x_trt != ref_trt) summary_dt[.trt == x_trt, c(grp,"iter","value"), with=FALSE]
  else data.table::copy(baseline_dt[, c(grp,"iter","value"), with=FALSE])
  dy <- if (y_trt != ref_trt) summary_dt[.trt == y_trt, c(grp,"iter","value"), with=FALSE]
  else data.table::copy(baseline_dt[, c(grp,"iter","value"), with=FALSE])
  
  data.table::setnames(dx, "value", "value.x")
  data.table::setnames(dy, "value", "value.y")
  
  keycols <- c("iter", grp)
  data.table:: setkeyv(dx, keycols)
  data.table:: setkeyv(dy, keycols)
  
  out <- dy[dx, nomatch = 0L]  # inner join on keys
  out[, value := value.x - value.y]
  out[, .trt := paste0("[", x_trt, "] vs. [", y_trt, "]")]
  
  out[, c("value.x","value.y") := NULL]
  
  data.table::setattr(out, "class", class(summary_dt))
  data.table::setattr(out, "reference_treatment", ref_trt)
  data.table::setattr(out, "ndraws", attr(summary_dt, "ndraws"))
  data.table::setattr(out, "ntreatments", attr(summary_dt, "ntreatments"))
  data.table::setattr(out, "treatments", attr(summary_dt, "treatments"))
  out
}

marginalize <- function(x, data = NULL, block = NULL, margins, ...) {
  UseMethod("marginalize")
}

marginalize.summary_brms_nma <- function(x, data, block = NULL, margins, ...) {
  
  block <- setdiff(block, ".trt")
  if(missing(margins) || is.null(margins)) margins <- setdiff(colnames(data), c(".trt", block, ".sample_size"))
  # get cell probabilities for all combinations of margins
  # that exist in data by getting mean values
  # ie, if bias = low, depressed = yes, def = no, 
  # what is the proportion of data matching this
  # in theory should be an expectation
  newdata <- data %>% 
    select(.trt,
           dplyr::all_of(block), 
           dplyr::all_of(margins),
           .sample_size) %>% 
    group_by(.trt,
             dplyr::across(dplyr::all_of(block)), 
             dplyr::across(dplyr::all_of(margins))) %>% 
    # summarize(
    #   n = sum(.sample_size, na.rm = TRUE),
    #   .groups = "drop"
    # ) %>%
    summarize(
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    group_by(.trt,
             dplyr::across(dplyr::all_of(block))
             ) %>% 
    mutate(prop = n / sum(n)) %>%
    select(-n) %>% 
    ungroup() %>% 
    mutate(.trt = .trt %>% stringr::str_to_title())
  
  block <- stats::setNames(block, block)
  margins <- stats::setNames(margins, margins)
  
  d_marg <- x %>%
    filter(.observed == TRUE) %>% 
    dplyr::left_join(
      newdata,
      by = c(".trt" = ".trt",
             block, margins),
      relationship = "many-to-one") %>% 
    group_by(.trt, across(all_of(block)), iter) %>% 
    summarize(
      value = sum(value * prop),
      .groups = "drop"
    )
  
  
  
  if(inherits(d_marg, "data.table")) {
    data.table::setattr(d_marg, "class", c("summary_brms_nma", class(d_marg)))
    data.table::setattr(d_marg, "reference_treatment", get_reference_treatment(x))
    data.table::setattr(d_marg, "ndraws", attr(x, "ndraws"))
    data.table::setattr(d_marg, "ntreatments", attr(x, "ntreatments"))
    data.table::setattr(d_marg, "treatments", attr(x, "treatments"))
  } else if (inherits(d_marg, "data.frame")) {
    class(d_marg) <- c("summary_brms_nma",class(d_marg))
    attr(d_marg, "reference_treatment") <- get_reference_treatment(x)
    attr(d_marg, "ndraws") <- attr(x, "ndraws")
    attr(d_marg, "ntreatments") <- attr(x, "ntreatments")
    attr(d_marg, "treatments") <- attr(x, "treatments")
  }
  
  return(d_marg)
  
}

marginalize.vitfit <- function(x, data = NULL, block = NULL, margins, ...) {
  check_vitfit(x)
  
  d_sum <- summary_brms_nma(x, keep = c(block,margins), ...)
  
  d_marg <- d_sum %>% marginalize.summary_brms_nma(
    data = x$prep$network$agd_contrast %>% 
      mutate(.trt = intervention),
    block = block,
    margins = margins
  )
  
  return(d_marg)
}

#' plot brms nma fit
plot_brms_nma <- function(fit, prep, ref_trt = "placebo", ...) {
  re <- brms::posterior_samples(fit, pars = "^r_mu\\[") %>%
    as.data.frame() %>%
    tidyr::pivot_longer(everything(), names_to = "term", values_to = "value") %>%
    mutate(term = term %>% str_replace_all("r_mu\\[","") %>% str_replace_all("\\]","")) %>%
    separate(term, into = c("trt","param"), sep = ",") %>%
    mutate(trt = prep$X_nma[,grep("^.trt", colnames(prep$X_nma))][1,as.numeric(trt)]) %>%
    mutate(param = ifelse(param == "1", "intercept", param)) %>%
    group_by(trt, param) %>%
    summarize(median = median(value),
              lower = quantile(value, 0.025),
              upper = quantile(value, 0.975),
              .groups = "drop") %>%
    mutate(trt = factor(trt, levels = c(ref_trt, setdiff(unique(trt), ref_trt)))) 
  
  p <- ggplot2::ggplot(re, ggplot2::aes(x = trt, y = median)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.2) +
    ggplot2::facet_wrap(~param, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::xlab("Treatment") +
    ggplot2::ylab("Relative effect (SMD)") +
    ggplot2::coord_flip() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red")
  
  return(p)
}

multinma_loo <- function(fit, seed, cores = NULL) {
  call <- match.call()
  if(!is.null(cores) && cores > 1) {
    cat("Registering parallel backend with", cores, "cores\n")
    cl <- parallel::makeCluster(cores)
    
    on.exit({
      parallel::stopCluster(cl)
    }, add = TRUE)
    
    doParallel::registerDoParallel(cl)
    
    doRNG::registerDoRNG(seed)
  } else {
    cat("Registering sequential backend\n")
    foreach::registerDoSEQ()
    
    set.seed(seed)
  }
  
  net <- fit_RE$network
  data<- net$agd_contrast
  
  studies <- data$.study
  unique_studies <- unique(studies)
  nstudy <- length(unique_studies)
  
  cat(sprintf("Running leave-one-study-out analyses on %i studies based on model `%s`...\n", 
              nstudy, prettify_name(call$fit)))
  foreach::foreach(i = 1:nstudy, 
                   .packages = c("dplyr","brms","multinma"),
                   .export = c("unique_studies",
                               "control", "iter", "warmup", "chains",
                               "main","interactive","confounders", "addl_form",
                               "family",
                               "fit","mod_sc","sm","studies","stand",
                               "data")) %dopar% {
    study_i <- unique_studies[i]
    cat("  Fitting without study", i, ":", as.character(study_i), "\n")
    
    loo_dat <- data %>% filter(.study != study_i)
    
    loo_net <- construct_nma_network(
      loo_dat,
      trt_ref = as.character(net$treatments[1])
    )
    
    fit_loo <- multinma::nma(loo_net,
                             trt_effects = fit_RE$trt_effects,
                             prior_trt = fit_RE$prior_trt,
                             prior_het = fit_RE$prior_het,
                             prior_intercept = fit_RE$prior_intercept,
                             prior_covariates = fit_RE$prior_covariates,
                             control = fit_RE$control,
                             refresh = refresh,
                             open_progress = open_progress,
                             show_messages = show_messages)
    
    
    
    # fit_loo <- prep_brms_nma(data = data %>%
    #                            dplyr::filter(study != study_i)
    # ) %>% fit_brms_nma(
    #   cores = 1L,
    #   iter = iter, warmup = warmup,
    #   chains = chains,
    #   control = control,
    #   silent = 2L, refresh = 0)
    
  } -> sum_loo
  
  return(sum_loo)
  
}

extract_metric_brms <- function(fit) {
  
  if(!inherits(fit, "brmsfit")) {
    rlang::abort("fit must be of class 'brmsfit'")
  }
  sa <- fit$fit@stan_args
  if(!fit$backend == "cmdstanr") {
    rlang::warn("fit does not use cmdstanr backend; adaptation not available")
    lapply(1:length(sa), list(step_size = NULL, inv_metric = NULL, metric = NULL))
  }
  
  
  # function to parse single chain's adaptation_info
  parse_one_chain <- function(sa_chain) {
    ai     <- sa_chain$adaptation_info
    metric <- sa_chain$metric   # "diag_e" or "dense_e"
    
    # split adaptation_info into lines
    lines <- unlist(strsplit(ai, "\n"))
    
    ## --- Step size ---
    step_line <- grep("Step size", lines, value = TRUE)
    step_size <- as.numeric(sub(".*= *", "", step_line))
    
    ## --- Which header to look for? ---
    header_pattern <- if (metric == "diag_e") {
      "Diagonal elements of inverse mass matrix"
    } else {
      "Elements of inverse mass matrix"
    }
    
    idx <- grep(header_pattern, lines)
    if (!length(idx)) {
      stop("Could not find inverse mass matrix header in adaptation_info")
    }
    
    ## --- Grab all following '# ...' lines with numbers ---
    num_lines_idx <- integer(0)
    for (j in (idx + 1):length(lines)) {
      line <- lines[j]
      if (!grepl("^#", line)) break          # stop when comments end
      if (grepl("^#\\s*$", line)) break      # stop at bare '#'
      num_lines_idx <- c(num_lines_idx, j)
    }
    
    if (!length(num_lines_idx)) {
      stop("No numeric lines found after inverse mass matrix header")
    }
    
    # strip "#", glue, and split into numbers
    numeric_lines <- gsub("^#\\s*", "", lines[num_lines_idx])
    numeric_lines <- gsub("#\\s*$", "", numeric_lines)
    
    nums_chr <- strsplit(numeric_lines, ",\\s*") %>% lapply(trimws)
    nums_val <- nums_chr %>% lapply(as.numeric) 
    vals     <- do.call("rbind", nums_val)
    
    if (any(is.na(vals))) {
      warning("Some inverse-metric entries could not be parsed as numeric.")
    }
    
    ## --- Shape based on metric type ---
    inv_metric <- if (metric == "diag_e") {
      c(vals)
    } else {
      # full matrix; should be K^2 entries
      if (nrow(vals) != ncol(vals)) {
        warning("Length of dense inverse mass matrix is not a perfect square.")
      }
      vals
    }
    
    list(
      step_size  = step_size,
      inv_metric = inv_metric,
      metric     = metric
    )
  }
  
  # apply to all chains
  lapply(sa, parse_one_chain)
}

average_metric_info <- function(list) {
  n_chains <- length(list)
  if (n_chains == 0) {
    stop("Input list is empty.")
  }
  
  # Check that all chains have the same metric type
  metric_types <- sapply(list, function(x) x$metric)
  if (length(unique(metric_types)) != 1) {
    stop("All chains must have the same metric type.")
  }
  
  metric_type <- metric_types[1]
  
  # Average step sizes
  step_sizes <- sapply(list, function(x) x$step_size)
  avg_step_size <- mean(step_sizes)
  
  # Average inverse metrics
  if (metric_type == "diag_e") {
    inv_metrics <- sapply(list, function(x) x$inv_metric)
    avg_inv_metric <- rowMeans(inv_metrics)
  } else if (metric_type == "dense_e") {
    inv_metrics <- lapply(list, function(x) x$inv_metric)
    sum_matrix <- Reduce("+", inv_metrics)
    avg_inv_metric <- sum_matrix / n_chains
  } else {
    stop("Unknown metric type: ", metric_type)
  }
  
  list(
    step_size  = avg_step_size,
    inv_metric = avg_inv_metric,
    metric     = metric_type
  )
}

get_parameters_block_names <- function(brms_fit) {
  code  <- brms::stancode(brms_fit)
  lines <- strsplit(code, "\n")[[1]]
  
  # locate the parameters block
  start <- grep("^parameters\\s*\\{", lines)
  if (!length(start)) stop("No parameters block found")
  # first closing brace after 'parameters {'
  end <- which(grepl("^\\}", lines) & seq_along(lines) > start)[1]
  if (is.na(end)) stop("Could not find end of parameters block")
  
  block <- lines[(start + 1):(end - 1)]
  
  # clean up: trim, drop empty lines and full-line comments
  block <- trimws(block)
  block <- block[nzchar(block) & !startsWith(block, "//")]
  # strip trailing comments
  block <- sub("//.*$", "", block)
  block <- trimws(block)
  block <- block[nzchar(block)]
  
  # keep only lines with a ';'
  block <- block[grepl(";", block)]
  
  # extract the last identifier before ';'
  # e.g. "  vector[J] z_u;  " -> "z_u"
  names <- sub(".*\\b([A-Za-z_][A-Za-z0-9_]*)\\s*;\\s*$", "\\1", block)
  unique(names)
}

get_parameter_dimlength <- function(brms_fit) {
  if(!inherits(brms_fit, "brmsfit")) {
    rlang::abort("brms_fit must be of class 'brmsfit'")
  }
  parnames <- get_parameters_block_names(brms_fit)
  
  # use par_dims and expand to get dimension lengths
  pdims <- brms::par_dims(brms_fit)[parnames]
  dim_lengths <- list()
  for (p in names(pdims)) {
    dims <- pdims[[p]]
    if (length(dims) == 0) {
      dim_lengths[[p]] <- parnames[p]  # scalar
    } else {
      dim_lengths[[p]] <- rep(parnames, dims)
    }
  }
  
  return(dim_lengths)
}

get_param_layout_true <- function(brms_fit) {
  sf <- brms_fit$fit
  
  params <- get_parameters_block_names(brms_fit)
  
  dims <- sapply(sf@par_dims[params], function(d)
    if (length(d) == 0L) 1L else prod(d)
  )
  
  offsets0 <- base::cumsum(c(0L, utils::head(dims, -1L)))  # 0-based offsets
  names(offsets0) <- params
  
  par_index <- function(name) {
    i0 <- offsets0[[name]]
    k  <- dims[[name]]
    seq.int(i0 + 1L, i0 + k)                  # 1-based indices in flat vector
  }
  
  list(
    params    = params,
    dims      = dims,
    offsets0  = offsets0,
    par_index = par_index
  )
}

drop_study_re_from_inv_metric <- function(studies, study_i, layout, inv_metric) {
  
  study_idx <- which(studies == study_i)
  drop.idx  <- layout$par_index("z_u")[study_idx]
  
  if(is.matrix( inv_metric )) {
    metric      <- solve(inv_metric)
    inv_metric2 <- solve(metric[-drop.idx, -drop.idx, drop = FALSE])
  } else {
    inv_metric2 <- inv_metric[-drop.idx ]
  }
  inv_metric2
}

make_init_from_brms <- function(fit, layout, studies, study_i, chain = NULL) {
  # fit: brmsfit
  # chain: optional chain number (if NULL, sample from all chains)
  
  dropidx <- which(studies == study_i)
  
  sf <- fit$fit
  
  return(posterior::as_draws_array(sf))
  
  # Extract draws in *constrained* space
  # draws <- posterior::as_draws_rvars(fit$fit)  # chains x iterations
  # 
  # # get one resample:
  # sample <- posterior::resample_draws(draws, ndraws = 1)
  # 
  # init_list <- list()
  # for (p in names(draws)) {
  #   if(p == "r_scale__contrastScaleBaselineOutcomeMean") browser()
  #   # if (p == "z_u") browser()
  #   init_list[[p]] <- sample[[p]] %>% posterior::draws_of()
  #   init_list[[p]] <- array(init_list[[p]], dim = dim(init_list[[p]])[-1])  # drop draws dim
  #   if (p == "z_u") {
  #     init_list[[p]] <- init_list[[p]][-dropidx, drop = FALSE]
  #   }
  # }
  # 
  # return(init_list)
  
  draws <- rstan::extract(sf, permuted = FALSE)  # chains x iterations
  
  # Pick chain / iter if not given
  nchains <- brms::nchains(fit)
  niters  <- brms::niterations(fit)
  
  draw_idx <- sample.int(nchains * niters, 1)
  
  # Helper: pull constrained draw for parameter 'p'
  get_constrained <- function(p) {
    idx <- grep(p, dimnames(draws)[[3]])
    x <- draws[iter, chain, idx]
    # scalars become numeric; vectors/matrices are preserved
    x
  }
  
  # We want only parameters that appear in the Stan `parameters` block  
  # (not transformed parameters or internals)
  param_names <- get_parameters_block_names(fit)
  
  # Build init list: each entry must match Stan-constrained form
  init_list <- list()
  for (p in param_names) {
    # if (p == "z_u") browser()
    init_list[[p]] <- unname(get_constrained(p))
    if (p == "z_u") {
      init_list[[p]] <- init_list[[p]][-dropidx, drop = FALSE]
    }
  }
  
  return(init_list)
}

has_mi_terms <- function(brms_model) {
  # Extract the brmsformula object
  bformula <- brms_model$formula
  
  # Function to check a single formula for the string "mi("
  check_formula_for_mi <- function(formula_part) {
    # Convert the formula to a character string and search for "mi("
    # The 'fixed = TRUE' ensures it searches for the literal string, not a regex pattern
    return(grepl("mi\\(", as.character(formula_part), fixed = FALSE))
  }
  
  # Check the primary formula
  if (check_formula_for_mi(bformula$formula)) {
    return(TRUE)
  }
  
  # If it's a multivariate model, check all secondary formulas
  if (!is.null(bformula$forms)) {
    for (form in bformula$forms) {
      if (check_formula_for_mi(form$formula)) {
        return(TRUE)
      }
    }
  }
  
  # Check non-linear formulas if they exist
  if (!is.null(bformula$nl)) {
    for (form in bformula$nl.forms) {
      if (check_formula_for_mi(form$formula)) {
        return(TRUE)
      }
    }
  }
  
  # If no "mi(" was found in any part of the formula
  return(FALSE)
}

get_object_name_rlang <- function(x) {
  # name_sym <- rlang::enexpr(x)
  name_sym <- x
  name_str <- as.character(name_sym)
  if (length(name_str) > 1) {
    name_str <- name_str[1]
  }
  if (name_str == "") {
    name_str <- "<unknown>"
  }
  # format pretty so if longer than 30 chars will truncate
  if (nchar(name_str) > 30) {
    name_str <- paste0(substr(name_str, 1, 27), "...")
  }
  
  return(name_str)
}

prettify_name <- function(x) {
  # name_sym <- rlang::enexpr(x)
  name_str <- as.character(x)
  if (length(name_str) > 1) {
    name_str <- name_str[1]
  }
  if (name_str == "") {
    name_str <- "<unknown>"
  }
  # format pretty so if longer than 30 chars will truncate
  if (nchar(name_str) > 30) {
    name_str <- paste0(substr(name_str, 1, 27), "...")
  }
  
  return(name_str)
}

loo_estimates <- function(fit, seed = NULL, 
                          cores = NULL, verbose = FALSE, resp = "y",
                          keep = NULL, get_estimates = TRUE,
                          reuse_metric = TRUE,
                          backend = "cmdstanr") {
  call <- match.call()
  
  check_vitfit(fit)
  
  if(!is.null(cores) && cores > 1) {
    cat("Registering parallel backend with", cores, "cores\n")
    # cl <- parallel::makeCluster(cores)
    future::plan(future::multisession, workers = cores)
    oopts <- options(future.globals.maxSize = 1e10) # 10 GB
    doFuture::registerDoFuture()
    
    on.exit({
      future::plan(future::sequential)
      doParallel::stopImplicitCluster()
      foreach::registerDoSEQ()
      options(oopts)
    }, add = TRUE)
    
    doRNG::registerDoRNG(seed)
  } else {
    cat("Registering sequential backend\n")
    future::plan(future::sequential)
    doFuture::registerDoFuture()
    foreach::registerDoSEQ()
    set.seed(seed)
  }
  
  data    <- fit$prep$data
  studies <- data$.study
  unique_studies <- unique(studies)
  nstudy  <- nstudies(fit)
  
  main <- fit$prep$orig_form$main
  interactions <- fit$prep$orig_form$interactions
  confounders <- fit$prep$orig_form$confounders
  addl_form <- fit$prep$orig_form$addl_form
  family <- fit$prep$orig_family
  
  control <- fit$stan_args$control
  iter    <- fit$fit@stan_args[[1]]$iter
  warmup  <- fit$fit@stan_args[[1]]$warmup
  chains  <- brms::nchains(fit)
  
  mod_sc <- sm <- NULL
  metric <- layout <- init_fun <- NULL
  
  if ( isTRUE(verbose) ) {
    refresh <- 100L
    silent  <- 0L
    open_progress <- TRUE
    show_messages <- TRUE
  } else {
    refresh <- 0L
    silent  <- 2L
    open_progress = FALSE 
    show_messages = FALSE
  }
  
  stanc <- fit %>% brms::make_stancode()
  stand <- fit %>% brms::make_standata()
  if (backend == "cmdstanr") {
    cat("Compiling base model for reuse with cmdstanr backend...\n")
    sample_iter   <- iter - warmup
    metric_info           <- extract_metric_brms(fit)
    metric                <- average_metric_info(metric_info)
    layout                <- get_param_layout_true(fit)
    adapt_engaged <- TRUE
    step_size     <- metric$step_size
    inv_metric    <- metric$inv_metric
    inits <- NULL
    cm  <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stanc))
    exe <- cm$exe_file()
    if (reuse_metric) {
      adapt_engaged <- FALSE
      orig  <- cm$sample(data = stand, 
                         metric = metric$metric,
                         inv_metric = inv_metric, 
                         adapt_engaged = adapt_engaged, 
                         step_size = step_size, 
                         iter_sampling = iter, iter_warmup = floor(warmup/2L), 
                         parallel_chains = min(cores, chains), chains = chains,
                         show_messages = show_messages, 
                         show_exceptions = FALSE,
                         diagnostics = NULL, refresh = refresh, 
                         max_treedepth = control$max_treedepth)
      inits <- posterior::as_draws_list(orig)
      warmup <- 1L
    }
  } else if (backend == "rstan") {
    # # edit stancode so it doesn't need to recompile each time
    cat("Compiling base model for reuse with rstan backend...\n")
    sm <- rstan::stan_model(model_code = stanc)
    if(!is.null(control$init_buffer)) {
      control$adapt_init_buffer <- control$init_buffer
      control$init_buffer <- NULL
    }
    if (! is.null(control$term_buffer)) {
      control$adapt_term_buffer <- control$term_buffer
      control$term_buffer <- NULL
    }
    if (! is.null(control$window)) {
      control$adapt_window <- control$window
      control$window <- NULL
    }
    # mod_sc <- rewrite_stan_efficient_data(stanc)
    # sm <- rstan::stan_model(model_code = mod_sc)
    # fit <- stats::update(fit, backend = "rstan", chains = 0, recompile = FALSE)
    # fit$prep$data <- data
  } else {
    rlang::abort("Unsupported backend: ", fit$backend, call = call)
  }
  
  cat(sprintf("Running leave-one-study-out analyses on %i studies based on model `%s`...\n", 
              nstudy, prettify_name(call$fit)))
  foreach::foreach(i = 1:nstudy, 
                   .packages = c("dplyr","brms","multinma",
                                 "rstan", "stringr"),
                   .export = c("unique_studies",
                               "control", "iter", "warmup", "chains",
                               "main","interactions","confounders", "addl_form",
                               "family","sample_iter",
                               "fit","mod_sc","sm","studies","stand",
                               "metric","layout","init_fun",
                               "data","cm","exe","adapt_engaged","step_size","inv_metric"
                               )) %dopar% {
                                 
      # pullout study to drop
      study_i       <- unique_studies[i]
      
      # which interventions are targeted in this study?
      flagged_intervention <- data %>% filter(.study == study_i) %>%
        select(.trt) %>% dplyr::pull() %>% unique() %>% stringr::str_to_title()
      
      # get data for the dropped study
      nd            <- get_desired_study_raw_data(fit, studies, study_i)
      
      # calculate original log-likelihood on the full model
      # do it somewhat inefficiently to account for MI models
      ll_orig       <- log_lik(fit,
                               allow_new_levels = TRUE, 
                               sample_new_levels = "uncertainty", 
                               resp = resp)[, which(unique_studies == study_i)]
      
      cat("  Fitting without study", i, ":", as.character(study_i), "\n")
      
      # get random seed
      seed_i <- sample.int(.Machine$integer.max, 1L)
      
      #two optional fitting methods depending on backend
      fit_loo     <- fit
      mod_stand   <- drop_desired_study(fit_loo, 
                                        study_i)
      if(fit$backend == "cmdstanr") { # reuses step size
        # init_fun <- function() make_init_from_brms(fit, layout, studies, study_i)
        m <- cmdstanr::cmdstan_model(cm$stan_file(), exe_file = exe, compile = FALSE)
        outs <- m$sample(data = mod_stand, 
                          init = inits,
                          metric = metric$metric,
                          inv_metric = inv_metric, 
                          adapt_engaged = adapt_engaged, 
                          step_size = step_size, 
                          iter_sampling = sample_iter, iter_warmup = warmup, 
                          parallel_chains = 1, chains = chains,
                          show_messages = show_messages, 
                          show_exceptions = FALSE,
                          diagnostics = NULL, refresh = refresh, 
                          seed = seed_i,
                          max_treedepth = control$max_treedepth)
        fit_loo$fit <- brms::read_csv_as_stanfit(outs$output_files(), model = m) 
      } else if(fit$backend == "rstan") {
      
        fit_loo$fit <- rstan::sampling(sm, data = mod_stand, chains = chains,
                                       iter = iter,
                                       warmup = warmup,
                                       cores = 1L,
                                       seed = seed_i,
                                       control = control, 
                                       open_progress = open_progress, 
                                       show_messages = show_messages,
                                       refresh = refresh,
                                       verbose = FALSE
        )
      } else {
        rlang::abort("Unsupported backend: ", fit$backend)
      }
      fit_loo$model <- stanc
      # fit_loo$data  <- fit$data %>% filter(.study != study_i)
      fit_loo       <- brms::rename_pars(fit_loo)
      if (!inherits(fit_loo,"vitfit")) {
        class(fit_loo) <- c("vitfit",class(fit_loo))
        fit_loo$prep <- fit$prep
      }
      
      # calculate log-likelihood on the dropped obs from the LOO model
      # tricky right now because MI models need to be handled carefully
      # right now just will return NA if any covariates are missing,
      # which is not ideal
      ll_newmod     <- log_lik(fit_loo, 
                           allow_new_levels = TRUE, 
                           sample_new_levels = "uncertainty",
                           resp = resp)
      ll_new        <- ll_newmod[,i]
      ll_est        <- data.frame(
        loo = ll_new,
        ll  = ll_orig
      )
      ll_est$iter   <- 1:nrow(ll_est)
      
      # obs log_lik of new model on observed data
      
      ll_avg <- rep(NA_real_, nstudy) 
      ll_avg[-i] <- ll_newmod[,-i] %>% colMeans()
      
      if (get_estimates) {
        sum_loo_i     <- summary_brms_nma(fit_loo, keep = keep)
      } else {
        sum_loo_i <- dplyr::tibble(
          .trt = NA_character_,
          estimate = NA_real_,
          iter = 1:brms::ndraws(fit)
        )
        class(sum_loo_i) <- c("summary_brms_nma", class(sum_loo_i))
        attr(sum_loo_i, "reference_treatment") <- get_reference_treatment(fit)
        attr(sum_loo_i, "ndraws") <- brms::ndraws(fit)
      }
      sum_loo_i$left_out_study <- study_i
      sum_loo_i     <- sum_loo_i %>% mutate(
        interventions_affected = ifelse(.trt %in% flagged_intervention, TRUE, FALSE)
      )
      sum_loo_i$loo_iter <- i
      sum_loo_i     <- sum_loo_i %>%
        dplyr::left_join(
          ll_est,
          by = c("iter" = "iter")
        )
      
      return(list(sum = sum_loo_i, ll_avg = ll_avg))
    
  } -> sum_loo
  
  cat("Combining results...\n")
  if (inherits (sum_loo[[1]]$sum, "data.table")) {
    out <- data.table::rbindlist(lapply(sum_loo, function(s) s$sum))
  } else {
    out <- lapply(sum_loo,function(s) s$sum) %>% dplyr::bind_rows()
  }
   ll_avg_mat <- do.call(rbind, lapply(sum_loo, function(x) x$ll_avg))
  
  if (inherits (out, "data.table")) {
    data.table::setattr(out, "class", c("vit_loo", class(out)) )
    data.table::setattr(out, "nchains", chains)
    data.table::setattr(out, "reference_treatment", get_reference_treatment(fit))
    data.table::setattr(out, "ndraws", brms::ndraws(fit))
    data.table::setattr(out, "ntreatments", ntreatments(fit))
    data.table::setattr(out, "treatments", fit$prep$network$treatments %>% levels())
    data.table::setattr(out, "ll_delta_avg", ll_avg_mat)
    data.table::setattr(out, "ll_avg", log_lik(fit) %>% colMeans())
  } else {
    class(out) <- c("vit_loo",class(out))
    attr(out, "reference_treatment") <- get_reference_treatment(fit)
    attr(out, "nchains") <- chains
    attr(out, "ndraws")  <- brms::ndraws(fit)
     attr(out, "ntreatments") <- ntreatments(fit)
     attr(out, "treatments") <- fit$prep$network$treatments %>% levels()
    attr(out, "ll_delta_avg")  <- ll_avg_mat
    attr(out, "ll_avg") <- log_lik(fit) %>% colMeans()
  }
  
  
  return(out)
}

loo.vit_loo <- function(x, approx = FALSE, ...) {
  if(! inherits(x, "vit_loo")) rlang::abort("Input must be of class 'vit_loo'")
  
  logMeanExp <- function (x) {
    # from Loo
    logS <- log(length(x))
    matrixStats::logSumExp(x) - logS
  }
  # compute LOO from stored log-lik values
  S <- attr(x, "ndraws")
  
    if(inherits(x, "data.table")) {
      x[, loo := data.table::fifelse(is.na(loo), ll, loo)]
    } else {
      x <- x %>% mutate(loo = ifelse(is.na(loo), ll, loo))
    }
  
  if(inherits(x, "data.table")) {
    DT <- x[, .(loo_iter, iter, loo)]
    
    # take first S rows per loo_iter (like slice_head)
    DTs <- DT[, head(.SD, S), by = loo_iter]
    
    # widen (pivot_wider)
    wide <- data.table::dcast(DTs, iter ~ loo_iter, value.var = "loo")
    
    # drop iter and convert to matrix
    loo_matrix <- as.matrix(wide[, -"iter"])
  }  else {
    loo_matrix <- x %>%
      select(loo_iter, iter, loo) %>%
      group_by(loo_iter) %>%
      dplyr::slice_head(n = S) %>% 
      dplyr::ungroup() %>% 
      tidyr::pivot_wider(names_from = loo_iter, values_from = loo) %>%
      select(-iter) %>%
      as.matrix()
  }
  
  
  ll_matrix <- x %>%
    select(loo_iter, iter, ll) %>%
    group_by(loo_iter) %>%
    dplyr::slice_head(n = S) %>% 
    dplyr::ungroup() %>% 
    tidyr::pivot_wider(names_from = loo_iter, values_from = ll) %>%
    select(-iter) %>%
    as.matrix()
  
  if(isTRUE(approx)) {
    nchains <- attr(x, "nchains")
    chain_id <- rep(seq_len(nchains), each = S / nchains)
    r_eff <- loo::relative_eff(exp(ll_matrix), chain_id = chain_id, ...)
    return(loo::loo(ll_matrix, r_eff = r_eff, ...))
  }
  
  N <- ncol(ll_matrix)
  
  e_i   <- apply(loo_matrix, 2, logMeanExp)
  p_loo <- apply(ll_matrix,  2, logMeanExp) - e_i
  looic <- -2 * e_i
    
  mat_est <- cbind(elpd_loo = e_i, p_loo = p_loo,
                   looic = looic)
  s_i <- sqrt(log(1 + matrixStats::colMeans2(sweep(exp(loo_matrix), 2, exp(e_i))^2)/e_i^2))
  est <- cbind(Estimate = matrixStats::colSums2(mat_est),
               SE = sqrt(N * matrixStats::colVars(mat_est)))
  loo_holder <- purrr::quietly(loo::loo)(loo_matrix, ...)$result
  
  loo_holder$estimates <- est
  loo_holder$pointwise[,c("elpd_loo","p_loo","looic")] <- mat_est
  loo_holder$pointwise[,"mcse_elpd_loo"] <- s_i
  
  loo_holder$diagnostics$pareto_k <- rep(0, N)
  loo_holder$diagnostics$neff     <- rep(S, N)
  
  return(loo_holder)
}

check_vitfit <- function(x) {
  if (!inherits(x,"vitfit")) {
    cli::cli_abort(message = c(i = "{.var x}} must be of class {.code vitfit}.",
    x = "{.var x} is of class {.cls {class(x)}}."))
  }
  return(TRUE)
}

inherits_check <- function(x, class_name) {
  if (!inherits(x, class_name)) {
    cli::cli_abort(message = c(i = "{.var x} must be of class {.code {class_name}}.",
    x = "{.var x} is of class {.cls {class(x)}}."))
  }
  return(TRUE)
}

check_vitsum <- function(x) {
  inherits_check(x, "summary_brms_nma")
}

nstudies.vitfit <- function(x, ...) {
  check_vitfit(x)
  if(!is.null(x$prep$network))  {
    return(nstudies(x$prep$network))
  } else {
    return(nstudies(x$prep$data %>% dplyr::filter(!is.na(y)), study = ".study"))
  }
}

nstudies.nma_data <- function(x, ...) {
  inherits_check(x, "nma_data")
  nlevels(x$studies)
}

nstudies.data.frame <- function(x, study = NULL) {
  x[[study]] %>% unique() %>% length()
}

nstudies <- function(x, ...) {
  UseMethod("nstudies", x)
}

ntreatments.vitfit <- function(x) {
  check_vitfit(x)
  ntr <- if(is.null(x$prep$network)) {
    length(unique(x$prep$data$intervention))  
  } else {
    ntreatments(x$prep$network)
  }
  ntr
}

ntreatments.summary_brms_nma <- function(x) {
  check_vitsum(x)
  ntrt <- attr(x, "ntreatments")
  if(is.null(ntrt)) {
    ntrt <- length(unique(x$.trt)) + 1
  }
  ntrt
}

ntreatments.nma_data <- function(x) {
  inherits_check(x, "nma_data")
  nlevels(x$treatments)
}

ntreatments <- function(x) {
  UseMethod("ntreatments", x)
}

get_reference_treatment.vitfit <- function(x) {
  check_vitfit(x)
  ref <-as.character(x$prep$reference_treatment)
  if(length(ref) == 0) {
    ref <- levels(x$prep$network$treatments)[1]
  }
  if(length(ref) > 1) {
    rlang::abort("Multiple reference treatments found.")
  }
  if(is.na(ref) || is.null(ref)) {
    rlang::abort("No reference treatment found.")
  }
  ref
}

get_reference_treatment.summary_brms_nma <- function(x) {
  check_vitsum(x)
  ref <-attr(x, "reference_treatment")
  if(length(ref) == 0) {
    rlang::abort("No reference treatment found in summary object.")
    ref <- unique(x$.trt)[1]
  }
  if(length(ref) > 1) {
    rlang::abort("Multiple reference treatments found.")
  }
  if(is.na(ref) || is.null(ref)) {
    rlang::abort("No reference treatment found.")
  }
  ref
}

get_reference_treatment.vit_loo <- function(x) {
  
  ref <-attr(x, "reference_treatment")
  if(length(ref) == 0) {
    rlang::abort("No reference treatment found in summary object.")
    ref <- unique(x$.trt)[1]
  }
  if(length(ref) > 1) {
    rlang::abort("Multiple reference treatments found.")
  }
  if(is.na(ref) || is.null(ref)) {
    rlang::abort("No reference treatment found.")
  }
  ref
}


get_reference_treatment.nma_data <- function(x) {
  inherits_check(x, "nma_data")
  ref <- as.character(levels(x$treatments[1]))
  if(length(ref) > 1) {
    rlang::abort("Multiple reference treatments found.")
  }
  if(is.na(ref) || is.null(ref)) {
    rlang::abort("No reference treatment found.")
  }
  ref
}

get_reference_treatment <- function(x) {
  UseMethod("get_reference_treatment", x)
}

traffic_light <- function(x, size = 2) {
  highlight_overall <- x %>%
    dplyr::distinct(study, domain) %>%
    filter(domain == "Overall") %>%
    mutate(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
  
  x %>% 
    ggplot2::ggplot(ggplot2::aes(y = 0, x = 0, color = judgment
               # , fill = judgment
               # , shape = judgment
    )) +
    ggplot2::geom_rect(
      data = highlight_overall,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "grey93", color = NA#"grey60"
    ) +
    ggplot2::geom_point(size = size) +
    ggplot2::facet_grid(study ~ 
                          factor(domain), switch = "y", space = "free") +
    ggplot2::geom_hline(
      data = x %>%
        dplyr::distinct(study) %>%
        mutate(y = -0.5),
      ggplot2::aes(yintercept = y),
      linewidth = 0.25,
      color = "grey60"
    ) +
    ggplot2::scale_color_manual(values = #c("low" = "#2ca02c", "some concerns" = "#ff7f0e", "high" = "#d62728")
                         bias_study_colors,
                       name = "Bias judgment"
    ) +
    ggplot2::scale_y_continuous(limits = c(-0.5, 0.5), expand = c(0, 0)) +
    theme_vit() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(), #element_line(linewidth = 0.2, colour = "grey75"),
      panel.grid.minor.y = ggplot2::element_blank(), # remove y axis grid lines
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.spacing.x = ggplot2::unit(0.08, "lines"),
      # strip.text.y.left = element_text(angle = 0),
      legend.position = "bottom",
    )
}



create_baseline_rows <- function(object) {
  check_vitsum(object)
  ref_trt <- get_reference_treatment(object) %>% stringr::str_to_title()
  grp <- setdiff(names(object), c(".trt","iter","value"))
  
  if(inherits(object, "data.table")) {
    baseline_rows <- unique(object[, c(grp, "iter"), with = FALSE])
    baseline_rows[, `:=`(.trt = ref_trt, value = 0)]
  } else {
    baseline_rows <-
      object %>%
      distinct(
        iter,
        across(-c(.trt, value))
      ) %>%
      mutate(
        .trt  = ref_trt,
        value = 0
      )
  }
  
  
  return(baseline_rows)
}

rank_vit <- function(x, newdata = NULL,  
         lower_better = TRUE,
         probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
         sucra = FALSE,
         summary = TRUE, 
         filter = NULL, ...) {
  # Checks
  if (!rlang::is_bool(lower_better)) rlang::abort("`lower_better` should be TRUE or FALSE.")
  
  if (!rlang::is_bool(summary)) rlang::abort("`summary` should be TRUE or FALSE.")
  
  if (!rlang::is_bool(sucra)) rlang::abort("`sucra` should be TRUE or FALSE.")
  
  if (!rlang::is_double(probs, finite = TRUE) || any(probs < 0) || any(probs > 1)) rlang::abort("`probs` must be a numeric vector of probabilities between 0 and 1.")
  
  if(!inherits(x, "summary_brms_nma")) {
    check_vitfit(x)
    
    rel_eff <- summary_brms_nma(fit = x, newdata = newdata, ...) 
    
  } else {
    check_vitsum(x)
    rel_eff <- x
 } 
  
  
  # set rank fun
  rank_fun <-if (rlang::is_true(lower_better) ) {
     dplyr::min_rank
  } else {
    function(x) dplyr::min_rank(-x)
  }
  
  # Get reference treatment, number of treatments
  trt_ref <- get_reference_treatment(x) %>% stringr::str_to_title()
  ntrt    <- ntreatments(x)
               
  # All other checks handled by relative_effects()
  trt_lvs <- c(trt_ref, levels(rel_eff$.trt))
  
  if(ntrt != length(trt_lvs)) {
    # say number found and number expected
    rlang::abort(message = c(i = "Number of treatments found in summary object does not match expected number of treatments.",
                     i = "Expected: {ntrt}, Found: {length(trt_lvs)}"))
  }
  
  if(inherits(rel_eff, "data.table")) {
    rel_eff <- rel_eff[.observed == TRUE][, .observed := NULL]
    
    if (!is.null(filter)) {
      rel_eff <- rel_eff[ rlang::eval_tidy(filter, data = rel_eff) ]
    }
    
    rel_eff <- rel_eff %>% 
      rbind(create_baseline_rows(rel_eff)) %>% 
      .[, `:=`(.trt = factor(.trt, levels = trt_lvs))] 
  } else {
    rel_eff <- rel_eff %>% 
      dplyr::filter(.observed == TRUE) %>% 
      dplyr::select(-.observed) %>% 
      dplyr::filter(!!filter)
    
    rel_eff <- rel_eff %>% 
      dplyr::bind_rows(create_baseline_rows(rel_eff)) %>% 
      mutate(.trt = factor(.trt, levels = trt_lvs)) 
  }
  
  
  # get ranks of treatment effects at each iteration, by group
  ranks <- rel_eff %>%
    group_by(iter, across(-c(.trt, value))) %>%
    mutate(
      rank = rank_fun(value)    
    ) %>%
    ungroup()
  
  if (summary) {
    rk_summary <-  ranks %>%
      group_by(.trt, across(-c(iter, value, rank))) %>%
      dplyr::summarize(
        mean = mean(rank),
        sd   = stats::sd(rank),
        q = list(stats::quantile(rank, probs = probs, type = 8, na.rm = TRUE)),
        .groups = "drop"
      ) %>% 
      tidyr::unnest_wider(q, names_sep = "")
    
    if (sucra) {
      # Calculate SUCRA using scaled mean rank relation of Rucker and Schwarzer (2015)
      sucras <- unname((ntrt - rk_summary$mean) / (ntrt - 1))
      rk_summary <- tibble::add_column(rk_summary, sucra = sucras, .after = "sd")
    }
    
    out <- list(summary = rk_summary, sims = ranks)
  } else {
    out <- list(sims = ranks)
  }
  
  if (summary) {
    if(inherits(out, "data.table")) {
      data.table::setattr(out, "class", c("vit_ranks","summary_brms_nma"))
      data.table::setattr(out, "xlab", "Treatment")
      data.table::setattr(out, "ylab", "Posterior Rank")
    } else {
      class(out) <- c("vit_ranks", "summary_brms_nma")
      attr(out, "xlab") <- "Treatment"
      attr(out, "ylab") <- "Posterior Rank"
    }
  }
  return(out)
}

rank_probs_vit <- function(x, newdata = NULL, lower_better = TRUE,
                          cumulative = FALSE, sucra = FALSE, 
                          filter = NULL, ...) {
  # Checks
  if (!rlang::is_bool(cumulative)) cli::cli_abort(c(i = "{.var cumulative} should be TRUE or FALSE.", x = "{.var cumulative} is of class {.cls {class(cumulative)}}."))
  
  if (!rlang::is_bool(sucra)) cli::cli_abort(c(i = "{.var sucra} should be TRUE or FALSE.", x = "{.var sucra} is of class {.cls {class(sucra)}}."))
  
  if (!rlang::is_bool(lower_better)) cli::cli_abort(c(i = "{.var lower_better} should be TRUE or FALSE.", x = "{.var lower_better} is of class {.cls {class(lower_better)}}."))
  
  if (!is.null(filter)) {
    if(!rlang::is_quosure(filter)) {
      cli::cli_abort(c(i = "{.var filter} must be a quosure created with {.code rlang::quo()}.",
                       x = "{.var filter} is of class {.cls {class(filter)}}."))
    }
  }
  
  if (!inherits(x, "summary_brms_nma") && !inherits(x, "vitfit")) {
    cli::cli_abort(c(i = "{.var x} must be of class {.code summary_brms_nma} or {.code vitfit}.",
                     x = "{.var x} is of class {.cls {class(x)}}."))
  }
  
  check_total_fun <- function(x, call) {
    if( any(abs(x$total - 1) > 1e-6) ) {
      problem_trt <- x %>%
        filter( abs(total - 1) > 1e-6 ) %>%
        dplyr::pull(.trt)
      cli::cli_warn(c("!" = "Rank probabilities do not sum to 1 for all treatments. Check treatment(s): {paste0(problem_trt, collapse = ', ')}",
                       i = "This occurred in {call}")
                      )
    }
  }
  
  sucra <- rlang::is_true(sucra)
  cumulative <- rlang::is_true(cumulative)
  lower_better <- rlang::is_true(lower_better)
  
  # All other checks handled by posterior_ranks()
  rk <- rank_vit(x = x, newdata = newdata,  
                        lower_better = lower_better, summary = FALSE, filter = filter, ...)
  
  ntrt <- if(is.null(filter)) {
    nlevels(x$prep$network$treatments)
  } else {
    rk$sims %>%
      dplyr::distinct(.trt) %>%
      nrow()  
  }
  nsim <- brms::niterations(x) * brms::nchains(x)
  rank_probs <-
    rk$sims %>%
    { if (is.null(filter)) . else dplyr::filter(., !!filter) } %>% 
    mutate(rank = factor(rank, levels = 1:ntrt)) %>%
    droplevels() %>% 
    dplyr::count(
      .trt,
      rank,
      across(-c(iter, value)),
      .drop = FALSE
    ) %>%
    group_by(
      .trt,
      across(-c(rank, n))
    ) %>%
    mutate(
      prob = n / nsim,
      prob = tidyr::replace_na(prob, 0)
    ) 
  
  check_total_fun(rank_probs %>%
    summarize(total = sum(prob)), call = "rank_probs_vit() when calculating rank probabilities")
  
  if (cumulative) {
    rank_probs <- rank_probs %>%
      arrange(.trt, rank) %>%
      group_by(.trt, across(-c(rank, n, prob))) %>%
      mutate(
        prob = cumsum(prob)
      )
    
    rank_probs %>%
      filter(rank == ntrt) %>%
      summarize(total = prob) %>% 
    check_total_fun(call = "rank_probs_vit() when calculating cumulative probabilities")
  }
  
  if (sucra) {
    sucras <- rank_probs %>%
      { if (is.null(filter)) . else dplyr::filter(., !!filter) } %>% 
      ungroup() %>% 
      group_by(
        .trt,
        across(-c(rank, n, prob))
      ) %>%
      summarize(
        sucra = if (cumulative) {
          mean(prob[-ntrt])
        } else {
          mean((cumsum(prob)[-ntrt]))
        },
        .groups = "drop"
      ) %>% ungroup()
  }
  
  out <- rank_probs %>%
    { if (is.null(filter)) . else dplyr::filter(., !!filter) } %>% 
    ungroup() %>% 
    select(-n) %>% 
    tidyr::pivot_wider(
      names_from = rank,
      names_glue = "p_rank[{rank}]",
      values_from = prob
    )  
  
  if (sucra) {
    cn <- colnames(sucras)  %>% setdiff("sucra") %>% stats::setNames(.,.)
    out <- out %>%
      dplyr::left_join(sucras, by = cn) %>% 
      select(dplyr::all_of(cn), sucra, `p_rank[1]`,
             dplyr::everything())
  }
  
  if (inherits(out,"data.table")) {
    data.table::setattr(out, "class", c("vit_rank_probs", class(out)) )
    data.table::setattr(out, "xlab", "Rank")
    data.table::setattr(out, "ylab", if (cumulative) "Cumulative Rank Probability" else "Rank Probability")
  } else {
    class(out) <- c("vit_rank_probs", class(out))
    attr(out, "xlab") <- "Rank"
    attr(out, "ylab") <- if (cumulative) "Cumulative Rank Probability" else "Rank Probability"
  }
  
  return(out)
}

plot.vit_rank_probs <- function(x, ...) {
  xlab <- attr(x, "xlab")
  if (is.null(xlab)) xlab <- ""
  
  ylab <- attr(x, "ylab")
  if (is.null(ylab)) ylab <- ""
  
  ntrt <- nrow(x)
  
  dat <- tidyr::pivot_longer(
    x,
    cols = dplyr::starts_with("p_rank"),
    names_to = "rank",
    names_pattern = "^p_rank\\[(\\d+)\\]$",
    values_to = "probability",
    names_transform = list(rank = as.integer)) %>% 
    dplyr::rename(Treatment = .trt)
  
  extra_cols <- setdiff(colnames(dat), c("Treatment", "rank", "sucra","probability"))
  
  p <- dat %>% 
    ggplot2::ggplot(ggplot2::aes(x = rank, y = probability)) +
    ggplot2::geom_line(...) +
    ggplot2::ylab(ylab) +
    ggplot2::scale_x_continuous(xlab, 
                                breaks = 1:ntrt, minor_breaks = NULL) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    theme_vit() 
  
  if (length(extra_cols) > 0) {
    p <- p + ggh4x::facet_wrap2(
      ggplot2::vars(!!!rlang::syms(extra_cols)) ~ Treatment, axes = "all"
    )
  } else {
    p <- p + ggh4x::facet_wrap2(~Treatment, axes = "all", remove_labels = "all")
  }
  
  return(p)
}

# approximation for vitfit objects
contribution_matrix_approx_vitfit <- function(x, 
                                              approx = TRUE,
                                              ...) {
  check_vitfit(x)
  
  # sm <- summary_brms_nma(x, ...)
  # 
  # check_vitsum(sm)
  # trt_var  <- set_trt(x)
  # trt_cols <- get_trt_cols(trt_var, x)
  # trt_vars <- grep(paste(trt_cols, collapse = "|"), brms::variables(x), value = TRUE)
  # 
  # studies <- x$data$.study %>% as.character() %>% unique() %>% sort()
  # V <- x$prep$Sigma
  # 
  # 
  # fe <- brms::fixef(x, summary = FALSE)
  # sdat <- brms::standata(x)
  # # Sigma <- x %>% brms::as_draws_matrix(variable = trt_vars)  %>% stats::cov()
  # Sigma <- brms::fixef(x, summary = FALSE) %>% stats::cov()
  # X    <- cbind(sdat$X_main,sdat$X_interactions,sdat$X_confounders)
  # H <- X %*% Sigma %*% crossprod(X, solve(V))
  # 
  # obs_contrasts <- x$prep$data %>% distinct(intervention, control) %>% mutate(control = as.character(control) %>%
  #            gsub(",z\\.",", ",.) %>% gsub("z\\.","",.) %>%
  #              factor(.,levels = levels(x$prep$network$treatments))) %>%
  #   mutate(intervention = factor(intervention, levels = levels(x$prep$network$treatments)))
  # ref  <- get_reference_treatment(x)
  # interv_mat     <- stats::model.matrix(~intervention, data = obs_contrasts)
  # colnames(interv_mat) <- gsub("intervention", ".trt", colnames(interv_mat)) %>% special_clean_names()
  # contr_mat      <- stats::model.matrix(~control, data = obs_contrasts)
  # contrast_matrix <- (interv_mat - contr_mat)[,-1]
  # 
  # obs_cont_data <- 
  # 
  # prep_brms_nma()
  # 
  # stopifnot(colnames(contrast_matrix) %in% colnames(H))
  
  
  if(any(x$prep$data$.study %>% table() > 1)) {
    # hacky introduce missing contrast
    dat <- x$prep$network$agd_contrast  %>% 
      group_by(study) %>% distinct(intervention,control, .keep_all = TRUE) %>% 
      mutate(control = as.character(control) %>%
               gsub(",z\\.",", ",.) %>% gsub("z\\.","",.) %>%
               factor(.,levels = levels(x$prep$network$treatments))) %>%
      mutate(y = ifelse(is.na(y),0,y))
    nm_dat <- with(dat,
      meta::pairwise(treat = intervention,
                     n = final.N,
                     mean = y,
                     sd = .se,
                     studlab = .study))
  } else {
    dat <- x$prep$data %>% 
      mutate(control = as.character(control) %>%
               gsub(",z\\.",", ",.) %>% gsub("z\\.","",.) %>%
               factor(.,levels = levels(x$prep$network$treatments)))
    
    
    nm_dat <- within(dat, {
      studlab <- .study      # or your study id
      treat1  <- intervention     # or whatever your baseline arm is called
      treat2  <- control     # the other arm
      TE      <- y           # observed contrast estimate
      seTE    <- .se         # sampling SE used in vreal()
    })
  }
  
  # You need these columns (rename as needed):
  # studlab, treat1, treat2, TE, seTE
  
  
  # Keep only what netmeta needs
  nm_dat <- nm_dat[, c("studlab", "treat1", "treat2", "TE", "seTE")]
  
  nm <- with(nm_dat,
             netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    sm = "SMD",                  # change if needed
    reference.group = "placebo"  # optional
  ))
  cmat <- netmeta::netcontrib(nm, common = FALSE, random = TRUE, study = TRUE)   
  cmat
}


contribution_matrix_vit_loo <- function(x, fit, ...) {
  
  # check types
  if(!inherits(fit, "summary_brms_nma")) {
    check_vitfit(fit)
    
    # get summary object when type is vitfit
    sm <- summary_brms_nma(fit, ...) %>% 
      .[, .observed := NULL] %>% 
      upper_tri_contrasts(trt_var = set_trt(fit)) %>% 
      .[, placebo_in_pair := a == "Placebo" | b == "Placebo"] %>% 
    .[, `:=`(
        trt_lo = dplyr::case_when(
          placebo_in_pair & a == "Placebo" ~ b,
          placebo_in_pair & b == "Placebo" ~ a,
          TRUE ~ pmin(a, b)
        ),
        
        trt_hi = dplyr::case_when(
          placebo_in_pair ~ "Placebo",
          TRUE ~ pmax(a, b) )
        )] %>% 
    .[, `:=`(
        flipped = a != trt_lo,
        
        trt_key = paste0("[", a, "] vs. [", b, "]")
    
        )] %>% 
      .[, `:=`(
        value = dplyr::if_else(flipped, -value, value),
        a = trt_lo,
        b = trt_hi,
        .trt = trt_key
      )] %>% 
   .[, c("trt_key", "trt_lo", "trt_hi", "flipped", "placebo_in_pair") := NULL]
  } else {
    sm <- fit
  }
  if(! inherits(x, "vit_loo")) rlang::abort("`x` must be of class `vit_loo`.")
  
  process_one <- function(x, sm, col_nm){
    final_df <- x %>% 
      select(-c("interventions_affected",
                "loo_iter",
                "loo", "ll")) %>% 
      upper_tri_contrasts(get_trt_var(fit)) %>% 
      select(-.observed) %>% 
      dplyr::mutate(
        placebo_in_pair = a == "Placebo" | b == "Placebo",
        trt_lo = dplyr::case_when(
          placebo_in_pair & a == "Placebo" ~ b,
          placebo_in_pair & b == "Placebo" ~ a,
          TRUE ~ pmin(a, b)
        ),
        
        trt_hi = dplyr::case_when(
          placebo_in_pair ~ "Placebo",
          TRUE ~ pmax(a, b)
        ),
        
        flipped = a != trt_lo,
        
        trt_key = paste0("[", a, "] vs. [", b, "]"),
        value = dplyr::if_else(flipped, -value, value),
        a = trt_lo,
        b = trt_hi
      ) %>% 
      select(-c(trt_key, trt_lo, trt_hi, 
                flipped, placebo_in_pair)) %>% 
      dplyr::left_join(
        sm,
        by = col_nm %>% stats::setNames(col_nm)
      ) %>% 
      group_by(.trt, a, b, 
               left_out_study,
               dplyr::across(dplyr::any_of(col_nm %>% setdiff("iter"))))
    
    influence_measures <- final_df %>% 
      summarize(contribution = approxOT::wasserstein(X = value.x, Y = value.y, method = "univariate")^2)
    
    # summarize(
    #   distance = approxOT::wasserstein(X = value.x, Y = value.y, method = "sinkhorn")^2 - 0.5 * approxOT::wasserstein(X = value.x, Y = value.x, method = "sinkhorn")^2 - 0.5 * approxOT::wasserstein(X = value.y, Y = value.y, method = "sinkhorn")^2,
    #   .groups = "drop"
    # ) -> influence_measures
    
    influence_measures <- influence_measures %>% 
      ungroup() %>% 
      dplyr::rename(trt1 = a, trt2 = b) %>% 
      dplyr::rename(study = left_out_study) %>% 
      dplyr::mutate(comparison = paste(trt1, trt2, sep = ":")) %>%
      select(comparison, study, contribution, trt1, trt2)
    
    return(influence_measures)
  }

  process_one_dt <- function(x_dt, sm_dt, col_nm, trt_var) {
    # drop big cols early
    drop_cols <- c("interventions_affected","loo_iter","loo","ll")
    keep <- setdiff(names(x_dt), drop_cols)
    x_dt <- x_dt[, ..keep]
    
    # compute contrasts for this slice
    loo_contr <- upper_tri_contrasts_dt(x_dt, trt_var)
    
    # apply your placebo ordering / flipping logic in DT
    loo_contr[, placebo_in_pair := (a == "Placebo" | b == "Placebo")]
    loo_contr[, trt_lo := data.table::fifelse(placebo_in_pair & a == "Placebo", b,
                                              data.table::fifelse(placebo_in_pair & b == "Placebo", a,
                                          pmin(a, b)))]
    loo_contr[, trt_hi := data.table::fifelse(placebo_in_pair, "Placebo", pmax(a, b))]
    loo_contr[, flipped := (a != trt_lo)]
    loo_contr[flipped == TRUE, value := -value]
    loo_contr[, `:=`(a = trt_lo, b = trt_hi)]
    loo_contr[, c("placebo_in_pair","trt_lo","trt_hi","flipped") := NULL]
    
    # join to sm (this creates value.x/value.y like your dplyr left_join)
    # IMPORTANT: col_nm must include iter + any grouping cols used in sm
    bycols <- stats::setNames(col_nm, col_nm)
    final_dt <- data.table::merge.data.table(
      loo_contr,
      sm_dt,
      by = bycols,
      all.x = TRUE,
      suffixes = c(".x",".y"),
      allow.cartesian = TRUE
    )
    
    # group cols = your original group_by(.trt, a, b, left_out_study, across(any_of(col_nm %>% setdiff("iter"))))
    grp_cols <- unique(c(".trt","a","b","left_out_study", setdiff(col_nm, "iter")))
    
    # wasserstein per group
    # (approxOT::wasserstein is not vectorized; we call it once per group)
    out <- final_dt[, .(
      contribution = approxOT::wasserstein(X = value.x, Y = value.y, method = "univariate")^2
    ), by = grp_cols]
    
    out[, `:=`(
      trt1 = a,
      trt2 = b,
      study = left_out_study,
      comparison = paste(a, b, sep = ":")
    )]
    
    out[, .(comparison, study, contribution, trt1, trt2)]
  }
  
  # get additional column-names
   if (!all(c(".trt", "iter") %in% colnames(x))) {
     rlang::abort("`x` must contain columns `.trt` and `iter`.")
   }
  col_nm <- c(colnames(x),"a","b") %>% setdiff(c("value", ".observed",
                                      "left_out_study",
                                      "interventions_affected",
                                      "loo_iter",
                                      "loo", "ll"))

  # browser()
  # get loo estimates and get wasserstein distances per
  # cell. will loop through the left out studies
  # to save memory overhead
  ids   <- unique(x$loo_iter)
  n_id  <- length(ids)
  loo_i <- vector("list", n_id)
  ref_trt_orig <- tryCatch(get_reference_treatment(x),
                      error = function(e) {
                        rlang::warn("No reference treatment found in `x`. Defaulting to reference treatment in `fit`.")
                        get_reference_treatment(fit)
                      })
  
  pb <- cli::cli_progress_bar(
    "Contribution calculation",
    total = n_id,
    format = "{cli::pb_bar} {cli::pb_percent} | {study} | {cli::pb_elapsed} elapsed | {cli::pb_eta} remaining"
  )
  
  filt_dat <- study <- id <- NULL
  cli::cli_alert("Calculating contribution matrix from LOO estimates...")
  
  if(!inherits(x, "data.table")) {
    x_dt <- data.table::as.data.table(x) %>% 
      .[,.observed := NULL]
    data.table::setattr(x_dt, "class", c("summary_brms_nma", 
                                         class(x)[1],class(x_dt)))
  } else {
    x_dt <- data.table::copy(x) %>% 
      .[,.observed := NULL]
    data.table::setattr(x_dt, "class", c("summary_brms_nma",class(x_dt)))
  }
  
  data.table::setattr(x_dt, "reference_treatment", ref_trt_orig)
  data.table::setattr(x_dt, "ndraws", attr(x,"ndraws"))
  data.table::setattr(x_dt, "nchains", attr(x,"nchains"))
  data.table::setattr(x_dt, "ntreatments", attr(x,"ntreatments"))
  data.table::setattr(x_dt, "treatments", attr(x, "treatments"))
  data.table::setattr(x_dt, "ll_delta_avg", attr(x,"ll_delta_avg"))
  data.table::setattr(x_dt, "ll_avg", attr(x,"ll_avg"))
  # browser()
  data.table::setkey(x_dt, loo_iter)
  ref_trt <- ref_trt_orig %>% stringr::str_to_title()
  trt_var <- set_trt(fit)
  
  for (i in seq_along(ids)) {
    # message("Calculating contribution for left out study ", i, " of ", n_id, "...")
    # get loo estimates for left out study i
    id  <- ids[i]
    # filt_dat <- x %>% filter(loo_iter == id)
    # study    <- filt_dat$left_out_study[1]
    filt_dt <- x_dt[J(id), nomatch = 0L] 
    data.table::setattr(filt_dt, "reference_treatment", ref_trt)
    data.table::setattr(filt_dt, "ndraws", attr(x,"ndraws"))
    data.table::setattr(filt_dt, "ntreatments", attr(x,"ntreatments"))
    data.table::setattr(filt_dt, "treatments", attr(x, "treatments"))
    
    study <- filt_dt$left_out_study[1]
    # filt_dt[,left_out_study := NULL]
    loo_i[[i]] <- process_one_dt(filt_dt, sm, 
                                 col_nm, 
                                 trt_var)
    # loo_i[[i]][,left_out_study := study]
    # loo_i[[i]] <- process_one(filt_dat,
    #                           sm,
    #                           col_nm)
    if (i %% 5 == 0) gc(FALSE)
    cli::cli_progress_update()
  }
  cli::cli_progress_done()
  
  influence_measures_dt <- data.table::rbindlist(loo_i, use.names = TRUE, fill = TRUE)
  
  # normalize within trt1/trt2
  influence_measures_dt[, contribution := contribution / sum(contribution), by = .(trt1, trt2)]
  
  influence_measures <- as.data.frame(influence_measures_dt)
  class(influence_measures) <- c("vit_influence", class(influence_measures))
  
  return(influence_measures)
}

# generic function
contribution_matrix <- function(x, fit = NULL, approx = TRUE, ...) {
  
  if (inherits(x, "vitfit")) {
    if (rlang::is_true(approx)) {
      return(contribution_matrix_approx_vitfit(x, approx = TRUE, ...))
    } else {
      loo_estimates(x, ...)
    }
  } else if (inherits(x, "vit_loo")) {
    if(is.null(fit)) rlang::abort("`fit` must be provided when `x` is of class `vit_loo`.")
    check_vitfit(fit)
    
    return(contribution_matrix_vit_loo(x, fit, ...))
    
  } else {
    rlang::abort("No method for class {.cls {class(x)}}.")
  }
  
}

vit_nodesplit <- function(fit, seed, cores = 4L, ...) {
  check_vitfit(fit)
  network <- fit$prep$network
  lvls_trt <- levels(network$treatments)
  
  nodesplit <- multinma::get_nodesplits(network, include_consistency = FALSE)
  if(nrow(nodesplit) == 0) {
    rlang::warn("No comparisons with both direct and independent indirect evidence for node-splitting.")
    return(NULL)
  }
  ns_check <- dplyr::rowwise(nodesplit) %>%
    dplyr::mutate(direct = multinma::has_direct(network, .data$trt1, .data$trt2),
                  indirect = multinma::has_indirect(network, .data$trt1, .data$trt2),
                  valid = .data$direct && .data$indirect)
  
  if (any(!ns_check$valid)) {
    ns_valid <- dplyr::filter(ns_check, .data$valid) %>%
      dplyr::ungroup() %>%
      dplyr::select("trt1", "trt2")
    
    ns_invalid <- dplyr::filter(ns_check, !.data$valid) %>%
      dplyr::mutate(comparison = paste(.data$trt1, .data$trt2, sep = " vs. "))
    
    if (nrow(ns_valid)) {
      rlang::warn(glue::glue(
        "Ignoring node-split comparisons without both both direct and independent indirect evidence: ",
        glue::glue_collapse(ns_invalid$comparison, sep = ", ", width = 100), "."
      ))
      
      nodesplit <- ns_valid
    } else {
      abort("No valid comparisons for node-splitting given in `nodesplit`.\n Comparisons must have both direct and independent indirect evidence for node-splitting.")
    }
    
  }
  
  nodesplit$trt1 <- factor(nodesplit$trt1, levels = lvls_trt)
  nodesplit$trt2 <- factor(nodesplit$trt2, levels = lvls_trt)
  
  for (i in 1:nrow(nodesplit)) {
    if (as.numeric(nodesplit$trt1[i]) > as.numeric(nodesplit$trt2[i])) {
      nodesplit[i, ] <- rev(nodesplit[i, ])
    }
  }
  
  n_ns <- nrow(nodesplit) 
  ns_fits <- vector("list", n_ns)
  
  set.seed(seed)
  for (i in 1:nrow(nodesplit)) {
    
    rlang::inform(glue::glue("Fitting model {i} of {n_ns}, node-split: ",
                      as.character(nodesplit$trt2[i]),
                      " vs. ",
                      as.character(nodesplit$trt1[i])))
    nt1 <- nodesplit$trt1[i] %>% as.character()
    nt2 <- as.character( nodesplit$trt2[i]) 
    nt1c <- nt1 %>% stringr::str_replace_all(", ","\\,") %>%
      stringr::str_replace_all(" ","z\\.")
    nt2c <- nt2 %>% stringr::str_replace_all(", ","\\,") %>%
      stringr::str_replace_all(" ","z\\.")
    direct_study <- network$agd_contrast %>% 
      dplyr::distinct(.study, .trt) |>
      dplyr::summarise(
        hasA = any(.trt == nt1),
        hasB = any(.trt == nt2),
        .by = .study
      ) |>
      dplyr::filter(hasA & hasB) |>
      dplyr::pull(.study) %>% unique()
    
    
    direct <- network$agd_contrast %>% filter(.study %in% direct_study) %>% 
      select(-c(.study,.trt, .y, .se, .sample_size))
    indirect <- network$agd_contrast %>% filter(!(study %in% direct_study)) %>% 
      select(-c(.study,.trt, .y, .se, .sample_size))
    
    ns_fits[[i]] <- list(
      network_fit = fit,
      direct_fit = prep_brms_nma(data = direct,
                                     family = fit$prep$orig_family,
                                     main = fit$prep$orig_form$main,
                                     interactions = fit$prep$orig_form$interactions,
                                     confounders = fit$prep$orig_form$confounders,
                                     addl_form = fit$prep$orig_form$addl_form
                                     # , trt_ref = get_reference_treatment(fit)
                                     ) %>% 
        fit_brms_nma(
          iter =  fit$fit@stan_args[[1]]$iter,
          warmup = fit$fit@stan_args[[1]]$warmup,
          control = fit$fit@stan_args[[1]]$control,
          chains = fit$fit@stan_args %>% length(),
          cores = cores,
          seed = sample.int(.Machine$integer.max, 1),
          silent = fit$stan_args$silent,
          refresh = fit$stan_args$refresh,
          threads = fit$threads,
          ...
        ),
      indirect_fit =  prep_brms_nma(data = indirect,
                                        family = fit$prep$orig_family,
                                        main = fit$prep$orig_form$main,
                                        interactions = fit$prep$orig_form$interactions,
                                        confounders = fit$prep$orig_form$confounders,
                                        addl_form = fit$prep$orig_form$addl_form
                                        # , trt_ref = get_reference_treatment(fit)
      ) %>% 
        fit_brms_nma(
          iter =  fit$fit@stan_args[[1]]$iter,
          warmup = fit$fit@stan_args[[1]]$warmup,
          control = fit$fit@stan_args[[1]]$control,
          chains = fit$fit@stan_args %>% length(),
          cores = cores,
          seed = sample.int(.Machine$integer.max, 1),
          silent = fit$stan_args$silent,
          refresh = fit$stan_args$refresh,
          threads = fit$threads,
          ...
        ),
      nodesplit = forcats::fct_c(nodesplit$trt1[i], nodesplit$trt2[i])
    )
    ns_fits[[i]]$tau <- list(
      network = posterior::as_draws_matrix(fit,variable = "tau"),
      direct = posterior::as_draws_matrix(ns_fits[[i]]$direct_fit,variable = "tau"),
      indirect = posterior::as_draws_matrix((ns_fits[[i]]$indirect_fit),variable = "tau")
    )
    comparison_list <- list(direct = summary_brms_nma(ns_fits[[i]]$direct_fit) %>%
    {
      if (inherits(., "data.table")) {
        .[, .observed := NULL]
      } else {
        dplyr::select(., -.observed)
      }
    } %>%
      contrasts(nodesplit$trt1[i] %>% stringr::str_to_title(),
                                       nodesplit$trt2[i] %>% stringr::str_to_title()),
      indirect = summary_brms_nma(ns_fits[[i]]$indirect_fit)  %>%
      {
        if (inherits(., "data.table")) {
          .[, .observed := NULL]
        } else {
          dplyr::select(., -.observed)
        }
      } %>%
        contrasts(nodesplit$trt1[i] %>% stringr::str_to_title(),
                  nodesplit$trt2[i] %>% stringr::str_to_title()),
      network = summary_brms_nma(fit)  %>%
        {
          if (inherits(., "data.table")) {
            .[, .observed := NULL]
          } else {
            dplyr::select(., -.observed)
          }
        } %>%
        contrasts(nodesplit$trt1[i] %>% stringr::str_to_title(),
                                        nodesplit$trt2[i] %>% stringr::str_to_title())
    )
    ns_fits[[i]]$comparison <- comparison_list$indirect %>% rename(indirect = value)
    ns_fits[[i]]$comparison$direct <- comparison_list$direct$value
    ns_fits[[i]]$comparison$network <- comparison_list$network$value
    ns_fits[[i]]$comparison$value <- ns_fits[[i]]$comparison$indirect - ns_fits[[i]]$comparison$direct
    ns_fits[[i]]$summary <- ns_fits[[i]]$comparison %>% select(-c(direct,indirect,network)) %>% summary()
  }
  
  class(ns_fits) <- c("vit_nodesplit", class(ns_fits))
  return(ns_fits)
}

alter_netcontrib <- function(cm, fit) {
  if(!inherits(cm, "netcontrib")) {
    rlang::abort("Expected contribution.matrix to be of class 'netcontrib' from netmeta::netcontrib().")
  }
  trt_levels <- levels(fit$prep$network$treatments) %>% 
    stringr::str_to_title()
  cm  <- cm$study.random
  trts <- cm$comparison %>% stringr::str_split(":")
  cm <- cm %>% 
    mutate(trt1 = purrr::map_chr(trts, 1) %>% 
             stringr::str_to_title() %>% 
             factor(levels = trt_levels ) %>% 
             as.character(),
           trt2 = purrr::map_chr(trts, 2) %>% 
             stringr::str_to_title() %>% 
             factor(levels = trt_levels) %>% 
             as.character()) %>% 
    dplyr::mutate(
      placebo_in_pair = trt1 == "Placebo" | trt2 == "Placebo",
      trt_lo = dplyr::case_when(
        placebo_in_pair & trt1 == "Placebo" ~ trt2,
        placebo_in_pair & trt2 == "Placebo" ~ trt1,
        TRUE ~ pmin(trt1, trt2)
      ),
      
      trt_hi = dplyr::case_when(
        placebo_in_pair ~ "Placebo",
        TRUE ~ pmax(trt1, trt2)
      )
    ) %>% 
    dplyr::mutate(
      flipped = trt1 != trt_lo,
      
      comparison = paste0(trt_lo, ":", trt_hi),
      
    ) %>% 
    dplyr::mutate(
      trt1 = trt_lo %>% factor(levels = trt_levels ),
      trt2 = trt_hi %>% factor(levels = trt_levels )
    ) %>% 
    dplyr::filter(trt1 != trt2) %>% 
    dplyr::distinct(study, trt1, trt2, .keep_all = TRUE) %>% 
    dplyr::select(-c(trt_lo, trt_hi, flipped, placebo_in_pair))
  return(cm)
}

cinema <- function(fit, mcid = 0.5, 
                         contribution.matrix = NULL, 
                         comparison = NULL,
                         incoherence = NULL,
                           reporting_bias = "Some concerns", # placeholder for now; would need to implement publication bias assessment
                           indirectness = "Low risk", # placeholder for now; would need to implement indirectness
                         ...) {
  
  cm <- contribution.matrix
  if(is.null(incoherence)) incoherence <- "No concerns"
  if(!is.null(comparison)) comparison <- stringr::str_to_title(comparison)
  
  interval_region <- function(lo, hi, mcid) {
    dplyr::case_when(
      lo >  mcid              ~ "pos",
      hi < -mcid              ~ "neg",
      hi >= -mcid & hi < 0 ~ "neg_trivial",
      lo > 0 & lo <= mcid  ~ "pos_trivial",
      TRUE                    ~ "mixed"
    )
  }
  
  # Imprecision grading per CINeMA logic: CI/CrI vs MCID thresholds
  grade_imprecision <- function(lo, hi, mcid) {
    dplyr::case_when(
      lo >  mcid | hi < -mcid                 ~ "No concerns",
      lo > 0 & hi > 0                         ~ "No concerns",
      lo < 0 & hi < 0                         ~ "No concerns",
      lo > -mcid & hi < mcid                  ~ "No concerns",
      lo < -mcid & hi > 0 & hi < mcid         ~ "Some concerns",
      lo > -mcid & lo < 0 & hi > mcid         ~ "Some concerns",
      lo < -mcid & hi >  mcid                 ~ "Major concerns",
      TRUE                                    ~ "Some concerns"
    )
  }
  
  grade_heterogeneity <- function(ci_lo, ci_hi, pi_lo, pi_hi, mcid) {
    ci_reg <- grade_imprecision(ci_lo, ci_hi, mcid)
    pi_reg <- grade_imprecision(pi_lo, pi_hi, mcid)
    
    # "major" if CI spans clinically meaningful benefit and harm
    major <- (pi_reg == "Major concerns" &
                      ci_reg == "No concerns") 
    
    some  <- ((ci_reg == "No concerns" & 
                    pi_reg == "Some concerns" )|
                      (ci_reg == "Some concerns" & 
                        pi_reg == "Major concerns") )
    
    none <- (ci_reg == pi_reg)
    stopifnot(any((major + some + none) == 1))
    
    dplyr::case_when(
      major ~ "Major concerns",
      some  ~ "Some concerns",
      TRUE  ~ "No concerns"
    )
  }
  
  penalty <- c(
    "No concerns"    = 0L,
    "Low risk"       = 0L,
    "Some concerns"  = 1L,
    "Major concerns" = 2L,
    "High risk"      = 2L
  )
  
  score1 <- function(x) {
    out <- unname(penalty[x])
    if (is.na(out)) 1L else out
  }
  
  overall_from_domains <- function(within_study_bias,
                                   reporting_bias,
                                   indirecteness,
                                   imprecision,
                                   heterogeneity,
                                   incoherence,
                                   combine_mode = c("max",
                                                    "max_plus_one")
                                   ) {
    combine_mode <- match.arg(combine_mode)
    
    imp <- score1(imprecision)
    het <- score1(heterogeneity)
    
    variability <- max(imp, het)
    if (combine_mode == "max_plus_one") {
      variability <- min(variability + as.integer(imp >= 1L & het >= 1L), 2L)
    }
    
    ind <- score1(indirecteness)
    inc <- score1(incoherence)
    
    transitivity_penalty = max(ind, inc)
    if (combine_mode == "max_plus_one") {
      transitivity_penalty <- min(transitivity_penalty + as.integer(ind >= 1L & inc >= 1L), 2L)
    }
    
    
    wib <- score1(within_study_bias)
    rep <- score1(reporting_bias)
    
    bias_penalty = max(wib, rep)
    if (combine_mode == "max_plus_one") {
      bias_penalty <- min(bias_penalty + as.integer(wib >= 1L & rep >= 1L), 2L)
    }
    
    total <- variability +
      transitivity_penalty +
      bias_penalty
    
    
    dplyr::case_when(
      total <= 0 ~ "High",
      total == 1 ~ "Moderate",
      total == 2 ~ "Low",
      TRUE       ~ "Very low"
    )
  }
  
  reasons_from_domains <- function(within_study_bias,
                                   reporting_bias,
                                   indirecteness,
                                   imprecision,
                                   heterogeneity,
                                   incoherence,
                                   combine_mode = c("max",
                                                    "max_plus_one")
  ) {
    combine_mode <- match.arg(combine_mode)
    
    imp <- score1(imprecision)
    het <- score1(heterogeneity)
    
    variability <- max(imp, het)
    if (combine_mode == "max_plus_one") {
      variability <- min(variability + as.integer(imp >= 1L & het >= 1L), 2L)
    }
    
    ind <- score1(indirecteness)
    inc <- score1(incoherence)
    
    transitivity_penalty = max(ind, inc)
    if (combine_mode == "max_plus_one") {
      transitivity_penalty <- min(transitivity_penalty + as.integer(ind >= 1L & inc >= 1L), 2L)
    }
    
    
    wib <- score1(within_study_bias)
    rep <- score1(reporting_bias)
    
    bias_penalty = max(wib, rep)
    if (combine_mode == "max_plus_one") {
      bias_penalty <- min(bias_penalty + as.integer(wib >= 1L & rep >= 1L), 2L)
    }
    
    total <- variability +
      transitivity_penalty +
      bias_penalty
    
    
    if (total <= 0) return ("")
    
    out <- ""
    
    if(ind > 0) {
      out <- paste0(out, "Indirectness")
    }
    
    if(inc >0) {
      if (out != "") out <- paste0(out, "; ")
      out <- paste0(out, "Incoherence")
    }
    
    if (imp > 0) {
      if (out != "") out <- paste0(out, "; ")
      out <- paste0(out, "Imprecision")
    }
    
    if (het > 0) {
      if (out != "") out <- paste0(out, "; ")
      out <- paste0(out, "Heterogeneity")
    }
    
    if(wib > 0) {
      if (out != "") out <- paste0(out, "; ")
      out <- paste0(out, "Within-study bias")
    }
    
    if(rep > 0) {
      if (out != "") out <- paste0(out, "; ")
      out <- paste0(out, "Reporting bias")
    }
    
    return(out)
  }
    
  count_direct_studies <- function(fit, ...) {
    extra_vars <- c(...)
    fit$prep$data %>% mutate(
      trt1 = fit$prep$data$intervention,
      trt2 = fit$prep$data$control %>% 
        as.character() %>% 
         gsub(",z\\.",", ",.) %>%
        gsub("z\\.","",.) %>% factor(levels = levels(trt1))) %>% 
      group_by(trt1, trt2, 
               dplyr::across(dplyr::any_of(extra_vars))
               ) %>%
      summarize(n = dplyr::n())
  }
  create_two_trts <- function(x) {
    x %>% 
    mutate(trt1 = stringr::str_split(.trt, "\\] vs. \\[") %>% 
             purrr::map_chr(1) %>%
             gsub("\\[","",.) %>% 
             factor(levels = levels(fit$prep$network$treatments)),
           trt2 = stringr::str_split(.trt, "\\] vs. \\[") %>% 
             purrr::map_chr(2) %>% 
             gsub("\\]","",.) %>% 
             factor(levels = levels(fit$prep$network$treatments))) %>% 
      filter(trt1 != trt2)
  }
  
  grade_within_study_bias <- function(
    contrib_df,
    fit,
    trt1 ="trt1",
    trt2 ="trt2",
    study_col        = "study",
    contrib_col      = "contribution",
    rob_col          = "bias",
    rob_map = c(
      "low" = 0,
      "some concerns" = 1,
      "high" = 2
    ),
    # cutoffs for weighted average
    cutoffs = c(
      no_concerns = 0.5,
      some_concerns = 1.5
    )
  ) {
    stopifnot(
      is.data.frame(contrib_df)
    )
    check_vitfit(fit)
    
    dat <- contrib_df |>
      dplyr::left_join(
        fit$prep$data,
        by = stats::setNames(study_col, study_col)
      )
    
    if (any(is.na(dat[[rob_col]]))) {
      stop("Missing RoB judgments for one or more studies.")
    }
    
    # map RoB to numeric scores
    dat <- dat |>
      dplyr::mutate(
        rob_score = rob_map[.data[[rob_col]]]
      )
    
    if (any(is.na(dat$rob_score))) {
      stop("RoB values not found in rob_map.")
    }
    
    # compute contribution-weighted mean RoB per comparison
    summary <- dat |>
      dplyr::group_by(.data[[trt1]], .data[[trt2]]) |>
      dplyr::summarise(
        weighted_rob = sum(.data[[contrib_col]] * rob_score),
        .groups = "drop"
      )
    
    # map weighted score to CINeMA categories
    summary <- summary |>
      dplyr::mutate(
        within_study_bias = dplyr::case_when(
          weighted_rob < cutoffs["no_concerns"] ~ "No concerns",
          weighted_rob < cutoffs["some_concerns"] ~ "Some concerns",
          TRUE ~ "Major concerns"
        )
      )
    
    summary
  }
  assess_incoherence <- function(incoherence) {
    
  }
  # cinema_tbl <- tibble(
  #   comparison = sm$.trt,
  #   within_study_bias = wsb_grade,
  #   reporting_bias    = pubbias_grade,
  #   indirectness      = indirectness_grade,
  #   imprecision       = imprecision_grade,
  #   heterogeneity     = heterogeneity_grade,
  #   incoherence       = incoherence_grade,
  #   overall_confidence = overall_grade
  # )
  
  check_vitfit(fit)
  
  sm <- summary_brms_nma(fit, ...) %>% 
    # filter(.observed) %>% 
    .[, .observed := NULL] %>% 
    arrange(.trt) %>% 
    upper_tri_contrasts(".trt") %>% 
    rename(trt1 = a, trt2 = b) %>% 
    filter(trt1 != trt2)
  
  keep_vars <- setdiff(colnames(sm), c("iter","value","Estimate", "S.E.", "q2.5%", "q50%","q97.5%"))
  keep_vars <- stats::setNames(keep_vars, keep_vars)
  pop_vars  <- setdiff(keep_vars, c(".trt","trt1","trt2"))
  
  p_i <- prediction_intervals(fit, only_observed = FALSE, ...) 
  p_i %>% data.table::setnames(old = c("a", "b"), new = c("trt1", "trt2")) 
  p_i <- p_i %>% 
   .[trt1 != trt2] %>% 
    summary()
  p_i %>% 
     data.table::setnames(old = c("q2.5%", "q97.5%"), new = c(".pi_lci", ".pi_uci"))
  p_i[, c("Estimate", "S.E.", "q50%") := NULL]
  
  n_study <- count_direct_studies(fit, pop_vars) %>% 
    mutate(trt1 = trt1 %>% stringr::str_to_title(),
           trt2 = trt2 %>% stringr::str_to_title())
  
  n_study <- dplyr::bind_rows(n_study, n_study %>% rename(.trt1 = trt2, trt2 = trt1) %>% rename(trt1 = .trt1))
  
  if (inherits(cm, "vit_loo")) {
    cm <- contribution_matrix(cm, fit, approx = FALSE, ...)
  } else if (inherits(cm, "matrix")) {
    cm <- cm
  } else if (inherits(cm, "netcontrib")) {
    cm <- alter_netcontrib(cm, fit)
  } else if (is.null(cm)) {
    cm <- contribution_matrix(fit, approx = TRUE) %>% alter_netcontrib(fit)
  }
  
  wi_study_bias <- grade_within_study_bias(cm, fit)
  wi_study_bias <- dplyr::bind_rows(wi_study_bias ,
                                    wi_study_bias %>% rename(.trt1 = trt2, 
                                                             trt2 = trt1) %>% 
                                                      rename(trt1 = .trt1)
                                    )
  
  indirectness_df <- data.frame(
    trt1 = p_i$trt1,
    trt2 = p_i$trt2,
    indirectness = "Low risk"
  ) 
  
  if (indirectness %>% is.data.frame()) {
    stopifnot(c(".trt", "rating") %in% colnames(indirectness))
    
    for(i in 1:nrow(indirectness_df)) {
      trt1_i <- indirectness_df$trt1[i]
      trt2_i <- indirectness_df$trt2[i]
      
      rating_i <- indirectness %>%
        filter(.trt == trt1_i) %>%
        dplyr::pull(rating)
      
      if (length(rating_i) == 0) {
        warning(glue::glue("No indirectness rating found for treatment '{trt1_i}'. Defaulting to 'Low risk'."))
        rating_i <- "Low risk"
      }
      
      rating_j <- indirectness %>%
        filter(.trt == trt2_i) %>%
        dplyr::pull(rating)
      
      if (length(rating_j) == 0) {
        warning(glue::glue("No indirectness rating found for treatment '{trt2_i}'. Defaulting to 'Low risk'."))
        rating_j <- "Low risk"
      }
      
      indirectness_df$indirectness[i] <- score1(rating_i) + score1(rating_j)
      indirectness_df$indirectness[i] <- dplyr::case_when(
        indirectness_df$indirectness[i] <= 1 ~ "Low risk",
        indirectness_df$indirectness[i] == 2 ~ "Some concerns",
        indirectness_df$indirectness[i] >= 3 ~ "Major concerns",
        TRUE ~ "Low risk"
      )
    }
    
  }
  indirectness_df <- indirectness_df %>% distinct(trt1, trt2, indirectness)
  
  out <- sm |> 
    summary() |>
    dplyr::left_join(p_i, by = keep_vars) %>% 
    dplyr::left_join(wi_study_bias, by = c("trt1", "trt2")) %>%
    dplyr::left_join(n_study, by = c("trt1", "trt2", pop_vars)) %>% 
    dplyr::left_join(indirectness_df, by = c("trt1", "trt2")) %>% 
    dplyr::mutate(
        placebo_in_pair = trt1 == "Placebo" | trt2 == "Placebo",
        trt_lo = dplyr::case_when(
          placebo_in_pair & trt1 == "Placebo" ~ trt2,
          placebo_in_pair & trt2 == "Placebo" ~ trt1,
          TRUE ~ pmin(trt1, trt2)
        ),
      
        trt_hi = dplyr::case_when(
          placebo_in_pair ~ "Placebo",
          TRUE ~ pmax(trt1, trt2)
        )
    ) %>% 
    dplyr::mutate(
        flipped = trt1 != trt_lo,
        
        trt_key = paste0("[", trt_lo, "] vs. [", trt_hi, "]"),
    ) %>% 
    dplyr::mutate(
        Estimate = dplyr::if_else(flipped, -Estimate, Estimate),
        
        lci  = dplyr::if_else(flipped, -`q97.5%`, `q2.5%`),
        uci = dplyr::if_else(flipped, -`q2.5%`, `q97.5%`),
        
        .pi_lci2  = dplyr::if_else(flipped, -.pi_uci, .pi_lci),
        .pi_uci2 = dplyr::if_else(flipped, -.pi_lci, .pi_uci),
        trt1 = trt_lo,
        trt2 = trt_hi,
        `q2.5%` = lci,
        `q97.5%` = uci,
        .pi_lci = .pi_lci2,
        .pi_uci = .pi_uci2
    ) %>% 
    # filter(trt1 == trt_lo & trt2 == trt_hi) %>% 
    dplyr::mutate(
      comparison = trt_key,
      reporting_bias = reporting_bias, # placeholder for now; would need to implement publication bias assessment
      # indirectness = indirectness, # placeholder for now; would need to implement indirectness
      imprecision = grade_imprecision(`q2.5%`, `q97.5%`, mcid),
      heterogeneity = grade_heterogeneity(`q2.5%`, `q97.5%`, .pi_lci, .pi_uci, mcid),
      incoherence = incoherence,
      overall_confidence = purrr::pmap_chr(
        list(within_study_bias, reporting_bias, indirectness, imprecision, heterogeneity, incoherence),
        overall_from_domains
      ),
      reasons = purrr::pmap_chr(
        list(within_study_bias, reporting_bias, indirectness, imprecision, heterogeneity, incoherence),
        reasons_from_domains
      )
    ) |>
    dplyr::select(
      comparison,
      trt1, trt2,
      dplyr::any_of(keep_vars %>% setdiff(c("trt1", "trt2",".trt"))),
      n_study = n,
      Estimate,
      "prediction_q2.5%" = .pi_lci,
      `q2.5%`, `q97.5%`,
      "prediction_q97.5%" = .pi_uci,
      within_study_bias,
      reporting_bias,
      indirectness,
      imprecision,
      heterogeneity,
      incoherence,
      overall_confidence,
      reasons
    ) %>% arrange(is.na(n_study), trt1, trt2) 
  
  if(!is.null(comparison) && !(comparison %in% out$trt2)) {
     rlang::warn(glue::glue("Specified comparison '{comparison}' not found in results. Available comparisons: {glue::glue_collapse(unique(out$trt2), sep = ', ')}. Returning basic table for now."))
  } else if(!is.null(comparison)) {
    out <- out %>% 
      { 
        if (is.null(comparison)) . 
        else filter(., trt2 == !!comparison)
      }
  }
  
  if(inherits(out, "data.table")) {
    data.table::setattr(out, "class", c("vit_cinema", class(out)))
  } else {
    class(out) <- c("vit_cinema", class(out))
  }
  
  out
  
}

cinema_table <- function(x, evidence = c("all","mixed", "indirect"), ...) {
  stopifnot(inherits(x, "vit_cinema"))
  
  gt_cinema_colors <- function(gt_tbl,
                               cols,
                               palette = c(
                                 "High"      = "#1a9850",
                                 "Moderate"  = "#4575b4",
                                 "Low"       = "#fdae61",
                                 "Very low"  = "#d73027",
                                 "No concerns"    = "#1a9850",
                                 "Some concerns"  = "#fdae61",
                                 "Major concerns" = "#d73027",
                                 "Low risk"  = "#1a9850",
                                 "High risk" = "#d73027"
                               )) {
    
    # apply one tab_style per level, across all requested columns
    for (col in cols) {
      col_sym <- rlang::sym(col)
      
      for (lvl in names(palette)) {
        gt_tbl <- gt_tbl %>%
          gt::tab_style(
            style = gt::cell_fill(color = palette[[lvl]]),
            locations = gt::cells_body(
              columns = gt::all_of(col),
              rows = !!rlang::expr(!!col_sym == !!lvl)
            )
          )
      }
    }
    
    gt_tbl
  }
  
  cinema_cols <- c(
    "within_study_bias", "reporting_bias", "indirectness",
    "imprecision", "heterogeneity", "incoherence",
    "overall_confidence", "reasons"
  )
  
  evidence <- match.arg(evidence)
  if(evidence == "all") evidence <- c("mixed","indirect")
  
  x %>%
    select(-c(trt1, trt2, 
              # Estimate, `q2.5%`, `q97.5%`, `prediction_q2.5%`, `prediction_q97.5%`
              )) %>%
    mutate(comparison = stringr::str_replace_all(as.character(comparison), "\\] vs\\. \\[", " vs. ") %>% 
             stringr::str_replace_all("\\[|\\]", "")) %>% 
    mutate(comparison = comparison %>%  factor(levels = unique(comparison))) %>%
    mutate(evidence = ifelse(!is.na(n_study) & n_study > 0, "mixed", "indirect")) %>%
    filter(evidence %in% !!evidence) %>%
    mutate(evidence = factor(evidence, levels = c("mixed", "indirect"), 
                             labels = c("Mixed Estimates", "Indirect Estimates"))) %>%
    gt::gt(groupname_col = "evidence") %>%
    gt::fmt_number() %>% 
    gt::fmt_number(columns = n_study, decimals = 0) %>% 
    gt::cols_merge(
      columns = c(Estimate, `q2.5%`, `q97.5%`),
      pattern = "{1} ({2}, {3})"
    ) %>%
    gt::cols_merge(
      columns = c(`prediction_q2.5%`, `prediction_q97.5%`),
      pattern = "({1}, {2})"
    ) %>%
    gt::tab_header(title = "Confidence in Network Meta-Analysis (CINeMA)") %>%
    gt_cinema_colors(cols = cinema_cols) %>% 
    gt::cols_label(
      evidence = "Evidence",
      comparison = "Comparison",
      n_study = "No.of Studies",
      Estimate = "Effect (95% CrI)",
      `prediction_q2.5%` = "95% PI",
      within_study_bias = "Within-study bias",
      reporting_bias = "Reporting bias",
      indirectness = "Indirectness",
      imprecision = "Imprecision",
      heterogeneity = "Heterogeneity",
      incoherence = "Incoherence",
      overall_confidence = "Overall confidence",
      reasons = "Reasons for downgrading"
    ) %>% 
    gt::sub_missing()
  
  # x %>%
  #   gt() %>%
  #   tab_header(
  #     title = "Confidence in Network Meta-Analysis (CINeMA)"
  #   ) %>%
  #   fmt_markdown(columns = everything()) %>%
  #   tab_style(
  #     style = cell_fill(color = "#d73027"),
  #     locations = cells_body(
  #       rows = within_study_bias == "Very low"
  #     )
  #   ) %>%
  #   tab_style(
  #     style = cell_fill(color = "#d73027"),
  #     locations = cells_body(
  #       rows = overall_confidence == "Very low"
  #     )
  #   ) %>% 
  #   tab_style(
  #     style = cell_fill(color = "#fdae61"),
  #     locations = cells_body(
  #       rows = overall_confidence == "Low"
  #     )
  #   ) %>% 
  #   tab_style(
  #     style = cell_fill(color = "#fdae61"),
  #     locations = cells_body(
  #       rows = overall_confidence == "Moderate"
  #     )
  #   ) %>% 
  #   tab_style(
  #     style = cell_fill(color = "#1a9850"),
  #     locations = cells_body(
  #       rows = overall_confidence == "High"
  #     )
  #   )
}

dodge_pointinterval_layers <- function(p, width = 0.6) {
  for (i in seq_along(p$layers)) {
    lyr <- p$layers[[i]]
    stat_class <- class(lyr$stat)
    if (any(grepl("^StatPointinterval", stat_class))) {
      # p$layers[[i]]$position <- ggplot2::position_dodge2(
      #   width = width, padding = 0.35, preserve = "single"
      # )
      p$layers[[i]]$position <- ggplot2::position_dodge(
        width = width, preserve = "single"
      )
      # optional: tweak aesthetics/params of the existing layer
      # p$layers[[i]]$aes_params$size <- 3
      # if you mapped xmin/xmax yourself earlier, ensure they’re in p$mapping
      # and that subset is used for grouping:
      if (is.null(p$mapping$group)) p$mapping$group <- rlang::sym("bias")
    }
  }
  p
}

prediction_intervals  <- function(x, tau_name = "tau", only_observed = TRUE, all.contrasts = TRUE, ...) {
  check_vitfit(x)
  
  sm <- summary_brms_nma(x, ...) 
  tau_dt <- dplyr::tibble(
    tau = brms::as_draws_matrix(x, variable = tau_name) %>% 
    as.numeric()) %>% 
    mutate(iter = dplyr::row_number()) %>% 
    data.table::as.data.table()
  
  if (rlang::is_true(only_observed) ) {
    sm <- sm %>% filter(.observed)
  }
  data.table::set(sm, j = ".observed", value = NULL)
  
  smcnt <-if (rlang::is_true(all.contrasts)) {
     sm %>% upper_tri_contrasts(set_trt(x)) %>% .[a != b]
  } else {
     sm %>% mutate(a = .trt, b = "Placebo")
  }
  
  # sm <- sm %>%
  #   dplyr::left_join(tau, by = "iter") %>%
  #   dplyr::mutate(
  #     u = stats::rnorm(dplyr::n(), mean = 0, sd = .data[[tau_name]]),
  #     value = value + u
  #   ) %>%
  #   dplyr::select(-u, -dplyr::all_of(tau_name))
  # If tau has only one row per iter, this is a normal left join
  smcnt[tau_dt, on = "iter", tau := i.tau]
  # draw u and update value
  smcnt[, u := stats::rnorm(.N, mean = 0, sd = tau)]
  smcnt[, value := value + u]
  
  # drop u and tau column
  smcnt[, c("u", "tau") := NULL]
  
  return(smcnt)
}

change_and_combine_forest <- function(data, bias, width) {
  data %>% 
    mutate(bias = bias) %>%
    mutate(Treatment = stringr::str_wrap(Treatment, width = width)) 
}

re_at <- function(level) {
  levels_bias <- c("low","some concerns","high")
  nd <- data.frame(bias = factor(level, levels = levels_bias),
                   study = "1")
  re <- relative_effects(d_fit_RE_all_bias,
                         newdata = nd,
                         study = "study",
                         summary = FALSE)
  # convert to tidy draws (one row per draw per treatment comparison)
  # multinma provides as.data.frame() for relative_effects; set summary = FALSE to get draws
  re_df <- as.data.frame(re, summary = FALSE)
  # Expect columns like d[<trt>] or similar per comparison; pivot long:
  re_long <- re_df %>%
    mutate(.draw = row_number()) %>%
    pivot_longer(-.draw, names_to = "contrast", values_to = "effect") %>%
    mutate(bias = level)
  re_long
}

bias_colors <- c("low" = "#1a9850", 
                 "low/some concerns" = "#fdae61", 
                 "low/some concerns/high" = "#d73027")

bias_colors <- c(
  "low" = "#3A8F6A",        # muted green
  "low/unclear" = "#D9903D",  # warm amber
  "low/unclear/high" = "#C43C3C"       # softer red
)

def_colors <- c(
  "yes" = "#6a51a3",  # purple
  "no"  = "#5383c7",  # blue
  "unknown"  = "grey63"  # grey
)

width <- 18
pres.size <- 16
pdf.size <- 11


bias_study_colors <- bias_colors
names(bias_study_colors) <- c("low","unclear","high")

box::register_S3_method("print", "summary_brms_nma", print.summary_brms_nma)
box::register_S3_method("plot", "summary_brms_nma", plot.summary_brms_nma)
box::register_S3_method("plot", "vit_rank_probs", plot.vit_rank_probs)
box::register_S3_method("log_lik", "vitfit", log_lik.vitfit)
box::register_S3_method("summary", "vitfit", summary.vitfit)
box::register_S3_method("summary", "summary_brms_nma", summary.summary_brms_nma)
box::register_S3_method("marginalize", "summary_brms_nma", marginalize.summary_brms_nma)
box::register_S3_method("marginalize", "vitfit", marginalize.vitfit)
