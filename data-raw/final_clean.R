library(dplyr)
library(tidyr)
library(tibble)
library(forcats)
library(stringr)

#### Setup ####
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

#### Read in CSV ####
data <- read.csv("data-raw/raw_clean.csv") %>% 
  filter(Study.Year.Country != "")

#### Clean data ####
labels <- read.csv("data-raw/raw_clean.csv",nrows = 1, header = FALSE)[1,] %>% 
  as.character() 

# R useful column names
c("study", "title", "design",
  "blinded", "target",
  "population",
  "baseline.vitamin",
  "end.vitamin",
  "selection.criteria",
  "psych.dx","vit.def",
  "pop.dx","total.N", 
  "intervention", "dose",
  "frequency", "tdd",
  "number.treatment.arms",
  "assigned.N", "duration", 
  "scale","baseline.N","baseline.outcome.mean",
  "baseline.outcome.sd", 
  "final.N",
  "final.outcome.mean",
  "final.outcome.sd","final.outcome.se",
  "age","male",
  "findings","adverse.events",
  "bias.sequence", "bias.concealment",
  "bias.blinding.participants", "bias.blinding.outcome",
  "bias.incomplete.outcome", "bias.selective.outcome",
  "bias.other","bias.overall","notes","flag",
  "gain.score.mean","gain.score.sd","gain.score.se","p",
  "secondary.flag","secondary.gm.mean","secondary.gm.sd","adverse.rate") -> cn

colnames(data) <- cn

# get and separate vitamin info
vb <- data$baseline.vitamin 
ve <- data$end.vitamin

vitamin.baseline <- data.frame(
                      id = data$study,
                      intervention = data$intervention,
                      dose = data$tdd,
                      name = vb) %>% 
  distinct(id, intervention, dose, .keep_all = TRUE) %>% 
  separate_rows(name, sep = ",\\s*") %>%
  separate_wider_delim(name, delim = ": ", names = c("name", "value"),
                       too_few = "align_start"
                       ) %>%
  mutate(name = ifelse(str_starts(name, "normal|insufficient|deficient|deficiency"),
                       "",name)) %>% 
  mutate(value = ifelse(str_starts(value, "deficiency"),
                       "",value)) %>% 
  mutate(name = fct_recode(as.character(name), "vitamin D" = "25(OH)D",
                           "ferritin" = "Ferritin",
                           "vitamin B9" ="folate",
                           "vitamin D" = "vitamin d",
                           "vitamin B12" = "vitamin b12",
                           "vitamin B2" = "Vitamin B2",
                           "selenoproteinP" = "selenoprotein P",
                           "rbc.B9" = "RBC vitamin B9"),
         name = fct_na_level_to_value(name, extra_levels = "")) %>% 
  mutate(value = ifelse(name %>% is.na(), NA, value)) %>% 
  separate_wider_regex(
    value,
    patterns = c(
      "\\s*",
      value = "\\S+",   # first non-space chunk
      "\\s+",           # one or more spaces (dropped; unnamed)
      units = ".*"      # the rest
    ),
    too_few = "align_start"
  ) %>% 
  mutate(value = as.numeric(value)) %>% 
  mutate(
    std_value = ifelse(grepl("g",units), 
           vit$convert_conc(value, vit$molecular_mass[as.character(name)], units),
           vit$convert_to_molar(value, units))) %>% 
  distinct() %>% 
  mutate(cutoff = vit$cutoff_vitamins[as.character(name)]) %>% 
  mutate(pct_deficient = ifelse(name %>% is.na(), NA, 
                                vit$log_normal_prob_from_mom(std_value,
                                                             0.3 * std_value,
                                                             cutoff,
                                                             lower = TRUE
                                ))
  ) %>% 
  mutate(name = str_remove(as.character(name), "^vitamin\\s+") %>% as.factor()) 
vitamin.baseline[vitamin.baseline$id == "Nguyen/2009/Guatemala", "pct_deficient"] <- c(.155,.175,.287,.171,.231,.167, .167,.140)


  
vitamin.end <- data.frame(
                  id = data$study,
                  intervention = data$intervention,
                  dose = data$tdd,
                  name = ve) %>% 
  distinct(id, intervention, dose, .keep_all = TRUE) %>% 
  separate_rows(name, sep = ",\\s*") %>%
  separate_wider_delim(name, delim = ": ", 
                       names = c("name", "value"),
                       too_few = "align_start") %>% 
  mutate(name = ifelse(str_starts(name, "normal|insufficient|deficient|deficiency"),
                       "",name)) %>% 
  mutate(value = ifelse(str_starts(value, "deficiency"),
                        "",value)) %>% 
  mutate(name = fct_recode(as.character(name), "vitamin D" = "25(OH)D",
                           "ferritin" = "Ferritin",
                           "vitamin B9" ="folate",
                           "vitamin D" = "vitamin d",
                           "vitamin B12" = "vitamin b12",
                           "vitamin B12" = "vitaminb12",
                           "vitamin B2" = "Vitamin B2",
                           "selenoproteinP" = "selenoprotein P",
                           "rbc.B9" = "RBC vitamin B9"),
         name = fct_na_level_to_value(name, extra_levels = "")) %>% 
  mutate(value = ifelse(name %>% is.na(), NA, value)) %>% 
  separate_wider_regex(
    value,
    patterns = c(
      value = "\\S+",   # first non-space chunk
      "\\s+",           # one or more spaces (dropped; unnamed)
      units = ".*"      # the rest
    ),
    too_few = "align_start"
  ) %>% 
  mutate(value = as.numeric(value)) %>% 
  mutate(std_value = 
           ifelse(
             # units %in% c("ng/mL", "md/dL","ug/L","pg/mL","mg/dL", "ng/mL", "µg/L","µg/l"), 
             grepl("g",units),
                            vit$convert_conc(value, vit$molecular_mass[as.character(name)], units),
                            vit$convert_to_molar(value, units))) %>% 
  distinct() %>% 
  mutate(name = str_remove(as.character(name), "^vitamin\\s+") %>% as.factor())
  
# combine vitamin data
vitamins <- vitamin.baseline %>% 
  select(id, intervention, name, dose, std_value, pct_deficient) %>%
  pivot_wider(id_cols = c("id","intervention","dose"),
              names_from = "name",
              names_prefix = "baseline.",
              values_from = "std_value") %>% 
  select(-`baseline.NA`) %>% 
  full_join(vitamin.end %>%  
              select(id, intervention, name, dose, std_value) %>%
              pivot_wider(id_cols = c("id","intervention","dose"),
                          names_from = "name",
                          names_prefix = "end.",
                          values_from = "std_value") %>% 
              select(- `end.NA`) 
            , by = c("id", "intervention", "dose")
            )

vit_tbl <- tribble(
  ~var,                  ~cutoff_name,       ~direction, ~vitamin_label,
  "baseline.zinc",        "zinc",             "<=",       "zinc",
  "baseline.B9",          "vitamin B9",       "<=",       "vitamin B9",
  "baseline.homocysteine","homocysteine",     ">=",       "vitamin B12",   # functional marker
  "baseline.B1",          "vitamin B1",       "<=",       "vitamin B1",
  "baseline.B2",          "vitamin B2",       "<=",       "vitamin B2",
  "baseline.B6",          "vitamin B6",       "<=",       "vitamin B6",
  "baseline.D",           "vitamin D",        "<=",       "vitamin D",
  "baseline.B12",         "vitamin B12",      "<=",       "vitamin B12",
  "baseline.magnesium",   "magnesium",        "<=",       "magnesium",
  "baseline.ferritin",    "ferritin",         "<=",       "ferritin",
  "baseline.C",           "vitamin C",        "<=",       "vitamin C",
  "baseline.selenium",    "selenium",         "<=",       "selenium",
  "baseline.selenoproteinP","selenoproteinP", "<=",       "selenium",      # proxy
  "baseline.methylmalonate","methylmalonate", ">=",       "vitamin B12",    # functional marker
  "baseline.rbc.B9",      "rbc.folate",       "<=",       "vitamin B9"
)
cut_tbl <- enframe(vit$cutoff_vitamins, name = "cutoff_name", value = "cutoff")

# add deficiency flags
for(j in 1:nrow(vit_tbl)) vitamins[[paste0("def.", vit_tbl$vitamin_label[j])]] <- NA_real_
vn <- gsub("def.", "", names(vitamins)[which(startsWith(names(vitamins), "def."))])
for(i in 1:nrow(vitamins)) {
  for(j in 1:nrow(vit_tbl)) {
    var <- vit_tbl$var[j]
    cut_name <- vit_tbl$cutoff_name[j]
    direction <- vit_tbl$direction[j]
    vit_label <- vit_tbl$vitamin_label[j]
    cutoff <- cut_tbl %>% filter(cutoff_name == cut_name) %>% pull(cutoff)
    col_label <- paste0("def.", vit_label)
    if(vitamins[i, var] %>% is.na()) next
    vitamins[i, col_label] <- ifelse(is.na(vitamins[i, var]), NA_integer_, 0L)
    if(direction == "<=") {
      vitamins[i, col_label] <- ifelse(!is.na(vitamins[i, var]) & 
                                                        vitamins[i, var] <= cutoff,
                                                      1L, vitamins[i, col_label])
    } else {
      vitamins[i, col_label] <- ifelse(!is.na(vitamins[i, var]) & 
                                                        vitamins[i, var] >= cutoff,
                                                      1L, vitamins[i, col_label])
    }
  }
  check_row <- vitamins[i, which(startsWith(names(vitamins), "def."))]
  if (any(!is.na(check_row)) && any (check_row[!is.na(check_row)] == 1)) {
    def.ind <- which(vitamins[i, which(startsWith(names(vitamins), "def."))] == 1)
    vitamins[i,"calc.vit.def"] <- paste(
      vn[def.ind],
      collapse = ", "
    )
  } else if(any(!is.na(check_row)) && any(check_row[!is.na(check_row)] == 0)){
    vitamins[i,"calc.vit.def"] <- "none"
  } else {
    vitamins[i,"calc.vit.def"] <- NA_character_
  }
  
  if(!is.na(vitamins[i,"calc.vit.def"]) && vitamins[i,"calc.vit.def"] == "") vitamins[i,"calc.vit.def"] <- NA_character_
}

vitamins <- vitamins %>% 
  group_by(id) %>% 
  mutate(cc = strsplit(calc.vit.def, ", ")) %>%
  mutate(calc.vit.def.total = ifelse(!is.na(calc.vit.def),
                               paste(unique(unlist(cc[!is.na(calc.vit.def)])),collapse = ","),
                               NA_character_)
  ) %>% 
    select(-cc)

#### dose data ####
dose <- data %>% 
  select(study, intervention, tdd) %>% 
  distinct() %>% 
  mutate(tdd = as.character(tdd),
         group = intervention,
         orig.tdd = tdd,
         dummy = 1L) %>% 
  separate_rows(intervention, tdd, sep = ",\\s*") %>% 
  mutate(intervention = as.factor(intervention)) %>%
  mutate(intervention = fct_recode(intervention,
                                   "calcium"    = "Calcium",
                                   "iron"       = "Iron",
                                   "placebo"    = "Placebo",
                                   "antidepressant"   = "TCA/SSRI",
                                   "antidepressant"   = "SSRI/TCA",
                                   "magnesium"  = "Magnesium",
                                   "magnesium"  = "magnesium ",
                                   "multinutrient" = "multivitamin",
                                   "multinutrient" = "multi-nutrient",
                                   "vitamin B3" = "Vitamin B 3",
                                   "vitamin B3" = "Vitamin B3",
                                   "vitamin B1" = "Vitamin B1",
                                   "vitamin B1" = "vitamin b1",
                                   "vitamin B2" = "Vitamin B2",
                                   "vitamin B2" = "vitamin b2",
                                   "vitamin B12"= "vitamin b12",
                                   "vitamin B3" = "vitamin b3",
                                   "vitamin B5" = "vitamin b5",
                                   "vitamin B5" = "Vitamin B5",
                                   "vitamin B6" = "vitamin b6",
                                   "vitamin B6" = "Vitamin B6",
                                   "vitamin B7" = "vitamin b7",
                                   "vitamin B7" = "Vitamin B7",
                                   "vitamin B9" = "vitamin b9",
                                   "vitamin B9" = "Vitamin B9",
                                   "vitamin B12"= "Vitamin B12",
                                   "vitamin C"  = "Vitamin C",
                                   "vitamin D"  = "vitamin d",
                                   "vitamin D"  = "Vitamin D high",
                                   "vitamin D"  = "vitamin D low",
                                   "vitamin D"  = "Vitamin D",
                                   "zinc"       = "Zinc",
                                   "fluoxetine" = "Fluoxetine",
                                   "citalopram" = "Citalopram",
                                   "escitalopram" = "Escitalopram",
                                   "none" = "None",
                                   "aerobic_exercise" ="aerobic exercise",
  )) %>% 
  mutate(intervention.names = str_remove(as.character(intervention), 
                                         "^vitamin\\s+") ) %>%
  mutate(tdd.char = tdd) %>% 
  mutate(unit = stringr::str_extract(tdd.char, "[a-zA-Z/μµ]+"),
         tdd = stringr::str_extract(tdd.char, "\\d+\\.?\\d*")) %>%
  mutate(tdd = as.numeric(tdd)) %>%
  mutate(tdd = ifelse(is.na(tdd), 
                      0,
                      ifelse(unit %in% c("IU","iu") & intervention == "vitamin D",
                             tdd/40,
                             tdd))) %>%
  mutate(tdd = ifelse(unit %in% c("IU","iu") & intervention == "vitamin A",
                      tdd * 0.3,
                      tdd)) %>%
  mutate(unit = ifelse(unit %in% c("IU","iu"),"ug", unit)) %>%
  mutate(unit = ifelse(is.na(unit), "g", unit)) %>%
  mutate(orig.unit = unit) %>% 
  mutate(tdd = vit$convert_mass_to_g(tdd, unit)) %>%
  mutate(unit = "g") %>% 
  mutate(
    dose_band = case_when(
      # --- Vitamin A (IU) ---
      intervention == "vitamin A"    & tdd <  vit$convert_mass_to_g(900 * 0.3, "ug")   ~ "< 900 IU",
      intervention == "vitamin A"    & tdd <= vit$convert_mass_to_g(3000 * 0.3, "ug")  ~ "900 - 3000 IU",
      intervention == "vitamin A"    & tdd >  vit$convert_mass_to_g(3000 * 0.3, "ug")  ~ "> 3000 IU",
      
      # --- Vitamin D (IU) --- (optional, remove if not needed)
      # 1 IU = 0.025 µg
      intervention == "vitamin D"    & tdd <  vit$convert_mass_to_g(1000 * 0.025, "ug")   ~ "< 1000 IU",
      intervention == "vitamin D"    & tdd <= vit$convert_mass_to_g(3000 * 0.025, "ug")  ~ "1000 - 3000 IU",
      intervention == "vitamin D"    & tdd >  vit$convert_mass_to_g(3000 * 0.025, "ug")  ~ "> 3000 IU",
      
      # --- Vitamin B1 (mg) ---
      intervention == "vitamin B1"   & tdd <  vit$convert_mass_to_g(10,  "mg")        ~ "< 10 mg",
      intervention == "vitamin B1"   & tdd <= vit$convert_mass_to_g(100, "mg")        ~ "10 - 100 mg",
      intervention == "vitamin B1"   & tdd >  vit$convert_mass_to_g(100, "mg")        ~ "> 100 mg",
      
      # --- Vitamin B2 (mg) ---
      intervention == "vitamin B2"   & tdd <  vit$convert_mass_to_g(10, "mg")         ~ "< 10 mg",
      intervention == "vitamin B2"   & tdd <= vit$convert_mass_to_g(25, "mg")         ~ "10 - 25 mg",
      intervention == "vitamin B2"   & tdd >  vit$convert_mass_to_g(25, "mg")         ~ "> 25 mg",
      
      # --- Vitamin B3 (mg) ---
      intervention == "vitamin B3"   & tdd <  vit$convert_mass_to_g(50,  "mg")        ~ "< 50 mg",
      intervention == "vitamin B3"   & tdd <= vit$convert_mass_to_g(500, "mg")        ~ "50 - 500 mg",
      intervention == "vitamin B3"   & tdd >  vit$convert_mass_to_g(500, "mg")        ~ "> 500 mg",
      
      # --- Vitamin B5 (mg) ---
      intervention == "vitamin B5"   & tdd <  vit$convert_mass_to_g(10,  "mg")        ~ "< 10 mg",
      intervention == "vitamin B5"   & tdd <= vit$convert_mass_to_g(100, "mg")        ~ "10 - 100 mg",
      intervention == "vitamin B5"   & tdd >  vit$convert_mass_to_g(100, "mg")        ~ "> 100 mg",
      
      # --- Vitamin B6 (mg) ---
      intervention == "vitamin B6"   & tdd <  vit$convert_mass_to_g(10, "mg")         ~ "< 10 mg",
      intervention == "vitamin B6"   & tdd <= vit$convert_mass_to_g(50, "mg")         ~ "10 - 50 mg",
      intervention == "vitamin B6"   & tdd >  vit$convert_mass_to_g(50, "mg")         ~ "> 50 mg",
      
      # --- Vitamin B7 (biotin, µg) ---
      intervention == "vitamin B7"   & tdd <  vit$convert_mass_to_g(30,  "ug")        ~ "< 30 µg",
      intervention == "vitamin B7"   & tdd <= vit$convert_mass_to_g(300, "ug")        ~ "30 - 300 µg",
      intervention == "vitamin B7"   & tdd >  vit$convert_mass_to_g(300, "ug")        ~ "> 300 µg",
      
      # --- Vitamin B9 (folic acid, mg) ---
      intervention == "vitamin B9"   & tdd <  vit$convert_mass_to_g(0.4, "mg")        ~ "< 0.5 mg",
      intervention == "vitamin B9"   & tdd <= vit$convert_mass_to_g(5,   "mg")        ~ "0.5 - 5 mg",
      intervention == "vitamin B9"   & tdd >  vit$convert_mass_to_g(5,   "mg")        ~ "> 5 mg",
      
      # --- Vitamin B12 (µg) ---
      intervention == "vitamin B12"  & tdd <  vit$convert_mass_to_g(250,  "ug")       ~ "< 250 µg",
      intervention == "vitamin B12"  & tdd <= vit$convert_mass_to_g(1000, "ug")       ~ "250 - 1000 µg",
      intervention == "vitamin B12"  & tdd >  vit$convert_mass_to_g(1000, "ug")       ~ "> 1000 µg",
      
      # --- Vitamin C (mg) ---
      intervention == "vitamin C"    & tdd <  vit$convert_mass_to_g(200,  "mg")       ~ "< 200 mg",
      intervention == "vitamin C"    & tdd <= vit$convert_mass_to_g(1000, "mg")       ~ "200 - 1000 mg",
      intervention == "vitamin C"    & tdd >  vit$convert_mass_to_g(1000, "mg")       ~ "> 1000 mg",
      
      # --- Vitamin E (IU) ---
      # assume 1 IU ≈ 0.67 mg alpha-tocopherol
      intervention == "vitamin E"    & tdd <  vit$convert_mass_to_g(100 * 0.67, "mg") ~ "< 100 IU",
      intervention == "vitamin E"    & tdd <= vit$convert_mass_to_g(400 * 0.67, "mg") ~ "100 - 400 IU",
      intervention == "vitamin E"    & tdd >  vit$convert_mass_to_g(400 * 0.67, "mg") ~ "> 400 IU",
      
      # --- Iron (mg elemental) ---
      intervention == "iron"         & tdd <  vit$convert_mass_to_g(18, "mg")         ~ "< 18 mg",
      intervention == "iron"         & tdd <= vit$convert_mass_to_g(45, "mg")         ~ "18 - 45 mg",
      intervention == "iron"         & tdd >  vit$convert_mass_to_g(45, "mg")         ~ "> 45 mg",
      
      # --- Magnesium (mg) ---
      intervention == "magnesium"    & tdd <  vit$convert_mass_to_g(200, "mg")        ~ "< 200 mg",
      intervention == "magnesium"    & tdd <= vit$convert_mass_to_g(400, "mg")        ~ "200 - 400 mg",
      intervention == "magnesium"    & tdd >  vit$convert_mass_to_g(400, "mg")        ~ "> 400 mg",
      
      # --- Selenium (µg) ---
      intervention == "selenium"     & tdd <  vit$convert_mass_to_g(100, "ug")        ~ "< 100 µg",
      intervention == "selenium"     & tdd <= vit$convert_mass_to_g(200, "ug")        ~ "100 - 200 µg",
      intervention == "selenium"     & tdd >  vit$convert_mass_to_g(200, "ug")        ~ "> 200 µg",
      
      # --- Zinc (mg) ---
      intervention == "zinc"         & tdd <  vit$convert_mass_to_g(15, "mg")         ~ "< 15 mg",
      intervention == "zinc"         & tdd <= vit$convert_mass_to_g(30, "mg")         ~ "15 - 30 mg",
      intervention == "zinc"         & tdd >  vit$convert_mass_to_g(30, "mg")         ~ "> 30 mg",
      
      intervention == "placebo"  ~ "placebo",
      intervention == "none"     ~ "none",
      
      TRUE ~ NA_character_
    )
  ) %>%
  pivot_wider(id_cols = c("study","group", "orig.tdd"),
              names_from = c("intervention.names"),
              names_sort = TRUE,
              values_from = c("dummy","tdd", "dose_band"),
              values_fill = list(dummy = 0, tdd = 0, dose_band = NA_character_)
              # names_prefix = c("z.")
  ) %>% 
  rename_with(~ sub("^dummy_", "z\\.", .x), starts_with("dummy_")) %>%
  rename_with(~ sub("^tdd_",   "d\\.", .x), starts_with("tdd_")) %>%
  rename_with(~ sub("^dose_band_",   "cd\\.", .x), starts_with("dose_band_")) %>%
  mutate(z.antidepressant = case_when(z.citalopram == 1 ~ 1,
                                      z.escitalopram == 1 ~ 1,
                                      z.imipramine == 1 ~ 1,
                                      z.fluoxetine == 1 ~ 1,
                                      z.nortriptyline == 1 ~ 1,
                                      z.amitriptyline == 1 ~ 1,
                                      z.SSRI == 1 ~ 1,
                                      # z.SSRI_TCA == 1 ~ 1,
                                      z.antidepressant == 1 ~ 1,
                                      .default = 0)) %>%
  select(-c(z.amitriptyline,
            z.escitalopram,
            z.citalopram,
            z.imipramine,
            z.fluoxetine,
            z.nortriptyline,
            z.SSRI)) %>% 
  select(-c(d.amitriptyline,
            d.escitalopram,
            d.citalopram,
            d.imipramine,
            d.fluoxetine,
            d.nortriptyline,
            d.SSRI)) %>% 
  mutate(d.antidepressant = ifelse(z.antidepressant == 1,
                                   1, d.antidepressant)) %>%
  mutate(d.aerobic_exercise = ifelse(z.aerobic_exercise == 1,
                                   1, d.aerobic_exercise)) %>%
  mutate(cd.antidepressant = ifelse(z.antidepressant == 1,
                                   "various", cd.antidepressant)) %>%
  mutate(cd.aerobic_exercise = ifelse(z.aerobic_exercise == 1,
                                     "N/A", cd.aerobic_exercise)) %>%
  mutate(
    intervention = apply(select(., z.A:z.zinc) == 1, 1, function(x) {
      paste(names(select(., z.A:z.zinc))[x], collapse = ",")
    }),
    intervention = as.factor(intervention),
    tdd = apply(select(., d.A:d.zinc) , 1, function(x) {
      paste0(sub("d\\.","",names(select(., d.A:d.zinc))[which(x > 0)]),
             " (", x[which(x > 0)], " g)", collapse = ", ")
    }),
    tdd = ifelse(intervention == "z.placebo" | intervention == "z.none", NA_character_,tdd),
    tdd = gsub("antidepressant (1 g)", "antidepressant (various)",tdd),
    tdd = gsub("aerobic_exercise \\(1 g\\)", 
                  "aerobic exercise (N/A)",tdd) %>% 
      as.factor(),
    tdd.cat = apply(select(., cd.A:cd.zinc) , 1, function(x) {
      paste0(sub("cd\\.","",names(select(., cd.A:cd.zinc))[which(x != 0)]),
             " (", x[which(!is.na(x))], ")", collapse = ", ")
    }),
    tdd.cat = gsub("antidepressant \\(1 g\\)", "antidepressant (various)",tdd.cat),
    tdd.cat = ifelse(intervention == "z.placebo",
                     "placebo (N/A)", tdd.cat),
    tdd.cat = ifelse(intervention == "z.none",
                     "none (N/A)",tdd.cat),
    tdd.cat = gsub("aerobic_exercise \\(N/A\\)", 
                  "aerobic exercise (N/A)",tdd.cat)  
  ) %>% 
  select(-c(starts_with("cd."))) %>% 
  mutate(d.antidepressant = ifelse(z.antidepressant == 1,
                                   0, d.antidepressant)) %>%
  mutate(d.aerobic_exercise = ifelse(!is.na(d.aerobic_exercise),
                                     0, d.aerobic_exercise)) 
  

#### add vitamin data back to full data ####
# clean other parts of data
full_data <- data %>% 
  select(-starts_with("bias"), bias.overall) %>% 
  select(-c(baseline.vitamin, end.vitamin)) %>% 
  select(-c(notes)) %>% 
  select(-c(population, selection.criteria, number.treatment.arms)) %>% 
  select(-c(adverse.events, findings)) %>% 
  left_join(vitamins, by = c("study" = "id", "intervention","tdd"="dose")) %>% 
  mutate(age = as.numeric(age),
         male = as.numeric(male) / 100,
         total.N = as.numeric(total.N),
         assigned.N = as.numeric(assigned.N),
         baseline.N = as.numeric(baseline.N),
         final.N = as.numeric(final.N),
         duration = as.numeric(duration),
         blinded = as.factor(blinded),
         target = as.factor(target),
         psych.dx = as.factor(psych.dx) %>% fct_na_level_to_value(""),
         vit.def = as.factor(vit.def) %>% fct_na_level_to_value(""),
         report.vit.def = vit.def,
         calc.vit.def = as.factor(calc.vit.def),
         calc.vit.def.total = as.factor(calc.vit.def.total),
         pop.dx = as.factor(pop.dx) %>% fct_na_level_to_value(""),
         final.outcome.sd = ifelse(is.na(final.outcome.sd),
                       as.numeric(final.outcome.se) * sqrt(final.N),
                       as.numeric(final.outcome.sd)
                       ),
         gain.score.sd = ifelse(is.na(gain.score.sd) & flag == "gain scores",
                                as.numeric(gain.score.se) * sqrt(final.N),
                                as.numeric(gain.score.sd)
         ),
         final.outcome.mean = as.numeric(final.outcome.mean),
         baseline.outcome.mean = as.numeric(baseline.outcome.mean),
         basline.outcome.sd = as.numeric(baseline.outcome.sd),
         gain.score.mean = as.numeric(gain.score.mean),
         # dummy = 1L,
         group = intervention,
         orig.dose = dose,
         orig.tdd = tdd,
         # dose = tdd,
         frequency = as.factor(frequency)
         ) %>% 
  select(-c(intervention, tdd)) %>% 
  left_join(dose, by = c("study", "group", "orig.tdd")) %>%
  mutate(bias = factor(bias.overall, levels = c("low", "some concerns", "high"))) %>%
  filter(flag != "duplicate" & flag != "unusable") %>% 
  filter(!grepl("z.omega-3|z.calcium|z.multinutrient|z.A|z.folinate|z.glycine|z.phosphatidylserine", intervention))  %>% 
  group_by(study,scale) %>% 
  filter(n()>1) %>% 
  ungroup() %>% 
  mutate(
  # psych.dx.temp = case_when(
  #   scale == "Beck Depression Inventory" & baseline.outcome.mean >= 14 ~ "depression",
  #   scale == "von Zersen mood scale (ZMS)" & baseline.outcome.mean >= 27 ~ "depression",
  #   scale == "PHQ-9" & baseline.outcome.mean >= 10 ~ "depression",
  #   scale == "MADRS" & baseline.outcome.mean >= 7 ~ "depression",
  #   scale == "BDI-II" & baseline.outcome.mean >= 14 ~ "depression",
  #   scale == "HAM-D" & baseline.outcome.mean >= 8 ~ "depression",
  #   scale == "HADS-D" & baseline.outcome.mean >= 8 ~ "depression",
  #   scale == "DASS-21, depression" & baseline.outcome.mean >= 5 ~ "depression",
  #   scale == "DASS-21, anxiety" & baseline.outcome.mean >= 4 ~ "anxiety",
  #   scale == "CES-D" & baseline.outcome.mean >= 16 ~ "depression",
  #   scale == "HAM-A" & baseline.outcome.mean >= 8 ~ "anxiety",
  #   scale == "SIGH-SAD" & baseline.outcome.mean >= 20 ~ "depression",
  #   scale == "BAI" & baseline.outcome.mean >= 8 ~ "anxiety",
  #   scale == "Log BAI" & baseline.outcome.mean >= log(8) ~ "anxiety",
  #   scale == "EPDS" & baseline.outcome.mean >= 10 ~ "depression",
  #   scale == "HADS-A" & baseline.outcome.mean >= 8 ~ "anxiety",
  #   scale == "GDS" & baseline.outcome.mean >= 10 ~ "depression",
  #   scale == "STAI-S" & baseline.outcome.mean >= 40 ~ "anxiety",
  #   scale == "depression per HADS-D" & baseline.outcome.mean >= 8 ~ "depression",
  #   scale == "SCAARED" & baseline.outcome.mean >= 25 ~ "anxiety",
  #   scale == "MFQ" & baseline.outcome.mean >= 27 ~ "depression",
  #   scale == "GAD-7" & baseline.outcome.mean >= 10 ~ "anxiety",
  #   scale == "AUC BDI-II" & baseline.outcome.mean >= 14 ~ "depression",
  #   scale == "AUC MADRS" & baseline.outcome.mean >= 7 ~ "depression",
  #   scale == "HRS" & baseline.outcome.mean >= 18 ~ "anxiety",   # assumes HRS = HAM-A
  #   scale == "Yasavage and Brink Depression Scale" & baseline.outcome.mean >= 10 ~ "depression",
  #   scale == "BDI" & baseline.outcome.mean >= 10 ~ "depression",
  #   scale == "GDS-15" & baseline.outcome.mean >= 5 ~ "depression",
  #   scale == "MDI" & baseline.outcome.mean >= 20 ~ "depression",
  #   scale == "STAI-T" & baseline.outcome.mean >= 40 ~ "anxiety",
  #   scale == "SF-36 MH" & baseline.outcome.mean >= -42 ~ "depression",
  #   scale == "DASS-42, depression" & baseline.outcome.mean >= 10 ~ "depression",
  #   scale == "DASS-42, anxiety" & baseline.outcome.mean >= 8 ~ "anxiety",
  #   scale == "ΗDRS" & baseline.outcome.mean >= 8 ~ "depression",
  #   scale == "HSCL-25" & baseline.outcome.mean >= 1.75 ~ "depression",
  #   scale == "HADS-D >= 8" & baseline.outcome.mean >= 8 ~ "depression",
  #   scale == "HADS-A >=8" & baseline.outcome.mean >= 8 ~ "anxiety",
  #   TRUE ~ as.character(psych.dx)
  # ),
  psych.dx.temp = vit$psych_dx_scale(scale, baseline.outcome.mean),
  psych.dx.depression = ifelse(psych.dx.temp == "depression",1,0),
  psych.dx.anxiety = ifelse(psych.dx.temp == "anxiety",1,0)
  ) %>%
  group_by(study) %>% 
  mutate(psych.dx.depression = any(psych.dx.depression == 1),
         psych.dx.anxiety = any(psych.dx.anxiety == 1),
         psych.dx = case_when(psych.dx.depression & psych.dx.anxiety ~ "depression, anxiety",
                             psych.dx.depression ~ "depression",
                             psych.dx.anxiety ~ "anxiety",
                             TRUE ~ as.character(psych.dx))
         ) %>%
  select(-c(psych.dx.temp, psych.dx.depression, psych.dx.anxiety,
            orig.tdd, group, bias.overall)) %>%
  ungroup() %>% 
  group_by(study) %>% 
  mutate(vit.def = {
    parts <- c(vit.def, calc.vit.def.total) |>
      str_split(",\\s*") |> unlist()
    parts <- parts[!is.na(parts) & parts != ""] |> unique()
    if (length(parts) > 1) parts <- parts[parts != "none"]
    paste(parts, collapse = ", ")
  }) %>% 
  ungroup() %>% 
  mutate(vit.def = as.factor(vit.def) %>% fct_na_level_to_value("")) %>% 
  mutate(country = str_split(study, "\\/") %>% sapply(`[`,3) %>% as.factor(),
         year  = str_extract(study, "\\d{4}") %>% as.numeric(),
         author = str_extract(study, "[A-Z][a-z]+") %>% as.factor()
         )  %>% 
  mutate(country = fct_recode(country,
                              "USA" = "US",
                              "Tanzania" = "Tanzania-USA",
                              "Iran" = "Iran-USA",
                              "China" = "Australia-China",
                              "Australia" = "AUS",
                              "UK" = "Wales")) %>% 
  droplevels() 
# 3 warnings generated around age (some are NA), some male % not given,
# and some labels persisting in recode that were removed in original excel
dplyr::last_dplyr_warnings()



#### Create data for analyses ####
ranked_scales <- full_data %>% 
  count(scale, sort = TRUE) %>%       # count occurrences, most common first
  mutate(rank = row_number()) 
baseline.vit.names <-  colnames(vitamins)[startsWith(colnames(vitamins),"baseline.")]
#  Copy cutoffs
cutoffs_temp <- vit$cutoff_vitamins

#  Standardize cutoff names to match baseline style
names(cutoffs_temp) <- names(cutoffs_temp) |>
  str_replace("^vitamin\\s+", "") |>      # remove "vitamin "
  str_replace("cobalamin", "B12") |>      # harmonize synonyms
  str_replace("rbc\\.B9", "rbc.B9") |>    # keep as-is if already fine
  str_replace("selenium$", "selenium") |> # explicit (optional clarity)
  str_replace("^", "baseline.")           # add baseline prefix

dose_data <- full_data %>% 
  # mutate(intervention.class = case_when(z.citalopram == 1 ~ "SSRI",
  #                                       z.escitalopram == 1 ~ "SSRI",
  #                                       z.imipramine == 1 ~ "TCA",
  #                                       z.fluoxetine == 1 ~ "SSRI",
  #                                       z.nortriptyline == 1 ~ "TCA",
  #                                       z.amitriptyline == 1 ~ "TCA",
  #                                       z.SSRI == 1 ~ "SSRI",
  #                                       z.SSRI_TCA == 1 ~ "TCA",
  #                                       starts_with("z.B") == 1 ~ "B vitamin",
  #                                       z.D == 1 ~ "vitamin D/Calcium",
  #                                       z.calcium == 1 ~ "vitamin D/Calcium",
  #                                       .default = NA)) %>%
  # mutate(
  #   intervention = apply(select(., z.A:z.zinc) == 1, 1, function(x) {
  #     paste(names(select(., z.A:z.zinc))[x], collapse = ",")
  #   }),
  #   intervention = as.factor(intervention)
  # ) %>% 
  select(-c(`z.omega-3`, z.calcium, z.multinutrient, z.A, z.folinate, z.glycine, z.phosphatidylserine,
            `d.omega-3`, d.calcium, d.multinutrient, d.A, d.folinate, d.glycine, d.phosphatidylserine)) %>% 
  mutate(priority = case_when(
    intervention == "z.placebo" ~ 1,
    intervention == "z.antidepressant" ~ 2,
    intervention == "z.B12,z.B9,z.iron" ~ 3,
    intervention == "z.B9" ~ 4,
    intervention == "z.none" ~ 5,
    tdd.cat == "D (2.5e-05 g)" ~ 6,
    TRUE ~ 99  # default for all others
  ),
    final.outcome.mean = ifelse(is.na(final.outcome.mean) & !is.na(baseline.outcome.mean) & !is.na(gain.score.mean) & flag == "gain scores",
                                baseline.outcome.mean + gain.score.mean,
                                final.outcome.mean)
  # , final.outcome.sd = ifelse(is.na(final.outcome.mean) & !is.na(baseline.outcome.mean) & !is.na(gain.score.mean),
  #                           sqrt(gain.score.sd^2), # need other terms
  #                           final.outcome.sd
  # )
  ) %>% 
  group_by(study,target,scale) %>%
  reframe(
    idx = which.min(priority),
    contrast = (paste0(as.character(intervention), " - ",as.character(intervention)[idx])),
    treatment = (as.character(intervention)),
    control = (intervention[idx]),
    contrast.mean = final.outcome.mean[idx],
    contrast.gain = gain.score.mean[idx],
    contrast.sd = final.outcome.sd[idx],
    contrast.gain.sd = gain.score.sd[idx],
    contrast.N = final.N[idx],
    baseline.contrast.mean = baseline.outcome.mean[idx],
    baseline.contrast.sd = baseline.outcome.sd[idx],
    baseline.contrast.N = baseline.N[idx],
    is.comparator = 1:n() == idx,
    , pick(everything())# keeps all columns
    ) %>% 
  ungroup() %>% 
  mutate(
    pool.sd = sqrt(((final.N- 1) * final.outcome.sd^2 + (contrast.N-1) * contrast.sd^2)/(final.N + contrast.N - 2)),
    J = vit$J_fun(final.N + contrast.N - 2),
    baseline.J = vit$J_fun(baseline.N + baseline.contrast.N - 2),
    y = (final.outcome.mean - contrast.mean)/pool.sd * J ,
    sd.y = sqrt(y^2/(2 * (final.N + contrast.N - 2)) + J^2 * (final.N + contrast.N )/(final.N * contrast.N  )),
    baseline.y = (baseline.outcome.mean - baseline.contrast.mean)/pool.sd * baseline.J ,
    baseline.sd.y = sqrt(baseline.y^2/(2 * (baseline.N + baseline.contrast.N - 2)) + baseline.J^2 * (baseline.N + contrast.N )/(baseline.N * baseline.contrast.N  )),,
    pool.sd.gain = sqrt(((final.N- 1) * gain.score.sd^2 + (contrast.N-1) * contrast.gain.sd^2)/(final.N + contrast.N - 2)),
    y.gain = (gain.score.mean - contrast.gain)/pool.sd.gain  * J,
    sd.y.gain = sqrt(y.gain^2/(2 * (final.N + contrast.N - 2)) + J^2 * (final.N + contrast.N )/(final.N * contrast.N  ))
  ) %>% 
    mutate(
    # treatment = ifelse(treatment == control, NA_character_,
    #                         as.character(treatment)) %>% as.factor(),
         treatment = as.factor(treatment),
         control = as.factor(control),
         intervention = as.factor(intervention)
         ) %>%
  left_join(ranked_scales, by = "scale") %>%
  mutate(scale.preference = coalesce(rank, 99)) %>% 
  # mutate(depression.scale.preference = 
  #          case_when(
  #            scale == "Beck Depression Inventory" ~ 1,
  #            scale == "HAM-D" ~ 2,
  #            scale == "BDI-II" ~ 2
  #            scale == "AUC BDI-II" ~ 2,
  #            scale == "BAI" ~ 3,
  #            scale == "EPDS" ~ 4,
  #            scale == "PHQ-9" ~ 5,
  #            scale == "HADS-D" ~ 6,
  #            scale == "DASS-21, depression" ~ 7,
  #            scale == "CES-D" ~ 8,
  #            scale == "MADRS" ~ 9, 
  #            scale == "AUC MADRS" ~ 9,
  #            scale == "HADS-A" ~ 10,
  #            
  #            scale == "HAM-D" ~ 5,
  #            scale == "MADRS" ~ 6,
  #            scale == "Beck Depression Inventory II" ~ 7,
  #            scale == "AUC BDI-II" ~ 8,
  #            TRUE  ~ 99
  #          ) %>% as.numeric()) %>% 
  # mutate(anxiety.scale.preference = 
  #          case_when(
  #            scale == "BAI" ~ 1,
  #            scale == "HADS-A" ~ 2,
  #            scale == "STAI-S" ~ 3,
  #            scale == "SCAARED" ~ 5,
  #            scale == "log BAI" ~ 1,
  #            scale == "DASS-21, anxiety" ~ 4,
  #            TRUE  ~ 99
  #          ) %>% as.numeric()) %>% 
  # mutate(baseline.vit.of.z = ) %>% 
  mutate(y.gain = case_when(
    flag == "OLS" ~ vit$ols_to_g(b = gain.score.mean,
                             seb = gain.score.se,
                             sdy = gain.score.sd,
                             n1 = final.N,
                             n2 = contrast.N,
                             p = p)$g,
    flag == "LMM" ~ vit$lmer_to_g(b = gain.score.mean,
                              seb = gain.score.se,
                              sdy = gain.score.sd,
                              n1 = final.N,
                              n2 = contrast.N,
                              extra = p)$g,
    flag == "logistic" ~ vit$logistic_to_g(b = gain.score.mean,
                                       seb = gain.score.se,
                                       n1 = final.N,
                                       n2 = contrast.N,
                                       p = p)$g,
    flag == "hedge's g" ~ gain.score.mean,
    flag == "cohen's d" ~ vit$gfromd(gain.score.mean, final.N, contrast.N, p = 2),
    TRUE ~ y.gain
  ) %>% as.numeric()) %>% 
  mutate(
  sd.y.gain = case_when(
    flag == "OLS" ~ vit$ols_to_g(b = gain.score.mean,
                             seb = gain.score.se,
                             sdy = gain.score.sd,
                             n1 = final.N,
                             n2 = contrast.N,
                             p = p)$sd,
    flag == "logistic" ~ vit$logistic_to_g(b = gain.score.mean,
                                       seb = gain.score.se,
                                       n1 = final.N,
                                       n2 = contrast.N,
                                       p = p)$sd,
    flag == "hedge's g" ~ gain.score.sd,
    flag == "cohen's d" ~ vit$J_fun(final.N + contrast.N - 2) * gain.score.sd,
    TRUE ~ sd.y.gain
  ) %>% as.numeric()) %>% 
  mutate(treatment = ifelse(is.comparator,
                            NA_character_,
                            treatment %>% as.character()) %>% 
           as.factor(),
         contrast = ifelse(is.comparator,
                           NA_character_,
                           contrast %>% as.character) %>% 
           as.factor()
         ) %>% 
  group_by(study,target,scale) %>%
  reframe(
    y.gain = ifelse(1:n() == which.min(priority), NA_real_, y.gain),
    # sd.y.gain = ifelse(1:n() == which.min(priority), NA_real_, sd.y),
    y = ifelse(1:n() == which.min(priority), NA_real_, y),
    sd.y = ifelse(1:n() == which.min(priority), sd.y/sd.y * sqrt(1/contrast.N) * J, sd.y),
    baseline.y = ifelse(1:n() == which.min(priority), NA_real_, baseline.y),
    baseline.sd.y = ifelse(1:n() == which.min(priority), baseline.sd.y/baseline.sd.y * sqrt(1/contrast.N) * baseline.J, baseline.sd.y),
    sd.y.gain = ifelse(1:n() == which.min(priority), sd.y.gain/sd.y.gain * sqrt(1/contrast.N) * J, sd.y.gain),
    , pick(everything())# keeps all columns
  ) %>% 
  ungroup() %>% 
  filter(!(scale %in% c("MFQ","HSCL-25","von Zersen mood scale (ZMS)", "SF-36 MH", "SF-12 MCS"))) %>% 
  mutate(pct_dis = vit$log_normal_prob_dis(baseline.outcome.mean, baseline.outcome.sd, scale = scale, lower = FALSE)) %>%
  mutate(pct_dep = ifelse(target == "depression", pct_dis, NA_real_)) %>%
  mutate(pct_anx = ifelse(target == "anxiety", pct_dis, NA_real_)) %>% 
  # mutate(pct_dep = ifelse(grepl("depression", psych.dx)  & is.na(pct_dep) & target == "depression", 1.0, pct_dep)) %>%
  # mutate(pct_anx = ifelse(grepl("anxiety", psych.dx) & is.na(pct_anx) & target == "anxiety", 1.0, pct_anx)) %>%
  # mutate(pct_dep = ifelse(is.na(pct_dep) & target == "depression", 0.0, pct_dep)) %>%
  # mutate(pct_anx = ifelse(is.na(pct_anx) & target == "anxiety", 0.0, pct_anx)) %>%
  mutate(
  across(
    any_of(baseline.vit.names),
    ~ {
      mu <- .x
      # guardrails: need positive mean for log-normal moment parameterization
      ifelse(
        is.na(mu) | mu <= 0,
        NA_real_,
        vit$log_normal_prob_from_mom(
          mu,
          mu * 0.3,                 # SD = 0.3 * mean
          cutoff   = cutoffs_temp[[cur_column()]],
          lower    = TRUE
        )
      )
    },
    .names = "pdef_{.col}"
  ) ) %>%
  # rowwise() %>%
  mutate(
    pct_deficient = rowMeans(across(starts_with("pdef_baseline.")), na.rm = TRUE),
    n_measured    = rowSums(!is.na(across(starts_with("pdef_baseline."))))
  ) %>%
  # ungroup() %>% 
  mutate(
    pct_deficient = if_else(n_measured > 0, pct_deficient, NA_real_)) %>% 
  # mutate(pct_deficient = if_else(
  #   is.na(pct_deficient) & !is.na(vit.def), 1, pct_deficient
  # )) %>% 
  select(c(study, author, country, year, design, blinded,
           target, total.N, assigned.N, baseline.N, final.N,
           duration, scale, frequency, tdd, tdd.cat,
           age,
           male, 
           y, sd.y,
           treatment, control,
           contrast,
           y.gain, sd.y.gain,
           starts_with("baseline"),
           starts_with("end"),
           starts_with("z."),
           starts_with("d."),
           -baseline.J,
           bias,
           intervention,
           scale.preference,
           vit.def, calc.vit.def.total,
           psych.dx,
           flag, priority,
           pct_dep, pct_anx, pct_deficient
           ))

for(i in 1:nrow(dose_data)) {
  if(dose_data$control[i] != "z.placebo" && !is.na(dose_data$contrast[i])) {
    temp_control <- dose_data$control[i] %>% as.character()
    dose_data[[temp_control]][i] <- dose_data[[temp_control]][i]-1
  } else {
    dose_data$z.placebo[i] <- 0
  }
}

analysis_data <- dose_data %>% 
  filter(!(study %in% c(
    "Gugger/2019/Switzerland",
    "Narula/2017/Canada",
    "Nguyen/2009/Guatemala",
    "Venkatasubramanian/2013/India",
    "Mozaffari-Khosravi/2013/Iran"
  )))

#### get quality data ####
quality <- data %>% 
  select(c(study, starts_with("bias"))) %>% 
  distinct() %>%
  mutate(across(everything(), as.factor)) %>% 
  distinct()

full_data$bias <- quality$bias.overall[match(full_data$study, quality$study)]

# reset bias
highrisk <- quality %>% select(-c(bias.overall,study)) %>% `==`(.,"high") %>%
  rowSums(, na.rm = TRUE) %>% `>`(.,1)

somerisk <- quality %>% select(-c(bias.overall,study)) %>% `==`(.,"some concerns") %>%
  rowSums(, na.rm = TRUE) %>% `>`(.,3) |
  quality %>% select(-c(bias.overall,study)) %>% `==`(.,"high") %>%
  rowSums(, na.rm = TRUE) %>% `==`(.,1)

lowrisk <- quality %>% select(-c(bias.overall,study)) %>% `==`(.,"some concerns") %>%
  rowSums(, na.rm = TRUE) %>% `<=`(.,3) &
  quality %>% select(-c(bias.overall,study)) %>% `==`(.,"high") %>%
  rowSums(, na.rm = TRUE) %>% `==`(.,0)

quality2 <- quality %>%
  # mutate(bias.overall = case_when(
  #   highrisk ~ "high",
  #   somerisk & !highrisk ~ "some concerns",
  #   lowrisk ~ "low",
  #   TRUE ~ NA_character_
  # ) %>% factor(levels = c("low", "some concerns", "high"))
  # ) %>%
  filter(study %in% analysis_data$study) 


attr(quality2, "labels") <-
  c(labels[grep("Sequence",labels):grep("overall risk of bias",labels)])

attr(quality, "labels") <-
  c(labels[grep("Sequence",labels):grep("overall risk of bias",labels)])

analysis_data$bias <- quality2$bias.overall[match(analysis_data$study, quality2$study)]

dose_data$bias <- quality$bias.overall[match(dose_data$study, quality$study)]

study_rm <- quality$study[!(quality$study %in% analysis_data$study)]
if("Vyas/2023/US/2" %in% analysis_data$study) {
  study_rm <- study_rm[-which(study_rm == "Vyas/2023/US")]
}
print(study_rm)


study_rm_dose <- quality$study[!(quality$study %in% dose_data$study)]
if("Vyas/2023/US/2" %in% dose_data$study) {
  study_rm_dose <- study_rm_dose[-which(study_rm_dose == "Vyas/2023/US")]
}
print(study_rm_dose)

# check all studies have comparators
stopifnot("one or more studies don't have comparitors" = all(analysis_data$study %>% table() %>% `>`(1)))

# check all studies have comparators
stopifnot("one or more studies don't have comparitors" = all(dose_data$study %>% table() %>% `>`(1)))

#### Save Data ####
saveRDS(full_data, "data/vitamins_full.rds")
# saveRDS(analysis_data, "data/vitamins.rds")
saveRDS(dose_data, "data/vitamins.rds")
saveRDS(quality, "data/quality.rds")
