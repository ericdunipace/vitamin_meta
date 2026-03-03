library(dplyr)
library(ggplot2)
library(forcats)
library(tidyr)


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


#### Get Data ####
data_full <- vit$get_vitamin_data(outcome = NULL, simple_analysis = FALSE,
                                  include_full_bias = TRUE)
data <- vit$get_vitamin_data(outcome = NULL, simple_analysis = TRUE,
                             include_full_bias = TRUE)
qual <- readRDS("data/quality.rds") %>% 
  filter(study %in% data_full$study)

domain <- attr(qual, "labels") %>% stringr::str_to_title() 

rob <- qual %>% 
  pivot_longer(
    cols = -study,
    names_to = "domain",
    values_to = "judgment"
  ) %>% 
  mutate(domain = fct_recode(domain %>% as.factor(),
                        "Random sequence generation" = "bias.sequence",
                        "Allocation concealment" = "bias.concealment",
                        "Blinding of participants and personnel" = "bias.blinding.participants",
                        "Blinding of outcome assessment" = "bias.blinding.outcome",
                        "Incomplete outcome data" = "bias.incomplete.outcome",
                        "Selective outcome reporting" = "bias.selective.outcome",
                        "Other sources of bias" = "bias.other",
                        "Overall" = "bias.overall"
  )) %>%
  mutate(domain = fct_relevel(domain, 
                              "Random sequence generation",
                              "Allocation concealment"
                              ,"Blinding of participants and personnel",
                              "Blinding of outcome assessment",
                              "Incomplete outcome data",
                              "Selective outcome reporting",
                              "Other sources of bias",
                              "Overall"
  ) %>% fct_rev()) %>%
  mutate(judgment = fct_relevel(judgment %>% as.factor(),
                                "high",
                                "some concerns",
                                "low"
                                )) %>% 
  count(domain, judgment) %>%
  group_by(domain) %>%
  mutate(percent = n / sum(n) * 100)

rob_plot <-   rob %>% ungroup() %>% 
  mutate(domain = as.character(domain)) %>%
  # add_row(domain = "__gap__", judgment = "low", n = 0, percent = 0) %>%
  mutate(sep = case_when(
    domain == "Overall" ~ "under",
    TRUE ~ "over"
  )) %>%
  mutate(
    domain = fct_relevel(domain, 
                         "Random sequence generation",
                         "Allocation concealment"
                         ,"Blinding of participants and personnel",
                         "Blinding of outcome assessment",
                         "Incomplete outcome data",
                         "Selective outcome reporting",
                         "Other sources of bias",
                         # "__gap__",
                         "Overall") %>% fct_rev()
  ) %>% 
  mutate(judgment = fct_relevel(judgment %>% as.factor(),
                                 "high",
                                 "some concerns",
                                 "low"
  )) %>% 
 ggplot(aes(x = domain, y = percent, 
                        fill = judgment)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_fill_manual(values = #c("low" = "#2ca02c", "some concerns" = "#ff7f0e", "high" = "#d62728")
                    c(
                      "low" = "#3A8F6A",        # muted green
                      "some concerns" = "#D9903D",  # warm amber
                      "high" = "#C43C3C"       # softer red
                    )
                    ) +
  labs(
    x = "",
    y = "Percentage",
    fill = "Judgment",
    title = "Risk of Bias"
  ) +
  facet_grid(rows = vars(sep), scales = "free_y", space = "free_y") +
  vit$theme_vit() +
  theme(
    # axis.ticks.x = element_blank(),
    axis.text.y = element_text(margin = margin(r = -15)),
    panel.border = element_blank(),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.text.y.right = element_blank()
  ) 
  
  
  
  

pdf("outputs/pdfs/rob.pdf", width = 8, height = 5)
print(rob_plot)
dev.off()
jpeg("outputs/jpegs/rob.jpeg", width = 800, height = 500)
print(rob_plot)
dev.off()
