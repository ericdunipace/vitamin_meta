# This script will run all of the required analyses and generate figures and tables
# for the manuscript. It assumes that all required packages are installed and that
# the working directory is set to the root of the project.
# please see the readme file for the project structure if interested.
# to run this script in the background, you can do
# library(here); rstudioapi::jobRunScript(here("main.R"))

# load renv to manage package versions
library(renv)

# reload packages via renv
renv::restore()

# Load other required packages
suppressPackageStartupMessages({
library(here)
library(quarto)
library(dplyr)
library(rlang)
library(cli)
})

# Set working directory to the root of the project
setwd(here::here())

# update cmdstanr compilation options for faster fitting
# optimized to my computer, please modify as needed
optimize <- FALSE
if(optimize) {
  if (cmdstanr::cmdstan_make_local() != (cxx <- "CXXFLAGS+=-g -O3 -march=native -arch arm64 -ftemplate-depth-256 -fbracket-depth=1024 -Wno-deprecated-declarations -Wall -pipe -DSTAN_THREADS")) {
    cmdstanr::cmdstan_make_local(cpp_options = list(cxx), append = FALSE)
    cmdstanr::rebuild_cmdstan()
  }
}


run_function <- function() {
  # Set up logging to capture output and errors
  rlang::global_entrace(enable = TRUE, class = c("error", "warning", "message"))
  
  log_file <- file(here::here("job_output.txt"), open = "wt")
  sink(log_file, split = TRUE)      # 'split = TRUE' sends output to both file and console
  sink(log_file, type = "message") # Optional: Also capture warnings/errors
  
  # will run when function crashes or exists
  # clears file connections, sink, and prints errors etc
  on.exit({
    # rlang::last_messages() %>% print()
    rlang::last_warnings() %>% print()
    tryCatch(rlang::last_error() %>% print(), error = function(e) {NULL})
    sink(type = "message") # Stop capturing messages
    sink()                # Stop capturing output
    close(log_file)       # Close the log file connection
  }, add = TRUE)
  
  # clear out Stan files; can cause some difficulty if model code is changed
  # as BRMS won't always detect it; so better to remove them
  # if code and formula exactly the same, can leave the Stan files unchanged
  # because cause will speed up the code below by not needing to compile model
  # code again
  remove_stan_files <- TRUE # default is to remove Stan files; set to FALSE to keep them
  if(remove_stan_files) {
    for(i in c("depression","anxiety")) {
      stan_files <- c(
        list.files(here("outputs","saved_models",i, "stan"), pattern = "\\.rds$|\\.", full.names = TRUE),
        list.files(here("outputs","saved_models",i, "sensitivity","stan"), pattern = "\\.rds$|\\.", full.names = TRUE))
      file.remove(stan_files)
    }
  }
  
  # Run analyses and generate outputs
  cli::cli(cli::cli_h1("Starting analyses and output generation..."))
  rlang::with_interactive(
    {
      cli::cli(cli::cli_h2("Running main analysis"))
      source(here("R","analysis.R"))
      cli::cli(cli::cli_h2("Running sensitivity analysis"))
      source(here("R","sensitivity.R"))
      cli::cli(cli::cli_h2("Generating tables"))
      source(here("R","tables.R"))
      cli::cli(cli::cli_h2("Generating figures"))
      source(here("R","figures.R"))
      # cli::cli(cli::cli_h2("Generating sensitivity analysis markdown"))
      cli::cli(cli::cli_h2("Generating main document results text"))
      quarto::quarto_render(
        input =  here::here("outputs","documents","Results.qmd"),
        output_format = "all"
      )
    }, value = FALSE)
}

run_function()

