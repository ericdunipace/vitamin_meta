#### load libraries ####
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(multinma)
  library(brms)
  library(stringr)
  library(forcats)
  library(prismadiagramR)
  library(here)
  library(grDevices)
  library(magick)
  library(ggtext)
  library(patchwork)
  library(purrr)
  library(cli)
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

figure_save <- function(fig, outdir, outcome, filename,
                        height = 5, width = 7, res = 600) {
  
  # make sure base output directory exists
  if(!dir.exists(out_path <- here::here(outdir, outcome))) {
    dir.create(out_path, recursive = TRUE)
  }
  
  # create EPS files
  eps_path <- here::here(out_path, "eps", paste0(filename, ".eps"))
  if(!dir.exists(here::here(out_path, "eps"))) {
    dir.create(here::here(out_path, "eps"), recursive = TRUE)
  }
  
  postscript(
    file = eps_path,       # where to save
    width = width,         # width in inches
    height = height,       # height in inches
    onefile = FALSE,       # IMPORTANT for EPS
    horizontal = FALSE,    # portrait orientation
    paper = "special"      # avoids default page sizes
  )
  print(fig)
  dev.off()
  
  # create TIFF files
  tiff_path <- here::here(out_path, "tiff", paste0(filename, ".tiff"))
  if(!dir.exists(here::here(out_path, "tiff"))) {
    dir.create(here::here(out_path, "tiff"), recursive = TRUE)
  }
  tiff(
    filename = tiff_path, # where to save
    width = width * res,  # width in pixels
    height = height * res,# height in pixels
    res = res,            # resolution
    # compression = "lzw"   # compression type
  )
  print(fig)
  dev.off()
  
  # pdf files
  pdf_path <- here::here(out_path, "pdf", paste0(filename, ".pdf"))
  if(!dir.exists(here::here(out_path, "pdf"))) {
    dir.create(here::here(out_path, "pdf"), recursive = TRUE)
  }
  pdf(
    file = pdf_path,       # where to save
    width = width,         # width in inches
    height = height        # height in inches
  )
  print(fig)
  dev.off()
  
  # create jpeg files
  jpeg_path <- here::here(out_path, "jpeg", paste0(filename, ".jpeg"))
  if(!dir.exists(here::here(out_path, "jpeg"))) {
    dir.create(here::here(out_path, "jpeg"), recursive = TRUE)
  }
  jpeg_res <- res/2
  jpeg(
    filename = jpeg_path,       # where to save
    width = width * jpeg_res,   # width in pixels
    height = height * jpeg_res, # height in pixels
    res = jpeg_res              # resolution
  )
  print(fig)
  dev.off()
  
  return(invisible(NULL))
}

as_printable_gtable <- function(g) {
  structure(
    list(g = g),
    class = "printable_gtable"
  )
}

print.printable_gtable <- function(x, ...) {
  grid::grid.newpage()
  grid::grid.draw(x$g)
  invisible(x)
}

plot_cell_colored_network <- function(
    arms,
    study_col = ".study",
    trt_col   = ".trt",
    def_col   = "def",      # expected values like "yes"/"no" (or TRUE/FALSE)
    dep_col   = "depressed",      # expected values like "yes"/"no"
    palette = c(
      "Def=no, Dep=no"  = "#7a7a7a",
      "Def=yes, Dep=no" = "#1b9e77",
      "Def=no, Dep=yes" = "#7570b3",
      "Def=yes, Dep=yes"= "#d95f02",
      "Multiple"        = "#111111"
    ),
    edge_alpha = 0.9,
    edge_width_range = c(0.3, 2.2),
    node_size_range = c(3, 10),
    seed = 1
) {
  # --- normalize inputs / build cell label
  df <- arms %>%
    mutate(
      .study = .data[[study_col]],
      .trt   = .data[[trt_col]],
      .def   = as.character(.data[[def_col]]),
      .dep   = as.character(.data[[dep_col]])
    ) %>%
    mutate(
      .def = ifelse(.def %in% c("TRUE", "T", "1","yes"), "yes", .def),
      .def = ifelse(.def %in% c("FALSE", "F", "0","no"), "no", .def),
      .dep = ifelse(.dep %in% c("TRUE", "T", "1","yes"), "yes", .dep),
      .dep = ifelse(.dep %in% c("FALSE", "F", "0","no"), "no", .dep),
      cell = paste0("Def=", .def, ", Dep=", .dep)
    ) %>%
    filter(!is.na(.study), !is.na(.trt), !is.na(cell))
  
  # --- build edges within each study (all pairwise comparisons), per cell
  edges <- df %>%
    group_by(cell, .study) %>%
    summarise(trts = list(sort(unique(.trt))), .groups = "drop") %>%
    mutate(pairs = map(trts, ~ {
      if (length(.x) < 2) return(NULL)
      t(combn(.x, 2))
    })) %>%
    filter(!map_lgl(pairs, is.null)) %>%
    unnest(pairs) %>%
    transmute(
      cell,
      t1 = pairs[, 1] %>% as.character(),
      t2 = pairs[, 2] %>% as.character()
    ) %>%
    mutate(
      a = pmin(t1, t2),
      b = pmax(t1, t2)
    ) %>%
    select(cell, a, b) %>%
    count(cell, a, b, name = "k_studies")
  
  # --- collapse across cells to get "which cells support each edge?"
  edge_support <- edges %>%
    group_by(a, b) %>%
    summarise(
      cells = list(sort(unique(cell))),
      n_cells = n_distinct(cell),
      k_total = sum(k_studies),
      .groups = "drop"
    ) %>%
    mutate(
      edge_group = ifelse(n_cells == 1, map_chr(cells, 1), "Multiple")
    )
  
  # --- node weights (how often a treatment appears; optional but nice)
  node_counts <- df %>%
    count(.trt, name = "n_arms") %>%
    rename(name = .trt)
  
  # --- igraph object
  g <- igraph::graph_from_data_frame(
    d = edge_support %>% transmute(from = a, to = b, edge_group, k_total),
    directed = FALSE,
    vertices = node_counts
  )
  
  # --- fixed layout from full network
  set.seed(seed)
  lay <- create_layout(g, layout = "fr")
  
  # --- plot
  ggraph(lay) +
    geom_edge_link(
      aes(color = edge_group, width = k_total),
      alpha = edge_alpha,
      lineend = "round"
    ) +
    geom_node_point(
      aes(size = n_arms),
      shape = 21, fill = "white", color = "black", stroke = 0.6
    ) +
    geom_node_text(
      aes(label = name),
      repel = TRUE, size = 3.3
    ) +
    scale_edge_color_manual(values = palette, drop = FALSE) +
    scale_edge_width(range = edge_width_range) +
    scale_size(range = node_size_range) +
    guides(
      edge_width = guide_legend(title = "Total studies\nper comparison"),
      color      = guide_legend(title = "Supported by\npopulation cell"),
      size       = guide_legend(title = "Arms per\ntreatment")
    ) +
    theme_void() +
    theme(legend.position = "right")
}

net_plot_combined_def_dis <- function(dat, size = NULL, nudge = 0) {
  get_global_edge_max <- function(plots) {
    max(sapply(plots, function(p) {
      built <- ggplot_build(p)$data
      idx <- which(vapply(built, function(d) "edge_width" %in% names(d), logical(1)))[1]
      if (is.na(idx)) return(NA_real_)
      max(built[[idx]]$edge_width, na.rm = TRUE)
    }), na.rm = TRUE)
  }
  
  get_edge_nstudy_range_from_plot <- function(p) {
    gr <- attr(p$data, "graph")
    if (is.null(gr)) stop("No graph found in plot data attribute.")
    ed <- tidygraph::activate(gr, "edges") %>% as_tibble()
    range(ed$.nstudy, na.rm = TRUE)
  }
  
  get_global_node_max <- function(plots) {
    max(sapply(plots, function(p) {
      built <- ggplot_build(p)$data
      idx <- which(vapply(built, function(d) "size" %in% names(d), logical(1)))[1]
      if (is.na(idx)) return(NA_real_)
      max(built[[idx]]$size, na.rm = TRUE)
    }), na.rm = TRUE)
  }
  resize_node_labels <- function(p, size = 2.5) {
    # remove existing node text/label layers
    p$layers <- p$layers[!sapply(p$layers, function(l)
      inherits(l$geom, "GeomText") | inherits(l$geom, "GeomLabel")
    )]
    
    # add new label layer with desired size
    p + ggraph::geom_node_text(
      aes(label = name),
      size = size,
      hjust = "outward",
      vjust = "outward"
    )
  }
  
  get_edge_ranges <- function(x) {
    b <- ggplot_build(x)
    # find the layer containing edge_width
    idx <- which(vapply(b$data, function(d) "edge_width" %in% names(d), logical(1)))[1]
    range(b$data[[idx]]$edge_width, na.rm = TRUE)
  }
  
  get_node_ranges <- function(x) {
    b <- ggplot_build(x)
    # find the layer containing size
    idx <- which(vapply(b$data, function(d) "size" %in% names(d), logical(1)))[1]
    b$data[[idx]]$size
  }
  
  parallel_networks <- lapply(dat$def %>% levels(), function(def_level) {
    lapply(dat[[outcome_label]] %>% unique() %>% sort(), function(dep_level) {
      # browser()
      net <- tryCatch(vit$construct_nma_network(dat %>% 
                                  filter(def == def_level) %>% 
                                  filter(.data[[outcome_label]] == dep_level) 
      ), error = function(e) {
        cli::cli_alert_warning(glue::glue("Error constructing network for def={def_level}, dep={dep_level}: {e$message}"))
        return(NULL)
      })
      if(is.null(net)) return(NULL)
      
      
      net %>%
        plot(level = "treatment", weight_nodes = TRUE, nudge = nudge) +
        theme(legend.position = "none") +
        ggtitle(paste0("Deficiency: ", switch(def_level,
          "yes" = "yes",
          "no" = "no",
          "NA" = "unknown")
          , glue::glue(", {outcome %>% stringr::str_to_title()}: "), 
          switch(dep_level %>% as.character(),
          "1" = "yes",
          "0" = "no",
        NA))) +
        theme(
          plot.margin = margin(12, 50, 12, 50),
          legend.position = "bottom",
          text = element_text(family = "Helvetica", size = 8),
          plot.title = element_text(size = 8, face = "plain", hjust = 0.5,
                                    vjust = 2.5,
                                    margin = margin(t=6, b = 2)),
          legend.title = element_text(size = 8),
          legend.text  = element_text(size = 7),
          legend.key.height = unit(3.2, "mm"),
          legend.key.width  = unit(6, "mm"),
        ) + 
        guides(edge_colour = "none")
    }) 
  })
  
  plots_flat <- parallel_networks %>% purrr::list_flatten()
  if(any(sapply(plots_flat, is.null))) {
    cli::cli_alert_warning("Some networks failed to construct and will be skipped in the combined plot.")
    plots_flat <- plots_flat[!sapply(plots_flat, is.null)]
  }
  
  edge_range_study <- range(vapply(plots_flat, get_edge_nstudy_range_from_plot, numeric(2)), na.rm = TRUE)
  node_max <- max(sapply(plots_flat, \(p) max(p$data$.sample_size, na.rm = TRUE)), na.rm = TRUE)
  node_min <- max(sapply(plots_flat, \(p) min(p$data$.sample_size, na.rm = TRUE)), na.rm = TRUE)
  node_ranges <- sapply(plots_flat, get_node_ranges ) %>% unlist() %>% range()
  edge_ranges <- sapply(plots_flat, get_edge_ranges ) %>% unlist() %>% range()
  
  if(!is.null(size)) {
    plots_flat <- lapply(plots_flat, resize_node_labels, size = size)
  }
  
  pretty_plus_one <- function(x) {
    sort(unique(c(1, scales::breaks_pretty(3)(x))))
  }
  
  pretty_plus_nmm <- function(x) {
    sort(unique(c(node_min, scales::breaks_pretty(2)(x), node_max)))
  }
  
  plots_scaled <- lapply(plots_flat, function(p) {
    p +
      ggplot2::scale_size_continuous(
        "Total sample size",
        limits = c(1, node_max),
        range  = c(node_ranges[1], node_ranges[2]),             # match multinma max_size
        breaks = pretty_plus_nmm
      ) +
      ggraph::scale_edge_width_continuous(
        "Number of studies",
        limits = c(1, edge_range_study[2]),
        range  = c(.25, min(3,edge_ranges[2])),
        breaks = pretty_plus_one
      )
  })
  pw <- patchwork::wrap_plots(plots_scaled, nrow = 3, guides = "collect", widths = 0.5)  +
    plot_annotation(
      tag_levels = "A"  # A, B, C...
    ) &
    theme(legend.position = "bottom", legend.box = "horizontal",
          plot.tag = element_text(family = "Helvetica", size = 9, face = "bold"),
          plot.tag.position = c(-.01, 1.02)  # top-left inside each panel
    )
  
  pw
}

#### load data ####
# full data
data_full <- vit$get_vitamin_data(outcome = NULL, simple_analysis = FALSE,
                             include_full_bias = TRUE, additive_tx = TRUE)
data <- vit$get_vitamin_data(outcome = NULL, simple_analysis = TRUE,
                                    include_full_bias = TRUE)

#### set overall outdirs ####
if(!dir.exists(overall_outdir <- here::here("outputs","overall"))) {
  dir.create(overall_outdir, recursive = TRUE)
}
if(!dir.exists(overall_pdf <- here::here(overall_outdir,"pdf"))) {
  dir.create(overall_pdf, recursive = TRUE)
}
if(!dir.exists(overall_jpeg <- here::here(overall_outdir,"jpeg"))) {
  dir.create(overall_jpeg, recursive = TRUE)
}
if(!dir.exists(overall_eps <- here::here(overall_outdir,"eps"))) {
  dir.create(overall_eps, recursive = TRUE)
}
if(!dir.exists(overall_tiff <- here::here(overall_outdir,"tiff"))) {
  dir.create(overall_tiff, recursive = TRUE)
}

#### prisma flow diagram ####
{
  prisma_data <- list(
    N_main_analysis = vit$nstudies(data %>% filter(bias == "low"), "study"))
  prisma_data$N_sensitivity_only <- (vit$get_vitamin_data(
      outcome = NA_character_, 
      simple_analysis = FALSE,
      include_full_bias = TRUE,
      additive_tx = TRUE
    ) %>%  vit$nstudies("study")) - prisma_data$N_main_analysis
  prisma_data$N_initial_studies <- 5036L
  prisma_data$N_exclude = 1886L
  prisma_data$N_screened = prisma_data$N_initial_studies - prisma_data$N_exclude # 3156
  prisma_data$N_exclude_screen = 2849L
  prisma_data$N_full_text = prisma_data$N_screened - prisma_data$N_exclude_screen # 301
  prisma_data$N_exclude_full_text = 190L
  prisma_data$N_extracted = prisma_data$N_full_text - prisma_data$N_exclude_full_text # 111
  prisma_data$N_exclude_after_extract = prisma_data$N_extracted - prisma_data$N_sensitivity_only - prisma_data$N_main_analysis
  prisma_data$N_exclude_after_extract_duplicate <- 1
  prisma_data$N_exclude_after_extract_not_meeting_inclusion <- 3
  prisma_data$N_included = prisma_data$N_extracted - prisma_data$N_exclude_after_extract # 101
  saveRDS(prisma_data, here::here(overall_outdir, "prisma_data.rds"))
  
  
  # prismaFormat <-
  #   data.frame(
  #     prismaLvl = c(1,1,2,2,3,3,4,4,5,5,6
  #                   # ,6,7
  #                   ),
  #     nodeType =  c("Node", "Filter", "Node", "Filter", "Node", "Filter", "Node", "Filter",
  #                   # "Node", "Filter",
  #                   "Node", "Filter", "End"),
  #     prismaTxt = c(
  #       sprintf("Records identified from databases\n(n = %s)",prisma_data$N_initial_studies),
  #       sprintf("Duplicates excluded\n(n = %s)", prisma_data$N_exclude),
  #       sprintf("Records screened\n(n = %s)",prisma_data$N_screened),
  #       paste0("Records excluded as irrelevant\n(n = ", prisma_data$N_exclude_screen,")"),
  #       paste0("Full-text studies assessed for eligibility\n(n = ", prisma_data$N_full_text, ")"),
  #       paste0("Records excluded\n(n = ", prisma_data$N_exclude_full_text, ")"),
  #       paste0("Studies included for extraction\n(n = ", prisma_data$N_extracted, ")"),
  #       paste0("Studies excluded after extraction\nfor not meeting inclusion criteria,\nhaving duplicate data,\n or not having usable data\n(n = ", prisma_data$N_exclude_after_extract, ")"),
  #       # Kennedy/2011/UK multinutrient
  #       # Katta/2023/New Zealand multinutrient
  #       # Mech/2016/USA multinutrient intervention
  #       # Kaviani/2020/Iran  duplicate
  #       # Khajehnasiri/2015/Iran duplicate and excluded treatments
  #       # Calissendorff/2015/Sweden data not usable due to not reporting means in each group 
  #       # Muhihi/2022/Tanzania-USA not a usuable depression or anxiety scale
  #       paste0("Studies included in analyses\n(n = ", prisma_data$N_included, ")"),
  #       # paste0("Studies included in additive treatment sensitivity analysis\n(n = ", N - 1886- 2849 - 190-7, ")"),
  #       paste0("Studies only included\nin sensitivity analyses due to:\n• Only usable in additive\ntreatment models (n = ",1, ")","\n• Only having dose\ncomparisons (n = ", 4, ")"),
  #       # paste0("Studies included in dose sensitivity analysis\n(n = ", N - 1886- 2849  - 190-7-1, ")"),
  #       # paste0("Studies excluded due to only have dose comparison arms\n(n = ", 4, ")"),
  #       paste0("Studies included in main analysis\n(n = ", prisma_data$N_main, ")")
  #     )
  #   )
  # prismaFormat$fontSize <- 20
  # prismaFormat$width = 3
  # prismaFormat$height = 0.9
  # base_flow_plot <- prismadiagramR::getPrisma(studyStatus = NULL, prismaFormat = prismaFormat) %>% 
  #   DiagrammeR::grViz() %>% 
  #   DiagrammeRsvg::export_svg()  %>%
  #   charToRaw()
  
  base_flow_plot <- DiagrammeR::grViz(glue::glue("digraph prisma2020 {
    
    graph [
        layout  = dot,
        rankdir = TB,
        splines = ortho,
        nodesep = 0.6,
    ]
    
    # ---- Global defaults ----
    node [shape = box, fontname = 'Helvetica']
    
    # ---- Title bar ----
    title [
        label    = 'Identification of studies via databases',
        style    = filled,
        fillcolor= '#FFC000',
        fontsize = 20,
        width    = 6,
        height   = 0.5
    ]
    
    # ---- Phase bars on the left ----
    # node [
    #     shape    = box,
    #     style    = filled,
    #     fillcolor= '#D9E1F2',
    #     fontsize = 12,
    #     width    = 1.0,
    #     height   = 1.2
    # ]
    # 
    # id_phase        [label = 'Identification']
    # screening_phase     [label = 'Screening', height = 5]
    # included_phase  [label = 'Included']
    
    # ---- Main boxes (reset defaults) ----
    node [
        shape    = box,
        style    = filled,
        fillcolor= 'white',
        fontsize = 14,
        width    = 3.3,
        height   = 1.0
    ]
    
    # Row 1 (Identification)
    rec_id [
        label = 'Records identified from:\\nDatabases (n = {{prisma_data$N_initial_studies}})'
    ]
    
    rec_removed [
        label = 'Records removed before screening:\\nDuplicate records removed\\n(n = {{prisma_data$N_exclude}})'
    ]
    
    # Row 2 (Screening)
    rec_screened [
        label = 'Records screened\\n(n = {{prisma_data$N_screened}})'
    ]
    
    rec_excluded [
        label = 'Records excluded\\n(n = {{prisma_data$N_exclude_screen}})'
    ]
    
    # # Row 3 (Screening – retrieval)
    # rep_sought [
    #     label = 'Reports sought for retrieval\\n(n = )'
    # ]
    # 
    # rep_not_retr [
    #     label = 'Reports not retrieved\\n(n = )'
    # ]
    
    # Row 4 (Eligibility)
    rep_assessed [
        label = 'Reports assessed for eligibility\\n(n = {{prisma_data$N_full_text}})'
    ]
    
    rep_excluded [
        label = 'Reports excluded:\\n(n = {{prisma_data$N_exclude_full_text}})'
    ]
    
    # Row 5 (Extraction)
    rep_extracted [
        label = 'Studies included in extraction\\n(n ={{prisma_data$N_extracted}})'
    ]
    
    rep_exclude_from_extract [
        label = 'Reports exluded after extraction\\n(n ={{prisma_data$N_exclude_after_extract}})'
    ]
    
    # Row 6 (Included)
    studies_included [
        label = 'Studies included in analyses\\n(n = {{prisma_data$N_included}})\\nMain analysis\\n(n = {{prisma_data$N_main_analysis}} )\\nSensitivity analysis only\\n(n = {{prisma_data$N_sensitivity_only}})'
    ]
    
    # ---- Alignment (rows) ----
    { rank = same; title }
    
    { rank = same; rec_id rec_removed }
    
    { rank = same; rec_screened rec_excluded }
    
    { rank = same; rep_assessed rep_excluded }
    
    { rank = same; rep_extracted rep_exclude_from_extract }
    
    { rank = same; studies_included }
    
    # ---- Main flow (left column) ----
    rec_id       -> rec_screened
    # rec_screened -> rep_sought
    # rep_sought   -> rep_assessed
    rec_screened -> rep_assessed
    rep_assessed -> rep_extracted
    rep_extracted -> studies_included
    
    # ---- Right-hand boxes ----
    rec_id       -> rec_removed
    rec_screened -> rec_excluded
    # rep_sought   -> rep_not_retr
    rep_assessed -> rep_excluded
    rep_extracted -> rep_exclude_from_extract
    
    # ---- INVISIBLE edges to keep layout aligned ----
    edge [style = invis, arrowhead = none]
    # id_phase ->  screening_phase -> included_phase
    title -> rec_id
    title -> rec_removed      # helps center title
    # id_phase -> rec_id
    # screening_phase   -> rec_screened
    # screening_phase   -> rep_assessed
    # screening_phase-> rep_extracted
    # included_phase -> studies_included
}", .open = "{{", .close = "}}")) %>% 
    DiagrammeRsvg::export_svg()  %>%
    charToRaw()
  
  base_flow_plot %>% 
    rsvg::rsvg_pdf(here::here(overall_pdf,"prisma_flow_diagram.pdf"),
                   width  = 6.5 * 72,
                   height = 5  * 72)
  
  base_flow_plot %>% 
    rsvg::rsvg_ps(here::here(overall_eps,"prisma_flow_diagram.eps"),
                   width  = 6.5 * 72,
                   height = 5 * 72)
  
  base_flow_plot %>% 
    rsvg::rsvg_png(
                   width  = 6.5 * 300,
                   height = 5  * 300) %>% 
  magick::image_read() %>% 
  magick::image_write(
    path   = here::here(overall_jpeg,
                        "prisma_flow_diagram.jpeg"),
    format = "jpeg",
    compression = "none",  # uncompressed (best for journals)
    density = 300
  )
  
  magick::image_read(
    rsvg::rsvg_png(base_flow_plot, 
                   width  = 6.5 * 600,   # width (inches) * dpi
                   height = 5 * 600)   # height (inches) * dpi
  )  %>% 
    magick::image_write(
      path   = here::here(overall_tiff,
                          "prisma_flow_diagram.tiff"),
      format = "tiff",
      compression = "none",  # uncompressed (best for journals)
      density = 600
    )
}

#### plot ROB ####
{
    qual <- readRDS(here::here("data","quality.rds")) %>% 
    filter(study %in% data_full$study) %>% 
    mutate(study = forcats::fct_recode(study, "Vyas/2023/US" = "Vyas/2023/US/2"))
  
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
    ))  %>%
    mutate(judgment = fct_relevel(judgment %>% as.factor(),
                                  "high",
                                  "some concerns",
                                  "low"
    ) %>% fct_recode(
      "high" = "high",
      "unclear" = "some concerns",
      "low" = "low")
    ) %>%
    mutate(domain = fct_relevel(domain, 
                                "Random sequence generation",
                                "Allocation concealment"
                                ,"Blinding of participants and personnel",
                                "Blinding of outcome assessment",
                                "Incomplete outcome data",
                                "Selective outcome reporting",
                                "Other sources of bias",
                                "Overall"
    ) )
  
  rob_plot <-   rob %>%
    mutate(domain = domain %>% fct_rev()) %>% 
    count(domain, judgment) %>%
    group_by(domain) %>%
    mutate(percent = n / sum(n) * 100)%>% 
    ungroup() %>% 
    mutate(domain = as.character(domain)) %>%
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
                           "Overall") %>% fct_rev()
    ) %>% 
    mutate(judgment = fct_relevel(judgment %>% as.factor(),
                                  "high",
                                  "unclear",
                                  "low"
    )) %>% 
    ggplot(aes(x = domain, y = percent, 
               fill = judgment)) +
    geom_bar(stat = "identity", position = "stack") +
    coord_flip() +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    scale_fill_manual(values = #c("low" = "#2ca02c", "some concerns" = "#ff7f0e", "high" = "#d62728")
                        vit$bias_study_colors
    ) +
    labs(
      x = "",
      y = "Percentage",
      fill = "Judgment"
      # title = "Risk of Bias"
    ) +
    facet_grid(rows = vars(sep), scales = "free_y", space = "free_y") +
    vit$theme_vit(12) +
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
  
  rob2 <- rob %>% 
    mutate(
      domain = forcats::fct_relabel(
        domain,
        ~ stringr::str_wrap(.x, width = 12)
      )
    ) 
  
  
  full_traffic_light <- vit$traffic_light(rob2, size = 2) + 
    theme(strip.text.y.left = element_text(size = 6.5, lineheight = 0.95))
  
  half_levels <- rob2 %>%
    distinct(study) %>%
    arrange(study) %>%
    slice_head(n = ceiling(nrow(.) / 2)) %>%
    pull(study)
  
  first_half_tl <- vit$traffic_light(rob2 %>% 
    filter(study %in% half_levels),
    size = 5.5)
  
  second_half_tl <- vit$traffic_light(rob2 %>% 
                                       filter(!(study %in% half_levels)),
                                      size = 5.5
  )
  
  figure_save(rob_plot,
              outdir = overall_outdir,
              outcome = "",
              filename = "rob",
              height = 5,
              width = 8)
  
  figure_save(rob_plot,
              outdir = overall_outdir,
              outcome = "",
              filename = "rob",
              height = 5,
              width = 8)
  
  figure_save(full_traffic_light,
              outdir = overall_outdir,
              outcome = "",
              filename = "rob_full_traffic_light",
              height = 11,
              width = 8.5)
  
  figure_save(first_half_tl,
              outdir = overall_outdir,
              outcome = "",
              filename = "rob_first_half_traffic_light",
              height = 11,
              width = 8.5)
  
  figure_save(second_half_tl,
              outdir = overall_outdir,
              outcome = "",
              filename = "rob_second_half_traffic_light",
              height = 11,
              width = 8.5)
  
  # pdf(here::here(overall_outdir,"pdf","rob.pdf"), width = 8, height = 5)
  # print(rob_plot)
  # dev.off()
  # 
  # jpeg(here::here(overall_outdir,"jpeg","rob.jpeg"), width = 800, height = 500)
  # print(rob_plot + theme(text = element_text(size = 16)))
  # dev.off()
  # 
  # pdf(here::here(overall_outdir,"pdf","rob_full_traffic_light.pdf"), width = 8.5, height = 11)
  # print(full_traffic_light )
  # dev.off()
  # 
  # pdf(here::here(overall_outdir,"pdf","rob_first_half_traffic_light.pdf"), width = 8.5, height = 11)
  # print(first_half_tl )
  # dev.off()
  # 
  # pdf(here::here(overall_outdir,"pdf","rob_second_half_traffic_light.pdf"), width = 8.5, height = 11)
  # print(second_half_tl )
  # dev.off()
}

#### plot specific to outcomes ####

outcomes <- c("depression", "anxiety")
if (rlang::is_interactive()) {
  outcome <- "depression"
}

for(outcome in outcomes) {
  cli::cli_alert_info(glue::glue("Creating figures for {outcome}"))
  outcome_label <- switch(outcome,
                          depression = "depressed",
                          anxiety = "anxiety")
  
  resp_label    <- switch(outcome,
                          depression = "y",
                          anxiety = "yanx")
  
  resp <- switch(outcome,
                depression = "y",
                anxiety = "y_anx")

  # # create output directories
  # pdf_outdir <- here::here("outputs", outcome, "pdf")
  # jpeg_outdir <- here::here("outputs", outcome, "jpeg")
  # 
  # if(!dir.exists(pdf_outdir)) dir.create(pdf_outdir, recursive = TRUE)
  # if(!dir.exists(jpeg_outdir)) dir.create(jpeg_outdir, recursive = TRUE)
  
  dat <- data %>% filter(target == outcome)
  
  # load models
  {
    fit_low <- readRDS(
      here::here("outputs", "saved_models", outcome, 
                 glue::glue("model_", 
                            outcome, 
                            "_low_bias_data.rds"))
    ) 
    
    fit_some <- readRDS(
      here::here("outputs", "saved_models", outcome,
                 glue::glue("model_",
                            outcome,
                            "_some_and_low_bias_data.rds"))
    )

    fit_high <- readRDS(
      here::here("outputs", "saved_models", outcome,
                 glue::glue("model_",
                            outcome,
                            "_all_data.rds"))
    )
    
    # if(outcome == "depression") {
      fit_ssri_dep <- readRDS(here::here("outputs", "saved_models", "depression", 
                                         glue::glue("model_", 
                                                    "depression",
                                                    "_low_bias_ssri_data.rds")
      )
      )
    # }
    
    fit_some_only <- readRDS(
      here::here("outputs", "saved_models", outcome,
                 glue::glue("model_", 
                            outcome, 
                            "_some_bias_data.rds"))
    ) 
    
    fit_high_only <- readRDS(
      here::here("outputs", "saved_models", outcome,
                 glue::glue("model_",
                            outcome,
                            "_high_bias_data.rds"))
    )
    
    fit_def <- readRDS(
      here::here("outputs", "saved_models", outcome,
                 glue::glue("model_", 
                            outcome, 
                            "_all_data_by_deficiency_depression_and_bias.rds"))
    )
    
    fit_def_low <- readRDS(
      here::here("outputs", "saved_models", outcome,
                 glue::glue("model_", 
                            outcome, 
                            "_low_bias_data_by_deficiency_depression.rds"))
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
    
    fit_dose_low <- readRDS(
      here::here("outputs", "saved_models", outcome, "sensitivity",
                 glue::glue("model_", 
                            outcome, 
                            "_low_bias_data_by_dose.rds"))
    )
    
    fit_dose <- readRDS(
      here::here("outputs", "saved_models", outcome, "sensitivity",
                 glue::glue("model_", 
                            outcome, 
                            "_all_data_by_dose.rds"))
    )
    
    # bias_reg <- readRDS(here::here("outputs", 
    #                                "saved_models", 
    #                                outcome,  
    #                                "sensitivity",
    #                                "bias_regression_fit.rds"))
    bias_reg <- readRDS(here::here("outputs", "saved_models", "bias_regression_fit_overall.rds"))
    
    fit_low_nogain <- readRDS(
      here::here("outputs", 
                 "saved_models", 
                 outcome,  
                 "sensitivity",
                 glue::glue("model_",
                            outcome,
                            "_low_bias_data_nogain.rds"))
    )
    
    fit_some_nogain <- readRDS(
      file = here::here("outputs", 
                        "saved_models", 
                        outcome,  
                        "sensitivity",
                        glue::glue("model_",
                                   outcome,
                                   "_some_and_low_bias_data_nogain.rds"))
    )
    
    fit_high_nogain <- readRDS(
      file = here::here("outputs", 
                        "saved_models", 
                        outcome,  
                        "sensitivity",
                        glue::glue("model_",
                                   outcome,
                                   "_all_data_nogain.rds"))
    )
    
    fit_time <- readRDS(
      file = here::here("outputs", 
                        "saved_models", 
                        outcome,  
                        "sensitivity",
                        glue::glue("model_",
                                   outcome,
                                   "_all_data_by_time_cat.rds"))
    )
    
    fit_time_low <- readRDS(
      file = here::here("outputs", 
                        "saved_models", 
                        outcome,  
                        "sensitivity",
                        glue::glue("model_",
                                   outcome,
                                   "_low_bias_data_by_time_cat.rds"))
    )
    
    fit_pcts_approx <- readRDS(
      file = here::here("outputs", 
                        "saved_models", 
                        outcome,  
                        "sensitivity",
                                   glue::glue("model_", 
                                              outcome, 
                                              "_low_bias_data_by_pct_def_and_dep.rds"))
    )
    
  }
  
  #Simple NMA networks #
  {
    # full network
    netfull <- vit$construct_nma_network(dat) %>% 
      plot(level = "treatment", 
           weight_nodes = TRUE, nudge = 0.05) +
      # theme(legend.position = "bottom", legend.box = "vertical")
      theme(legend.position = "none")
    
    # some + low
    netsome <- vit$construct_nma_network(dat %>%
                                           filter(bias != "high") %>% 
                                           filter(control != "z.antidepressant")) %>% 
      
      plot(level = "treatment", weight_nodes = TRUE, nudge = 0.05) +
      # theme(legend.position = "bottom", legend.box = "vertical") +
      theme(legend.position = "none")
    
    # low bias only
    net <- vit$construct_nma_network(dat %>%
                                       filter(bias == "low")
                                     # %>% 
                                     #   filter(control != "z.antidepressant")
    ) %>% 
      plot(level = "treatment", weight_nodes = TRUE, nudge = 0.05) +
      theme(legend.position = "none")
    saveRDS(net, here::here("outputs",outcome, "nma_network_low_bias.rds"))
    
    # dose data network
    full_dose_net <- vit$construct_nma_network(data_full %>% filter(target == outcome)) %>%
      plot(level = "treatment", weight_nodes = TRUE,
           nudge = .05) +
      # theme(legend.position = "bottom", legend.box = "vertical")
      theme(legend.position = "none")
    
    # def and depression networks
    low_def_dis_net <- net_plot_combined_def_dis(dat %>% filter(bias == "low") %>% 
                                                   mutate(intervention = stringr::str_wrap(intervention %>% as.character(), width = 17)), nudge = 0.2, size = 2)
    
    
    figure_save(low_def_dis_net,
                outdir = "outputs",
                outcome = outcome,
                filename = "nma_def_and_dep_networks_low_bias",
                width = 7.5, height = 9.0)
    
    # save full network plots
    figure_save(netfull,
                outdir = "outputs",
                outcome = outcome,
                filename = "nma_network_full",
                height = 5.5,
                width = 8.5)
    
    figure_save(net,
                outdir = "outputs",
                outcome = outcome,
                filename = "nma_network",
                height = 4,
                width = 6)
    
    figure_save(full_dose_net,
                outdir = "outputs",
                outcome = outcome,
                filename = "nma_dose_network",
                height = 5.5,
                width = 8.5)
    
    # pdf(here::here(pdf_outdir, "nma_network_full.pdf"), 
    #     width = 7, height = 5)
    # print(netfull)
    # dev.off()
    # 
    # jpeg(here::here(jpeg_outdir, 
    #                 "nma_network_full.jpeg"), 
    #      width = 500, height = 500)
    # print(netfull)
    # dev.off()
    
    # save network plots
    # pdf(here::here(pdf_outdir,
    #                "nma_network.pdf"), width = 7, height = 5)
    # print(net)
    # dev.off()
    # 
    # jpeg(here::here(jpeg_outdir, "nma_network.jpeg"), 
    #      width = 500, height = 500)
    # print(net)
    # dev.off()
    
    # pdf(here::here(pdf_outdir,
    #                "nma_dose_network.pdf"), width = 7, height = 5)
    # print(full_dose_net)
    # dev.off()
    # 
    # jpeg(here::here(jpeg_outdir, "nma_dose_network.jpeg"), 
    #      width = 500, height = 500)
    # print(full_dose_net)
    # dev.off()
  }
  
  # get summary objects
  {
    low_sum  <- fit_low %>% vit$summary_brms_nma(keep = ".trt") %>% mutate(bias = "low")
    some_sum <- fit_some %>% vit$summary_brms_nma(keep = ".trt") %>% mutate(bias = "low/unclear")
    high_sum <- fit_high %>% vit$summary_brms_nma(keep = ".trt") %>% mutate(bias = "low/unclear/high")
    some_sum_only <- fit_some_only %>% vit$summary_brms_nma(keep = ".trt") %>% mutate(bias = "unclear")
    high_sum_only <- fit_high_only %>% vit$summary_brms_nma(keep = ".trt") %>% mutate(bias = "high")
    all_sum  <- bind_rows(low_sum, 
                          some_sum, 
                          high_sum)
    
    def_low_sum <- fit_def_low %>% 
      vit$summary_brms_nma(keep = c(".trt","def",outcome_label)) 
    
    
    bias_summary <- vit$summary_brms_nma(bias_reg, 
                         keep = c(".trt","def",outcome_label,"bias.overall",".idx_anx",".idx_dep"),
                         index = ".obs",
                         resp = resp_label) %>% 
      filter(.observed == TRUE)  %>% 
      select(-c(.idx_anx, .idx_dep)) %>% 
      mutate(!!outcome_label := factor(.data[[outcome_label]], levels = 1:0, labels = c("yes", "no"))) %>%
      mutate(!!outcome_label := forcats::fct_relevel(.data[[outcome_label]], "yes","no")) %>% 
      mutate(def = forcats::fct_recode(as.factor(def), "unknown" = "NA")) %>%
      mutate(def = forcats::fct_relevel(def, "unknown", "no","yes")) %>%
      rename(deficiency = def,
             bias = bias.overall)  %>% 
      mutate(bias = forcats::fct_recode(bias,
                                        "unclear" = "some concerns")) %>% 
      mutate(bias = forcats::fct_relevel(bias, "low", "unclear", "high"))
  }
  
  # plot relative effects
  {
    low_bias_plot <- fit_low %>% 
      vit$summary_brms_nma(keep = ".trt") %>% 
      filter(.observed) %>% select(-.observed ) %>%
      plot()
    low_bias_all_est_plot <- all_sum %>% 
      filter(.trt %in% low_sum$.trt) %>% 
      filter(.observed) %>% select(-.observed ) %>%
      plot(color = "bias") +
      scale_color_manual(values = vit$bias_colors,
                         name = "Risk of bias") +
      scale_y_discrete(expand = expansion(mult = c(0.02, 0.02))) +
      coord_cartesian(clip = "off") +
      theme(
        panel.spacing.y = unit(10, "pt"),
        strip.placement = "outside",
        strip.clip = "off"
      )
    combined_plot <-  all_sum %>%
      filter(.observed) %>% select(-.observed ) %>%
      plot(color = "bias") + 
      scale_color_manual(values = vit$bias_colors,
                         name = "Risk of bias") +
      scale_y_discrete(expand = expansion(mult = c(0.02, 0.02))) +
      coord_cartesian(clip = "off") +
      theme(
        panel.spacing.y = unit(10, "pt"),
        strip.placement = "outside",
        strip.clip = "off"
      )
    
    # if(outcome == "depression") {
     fit_ssri_dep %>% 
        vit$summary_brms_nma(keep = ".trt") %>%
        filter(.observed) %>% select(-.observed) %>% 
        plot() %>% 
        figure_save(outdir = "outputs",
                    outcome = outcome,
                    filename = "relative_effects_low_bias_ssri_only",
                    height = 5,
                    width = 7)
    # }
    
    figure_save(low_bias_plot,
                outdir = "outputs",
                outcome = outcome,
                filename = "relative_effects_low_bias",
                height = 5,
                width = 7)
    
    figure_save(low_bias_all_est_plot,
                outdir = "outputs",
                outcome = outcome,
                filename = "relative_effects_all_estimates_for_which_low_bias_exists",
                height = 5,
                width = 7)
    
    figure_save(combined_plot,
                outdir = "outputs",
                outcome = outcome,
                filename = "relative_effects_all",
                height = 11,
                width = 7)
    
    # pdf(here::here(pdf_outdir,
    #                "relative_effects_low_bias.pdf"), width = 7, height = 5)
    # print(low_bias_plot)
    # dev.off()
    # 
    # jpeg(here::here(jpeg_outdir, "relative_effects_low_bias.jpeg"), 
    #      width = 500, height = 500)
    # print(low_bias_plot)
    # dev.off()
    
    # pdf(here::here(pdf_outdir,
    #                "relative_effects_all_estimates_for_which_low_bias_exists.pdf"), width = 7, height = 5)
    # print(low_bias_plot)
    # dev.off()
    # 
    # jpeg(here::here(jpeg_outdir, "relative_effects_all_estimates_for_which_low_bias_exists.jpeg"), 
    #      width = 500, height = 500)
    # print(low_bias_plot)
    # dev.off()
    
    # pdf(here::here(pdf_outdir,
    #                "relative_effects_all.pdf"), width = 7, height = 11)
    # print(combined_plot ) 
    # dev.off()
  }
  
  # vitamin deficiency and disease status plots
  {
    
    def_low_sum %>% mutate(bias = "low") %>% mutate(deficiency = def) %>% filter(.observed) %>% vit$def_plot() %>% 
      figure_save(
        outdir = "outputs",
        outcome = outcome,
        filename = "deficiency_summary_low_bias",
        height = 5,
        width = 7
      )
    
    summary_def <- fit_def %>% 
      vit$summary_brms_nma(keep = c(".trt","def",outcome_label,"bias")) %>% 
      mutate(!!outcome_label := if(!is.factor(.data[[outcome_label]]) && is.numeric(.data[[outcome_label]])) {factor(.data[[outcome_label]], levels = 1:0, labels = c("yes", "no"))} else {.data[[outcome_label]]}) %>%
      mutate(!!outcome_label := forcats::fct_relevel(.data[[outcome_label]], "yes","no")) %>% 
      mutate(def = forcats::fct_recode(as.factor(def), "unknown" = "NA")) %>%
      mutate(def = forcats::fct_relevel(def, "unknown", "no","yes")) %>%
      rename(deficiency = def)  %>% 
      mutate(bias = forcats::fct_recode(bias,
                                        "unclear" = "some concerns")) %>% 
      mutate(bias = forcats::fct_relevel(bias, "low", "unclear", "high"))
    
    
    def_plot_D <- summary_def %>%
      filter(.trt == "D") %>% 
      vit$def_plot() +
      # geom_vline(xintercept = -0.5) +
      scale_x_continuous(breaks = seq(-1,.5, by = 0.5) %>% round(digits = 2),
                         limits = c(-1,.5),
                         labels = c(-1.0,-0.5,0.0,""),
                         expand = c(0,0)) +
      theme(strip.text.y.right = element_blank(),
            axis.text.y.left = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank()) 
    
     figure_save(def_plot_D,
                 outdir = "outputs",
                 outcome = outcome,
                 filename = "deficiency_plots_vitamin_D",
                 height = 5,
                 width = 7.5)
     
     sum_def_d_low <- summary_def %>%
       filter(.trt == "D",
              bias == "low")
     
     sum_def_d_low <- def_low_sum %>% 
       filter(.trt == "D") %>% 
       filter(.observed) %>% select(-.observed) %>%
       mutate(!!outcome_label := if(!is.factor(.data[[outcome_label]]) && is.numeric(.data[[outcome_label]])) {factor(.data[[outcome_label]], levels = 1:0, labels = c("yes", "no"))} else {.data[[outcome_label]]}) %>%
       mutate(!!outcome_label := forcats::fct_relevel(.data[[outcome_label]], "yes","no")) %>% 
       mutate(def = forcats::fct_recode(as.factor(def), "unknown" = "NA")) %>%
       mutate(def = forcats::fct_relevel(def, "unknown", "no","yes")) %>%
       rename(deficiency = def) 
     
     main_plot_data <- low_sum %>% mutate(
       deficiency = "all",
       !!outcome_label := ""
     ) %>%
       bind_rows(sum_def_d_low) %>%
       select(-bias) %>% 
       mutate(
         y_key = interaction(.trt, !!outcome_label, drop = TRUE),
         y_lab = if_else(.trt == "D", as.character(.data[[outcome_label]]), "")
       ) %>% 
       mutate(y_lab = ifelse(y_lab == "yes", glue::glue("{outcome_label}: yes"), y_lab)) %>% 
       mutate(y_lab = forcats::fct_relevel(y_lab, 
                                           "",
                                           "no",
                                      glue::glue("{outcome_label}: yes"))) %>% 
       mutate(deficiency = forcats::fct_relevel(deficiency, 
                                                "all",
                                                "unknown",
                                                "no",
                                                "yes"))
     
     xmin <- if(outcome_label == "anxiety") {
       -1.5
     } else {
       -2.5
     }
     xmax <- 1
     sx <- if(outcome_label == "anxiety") {
       scale_x_continuous(breaks = seq(xmin,xmax, by = 0.5) %>% round(digits = 2),
                          labels = c("",format(seq(xmin+0.5,0.5, by = 0.5), digits = 2),
                                     ""))
     } else {
       scale_x_continuous(breaks = seq(xmin,xmax, by = 0.5) %>% round(digits = 2),
                          labels = c("",format(seq(xmin + 0.5,1, by = 0.5), digits = 2)))
     }
     
     main_plot <- main_plot_data %>% 
       plot(color = "deficiency", 
            position = position_dodge(width = .75),
            add_left_wall = FALSE) +
       sx +
       theme(
         panel.grid.minor = element_blank()
       ) + 
       # scale_y_discrete(
       #   labels = setNames(main_plot$y_lab, main_plot$y_key)
       # ) + 
       ggplot2::scale_color_manual(values = c(
         "yes" = "#6a51a3",  # purple
         "no"  = "#5383c7",  # blue
         "unknown"  = "grey63",  # grey,
         "all" = "black"
         ),
         breaks = c("yes", "no", "unknown","all"),
         labels = c("yes" = "Deficient",
                    "no" = "Not Deficient",
                    "unknown" = "Unknown",
                    "all" = "Overall Estimate"),
         name = "Vitamin D Deficiency\nStatus"
       ) + aes(y = y_lab) + ylab(NULL) 
       
     txt_lbl <- if(outcome == "anxiety") {
       "Anxious"
     } else {
       "Depressed"
     }
     
     main_plot_final <- main_plot + 
       geom_text(
         data = data.frame(
           .trt = factor(c("D","D"), levels = levels(main_plot@data$.trt)),
           x = if(outcome == "depression") {
             rep(-.9,2)
             } else {
               rep(1.25,2)
             },              # choose x inside your plot range
           y = c(3.2,2.15),             # vertical position
           label = c(glue::glue("{txt_lbl}: Yes"), "No")
         ) %>% mutate(
           !!outcome_label := c(glue::glue("{txt_lbl %>% stringr::str_to_lower()}: yes"), "no")
         ),
         aes(x, y, label = label),
         inherit.aes = FALSE,
         size = 2.5, hjust = 1
       ) +  scale_y_discrete(position = "left") +
       theme(
         axis.text.y  = element_blank(),       # remove tick labels
         axis.ticks.y = element_blank()        # remove tick marks
         # panel.border   = ggplot2::element_rect(fill = NA, color = "grey88", linewidth = 0.25),
       ) 
       # labs(x = "Favors treatment <span style='font-size:20pt;'>\u2190</span>  <span style='font-size:20pt;'>\u2192</span> Favors placebo") +
       # theme(axis.title.x = ggtext::element_markdown())
     gb <- ggplot_build(main_plot_final)
     xm <- gb$layout$panel_params[[1]]$x.range
     th <- gb$plot$theme
     text_pt <- th$text$size %||% theme_get()$text$size
     text_size <- (text_pt - 2) / ggplot2::.pt
     # axis_pt <- th$axis.title.x$size %||% th$text$size %||% theme_get()$text$size
     # axis_size <- axis_pt / ggplot2::.pt
     
     arrow_length <- switch(outcome,
                            anxiety = 0.5,
                            depression = 0.75)
     arrow_min    <- switch(outcome,
                            anxiety = 0.25,
                            depression = 0.5)
     arrows <- ggplot() +
       # arrows starting at x=0 in *data* coordinates
       annotate("segment", x = -arrow_min, xend = -arrow_min-arrow_length, y = 0, yend = 0,
                arrow = arrow(type = "closed", length = unit(2.5, "mm"))) +
       annotate("segment", x = arrow_min, xend = arrow_min + arrow_length, y = 0, yend = 0,
                arrow = arrow(type = "closed", length = unit(2.5, "mm"))) +
       annotate("text", x = -(arrow_min + 0.5 *arrow_length), y = -.5, vjust = -.1,
                label = "Favors\ntreatment",
                size = text_size) +
       annotate("text", x =  (arrow_min + 0.5 *arrow_length), y = -.5, vjust = -.1,
                label = "Favors\nplacebo",
                size = text_size) +
       scale_x_continuous(limits = xm) +
       coord_cartesian(clip = "off") +
       theme_void() +
       theme(plot.margin = margin(t = 2, r = 5.5, b = 0, l = 5.5))
     main_arrow <- (main_plot_final + xlab(NULL) + theme(legend.position = "none") ) / (arrows) + plot_layout(heights = c(1, 0.1))
     
     
     if (outcome == "depression") {
       main_leg <- cowplot::get_legend(main_plot)
       leg_plot <- cowplot::ggdraw(main_leg)
       
       saveRDS(leg_plot, here::here("outputs", outcome, "main_plot_legend_only.rds"))
       figure_save(leg_plot,
                   outdir = "outputs",
                   outcome = outcome,
                   filename = "main_plot_legend_only",
                   height = 1,
                   width = 6.5)
     }
     
     saveRDS(main_arrow, here::here("outputs", outcome, "main_arrow.rds"))
     
     
     figure_save(main_plot_final,
                 outdir = "outputs",
                 outcome = outcome,
                 filename = "main_plot_low_bias_w_legend",
                 height = 5,
                 width = 7.5)
     
     figure_save(main_arrow,
                 outdir = "outputs",
                 outcome = outcome,
                 filename = "main_plot_low_bias_no_legend",
                 height = 6,
                 width = 6.5/2)
     
     
  }
  
  # cumulative rank plots
  # {
  #   low_rank_plot  <- fit_low %>% 
  #     vit$rank_probs_vit(cumulative = TRUE, sucra = TRUE) %>%
  #              plot()
  #   some_rank_plot <- fit_some %>% vit$rank_probs_vit(cumulative = TRUE) %>% plot()
  #   high_rank_plot <- fit_high %>% vit$rank_probs_vit(cumulative = TRUE) %>% plot()
  #   
  #   figure_save(low_rank_plot + 
  #               scale_x_continuous(breaks = seq(1, vit$ntreatments(fit_low), by = 2), name = "Rank" ),
  #               outdir = "outputs",
  #               outcome = outcome,
  #               filename = "cumulative_ranks_low_bias",
  #               height = 4.5,
  #               width = 6.5)
  #   
  #   figure_save(some_rank_plot+ 
  #                 scale_x_continuous(breaks = seq(1, vit$ntreatments(fit_some), by = 2), name = "Rank" ),
  #               outdir = "outputs",
  #               outcome = outcome,
  #               filename = "cumulative_ranks_some_and_low_bias",
  #               height = 4.5,
  #               width = 6.5)
  #   
  #   figure_save(high_rank_plot + 
  #                 scale_x_continuous(breaks = seq(1, vit$ntreatments(fit_high), by = 5), name = "Rank" ),
  #               outdir = "outputs",
  #               outcome = outcome,
  #               filename = "cumulative_ranks_all_bias",
  #               height = 4.5,
  #               width = 6.5)
  #   
  #   # pdf(here::here(pdf_outdir,
  #   #                "cumulative_ranks_low_bias.pdf"), width = 7, height = 5)
  #   # print(low_rank_plot + 
  #   #         scale_x_continuous(breaks = seq(1, vit$ntreatments(fit_low), by = 2), name = "Rank" ))
  #   # dev.off()
  #   # 
  #   # jpeg(here::here(jpeg_outdir, "cumulative_ranks_low_bias.jpeg"), 
  #   #      width = 500, height = 500)
  #   # print(low_rank_plot + 
  #   #         scale_x_continuous(breaks = seq(1, vit$ntreatments(fit_low), by = 2), name = "Rank" ))
  #   # dev.off()
  #   
  #   # pdf(here::here(pdf_outdir,
  #   #                "cumulative_ranks_some_and_low_bias.pdf"), width = 7, height = 5)
  #   # print(some_rank_plot+ 
  #   #         scale_x_continuous(breaks = seq(1, vit$ntreatments(fit_some), by = 2), name = "Rank" ))
  #   # dev.off()
  #   # 
  #   # jpeg(here::here(jpeg_outdir, "cumulative_ranks_some_and_low_bias.jpeg"), 
  #   #      width = 500, height = 500)
  #   # print(some_rank_plot+ 
  #   #         scale_x_continuous(breaks = seq(1, vit$ntreatments(fit_some), by = 2), name = "Rank" ))
  #   # dev.off()
  #   # 
  #   # pdf(here::here(pdf_outdir,
  #   #                "cumulative_ranks_all_bias.pdf"), width = 7, height = 5)
  #   # print(high_rank_plot + 
  #   #         scale_x_continuous(breaks = seq(1, vit$ntreatments(fit_high), by = 5), name = "Rank" ))
  #   # dev.off()
  #   # 
  #   # jpeg(here::here(jpeg_outdir, "cumulative_ranks_all_bias.jpeg"), 
  #   #      width = 500, height = 500)
  #   # print(high_rank_plot)
  #   # dev.off()
  # }
  
  # setup funnel plots
  {
    demeaned_data <- fit_low$data %>%
      mutate(E = colMeans(posterior_epred(fit_low,
                                          resp = resp))) %>%
      mutate(y.dm = y - E) %>% 
      mutate(intervention = as.character(fit_low$prep$data %>% filter(!is.na(.data[["y"]])) %>% pull(intervention)) %>%
               stringr::str_to_title() %>%
               factor())
    
    demeaned_data_high <- fit_high$data %>%
      mutate(E = colMeans(posterior_epred(fit_high,
                                          resp = resp))) %>%
      mutate(y.dm = y - E) %>% 
      mutate(intervention = as.character(fit_high$prep$data %>% filter(!is.na(.data[["y"]])) %>% pull(intervention)) %>%
               stringr::str_to_title() %>%
               factor())
    
    demeaned_data_def <- fit_def$data %>%
      mutate(E = colMeans(posterior_epred(fit_def,
                                          resp = resp))) %>%
      mutate(y.dm = y - E) %>% 
      mutate(!!outcome_label := factor(.data[[outcome_label]], levels = 1:0, labels = c("yes", "no"))) %>%
      mutate(def = forcats::fct_recode(as.factor(def), "unknown" = "NA")) %>%
      mutate(intervention = as.character(fit_def$prep$data %>% filter(!is.na(.data[["y"]])) %>% pull(intervention)) %>%
               stringr::str_to_title() %>%
               factor())
    
    demeaned_data_def_low <- fit_def_low$data %>%
      mutate(E = colMeans(posterior_epred(fit_def_low,
                                          resp = resp))) %>%
      mutate(y.dm = y - E) %>% 
      mutate(!!outcome_label := factor(.data[[outcome_label]], levels = 1:0, labels = c("yes", "no"))) %>%
      mutate(def = forcats::fct_recode(as.factor(def), "unknown" = "NA")) %>%
      mutate(intervention = as.character(fit_def_low$prep$data %>% filter(!is.na(.data[["y"]])) %>% pull(intervention)) %>%
               stringr::str_to_title() %>%
               factor())
    
    
    bias_demeaned_data <- bias_reg$data %>%
      filter(.data[[resp]] %>% is.na() %>% `!`()) %>%
      bind_cols(E = colMeans(posterior_epred(bias_reg, 
                                             resp = resp_label))) %>%
      mutate(y.dm = y - E) %>% 
      mutate(!!outcome_label := factor(.data[[outcome_label]], levels = 1:0, labels = c("yes", "no"))) %>%
      mutate(def = forcats::fct_recode(as.factor(def), 
                                       "unknown" = "NA")) %>%
      mutate(intervention = as.character(bias_reg$prep$data %>% filter(!is.na(.data[[resp]])) %>% pull(intervention)) %>%
               stringr::str_to_title() %>%
               factor()) %>%
      filter(!is.na(y))
    
    xlab_funnel <- latex2exp::TeX("SMD - \\widehat{\\mu}")
    pval <- 1-0.05/2
    p01<- "grey95"
    p05 <- "grey85"
    p10 <- "grey70"
    
    se_max <- 4  
    se_seq <- seq(0, se_max, length.out = se_max * 1000)
    
    contours <- tibble(
      se_lo = head(se_seq, -1),
      se_hi = tail(se_seq, -1)
    ) %>%
      mutate(
        x10 = qnorm(0.95)  * se_hi,   # two-sided p=.10
        x05 = qnorm(0.975) * se_hi,   # two-sided p=.05
        x01 = qnorm(0.995) * se_hi    # two-sided p=.01
      )
    
    contour_legend <- tibble(
      region = factor(
        c("p > 0.10",
          "0.05 < p ≤ 0.10",
          "0.01 < p ≤ 0.05",
          "p ≤ 0.01"),
        levels = c("p > 0.10",
                   "0.05 < p ≤ 0.10",
                   "0.01 < p ≤ 0.05",
                   "p ≤ 0.01")
      ),
      x = NA,
      y = NA
    )
    # Add these layers BEFORE geom_point()
    contour_layers <-
      # 0.05 < p <= 0.10  (between x10 and x05)
      list(
        geom_rect(data = contours,
                  aes(xmin = -x05, xmax = -x10, ymin = se_lo, ymax = se_hi),
                  inherit.aes = FALSE, fill = p10),
        geom_rect(data = contours,
                  aes(xmin =  x10, xmax =  x05, ymin = se_lo, ymax = se_hi),
                  inherit.aes = FALSE, fill = p10),
        
        # 0.01 < p <= 0.05  (between x05 and x01)
        geom_rect(data = contours,
                  aes(xmin = -x01, xmax = -x05, ymin = se_lo, ymax = se_hi),
                  inherit.aes = FALSE, fill = p05),
        geom_rect(data = contours,
                  aes(xmin =  x05, xmax =  x01, ymin = se_lo, ymax = se_hi),
                  inherit.aes = FALSE, fill = p05),
        
        # p <= 0.01  (beyond x01 to infinity)
        geom_rect(data = contours,
                  aes(xmin = -Inf, xmax = -x01, ymin = se_lo, ymax = se_hi),
                  inherit.aes = FALSE, fill = p01),
        geom_rect(data = contours,
                  aes(xmin =  x01, xmax =  Inf, ymin = se_lo, ymax = se_hi),
                  inherit.aes = FALSE, fill = p01),
        geom_point(
          data = contour_legend,
          aes(x = x, y = y, fill = region),
          shape = 22,
          size = 5,
          inherit.aes = FALSE
        ),
        scale_fill_manual(
          values = c(
            "p > 0.10" = "white",
            "0.05 < p ≤ 0.10" = p10,
            "0.01 < p ≤ 0.05" = p05,
            "p ≤ 0.01" = p01
          ),
          name = "",
        ),
        geom_abline(aes(slope = -1/qnorm(pval), 
                        intercept = 0), linetype = "dashed"),
        geom_abline(aes(slope =  1/qnorm(pval), 
                          intercept = 0), linetype = "dashed"),
        geom_vline(xintercept = 0, color = "gray50",
                   linetype = "dotted")
      )
    
    funnel_centered <- demeaned_data %>% 
      filter(!is.na(y.dm) & !is.na(y)) %>%
      ggplot(aes(x = y.dm, y = .se)) +
      contour_layers +
      geom_point(shape = 17) +
      ylab("Standard Error") + 
      xlab(xlab_funnel) + 
      vit$theme_vit() + 
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      coord_cartesian(
        xlim = c(-1,1) * 1.05 * max(abs(demeaned_data$y.dm)),
        ylim = c(0, max(demeaned_data$.se[!is.na(demeaned_data$y.dm)])*1.05)
      ) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    funnel_facet_centered <-  funnel_centered +
      facet_wrap(~intervention)
    
    funnel_centered_high <- demeaned_data_high %>% 
      filter(!is.na(y.dm) & !is.na(y)) %>%
      ggplot(aes(x = y.dm, y = .se)) +
      contour_layers +
      geom_point(shape = 17) +
      ylab("Standard Error") + 
      xlab(xlab_funnel) + 
      vit$theme_vit() + 
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      coord_cartesian(
        xlim = c(-1,1) * 1.05 * max(abs(demeaned_data_high$y.dm)),
        ylim = c(0, max(demeaned_data_high$.se[!is.na(demeaned_data_high$y.dm)])*1.05)
      ) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    funnel_facet_centered_high <-  funnel_centered_high +
      facet_wrap(~intervention)
    
    funnel_facet_centered_def <-  demeaned_data_def %>% 
      filter(!is.na(y.dm) & !is.na(y)) %>%
      mutate(bias = forcats::fct_recode(bias, 
                                        "low" = "low",
                                        "unclear" = "some concerns",
                                        "high" = "high")) %>%
      mutate(bias = factor(bias, levels = c("low", "unclear", "high"), ordered = TRUE)) %>%
      ggplot(aes(x = y.dm, y = .se, color = bias)) +
      contour_layers +
      geom_point(shape = 17) +
      scale_color_manual(values = vit$bias_study_colors) +
      ylab("Standard Error") + 
      xlab(xlab_funnel) + 
      vit$theme_vit() + 
      theme(legend.position = "bottom") +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      coord_cartesian(
        xlim = c(-1,1) * 1.05 * max(abs(demeaned_data_def$y.dm)),
        ylim = c(0, max(demeaned_data_def$.se[!is.na(demeaned_data_def$y.dm)])*1.05)
      ) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      facet_wrap(~intervention)
    
    funnel_overall_debias_plot <- bias_demeaned_data %>% 
      filter(!is.na(y.dm) & !is.na(y)) %>%
      mutate(bias = forcats::fct_recode(bias.overall, 
                                        "low" = "low",
                                        "unclear" = "some concerns",
                                        "high" = "high")) %>%
      mutate(bias = factor(bias, 
                           levels = c("low", "unclear", "high"),
                           ordered = TRUE)) %>%
      ggplot(aes(x = y.dm, 
                 y = .se, 
                 color = bias)) +
      contour_layers + 
      geom_point(shape = 17) +
      scale_color_manual(values = vit$bias_study_colors) +
      ylab("Standard Error") + 
      xlab(xlab_funnel) + 
      vit$theme_vit() + 
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      coord_cartesian(
        xlim = c(-1,1) * 1.05 * max(abs(bias_demeaned_data$y.dm)),
        ylim = c(0, max(bias_demeaned_data$.se[!is.na(bias_demeaned_data$y.dm)])*1.05)
      ) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "bottom")
    
    funnel_facet_debias <- funnel_overall_debias_plot +
      facet_wrap(~intervention) 
    
    funnel_def <- demeaned_data_def %>% 
      filter(!is.na(y.dm) & !is.na(y)) %>%
      mutate(bias = forcats::fct_recode(bias, 
                                        "low" = "low",
                                        "unclear" = "some concerns",
                                        "high" = "high")) %>%
      mutate(bias = factor(bias, levels = c("low", "unclear", "high"), ordered = TRUE)) %>%
      ggplot(aes(x = y.dm, y = .se, color = bias)) +
      contour_layers +
      geom_point(shape = 17) +
      scale_color_manual(values = vit$bias_study_colors) +
      ylab("Standard Error") + 
      xlab(xlab_funnel) + 
      vit$theme_vit() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      coord_cartesian(
        xlim = c(-1,1) * 1.05 * max(abs(demeaned_data_def$y.dm)),
        ylim = c(0, max(demeaned_data_def$.se[!is.na(demeaned_data_def$y.dm)])*1.05)
      )
    
    funnel_def_low <- demeaned_data_def_low %>% 
      filter(!is.na(y.dm) & !is.na(y)) %>%
      ggplot(aes(x = y.dm, y = .se)) +
      contour_layers +
      geom_abline(aes(slope = -1/qnorm(pval), intercept = 0), linetype = "dashed") +
      geom_abline(aes(slope =  1/qnorm(pval), intercept = 0), linetype = "dashed") +
      geom_point(shape = 17) +
      ylab("Standard Error") + 
      xlab(xlab_funnel) + 
      vit$theme_vit()  +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      coord_cartesian(
        xlim = c(-1,1) * 1.05 * max(abs(demeaned_data_def_low$y.dm)),
        ylim = c(0, max(demeaned_data_def_low$.se[!is.na(demeaned_data_def_low$y.dm)])*1.05)
      )
    
    funnel_facet_def_low <- demeaned_data_def_low %>% 
      filter(!is.na(y.dm) & !is.na(y)) %>%
      ggplot(aes(x = y.dm, y = .se)) +
      contour_layers +
      geom_abline(aes(slope = -1/qnorm(pval), intercept = 0), linetype = "dashed") +
      geom_abline(aes(slope =  1/qnorm(pval), intercept = 0), linetype = "dashed") +
      geom_point(shape = 17) +
      ylab("Standard Error") + 
      xlab(xlab_funnel) + 
      vit$theme_vit()  +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      coord_cartesian(
        xlim = c(-1,1) * 1.05 * max(abs(demeaned_data_def_low$y.dm)),
        ylim = c(0, max(demeaned_data_def_low$.se[!is.na(demeaned_data_def_low$y.dm)])*1.05)
      ) + 
      facet_wrap(~intervention)
    
    figure_save(funnel_centered,
                outdir = "outputs",
                outcome = outcome,
                filename = "funnel_plot",
                height = 6,
                width = 8)
    
    figure_save(funnel_facet_centered,
                outdir = "outputs",
                outcome = outcome,
                filename = "funnel_plot_by_intervention_centered",
                height = 6,
                width = 8)
    
    figure_save(funnel_centered_high,
                outdir = "outputs",
                outcome = outcome,
                filename = "funnel_plot_high",
                height = 6,
                width = 8)
    
    figure_save(funnel_facet_centered_high,
                outdir = "outputs",
                outcome = outcome,
                filename = "funnel_plot_by_intervention_centered_high",
                height = 6,
                width = 8)
    
    figure_save(funnel_facet_centered_def,
                outdir = "outputs",
                outcome = outcome,
                filename = "funnel_plot_by_intervention_and_def_centered",
                height = 6,
                width = 8)
    
    figure_save(funnel_overall_debias_plot,
                outdir = "outputs",
                outcome = outcome,
                filename = "funnel_plot_debiasing_regression_adjusted",
                height = 5,
                width = 7
    )
    
    figure_save(funnel_facet_debias,
                outdir = "outputs",
                outcome = outcome,
                filename = "funnel_plot_by_intervention_debiased",
                height = 6,
                width = 8)
    
    figure_save(funnel_def,
                outdir = "outputs",
                outcome = outcome,
                filename = "funnel_plot_by_def_centered_with_contours",
                height = 6,
                width = 8)
    
     figure_save(funnel_def_low,
                 outdir = "outputs",
                 outcome = outcome,
                 filename = "funnel_plot_by_def_centered_low_bias_only_with_contours",
                 height = 6,
                 width = 8)
     
     figure_save(funnel_facet_def_low,
                 outdir = "outputs",
                 outcome = outcome,
                 filename = "funnel_plot_by_def_centered_low_bias_only_with_contours_facet",
                 height = 6,
                 width = 8)
  }
  
  # country sd effects
  {
    
    sd_col_label <- switch(outcome,
                           anxiety = "yanx_confounders_s_se_anx",
                           depression = "y_interactions_s_se")
    sd_b_label <- switch(outcome,
                         anxiety = "b_yanx_confounders_s_se_anx",
                         depression = "b_y_confounders_s_se")
    
    count_intercept_label <- switch(outcome,
                                    anxiety = "yanx_confounders_Intercept",
                                   depression = "y_interactions_Intercept")
    
    sd_se <- switch(outcome,
                    anxiety = bias_reg$prep$data$sd_se_anx %>% max(na.rm=TRUE),
                    depression = bias_reg$prep$data$sd_se[1])
    m_se  <- switch(outcome,
                    anxiety = bias_reg$prep$data$m_se_anx %>% max(na.rm=TRUE),
                    depression = bias_reg$prep$data$m_se[1])
    re <- brms::ranef(bias_reg, summary = FALSE, probs = c(0.025, 0.17,0.83,.975))$country
    re[,,sd_col_label] <- sweep(re[,,sd_col_label],1,posterior::as_draws_matrix(bias_reg, variable = sd_b_label)[,1,drop = FALSE], FUN = "+")
    # re[,,"y_interactions_Intercept"] <- re[,,"y_interactions_Intercept"] - (re[,,"y_interactions_s_se"]  * m_se/sd_se)
    re_sum <- brms::posterior_summary(re, probs = c(0.025,0.17,0.83,0.975))
    re_int <- re_sum[,,count_intercept_label] %>% as.data.frame() %>% 
      mutate(country = rownames(re_sum),
             variable = "Intercept") 
    re_coef<- re_sum[,,sd_col_label] %>% as.data.frame() %>% 
      mutate(country = rownames(re_sum),
             variable = "S.E. Coefficient") %>% 
      mutate(sd.y.scale = sd_se) %>% 
      mutate(Estimate = Estimate/sd.y.scale,
             `Q2.5` = `Q2.5`/sd.y.scale,
             Q17 = Q17/sd.y.scale,
             Q83 = Q83/sd.y.scale,
             `Q97.5` = `Q97.5`/sd.y.scale)
    
    country_coef_plt <- bind_rows(re_int,
                                  re_coef) %>%
      ggplot(aes(
        x    = Estimate,
        y    = forcats::fct_rev(country)
      )) +
      geom_vline(xintercept = 0, linetype = "dashed",color = "gray70") +
      ggdist::geom_pointinterval(aes(xmin = Q2.5,
                                     xmax = Q97.5,),
                                 linewidth = 0.6) +
      ggdist::geom_pointinterval(aes(xmin = Q17,
                                     xmax = Q83),
                                 size = 4) +
      
      ggplot2::facet_grid(cols = vars(variable), rows = vars(country), scales = "free",
                          switch = "y")  +
      ylab("Country") +
      vit$theme_vit() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
      )
    
    figure_save(country_coef_plt,
                outdir = "outputs",
                outcome = outcome,
                filename = "bias_regression_country_effects",
                height = 5,
                width = 6)
  }
  
  # country specific funnel plot
  {
    iran.ests <- re[,"Iran",c(count_intercept_label,sd_col_label)]
    
    pval <- 0.05
    sdyseq <- seq(-0.9,9.5,.1)
    preds <- cbind(1,  sdyseq) %*% t(iran.ests )
    preds.df <- data.frame(
      y.dm = rowMeans(preds),
      low =  apply(preds, 1, quantile, .025),
      high = apply(preds, 1, quantile, .975),
      .se =  (sdyseq * sd_se + m_se)
    )
    
    coef_E <- colMeans( iran.ests )
    
    # dat %>% 
    #   filter(!is.na(y)) %>% 
    #   mutate(y.dm = y - fitted.values(fit_high)[,"Estimate"]) %>%
    #   select(y.dm, sd.y, country) %>% 
      
      
    iran_funnel <- demeaned_data_def %>%
      mutate(country = fit_def$prep$data$country) %>% 
      filter(!is.na(y)) %>%
      ggplot(aes(x = y.dm, y = .se)) +
      geom_abline(slope = -1/qnorm(pval), intercept = 0, linetype = "dashed") +
      geom_abline(slope =  1/qnorm(pval), intercept = 0, linetype = "dashed") +
      geom_vline(xintercept = 0, color = "gray") +
      geom_ribbon(data = preds.df,
                  aes(xmin = low, xmax = high), alpha = 0.2) +
      geom_point(data = ~ dplyr::filter(.x, country == "Iran"), 
                 color = "red") +
      geom_point(data = ~ dplyr::filter(.x, country != "Iran")) +
      
      geom_line(data = preds.df, color = "#e6550d", linewidth = 1) +
      scale_y_continuous(limits = c(0, 3)
                         , expand = c(0,0)
      ) +
      xlab(latex2exp::TeX("SMD - \\widehat{\\mu}")) +
      ylab("Standard error") +
      # ggtitle("Funnel plot with predicted effects for Iranian studies") +
      vit$theme_vit()
    
    figure_save(iran_funnel,
                outdir = "outputs",
                outcome = outcome,
                filename = "funnel_plot_iran_studies_highlighted",
                height = 5,
                width = 7)
    
  }
  
  # dose response plot
  {
    dose_summ_low <- fit_dose_low %>%
      vit$summary_brms_nma(keep = c(".trt",".dose")) 
    
    ld_p <- dose_summ_low %>% 
      filter(.trt == "D") %>%
      mutate(.dose = factor(.dose,
                            levels = c("D (< 1000 IU)",
                                       "D (1000 - 3000 IU)",
                                       "D (> 3000 IU)"),
                            ordered = TRUE)) %>% 
      # mutate(deficiency = forcats::fct_relevel(deficiency, "unknown", "no","yes")) %>%
      droplevels() %>% 
      filter(.observed) %>% select(-.observed) %>% 
      plot(group = ".dose", position = position_dodge(width = .5)) +
      aes(shape = .dose) +
      ggplot2::scale_shape_manual( #legend for doses
        name   = "Total daily dose",
        values = c("D (< 1000 IU)" = 15, "D (1000 - 3000 IU)" = 16, "D (> 3000 IU)" = 17)
      ) +
      theme(strip.text.y.right = element_blank(),
            # axis.text.y.left = element_blank(),
            strip.text.y.left = element_blank())
      
     figure_save(ld_p,
                 outdir = "outputs",
                 outcome = outcome,
                 filename = "vitamin_d_dose_response_low_bias",
                 height = 3.5,
                 width = 6)  
    
    # mutate(.dose = factor(.dose, levels = vit$restore_dose_names(fit_dose_low$prep$network$agd_contrast[[fit_dose_low$prep$class_ids$varname]])$varname,
      #                       labels = vit$restore_dose_names(fit_dose_low$prep$network$agd_contrast[[fit_dose_low$prep$class_ids$varname]])$label))
    # 
    # dose_summ_low %>% 
    #   arrange(.trt, .dose) %>% 
    #   filter(.trt == "D") %>%
    #   # filter(bias == "low") %>% 
    #   mutate(.dose = factor(.dose,
    #                         levels = c("D (< 1000 IU)",
    #                                    "D (1000 - 3000 IU)",
    #                                    "D (> 3000 IU)"),
    #                         ordered = TRUE))
    
    dose_summ <- fit_dose %>% 
      vit$summary_brms_nma(keep = c(".trt",".dose", "bias")) %>% 
      # mutate(.dose = factor(.dose, levels = vit$restore_dose_names(fit_dose$prep$network$agd_contrast[[fit_dose$prep$class_ids$varname]])$varname,
      #                       labels = vit$restore_dose_names(fit_dose$prep$network$agd_contrast[[fit_dose$prep$class_ids$varname]])$label)) %>% 
      mutate(bias = forcats::fct_recode(bias,
                           "unclear" = "some concerns"),
             # !!outcome_label := factor(!!sym(outcome_label), levels = 1:0, labels = c("yes", "no")),
             # def = forcats::fct_recode(def,
             #                                  "unknown" = "NA")
             ) %>% 
      select(bias, 
             # deficiency = def,
             # !!outcome_label,
             .trt, .dose, iter, value)
    
    pb <- (dose_summ %>% 
      arrange(.trt, bias) %>%
      filter(.trt == "D") %>%
      mutate(.dose = factor(.dose,
                            levels = c("D (< 1000 IU)",
                                       "D (1000 - 3000 IU)",
                                       "D (> 3000 IU)"),
                            ordered = TRUE)) %>% 
      # mutate(deficiency = forcats::fct_relevel(deficiency, "unknown", "no","yes")) %>%
      droplevels() %>% 
        mutate(bias = forcats::fct_recode(bias,
                                        "unclear" = "some concerns")) %>%
      mutate(bias = forcats::fct_relevel(bias, "low",
                                           "unclear",
                                           "high")) %>%
      # vit$def_plot(group = .dose, position = position_dodge(preserve = "single", orientation = "y")) +
        plot(group = ".dose", position = position_dodge(width = .5)) +
      aes(shape = .dose) +
      ggplot2::scale_shape_manual( #legend for doses
        name   = "Total daily dose",
        values = c("D (< 1000 IU)" = 15, "D (1000 - 3000 IU)" = 16, "D (> 3000 IU)" = 17)
      ) +
      theme(strip.text.y.right = element_blank(),
            # axis.text.y.left = element_blank(),
            strip.text.y.left = element_blank()) +
        scale_y_discrete(
          name = NULL,
          limits = c(
            "bias: low",
            "bias: unclear",
            "bias: high"
          ),
          labels = c(
            "low",
            "unclear",
            "bias: high"
          )
        )
      
      ) %>% 
    ggplot_build()
    
    # extract and adjust y positions
    y_pos <- round(pb$data[[2]]$y)
    shape <- pb$data[[2]]$shape
    color <- pb$data[[2]]$colour
    
    pb$data[[2]]$y <- y_pos + ifelse(shape == 15, -0.2, ifelse(shape == 16, 0, 0.2))
    
    #convert back to grob
    g <- pb %>% ggplot_gtable() %>% 
      as_printable_gtable()
    
    figure_save(g,
                outdir = "outputs",
                outcome = outcome,
                filename = "vitamin_d_dose_response_by_bias",
                height = 5,
                width = 7.5)
    
  }
  
  # loo plots
  {
    sensitivity_outdir <- here::here("outputs", "saved_models", outcome,
                                     "sensitivity")

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

    vit$loo_plot(fit_low, loo_low) %>%
      figure_save(
        outdir = "outputs",
        outcome = outcome,
        filename = "loo_plot_low_bias",
        height = 5,
        width = 7
      )
    vit$loo_plot(fit_some, loo_some) %>%
      figure_save(
        outdir = "outputs",
        outcome = outcome,
        filename = "loo_plot_some_and_low_bias",
        height = 5,
        width = 7
      )
    vit$loo_plot(fit_high, loo_high) %>%
      figure_save(
        outdir = "outputs",
        outcome = outcome,
        filename = "loo_plot_all_bias",
        height = 9,
        width = 7
      )
  }
  
  # gain vs non-gain score estimates
  {
    low_nogain <- bind_rows(
      fit_low_nogain %>% vit$summary_brms_nma(keep = ".trt") %>% mutate(gain = "endpoint"),
      low_sum %>% mutate(gain = "full") %>% select(-bias))
    
    some_nogain <- bind_rows(
      fit_some_nogain %>% vit$summary_brms_nma(keep = ".trt") %>% mutate(gain = "endpoint"),
      some_sum %>% mutate(gain = "full")  %>% select(-bias)
      )
    
    high_nogain <- bind_rows(
      fit_high_nogain %>% vit$summary_brms_nma(keep = ".trt") %>% mutate(gain = "endpoint"),
      high_sum %>% mutate(gain = "full") %>% select(-bias)
    )
    
    low_nogain_plot <- low_nogain %>% 
      plot(color = "gain", 
           position = position_dodge(preserve = "single", 
                                     width = .2, )) +
      aes(y = NULL, group = gain) +
      scale_color_manual(values = c("full" = "black",
                                    "endpoint"  = "gray63"
      ),      name = "Endpoint-only vs Full Estimates") +
      ylab("")
    
    low_nogain_plot %>% 
      figure_save(
                  outdir = "outputs",
                  outcome = outcome,
                  filename = "gain_vs_no_gain_low_bias",
                  height = 5,
                  width = 7.5)
      
    some_nogain_plot <- some_nogain %>% 
      plot(color = "gain", 
           position = position_dodge(preserve = "single", 
                                     width = .3, )) +
      aes(y = NULL, group = gain) +
      scale_color_manual(values = c("full" = "black",
                                    "endpoint"  = "gray63"
      ),      name = "Endpoint-only vs Full Estimates") +
      ylab("")
      
      
      some_nogain_plot %>% 
      figure_save(
        outdir = "outputs",
        outcome = outcome,
        filename = "gain_vs_no_gain_low_and_some_bias",
        height = 5,
        width = 7.5)
    
    high_nogain_plot <- high_nogain %>% 
      plot(color = "gain", 
           position = position_dodge(preserve = "single", 
                                     width = .5, )) +
      aes(y = NULL, group = gain) +
      scale_color_manual(values = c("full" = "black",
                                    "endpoint"  = "gray63"
      ),      name = "Endpoint-only vs Full Estimates") +
      ylab("")
    
    high_nogain_plot %>% 
      figure_save(
        outdir = "outputs",
        outcome = outcome,
        filename = "gain_vs_no_gain_all_bias",
        height = 7,
        width = 7.5)
    
  }
  
  # heterogeneity plot
  {
    bias_het_name <- switch(
      outcome,
      depression = "tau_y",
      anxiety = "tau_yanx"
    )
    heterogeneity_plot <- bind_rows(
      posterior::as_draws_df(bias_reg, variable = bias_het_name) %>%
        select(any_of(bias_het_name)) %>% mutate(tau = .data[[bias_het_name]],
                                 model = "meta\nregression"),
      fit_def_no_iran %>% 
        posterior::as_draws_df(variable = "tau") %>%
        select(tau) %>%
        mutate(model = "deficiency v.s.\ndepressed\nNo Iranian studies"),
      posterior::as_draws_df(fit_def, variable = "tau") %>%
        select(tau) %>%
        mutate(model = "deficiency v.s.\ndepressed"),
      posterior::as_draws_df(fit_high, variable = "tau") %>%
        select(tau) %>%
        mutate(model = "simple model") 
    ) %>% 
      mutate(model = forcats::fct_relevel(model,
                                    "meta\nregression",
                                    "deficiency v.s.\ndepressed\nNo Iranian studies",
                                    "deficiency v.s.\ndepressed",
                                    "simple model")) %>% 
      ggplot(aes(x = tau, fill = model)) +
      ggridges::geom_ridgeline(aes(y = model, 
                                   height = after_stat(density), 
                                   group = model),
                               stat = "density_ridges", scale = 0.25,
                               alpha = 0.6) +
      scale_fill_manual(values = c(
        "meta\nregression" = "#e6550d",
        "deficiency v.s.\ndepressed\nNo Iranian studies" = "#3182bd",
        "deficiency v.s.\ndepressed" = "#31a354",
        "simple model" = "#6a51a3"
      ),
      breaks = c()
      ) +
      xlab(latex2exp::TeX("$\\tau$")) +
      ggtitle(latex2exp::TeX("Posterior distribution of study heterogeneity ($\\tau$) for all data")) +
      vit$theme_vit()
    
    heterogeneity_plot_low <- bind_rows(
      posterior::as_draws_df(fit_def_low_no_iran, variable = "tau") %>%
        select(tau) %>%
        mutate(model = "deficiency v.s.\ndepressed\nNo Iranian studies"),
      posterior::as_draws_df(fit_def_low, variable = "tau") %>%
        select(tau) %>%
        mutate(model = "deficiency v.s.\ndepressed"),
      posterior::as_draws_df(fit_low, variable = "tau") %>%
        select(tau) %>%
        mutate(model = "simple model") 
    ) %>% 
      mutate(model = forcats::fct_relevel(model,
                                          "deficiency v.s.\ndepressed\nNo Iranian studies",
                                          "deficiency v.s.\ndepressed",
                                          "simple model")) %>% 
      ggplot(aes(x = tau, fill = model)) +
      ggridges::geom_ridgeline(aes(y = model, 
                                   height = after_stat(density), 
                                   group = model),
                               stat = "density_ridges", scale = 0.25,
                               alpha = 0.6) +
      scale_fill_manual(values = c(
        "deficiency v.s.\ndepressed\nNo Iranian studies" = "#3182bd",
        "deficiency v.s.\ndepressed" = "#31a354",
        "simple model" = "#6a51a3"
      ),
      breaks = c()
      ) +
      xlab(latex2exp::TeX("$\\tau$")) +
      ggtitle(latex2exp::TeX("Posterior distribution of study heterogeneity ($\\tau$) for all data")) +
      vit$theme_vit()
    
    figure_save(heterogeneity_plot,
                outdir = "outputs",
                outcome = outcome,
                filename = "heterogeneity_comparison_all_data",
                height = 5,
                width = 7)
    
    figure_save(heterogeneity_plot_low,
                outdir = "outputs",
                outcome = outcome,
                filename = "heterogeneity_comparison_low_bias",
                height = 5,
                width = 7)
  }
  
  
  # treatment ranks
  # only makes sense for non-deficient populations
  {
    low_ranks_dis <- tryCatch(
      fit_def_low %>% 
      vit$rank_probs_vit(sucra = TRUE,
                         cumulative = TRUE,
                         filter = rlang::quo(def == "no" & .data[[outcome_label]] == "yes") ) %>%
      select(-c(def, any_of(outcome_label))),
      error = function(e) {NULL})
    
    if(!is.null(low_ranks_dis)) {
      seq(1, 
          length(unique(low_ranks_dis$.trt)), 
          by = 2) -> breaks_seq
      if(length(breaks_seq) < 10) {
        breaks_seq <- seq(1, 
                          length(unique(low_ranks_dis$.trt)), 
                          by = 1)
      }
      (low_ranks_dis %>% 
        plot() +
        # ggtitle("Treatment ranks (low bias studies)") +
        scale_x_continuous(breaks =breaks_seq, name = "Rank" )) %>% 
        figure_save(
          outdir = "outputs",
          outcome = outcome,
          filename = "treatment_ranks_low_bias_data_no_def_dis",
          height = 5,
          width = 7
        )
    }
    
    low_ranks_no_dis <- tryCatch(
      fit_def_low %>% 
      vit$rank_probs_vit(sucra = TRUE,
                         cumulative = TRUE,
                         filter = rlang::quo(def == "no" & .data[[outcome_label]] == "no") ) %>%
      select(-c(def, any_of(outcome_label))) , error = function(e) {NULL})
    
    if(!is.null(low_ranks_no_dis)) {
      seq(1, 
          length(unique(low_ranks_no_dis$.trt)), 
          by = 2) -> breaks_seq
      if(length(breaks_seq) < 10) {
        breaks_seq <- seq(1, 
                          length(unique(low_ranks_no_dis$.trt)), 
                          by = 1)
      }
      (low_ranks_no_dis %>% 
        plot() +
        # ggtitle("Treatment ranks (low bias studies)") +
        scale_x_continuous(breaks = breaks_seq, name = "Rank" )) %>% 
        figure_save(
          outdir = "outputs",
          outcome = outcome,
          filename = "treatment_ranks_low_bias_data_no_def_no_dis",
          height = 5,
          width = 7
        )
    }
    
    low_ranks_dis_unknown_def <- tryCatch(fit_def_low %>% 
                                               vit$rank_probs_vit(sucra = TRUE,
                                                                  cumulative = TRUE,
                                                                  filter = rlang::quo(def == "NA" & .data[[outcome_label]] == "yes") ) %>%
                                               select(-c(def, any_of(outcome_label))) , error = function(e){NULL})
    
    if(!is.null(low_ranks_dis_unknown_def)) {
      seq(1, 
          length(unique(low_ranks_dis_unknown_def$.trt)), 
          by = 2) -> breaks_seq
      if(length(breaks_seq) < 10) {
        breaks_seq <- seq(1, 
                          length(unique(low_ranks_dis_unknown_def$.trt)), 
                          by = 1)
      }
      (low_ranks_dis_unknown_def %>% 
          mutate(`p_rank[1]` = ifelse(is.na(`p_rank[1]`), 0, `p_rank[1]`)) %>%
        plot() +
        # ggtitle("Treatment ranks (low bias studies)") +
        scale_x_continuous(breaks = breaks_seq, name = "Rank" )) %>% 
        figure_save(
          outdir = "outputs",
          outcome = outcome,
          filename = "treatment_ranks_low_bias_data_by_unknown_def_and_dis",
          height = 5,
          width = 7
        )
    }
    
    
    low_ranks_no_dis_unknown_def <- tryCatch(fit_def_low %>% 
                                            vit$rank_probs_vit(sucra = TRUE,
                                                               cumulative = TRUE,
                                                               filter = rlang::quo(def == "NA" & .data[[outcome_label]] == "no") ) %>%
                                            select(-c(def, any_of(outcome_label))) , error = function(e){NULL})
    
    if(!is.null(low_ranks_no_dis_unknown_def)) {
      seq(1, 
          length(unique(low_ranks_no_dis_unknown_def$.trt)), 
          by = 2) -> breaks_seq
      if(length(breaks_seq) < 10) {
        breaks_seq <- seq(1, 
                          length(unique(low_ranks_no_dis_unknown_def$.trt)), 
                          by = 1)
      }
      (low_ranks_no_dis_unknown_def %>% 
        plot() +
        # ggtitle("Treatment ranks (low bias studies)") +
        scale_x_continuous(breaks = breaks_seq, name = "Rank" ) ) %>% 
        figure_save(
          outdir = "outputs",
          outcome = outcome,
          filename = "treatment_ranks_low_bias_data_by_unknown_def_and_no_dis",
          height = 5,
          width = 7
        )
    }

    # all_ranks <- fit_high %>%
    #   vit$rank_probs_vit(sucra = TRUE,
    #                      cumulative = TRUE) %>%
    #   plot() +
    #   scale_x_continuous(breaks = seq(1, vit$ntreatments(fit_high), by = 2), name = "Rank" )

    # figure_save(low_ranks,
    #             outdir = "outputs",
    #             outcome = outcome,
    #             filename = "treatment_ranks_low_bias_data",
    #             height = 5,
    #             width = 7)
    # 
    # figure_save(all_ranks,
    #             outdir = "outputs",
    #             outcome = outcome,
    #             filename = "treatment_ranks_all_data",
    #             height = 5,
    #             width = 7)

  }
  
  # inconsistency plot
  {
    low_inconsistency <- here::here("outputs", 
                                    "saved_models", 
                                    outcome,  
                                    "sensitivity", "nodesplit_low.rds") %>%
      readRDS()
    
    some_inconsistency <- here::here("outputs", 
                                     "saved_models", 
                                     outcome,  
                                     "sensitivity", "nodesplit_some.rds") %>%
      readRDS()
    
    high_inconsistency <- here::here("outputs", 
                                     "saved_models", 
                                     outcome,  
                                     "sensitivity", "nodesplit_high.rds") %>%
      readRDS()
    
    if (!is.null(low_inconsistency)) {
      low_inc_plot <- low_inconsistency[[1]]$comparison %>% 
        select(-value) %>% 
        tidyr::pivot_longer(cols = c("direct","indirect","network"),
                            names_to = "evidence",
                            values_to = "value") %>% 
        {
          class(.) <- c("summary_brms_nma", class(.))
          attr(., "reference_treatment") <-  vit$get_reference_treatment(low_inconsistency[[1]]$comparison)
          attr(., "ndraws") <- attr(low_inconsistency[[1]]$comparison, "ndraws")
          .
        } %>% 
        mutate(evidence = factor(evidence,
                                 levels = c("network","indirect","direct"),
                                 labels = c("Network",
                                            "Indirect",
                                            "Direct"))) %>%
        mutate(.trt = gsub("\\[","",.trt) %>% 
                 gsub("\\]","",.)) %>%
        mutate(.trt = stringr::str_wrap(.trt, width = 20)) %>%
        ggplot2::ggplot(ggplot2::aes( 
          y = evidence,
          x = value
        )) + 
        ggplot2::geom_vline(xintercept = 0, lty = 2, col = "grey50") +
        ggdist::stat_pointinterval() +
        ggplot2::xlab("SMD") + 
        vit$theme_vit(base_size = 11) +
        facet_wrap(~.trt, ncol = 1, scales = "free_y", strip.position = "right") +
        ylab("Evidence")
      
      figure_save(low_inc_plot,
                  outdir = "outputs",
                  outcome = outcome,
                  filename = "inconsistency_low_bias",
                  height = 3,
                  width = 3.5)
    }
    
    if(! is.null(some_inconsistency) ) {
      some_inc_plot <- some_inconsistency[[1]]$comparison %>% 
        select(-value) %>% 
        tidyr::pivot_longer(cols = c("direct","indirect","network"),
                            names_to = "evidence",
                            values_to = "value") %>% 
        {
          class(.) <- c("summary_brms_nma", class(.))
          attr(., "reference_treatment") <-  vit$get_reference_treatment(low_inconsistency[[1]]$comparison)
          attr(., "ndraws") <- attr(low_inconsistency[[1]]$comparison, "ndraws")
          .
        } %>% 
        mutate(evidence = factor(evidence,
                                 levels = c("network","indirect","direct"),
                                 labels = c("Network",
                                            "Indirect",
                                            "Direct"))) %>%
        mutate(.trt = gsub("\\[","",.trt) %>% 
                 gsub("\\]","",.)) %>%
        mutate(.trt = stringr::str_wrap(.trt, width = 20)) %>%
        ggplot2::ggplot(ggplot2::aes( 
          y = evidence,
          x = value
        )) + 
        ggplot2::geom_vline(xintercept = 0, lty = 2, col = "grey50") +
        ggdist::stat_pointinterval() +
        ggplot2::xlab("SMD") + 
        vit$theme_vit(base_size = 11) +
        facet_wrap(~.trt, ncol = 1, scales = "free_y", strip.position = "right") +
        ylab("Evidence")
      
      figure_save(some_inc_plot,
                  outdir = "outputs",
                  outcome = outcome,
                  filename = "inconsistency_some_and_low_bias",
                  height = 3,
                  width = 3.5)
      
    }
    
    if(!is.null(high_inconsistency)) {
      high_inc_plot <- high_inconsistency %>% 
        lapply(function(x){
          x$comparison
          
        }) %>% 
        bind_rows() %>%
        mutate(.trt = gsub("\\[","",.trt) %>% 
                 gsub("\\]","",.)) %>%
        mutate(
          across(
            c(direct, indirect, network),
            ~ ifelse(sub(" .*", "", .trt) == "placebo", -.x, .x)
          )
        ) %>% 
        mutate(.trt = ifelse(sub(" .*", "", .trt) == "placebo", paste0(.trt, " vs. placebo"),.trt),
               .trt = gsub("placebo vs. ", "", .trt)
               ) %>%
        mutate(.trt = stringr::str_wrap(.trt, width = 20)) %>% 
        select(-value) %>% 
        tidyr::pivot_longer(cols = c("direct","indirect","network"),
                            names_to = "evidence",
                            values_to = "value") %>% 
        {
          class(.) <- c("summary_brms_nma", class(.))
          attr(., "reference_treatment") <-  vit$get_reference_treatment(high_inconsistency[[1]]$comparison)
          attr(., "ndraws") <- attr(high_inconsistency[[1]]$comparison, "ndraws")
          .
        } %>% 
        mutate(evidence = factor(evidence,
                                 levels = c("network","indirect","direct"),
                                 labels = c("Network",
                                            "Indirect",
                                            "Direct"))) %>%
        ggplot2::ggplot(ggplot2::aes( 
          y = evidence,
          x = value
        )) + 
        ggplot2::geom_vline(xintercept = 0, lty = 2, col = "grey50") +
        ggdist::stat_pointinterval() +
        ggplot2::xlab("SMD") + 
        vit$theme_vit(base_size = 11) +
        facet_wrap(~.trt, ncol = 1, scales = "free_y", strip.position = "right") +
        theme(panel.spacing.y = ggplot2::unit(1, "lines")) +
        ylab("Evidence")
      
      
      figure_save(high_inc_plot,
                  outdir = "outputs",
                  outcome = outcome,
                  filename = "inconsistency_all_bias",
                  height = 5,
                  width = 3.5)
      }
    
                
  }
  
  # prediction intervals 
  {
    
    low_pred <- fit_low %>% vit$summary_brms_nma() %>% 
      filter(.observed) %>% select(-.observed) %>% 
      plot() + 
      ggdist::stat_pointinterval(data = fit_low %>% 
            vit$prediction_intervals(keep = c(".trt"),only_observed = FALSE) %>% 
              filter(b == "Placebo") %>%
              mutate(.trt = a %>% stringr::str_wrap(width = vit$width)),
          .width = .95, linewidth = .1, 
          color = "red", size = 0.001) +
      ggdist::stat_pointinterval()
    
    some_pred <- fit_some %>% vit$summary_brms_nma() %>% 
      filter(.observed) %>% select(-.observed) %>% 
      plot() + ggdist::stat_pointinterval(data =fit_some %>% 
                                            vit$prediction_intervals(keep = c(".trt"), only_observed = TRUE)  %>% 
                                            filter(b == "Placebo") %>%
                                            mutate(.trt = a %>% stringr::str_wrap(width = vit$width)) ,
                                          .width = .95, linewidth = .1, 
                                          color = "red", size = 0.001) +
      ggdist::stat_pointinterval()
    
    high_pred <- fit_high %>% vit$summary_brms_nma() %>% 
      filter(.observed) %>% select(-.observed) %>% 
      plot() + ggdist::stat_pointinterval(data = fit_high %>% 
           vit$prediction_intervals(keep = c(".trt"), only_observed = TRUE)  %>% 
           filter(b == "Placebo") %>%
           mutate(.trt = a %>% stringr::str_wrap(width = vit$width)),
                                          .width = .95, linewidth = .1, 
                                          color = "red", size = 0.001) +
      ggdist::stat_pointinterval()
    
    def_low_pred <- fit_def_low %>% 
      vit$summary_brms_nma(keep = c(".trt","def",outcome_label)) %>% 
      mutate(def = forcats::fct_recode(def,
                                       "unknown" = "NA")) %>%
      rename(deficiency = def) %>%
      mutate(bias = "low") %>% 
      filter(.observed) %>% 
      vit$def_plot() +
      ggdist::stat_pointinterval(data = fit_def_low %>% 
                                                   vit$prediction_intervals(keep = c(".trt","def",outcome_label),
                                                                            only_observed = TRUE) %>%
                                   mutate(def = forcats::fct_recode(def,
                                                                    "unknown" = "NA")) %>%
                                   rename(deficiency = def) %>%
                                   mutate(bias = "low")  %>% 
                                   filter(b == "Placebo") %>%
                                   mutate(.trt = a %>% stringr::str_wrap(width = vit$width)),
                                                 .width = .95, linewidth = .1, 
                                                 color = "red", size = 0.001) +
      ggdist::stat_pointinterval()
    
    def_low_pred_no_iran <- fit_def_low_no_iran %>% 
      vit$summary_brms_nma(keep = c(".trt","def",outcome_label)) %>% 
      mutate(def = forcats::fct_recode(def,
                                       "unknown" = "NA")) %>%
      rename(deficiency = def) %>%
      mutate(bias = "low") %>% 
      filter(.observed) %>% 
      select(-.observed) %>% 
      vit$def_plot() +
      ggdist::stat_pointinterval(data = fit_def_low_no_iran %>% 
                                   vit$prediction_intervals(keep = c(".trt","def",outcome_label), only_observed = TRUE) %>%
                                   mutate(def = forcats::fct_recode(def,
                                                                    "unknown" = "NA")) %>%
                                   rename(deficiency = def) %>%
                                   mutate(bias = "low")  %>% 
                                   filter(b == "Placebo") %>%
                                   mutate(.trt = a %>% stringr::str_wrap(width = vit$width)),
                                 .width = .95, linewidth = .1, 
                                 color = "red", size = 0.001) +
      ggdist::stat_pointinterval()
    
    
    def_pred <- fit_def %>% 
      vit$summary_brms_nma(keep = c(".trt","def",outcome_label,"bias")) %>% 
      mutate(def = forcats::fct_recode(def,
                                           "unknown" = "NA")) %>%
      rename(deficiency = def) %>%
      filter(.observed) %>% select(-.observed) %>% 
      mutate(bias = forcats::fct_recode(bias,
                                   "unclear" = "some concerns")) %>%
      mutate(bias = forcats::fct_relevel(
        bias,
        "low",
        "unclear",
        "high"
      )) %>% 
      vit$def_plot() +
      ggdist::stat_pointinterval(data = fit_def %>% 
                                   vit$prediction_intervals(keep = c(".trt","def",outcome_label, "bias"), only_observed = TRUE) %>%
                                   mutate(def = forcats::fct_recode(def,
                                                                    "unknown" = "NA")) %>%
                                   rename(deficiency = def) %>% 
                                   mutate(bias = forcats::fct_recode(bias,
                                                                                             "unclear" = "some concerns"))  %>% 
                                   filter(b == "Placebo") %>%
                                   mutate(.trt = a %>% stringr::str_wrap(width = vit$width)) ,
                                 .width = .95, linewidth = .1, 
                                 color = "red", size = 0.001) +
      ggdist::stat_pointinterval()
    
    def_pred_no_iran <- fit_def_no_iran %>% 
      vit$summary_brms_nma(keep = c(".trt","def",outcome_label,"bias")) %>% 
      mutate(def = forcats::fct_recode(def,
                                       "unknown" = "NA")) %>%
      rename(deficiency = def) %>%
      filter(.observed) %>% select(-.observed) %>% 
      mutate(bias = forcats::fct_recode(bias,
                                        "unclear" = "some concerns")) %>%
      mutate(bias = forcats::fct_relevel(
        bias,
        "low",
        "unclear",
        "high"
      )) %>% 
      vit$def_plot() +
      ggdist::stat_pointinterval(data = fit_def_no_iran %>% 
                                   vit$prediction_intervals(keep = c(".trt","def",outcome_label, "bias"), only_observed = TRUE) %>%
                                   mutate(def = forcats::fct_recode(def,
                                                                    "unknown" = "NA")) %>%
                                   rename(deficiency = def) %>% 
                                   mutate(bias = forcats::fct_recode(bias,
                                                                     "unclear" = "some concerns"))  %>% 
                                   filter(b == "Placebo") %>%
                                   mutate(.trt = a %>% stringr::str_wrap(width = vit$width)) ,
                                 .width = .95, linewidth = .1, 
                                 color = "red", size = 0.001) +
      ggdist::stat_pointinterval() 
    
    metareg_pred <- bias_summary %>% 
      filter(.observed) %>% select(-.observed) %>% 
      mutate(bias = forcats::fct_relevel(
        bias,
        "low",
        "unclear",
        "high"
      )) %>% 
      vit$def_plot() +
      ggdist::stat_pointinterval(data = bias_reg %>% 
                                   vit$prediction_intervals(
                                     tau_name = "tau_y",
                                     keep = c(".trt","def",outcome_label,"bias.overall",".idx_anx",".idx_dep"),
                                                          index = ".obs",
                                                          resp = resp_label, only_observed = TRUE) %>%
                                   select(-c(.idx_anx, .idx_dep)) %>% 
                                   mutate(!!outcome_label := factor(.data[[outcome_label]], levels = 1:0, labels = c("yes", "no"))) %>%
                                   mutate(!!outcome_label := forcats::fct_relevel(.data[[outcome_label]], "yes","no")) %>% 
                                   mutate(def = forcats::fct_recode(as.factor(def), "unknown" = "NA")) %>%
                                   mutate(def = forcats::fct_relevel(def, "unknown", "no","yes")) %>%
                                   rename(deficiency = def,
                                          bias = bias.overall)  %>% 
                                   mutate(bias = forcats::fct_recode(bias,
                                                                     "unclear" = "some concerns")) %>% 
                                   mutate(bias = forcats::fct_relevel(bias, "low", "unclear", "high"))  %>% 
                                   filter(b == "Placebo") %>%
                                   mutate(.trt = a %>% stringr::str_wrap(width = vit$width) ) ,
                                 .width = .95, linewidth = .1, 
                                 color = "red", size = 0.001) +
      ggdist::stat_pointinterval()
    
    figure_save(def_low_pred,
                outdir = "outputs",
                outcome = outcome,
                filename = "deficiency_model_low_bias_with_prediction_intervals",
                height = 5,
                width = 7.5)
    
    figure_save(def_pred,
                outdir = "outputs",
                outcome = outcome,
                filename = "deficiency_model_all_data_with_prediction_intervals",
                height = 5,
                width = 7.5)
    
    
    figure_save(def_low_pred_no_iran,
                outdir = "outputs",
                outcome = outcome,
                filename = "deficiency_model_low_bias_no_iran_with_prediction_intervals",
                height = 5,
                width = 7.5)
    
    figure_save(def_pred_no_iran,
                outdir = "outputs",
                outcome = outcome,
                filename = "deficiency_model_all_data_no_iran_with_prediction_intervals",
                height = 5,
                width = 7.5)
    
    figure_save(low_pred,
                outdir = "outputs",
                outcome = outcome,
                filename = "simple_model_low_bias_with_prediction_intervals",
                height = 7.5,
                width = 7.5)
    
    figure_save(some_pred,
                outdir = "outputs",
                outcome = outcome,
                filename = "simple_model_some_and_low_bias_with_prediction_intervals",
                height = 5,
                width = 7.5)
    
    figure_save(high_pred,
                outdir = "outputs",
                outcome = outcome,
                filename = "simple_model_all_bias_with_prediction_intervals",
                height = 8,
                width = 7.5)
    
    figure_save(metareg_pred,
                outdir = "outputs",
                outcome = outcome,
                filename = "meta_reg_with_prediction_intervals",
                height = 10,
                width = 7.5)
  }
  
  # duration of follow-up
  {
    followup_plot <- fit_time_low %>% 
      vit$summary_brms_nma(keep = c(".trt","time_cat")) %>% 
      .[.observed == TRUE] %>% 
      mutate(time_cat = factor(time_cat,
                               levels =c("0-8 weeks", "8-12 weeks", "12-26 weeks", "26+ weeks"))) %>%
      filter(.trt == "D") %>% 
      ggplot(aes(x = time_cat, y = value)) +
      ggdist::stat_pointinterval() +
      geom_hline(yintercept = 0, lty = 2, col = "grey50") +
      xlab("Mean duration of follow-up") +
      ylab("SMD") +
      vit$theme_vit()
    
    figure_save(followup_plot,
                outdir = "outputs",
                outcome = outcome,
                filename = "follow_up_duration_distribution",
                height = 5,
                width = 7.5)
    
    followup_plot_all <- fit_time %>% 
      vit$summary_brms_nma(keep = c(".trt","time_cat")) %>% 
      .[.observed == TRUE] %>% 
      mutate(time_cat = factor(time_cat,
                               levels =c("0-8 weeks", "8-12 weeks", "12-26 weeks", "26+ weeks"))) %>%
      filter(.trt == "D") %>% 
      ggplot(aes(x = time_cat, y = value)) +
      ggdist::stat_pointinterval() +
      geom_hline(yintercept = 0, lty = 2, col = "grey50") +
      xlab("Mean duration of follow-up") +
      ylab("SMD") +
      vit$theme_vit()
    
    figure_save(followup_plot_all,
                outdir = "outputs",
                outcome = outcome,
                filename = "follow_up_duration_distribution_all_data",
                height = 5,
                width = 7.5)
  }
  
  
  # pctage plot
  {
    grid_apx <- tidyr::expand_grid(.trt = fit_pcts_approx$prep$network$treatments, p_dep = c(0,1), p_def = c(0,1))
    grid_apx %>% model.matrix(~ .trt : p_dep + .trt : p_def + .trt : p_dep : p_def + .trt, data = .) %>%  as.data.frame() %>% vit$special_clean_names() %>% mutate(.study = fit_pcts_approx$data$.study[1], .se = fit_pcts_approx$data$.se[1], .cov = fit_pcts_approx$data$.cov[1], p_def = grid_apx$p_def, p_dep = grid_apx$p_dep, .trt = grid_apx$.trt) -> newdata_apx
    pct_sum_apx <- vit$summary_brms_nma(fit_pcts_approx, newdata = newdata_apx, keep = c(".trt", 'p_dep', "p_def")) %>% 
      mutate(p_dep = ifelse(p_dep == 1, "yes", "no"),
             p_def = ifelse(p_def == 1, "yes", "no")) %>% 
      rename(!!outcome_label := p_dep, deficiency = p_def)
    pct_sum_apx %>% mutate(bias = "low") %>% 
      filter(.trt == "D") %>% 
      vit$def_plot() %>% 
      figure_save(
        outdir = "outputs",
        outcome = outcome,
        filename = "percentage_plot_approximation",
        height = 5,
        width = 7.5
      )
  }
}

#### adverse events ####
fit_adv <- here::here("outputs", "saved_models", "adverse_events.rds") %>% readRDS()
newdata <- vit$nma_newdata_for_summary(fit_adv)
newdata$final.N <- 1
newdata$.study <- fit_adv$data$.study[1]
newdata$.obs_re <- as.factor(newdata$.obs_re)
adverse_summ <- vit$summary_brms_nma(fit_adv, newdata = newdata)

(plot(adverse_summ %>% 
       filter(.observed) %>% 
       select(-.observed) %>% 
       mutate(value = value)) + 
  theme(plot.margin = margin(5,10,5,5)) + xlab("Odds Ratio") + 
  scale_x_continuous(
    breaks = log(10^(-3:3)),
    labels = scales::trans_format("exp", scales::math_format(.x)))) %>% 
  figure_save(
    outdir = "outputs",
    outcome = "overall",
    filename = "adverse_events_odds_ratios_overall",
    height = 6, width = 6.5)

fit_adv_low <- here::here("outputs", "saved_models", "adverse_events_low_bias.rds") %>% readRDS()
newdata_low <- vit$nma_newdata_for_summary(fit_adv_low)
newdata_low$final.N <- 1
newdata_low$.study <- fit_adv_low$data$.study[1]
newdata_low$.obs_re <- as.factor(newdata_low$.obs_re)
adverse_summ_low <- vit$summary_brms_nma(fit_adv_low, newdata = newdata_low)


(plot(adverse_summ_low %>% 
       filter(.observed) %>% 
       select(-.observed) %>% 
       mutate(value = value)) + 
  theme(plot.margin = margin(5,10,5,5)) + xlab("Odds Ratio") + 
  scale_x_continuous(
    breaks = log(10^(-3:3)),
    labels = scales::trans_format("exp", scales::math_format(.x)))
  ) %>% 
  figure_save(
    outdir = "outputs",
    outcome = "overall",
    filename = "adverse_events_odds_ratios_low_bias_overall",
    height = 4, width = 6.5)


fit_adv_low_ssri <- here::here("outputs", "saved_models", "adverse_events_low_bias_ssri.rds") %>% readRDS()
newdata_ssri <- vit$nma_newdata_for_summary(fit_adv_low_ssri)
newdata_ssri$final.N <- 1
newdata_ssri$.study <- fit_adv_low_ssri$data$.study[1]
newdata_ssri$.obs_re <- as.factor(newdata_ssri$.obs_re)
adverse_summ_ssri <- vit$summary_brms_nma(fit_adv_low_ssri, newdata = newdata_ssri)


(plot(adverse_summ_ssri %>% 
       filter(.observed) %>% 
       select(-.observed) %>% 
       mutate(value = value)) + 
  theme(plot.margin = margin(5,10,5,5)) + xlab("Odds Ratio") + 
  scale_x_continuous(
    breaks = log(10^(-3:3)),
    labels = scales::trans_format("exp", scales::math_format(.x))) )%>% 
  figure_save(
    outdir = "outputs",
    outcome = "overall",
    filename = "adverse_events_odds_ratios_ssri_overall",
    height = 4, width = 6.5)

#### placebo effects over time ####
{
  readRDS(here::here("data","vitamins_full.rds")) %>% 
    filter(intervention == "z.placebo") ->
    placebo_dat
  placebo_dat %>% lme4::lmer(gain.score.mean ~ year + (1|scale), data = .) %>% summary() %>% .$coefficients -> placebo.test
  
  placebo_dat %>% filter(final.outcome.mean >= 0) %>% 
    lme4::lmer(final.outcome.mean ~ year + (1|scale), data = .) %>% summary() %>% .$coefficients -> placebo.test.outcome
  
  (ggplot(placebo_dat, aes(x = year, y = gain.score.mean)) +
      geom_point() +
      geom_abline(intercept = placebo.test["(Intercept)","Estimate"], slope = placebo.test["year","Estimate"], color = "red") +
      # ggtitle("Placebo effects over time") +
      xlab("Year of publication") +
      ylab("Mean gain score in placebo group") +
      vit$theme_vit() +
      scale_x_continuous(limits = c(2000,2025))) %>% 
    figure_save(
      outdir = "outputs",
      outcome = outcome,
      filename = "placebo_gain_score_effects_over_time",
      height = 5,
      width = 7.5
    )
  
    (placebo_dat %>% filter(final.outcome.mean >= 0) %>% 
        ggplot(aes(x = year, y = final.outcome.mean)) +
         geom_point() +
         geom_abline(intercept = placebo.test.outcome["(Intercept)","Estimate"], slope = placebo.test.outcome["year","Estimate"], color = "red") +
         # ggtitle("Placebo effects over time") +
         xlab("Year of publication") +
         ylab("Mean end score in placebo group") +
         vit$theme_vit() ) %>% 
    figure_save(
      outdir = "outputs",
      outcome = outcome,
      filename = "placebo_effects_over_time",
      height = 5,
      width = 7.5
    )
}


#### Add figures for vitamin levels ####

level_low <- readRDS(here::here("outputs", "saved_models", "vitamin_levels_low_bias.rds"))
level_ssri_low <- readRDS(here::here("outputs", "saved_models", "vitamin_levels_low_bias_ssri.rds"))
level_ssri <- readRDS(here::here("outputs", "saved_models", "vitamin_levels_all_ssri.rds"))
level_all <- readRDS(here::here("outputs", "saved_models", "vitamin_levels_all.rds"))

level_plot <- function(lev, nm) {
  for(i in seq_along(lev)) {
    unit <- lev[[i]]$unit
    name <- names(lev)[i] %>% stringr::str_remove("end\\.")
    if(is.null(lev[[i]])) next
    llp <- lev[[i]] %>% 
      vit$summary_brms_nma(keep = c(".trt")) %>% 
      filter(.observed) %>% 
      select(-.observed) %>% 
      plot() +
      xlab(glue::glue("M.D. ({unit})"))
    
    figure_save(llp,
                outdir = "outputs",
                outcome = "overall",
                filename = glue::glue("vitamin_levels_{nm}_{name}"),
                height = 3.5,
                width = 3.25)
  }
  
}

level_plot(level_low, "low_bias")
level_plot(level_ssri_low, "low_bias_ssri")
level_plot(level_ssri, "all_data_ssri")
level_plot(level_all, "all_data")

#### make combined plot easier to read in ####
main_plot_dep <- readRDS(here::here("outputs", "depression", "main_arrow.rds"))
main_plot_anx <- readRDS(here::here("outputs", "anxiety", "main_arrow.rds"))
leg           <- readRDS(here::here("outputs", "depression", "main_plot_legend_only.rds"))

main_plot_cmb <- (wrap_plots(main_plot_dep, main_plot_anx, ncol = 2) / leg) +
  plot_layout(heights = c(10, 2))

figure_save(main_plot_cmb,
            outdir = "outputs",
            outcome = "overall",
            filename = "main_plot_combined",
            height = 8,
            width = 6.5)

net_dep <- readRDS(here::here("outputs","depression", "nma_network_low_bias.rds"))
net_anx <- readRDS(here::here("outputs","anxiety", "nma_network_low_bias.rds"))

main_net <- ((net_dep + theme(plot.margin = margin(r = 80))) + 
               (net_anx + theme(plot.margin = margin(l = 80))) +
                plot_annotation(tag_levels = "A"))

figure_save(main_net,
            outdir = "outputs",
            outcome = "overall",
            filename = "network_combined",
            height = 4,
            width = 13)
