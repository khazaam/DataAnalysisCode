# ==============================================================
#   1. Number of towers placed
#   2. Number of spatial tower clusters
#   3. DBSCAN noise fraction / isolated-placement proportion
#   4. Median interval between tower placements
#   5. Percentage of placements in first half of placement period
#
# Statistics:
#   - participant-level paired A/B observations
#   - delta = B - A
#   - median paired difference
#   - 10,000-resample bootstrap 95% CI
#   - paired Wilcoxon signed-rank test
#       exact = FALSE
#       correct = FALSE
#   - matched-pairs rank-biserial correlation
#   - direction counts
#   - BH/FDR across exactly five predefined Finding 3 tests
#
# Tower DBSCAN:
#   eps = 2.0 m
#   minPts = 2
#
# Figure symbols:
#   square   = ground tower
#   triangle = bird-nest tower (elevated)
#   cross    = enemy spawner
#   circle   = home base
#
# ==============================================================


# --------------------------------------------------------------
# 0. PACKAGES
# --------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(dbscan)
  library(ggplot2)
})


# --------------------------------------------------------------
# 1. CONFIGURATION
# --------------------------------------------------------------

set.seed(42)

DIR_TOWERS <- "data/towers" #Objects that players have placed
DIR_OUT    <- "out_data_analysis"
DIR_FIG    <- file.path(DIR_OUT, "figures")

dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)

DBSCAN_EPS  <- 2.0
DBSCAN_MINP <- 2

BOOT_N <- 10000

PARTICIPANTS <- sprintf("P%02d", 1:15)
CONDITIONS   <- c("A", "B")

# Known missing tower logs from the study.
#My data was missing little bit
KNOWN_MISSING <- tibble(
  participant = c("P01", "P02", "P03", "P14"),
  condition   = c("A",   "A",   "B",   "A"),
  reason      = c(
    "empty log",
    "empty log",
    "empty log",
    "not saved"
  )
)

# Previous tower colours.
COL_TOWER    <- "#4CAF50"
COL_BIRDNEST <- "#1E88E5"
COL_SPAWNER  <- "#D32F2F"
COL_HOMEBASE <- "#FFA000"

# Requested symbols.
SHAPE_TOWER    <- 22
SHAPE_BIRDNEST <- 24
SHAPE_SPAWNER  <- 4
SHAPE_HOMEBASE <- 21


msg <- function(...) {
  cat(
    "[",
    format(Sys.time(), "%H:%M:%S"),
    "]",
    ...,
    "\n"
  )
}


# --------------------------------------------------------------
# 2. EXPECTED RAW FILE PATH
# --------------------------------------------------------------
#
# IMPORTANT:
#
# Do not scan the directory and then try to infer which files
# correspond to which participant.
#
# Construct the expected raw filename directly:
#
#   Game-A-P1-Towers.csv
#   Game-B-P1-Towers.csv
#   ...
#
# This prevents unrelated/duplicate CSV files in data/towers
# from contaminating the analysis.
#
# --------------------------------------------------------------

get_tower_path <- function(participant, condition) {
  
  pnum <- as.integer(
    str_remove(
      participant,
      "^P"
    )
  )
  
  file.path(
    DIR_TOWERS,
    sprintf(
      "Game-%s-P%d-Towers.csv",
      condition,
      pnum
    )
  )
}


# --------------------------------------------------------------
# 3. TOWER FILE LOADER
# --------------------------------------------------------------

load_towers <- function(path) {
  
  # Missing file.
  if (!file.exists(path)) {
    
    return(
      list(
        data = NULL,
        status = "missing"
      )
    )
  }
  
  # Empty file.
  if (file.size(path) == 0) {
    
    return(
      list(
        data = NULL,
        status = "empty"
      )
    )
  }
  
  # Attempt read.
  df <- tryCatch(
    suppressMessages(
      read_csv(
        path,
        show_col_types = FALSE
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(df)) {
    
    return(
      list(
        data = NULL,
        status = "read_error"
      )
    )
  }
  
  if (nrow(df) == 0) {
    
    return(
      list(
        data = NULL,
        status = "empty"
      )
    )
  }
  
  required <- c(
    "TowerType",
    "PositionX",
    "PositionZ"
  )
  
  if (!all(required %in% names(df))) {
    
    return(
      list(
        data = NULL,
        status = "wrong_columns"
      )
    )
  }
  
  # Convert spatial coordinates safely.
  df <- df %>%
    mutate(
      PositionX = suppressWarnings(
        as.numeric(PositionX)
      ),
      PositionZ = suppressWarnings(
        as.numeric(PositionZ)
      )
    )
  
  # Do NOT silently drop the entire file if one coordinate is bad.
  valid_spatial <-
    is.finite(df$PositionX) &
    is.finite(df$PositionZ)
  
  if (!any(valid_spatial)) {
    
    return(
      list(
        data = NULL,
        status = "no_valid_coordinates"
      )
    )
  }
  
  df <- df[
    valid_spatial,
    ,
    drop = FALSE
  ]
  
  # ------------------------------------------------------------
  # Parse timestamp
  # ------------------------------------------------------------
  
  if ("Timestamp" %in% names(df)) {
    
    timestamp_raw <- as.character(
      df$Timestamp
    )
    
    parsed <- suppressWarnings(
      as.POSIXct(
        timestamp_raw,
        format = "%m/%d/%Y %H:%M:%S",
        tz = "UTC"
      )
    )
    
    # Alternative day/month format.
    missing_parse <- is.na(parsed)
    
    if (any(missing_parse)) {
      
      parsed_alt <- suppressWarnings(
        as.POSIXct(
          timestamp_raw[missing_parse],
          format = "%d/%m/%Y %H:%M:%S",
          tz = "UTC"
        )
      )
      
      parsed[missing_parse] <- parsed_alt
    }
    
    df$TimestampParsed <- parsed
    
  } else {
    
    df$TimestampParsed <- as.POSIXct(
      rep(
        NA_character_,
        nrow(df)
      ),
      tz = "UTC"
    )
  }
  
  list(
    data = df,
    status = "ok"
  )
}


# --------------------------------------------------------------
# 4. DBSCAN METRICS
# --------------------------------------------------------------

compute_tower_dbscan <- function(tw) {
  
  if (is.null(tw)) {
    
    return(
      tibble(
        n_clusters = NA_real_,
        noise_fraction = NA_real_
      )
    )
  }
  
  if (nrow(tw) < DBSCAN_MINP) {
    
    return(
      tibble(
        n_clusters = NA_real_,
        noise_fraction = NA_real_
      )
    )
  }
  
  coords <- as.matrix(
    tw[, c(
      "PositionX",
      "PositionZ"
    )]
  )
  
  fit <- dbscan(
    coords,
    eps = DBSCAN_EPS,
    minPts = DBSCAN_MINP
  )
  
  cluster_ids <- unique(
    fit$cluster[
      fit$cluster > 0
    ]
  )
  
  tibble(
    n_clusters =
      length(cluster_ids),
    
    noise_fraction =
      mean(
        fit$cluster == 0
      )
  )
}


# --------------------------------------------------------------
# 5. TIMING METRICS
# --------------------------------------------------------------

compute_timing_metrics <- function(tw) {
  
  if (is.null(tw)) {
    
    return(
      tibble(
        median_gap_s = NA_real_,
        pct_first_half = NA_real_,
        placement_period_s = NA_real_
      )
    )
  }
  
  timed <- tw %>%
    filter(
      !is.na(TimestampParsed)
    ) %>%
    arrange(
      TimestampParsed
    )
  
  # Need at least two valid timestamps.
  if (nrow(timed) < 2) {
    
    return(
      tibble(
        median_gap_s = NA_real_,
        pct_first_half = NA_real_,
        placement_period_s = NA_real_
      )
    )
  }
  
  t0 <- min(
    timed$TimestampParsed,
    na.rm = TRUE
  )
  
  t1 <- max(
    timed$TimestampParsed,
    na.rm = TRUE
  )
  
  duration_s <- as.numeric(
    difftime(
      t1,
      t0,
      units = "secs"
    )
  )
  
  if (
    !is.finite(duration_s) ||
    duration_s <= 0
  ) {
    
    return(
      tibble(
        median_gap_s = NA_real_,
        pct_first_half = NA_real_,
        placement_period_s = duration_s
      )
    )
  }
  
  elapsed_s <- as.numeric(
    difftime(
      timed$TimestampParsed,
      t0,
      units = "secs"
    )
  )
  
  gaps_s <- diff(
    elapsed_s
  )
  
  progress <- elapsed_s /
    duration_s
  
  tibble(
    median_gap_s =
      median(
        gaps_s,
        na.rm = TRUE
      ),
    
    pct_first_half =
      100 *
      mean(
        progress < 0.5,
        na.rm = TRUE
      ),
    
    placement_period_s =
      duration_s
  )
}


# --------------------------------------------------------------
# 6. PROCESS ONE SESSION
# --------------------------------------------------------------

process_session <- function(
    participant,
    condition
) {
  
  path <- get_tower_path(
    participant,
    condition
  )
  
  loaded <- load_towers(
    path
  )
  
  tw <- loaded$data
  status <- loaded$status
  
  known_row <- KNOWN_MISSING %>%
    filter(
      .data$participant == participant,
      .data$condition == condition
    )
  
  known_missing <-
    nrow(known_row) == 1
  
  known_reason <-
    if (known_missing) {
      known_row$reason[1]
    } else {
      NA_character_
    }
  
  if (status != "ok") {
    
    msg(
      participant,
      condition,
      "-",
      status,
      if (known_missing) {
        paste0(
          " [expected: ",
          known_reason,
          "]"
        )
      } else {
        " [UNEXPECTED]"
      }
    )
    
    return(
      tibble(
        participant = participant,
        condition = condition,
        filename = basename(path),
        file_status = status,
        known_missing = known_missing,
        known_missing_reason = known_reason,
        n_towers = NA_real_,
        n_clusters = NA_real_,
        noise_fraction = NA_real_,
        median_gap_s = NA_real_,
        pct_first_half = NA_real_,
        placement_period_s = NA_real_
      )
    )
  }
  
  # Valid file where we expected a missing file should be flagged.
  if (known_missing) {
    
    msg(
      participant,
      condition,
      "- valid file found although marked as known missing"
    )
  }
  
  db <- compute_tower_dbscan(
    tw
  )
  
  timing <- compute_timing_metrics(
    tw
  )
  
  msg(
    participant,
    condition,
    "- valid:",
    nrow(tw),
    "placements"
  )
  
  tibble(
    participant = participant,
    condition = condition,
    filename = basename(path),
    file_status = "ok",
    known_missing = known_missing,
    known_missing_reason = known_reason,
    n_towers = nrow(tw)
  ) %>%
    bind_cols(
      db,
      timing
    )
}


# --------------------------------------------------------------
# 7. PROCESS ALL 30 POSSIBLE SESSIONS
# --------------------------------------------------------------

session_metrics <- map_dfr(
  PARTICIPANTS,
  function(pid) {
    
    map_dfr(
      CONDITIONS,
      function(cond) {
        
        process_session(
          pid,
          cond
        )
      }
    )
  }
)


write_csv(
  session_metrics,
  file.path(
    DIR_OUT,
    "finding3_session_metrics.csv"
  )
)


# --------------------------------------------------------------
# 8. DATA INVENTORY / VALIDATION
# --------------------------------------------------------------

cat(
  "\n\n============================================\n"
)

cat(
  "TOWER DATA INVENTORY\n"
)

cat(
  "============================================\n\n"
)

print(
  session_metrics %>%
    select(
      participant,
      condition,
      filename,
      file_status,
      known_missing,
      n_towers
    ),
  n = Inf
)


n_valid_sessions <- sum(
  session_metrics$file_status == "ok"
)

n_invalid_sessions <- sum(
  session_metrics$file_status != "ok"
)

cat(
  "\nValid tower sessions:",
  n_valid_sessions,
  "\n"
)

cat(
  "Unavailable tower sessions:",
  n_invalid_sessions,
  "\n"
)


# --------------------------------------------------------------
# 9. CHECK OBSERVED MISSINGNESS AGAINST DOCUMENTED MISSINGNESS
# --------------------------------------------------------------

unexpected_invalid <- session_metrics %>%
  filter(
    file_status != "ok",
    !known_missing
  )

unexpected_valid <- session_metrics %>%
  filter(
    file_status == "ok",
    known_missing
  )

if (nrow(unexpected_invalid) > 0) {
  
  cat(
    "\nWARNING: Unexpected unavailable tower sessions:\n"
  )
  
  print(
    unexpected_invalid %>%
      select(
        participant,
        condition,
        file_status
      )
  )
}

if (nrow(unexpected_valid) > 0) {
  
  cat(
    "\nWARNING: Files documented as missing were found valid:\n"
  )
  
  print(
    unexpected_valid %>%
      select(
        participant,
        condition,
        file_status
      )
  )
}


# --------------------------------------------------------------
# 10. CREATE PAIRED DATA
# --------------------------------------------------------------

metric_names <- c(
  "n_towers",
  "n_clusters",
  "noise_fraction",
  "median_gap_s",
  "pct_first_half"
)

wide <- session_metrics %>%
  select(
    participant,
    condition,
    all_of(metric_names)
  ) %>%
  pivot_wider(
    names_from = condition,
    values_from = all_of(metric_names),
    names_glue = "{.value}_{condition}"
  )


write_csv(
  wide,
  file.path(
    DIR_OUT,
    "finding3_paired_metrics.csv"
  )
)


# --------------------------------------------------------------
# 11. REPORT COMPLETE PAIRS
# --------------------------------------------------------------

spatial_complete <- wide %>%
  filter(
    is.finite(n_towers_A),
    is.finite(n_towers_B)
  ) %>%
  pull(
    participant
  )

cat(
  "\nParticipants with valid A/B tower logs:",
  paste(
    spatial_complete,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Number of complete A/B tower pairs:",
  length(spatial_complete),
  "\n"
)

if (length(spatial_complete) != 11) {
  
  warning(
    paste0(
      "Expected 11 complete A/B tower pairs based on documented ",
      "missing logs, but found ",
      length(spatial_complete),
      ". Inspect finding3_session_metrics.csv before interpreting results."
    )
  )
}


# --------------------------------------------------------------
# 12. BOOTSTRAP MEDIAN CI
# --------------------------------------------------------------

bootstrap_median_ci <- function(
    d,
    B = BOOT_N,
    conf = 0.95
) {
  
  d <- d[
    is.finite(d)
  ]
  
  n <- length(d)
  
  if (n == 0) {
    
    return(
      c(
        low = NA_real_,
        high = NA_real_
      )
    )
  }
  
  boot_values <- replicate(
    B,
    median(
      sample(
        d,
        size = n,
        replace = TRUE
      ),
      na.rm = TRUE
    )
  )
  
  alpha <- (
    1 - conf
  ) / 2
  
  ci <- quantile(
    boot_values,
    probs = c(
      alpha,
      1 - alpha
    ),
    na.rm = TRUE,
    names = FALSE,
    type = 7
  )
  
  c(
    low = ci[1],
    high = ci[2]
  )
}


# --------------------------------------------------------------
# 13. MATCHED-PAIRS RANK-BISERIAL CORRELATION
# --------------------------------------------------------------

rank_biserial_paired <- function(d) {
  
  d <- d[
    is.finite(d)
  ]
  
  d_nonzero <- d[
    d != 0
  ]
  
  if (length(d_nonzero) == 0) {
    return(0)
  }
  
  ranks <- rank(
    abs(d_nonzero),
    ties.method = "average"
  )
  
  W_positive <- sum(
    ranks[
      d_nonzero > 0
    ]
  )
  
  W_negative <- sum(
    ranks[
      d_nonzero < 0
    ]
  )
  
  denominator <-
    W_positive +
    W_negative
  
  if (denominator == 0) {
    return(0)
  }
  
  (
    W_positive -
      W_negative
  ) /
    denominator
}


# --------------------------------------------------------------
# 14. ANALYSE ONE METRIC
# --------------------------------------------------------------

analyse_metric <- function(
    wide_data,
    metric
) {
  
  A_name <- paste0(
    metric,
    "_A"
  )
  
  B_name <- paste0(
    metric,
    "_B"
  )
  
  pair_data <- tibble(
    participant =
      wide_data$participant,
    
    A =
      wide_data[[A_name]],
    
    B =
      wide_data[[B_name]]
  ) %>%
    filter(
      is.finite(A),
      is.finite(B)
    ) %>%
    mutate(
      delta = B - A
    )
  
  n <- nrow(
    pair_data
  )
  
  if (n == 0) {
    
    return(
      list(
        summary = tibble(
          metric = metric,
          N = 0,
          median_A = NA_real_,
          median_B = NA_real_,
          median_delta = NA_real_,
          ci_low = NA_real_,
          ci_high = NA_real_,
          p_raw = NA_real_,
          r_rb = NA_real_,
          B_gt_A = 0,
          B_lt_A = 0,
          ties = 0
        ),
        
        deltas = pair_data
      )
    )
  }
  
  d <- pair_data$delta
  
  ci <- bootstrap_median_ci(
    d,
    B = BOOT_N
  )
  
  # If every paired difference is exactly zero,
  # wilcox.test cannot provide a meaningful signed-rank test.
  if (all(d == 0)) {
    
    p_value <- 1
    
  } else {
    
    p_value <- tryCatch(
      suppressWarnings(
        wilcox.test(
          pair_data$B,
          pair_data$A,
          paired = TRUE,
          exact = FALSE,
          correct = FALSE
        )$p.value
      ),
      error = function(e) NA_real_
    )
  }
  
  r_rb <- rank_biserial_paired(
    d
  )
  
  summary_row <- tibble(
    metric = metric,
    
    N = n,
    
    median_A =
      median(
        pair_data$A,
        na.rm = TRUE
      ),
    
    median_B =
      median(
        pair_data$B,
        na.rm = TRUE
      ),
    
    median_delta =
      median(
        d,
        na.rm = TRUE
      ),
    
    ci_low =
      unname(
        ci["low"]
      ),
    
    ci_high =
      unname(
        ci["high"]
      ),
    
    p_raw =
      p_value,
    
    r_rb =
      r_rb,
    
    B_gt_A =
      sum(
        d > 0
      ),
    
    B_lt_A =
      sum(
        d < 0
      ),
    
    ties =
      sum(
        d == 0
      )
  )
  
  delta_rows <- pair_data %>%
    mutate(
      metric = metric
    ) %>%
    select(
      metric,
      participant,
      A,
      B,
      delta
    )
  
  list(
    summary = summary_row,
    deltas = delta_rows
  )
}


# --------------------------------------------------------------
# 15. RUN EXACTLY FIVE PREDEFINED FINDING 3 TESTS
# --------------------------------------------------------------

analysis_list <- lapply(
  metric_names,
  function(m) {
    
    analyse_metric(
      wide,
      m
    )
  }
)

analysis_summary <- bind_rows(
  lapply(
    analysis_list,
    function(x) {
      x$summary
    }
  )
)

participant_deltas <- bind_rows(
  lapply(
    analysis_list,
    function(x) {
      x$deltas
    }
  )
)


# --------------------------------------------------------------
# 16. BH/FDR ACROSS EXACTLY FIVE TESTS
# --------------------------------------------------------------

analysis_summary <- analysis_summary %>%
  mutate(
    p_FDR =
      p.adjust(
        p_raw,
        method = "BH"
      )
  )


# --------------------------------------------------------------
# 17. HUMAN-READABLE LABELS
# --------------------------------------------------------------

metric_labels <- c(
  n_towers =
    "Number of towers",
  
  n_clusters =
    "Tower clusters (DBSCAN)",
  
  noise_fraction =
    "Isolated-placement fraction",
  
  median_gap_s =
    "Median placement interval (s)",
  
  pct_first_half =
    "Placements in first half of placement period (%)"
)

analysis_summary <- analysis_summary %>%
  mutate(
    label =
      unname(
        metric_labels[
          metric
        ]
      )
  ) %>%
  select(
    metric,
    label,
    everything()
  )


# --------------------------------------------------------------
# 18. SAVE ANALYSIS OUTPUTS
# --------------------------------------------------------------

write_csv(
  analysis_summary,
  file.path(
    DIR_OUT,
    "finding3_final_analysis.csv"
  )
)

write_csv(
  participant_deltas,
  file.path(
    DIR_OUT,
    "finding3_participant_deltas.csv"
  )
)


# --------------------------------------------------------------
# 19. PRINT FINAL ANALYSIS
# --------------------------------------------------------------

cat(
  "\n\n============================================\n"
)

cat(
  "FINDING 3 FINAL STATISTICAL ANALYSIS\n"
)

cat(
  "============================================\n\n"
)

print(
  analysis_summary,
  n = Inf
)


# --------------------------------------------------------------
# 20. TEXT SUMMARY
# --------------------------------------------------------------

summary_file <- file.path(
  DIR_OUT,
  "finding3_final_summary.txt"
)

sink(
  summary_file
)

cat(
  "Finding 3 / SQ3 final statistical analysis\n\n"
)

cat(
  "Expected missing logs: P01-A, P02-A, P03-B, P14-A\n"
)

cat(
  "Expected complete A/B tower pairs: 11\n"
)

cat(
  "Observed complete A/B tower pairs:",
  length(spatial_complete),
  "\n\n"
)

cat(
  "Analysis family: 5 predefined paired comparisons.\n"
)

cat(
  "Bootstrap resamples:",
  BOOT_N,
  "\n"
)

cat(
  "Wilcoxon convention: exact = FALSE, correct = FALSE.\n"
)

cat(
  "FDR correction: Benjamini-Hochberg across the 5 predefined Finding 3 comparisons only.\n"
)

cat(
  "DBSCAN eps:",
  DBSCAN_EPS,
  "m\n"
)

cat(
  "DBSCAN minPts:",
  DBSCAN_MINP,
  "\n"
)

cat(
  "Timing reference: tower-placement period from first to last recorded tower timestamp.\n\n"
)

for (
  i in seq_len(
    nrow(
      analysis_summary
    )
  )
) {
  
  r <- analysis_summary[
    i,
  ]
  
  cat(
    sprintf(
      paste0(
        "%s: N=%d, ",
        "median A=%.5f, ",
        "median B=%.5f, ",
        "median delta=%.5f, ",
        "95%% CI [%.5f, %.5f], ",
        "raw p=%.4f, ",
        "p_FDR=%.4f, ",
        "r_rb=%.3f, ",
        "B>A=%d, ",
        "B<A=%d, ",
        "ties=%d\n"
      ),
      
      r$label,
      r$N,
      r$median_A,
      r$median_B,
      r$median_delta,
      r$ci_low,
      r$ci_high,
      r$p_raw,
      r$p_FDR,
      r$r_rb,
      r$B_gt_A,
      r$B_lt_A,
      r$ties
    )
  )
}

sink()


# --------------------------------------------------------------
# 21. LATEX TABLE
# --------------------------------------------------------------

format_num <- function(
    x,
    digits = 3
) {
  
  if (
    length(x) == 0 ||
    is.na(x)
  ) {
    return("--")
  }
  
  formatC(
    x,
    format = "f",
    digits = digits
  )
}


latex_file <- file.path(
  DIR_OUT,
  "table_finding3_final.tex"
)

con <- file(
  latex_file,
  open = "w"
)

writeLines(
  "\\begin{table}[H]",
  con
)

writeLines(
  "\\centering",
  con
)

writeLines(
  paste0(
    "\\caption{Participant-level paired comparisons of tower-placement measures between Configurations A and B. ",
    "$\\Delta$ denotes the median paired difference (B $-$ A). ",
    "Confidence intervals are bootstrap 95\\% CIs based on 10{,}000 resamples of participant-level paired differences. ",
    "$p_{\\mathrm{FDR}}$ values are Benjamini--Hochberg adjusted across the five predefined Finding~3 comparisons. ",
    "$r_{rb}$ denotes the matched-pairs rank-biserial correlation.}"
  ),
  con
)

writeLines(
  "\\label{tab:finding3_towers}",
  con
)

writeLines(
  "\\small",
  con
)

writeLines(
  "\\begin{tabular}{lrrrrrr}",
  con
)

writeLines(
  "\\toprule",
  con
)

writeLines(
  paste0(
    "\\textbf{Metric} & ",
    "\\textbf{$N$} & ",
    "\\textbf{$\\Delta$} & ",
    "\\textbf{95\\% CI} & ",
    "\\textbf{$p_{\\mathrm{FDR}}$} & ",
    "\\textbf{$r_{rb}$} & ",
    "\\textbf{B$>$A / B$<$A} \\\\"
  ),
  con
)

writeLines(
  "\\midrule",
  con
)

for (
  i in seq_len(
    nrow(
      analysis_summary
    )
  )
) {
  
  r <- analysis_summary[
    i,
  ]
  
  direction_text <-
    if (r$ties > 0) {
      
      sprintf(
        "%d / %d / %d",
        r$B_gt_A,
        r$B_lt_A,
        r$ties
      )
      
    } else {
      
      sprintf(
        "%d / %d",
        r$B_gt_A,
        r$B_lt_A
      )
    }
  
  line <- paste0(
    r$label,
    " & ",
    r$N,
    " & ",
    format_num(
      r$median_delta
    ),
    " & [",
    format_num(
      r$ci_low
    ),
    ", ",
    format_num(
      r$ci_high
    ),
    "] & ",
    format_num(
      r$p_FDR
    ),
    " & ",
    format_num(
      r$r_rb,
      2
    ),
    " & ",
    direction_text,
    " \\\\"
  )
  
  writeLines(
    line,
    con
  )
}

writeLines(
  "\\bottomrule",
  con
)

writeLines(
  "\\end{tabular}",
  con
)

writeLines(
  "\\end{table}",
  con
)

close(
  con
)


# ==============================================================
# 22. BUILD FIGURE DATA FROM VALID FILES ONLY
# ==============================================================

figure_rows <- list()

for (pid in PARTICIPANTS) {
  
  for (cond in CONDITIONS) {
    
    path <- get_tower_path(
      pid,
      cond
    )
    
    loaded <- load_towers(
      path
    )
    
    if (loaded$status != "ok") {
      next
    }
    
    tw <- loaded$data %>%
      mutate(
        participant = pid,
        condition = cond
      )
    
    figure_rows[[
      length(figure_rows) + 1
    ]] <- tw
  }
}


tower_plot_data <- bind_rows(
  figure_rows
)


# --------------------------------------------------------------
# 23. STANDARDISE TOWER TYPES
# --------------------------------------------------------------

tower_plot_data <- tower_plot_data %>%
  mutate(
    tower_type_label =
      case_when(
        TowerType == "Tower" ~
          "Ground tower",
        
        TowerType == "BirdNest" ~
          "Bird-nest tower",
        
        TowerType == "EnemySpawner" ~
          "Enemy spawner",
        
        TowerType == "HomeBase" ~
          "Home base",
        
        TRUE ~
          as.character(TowerType)
      ),
    
    tower_type_label =
      factor(
        tower_type_label,
        levels = c(
          "Ground tower",
          "Bird-nest tower",
          "Enemy spawner",
          "Home base"
        )
      )
  )


# --------------------------------------------------------------
# 24. CENTRE FIGURE ON HOME BASE
# --------------------------------------------------------------

centre_towers <- function(df) {
  
  hb <- df %>%
    filter(
      TowerType == "HomeBase"
    )
  
  if (nrow(hb) > 0) {
    
    origin_x <- hb$PositionX[1]
    origin_z <- hb$PositionZ[1]
    
  } else {
    
    origin_x <- mean(
      df$PositionX,
      na.rm = TRUE
    )
    
    origin_z <- mean(
      df$PositionZ,
      na.rm = TRUE
    )
  }
  
  df %>%
    mutate(
      relative_x =
        PositionX -
        origin_x,
      
      relative_z =
        PositionZ -
        origin_z
    )
}


tower_plot_data <- tower_plot_data %>%
  group_by(
    participant,
    condition
  ) %>%
  group_modify(
    ~ centre_towers(.x)
  ) %>%
  ungroup()


# --------------------------------------------------------------
# 25. PANEL LABELS
# --------------------------------------------------------------

panel_levels <- as.vector(
  t(
    outer(
      PARTICIPANTS,
      CONDITIONS,
      paste,
      sep = "-"
    )
  )
)

tower_plot_data <- tower_plot_data %>%
  mutate(
    panel =
      paste(
        participant,
        condition,
        sep = "-"
      ),
    
    panel =
      factor(
        panel,
        levels = panel_levels
      )
  )


# --------------------------------------------------------------
# 26. FIGURE MAPPINGS
# --------------------------------------------------------------

shape_values <- c(
  "Ground tower" =
    SHAPE_TOWER,
  
  "Bird-nest tower" =
    SHAPE_BIRDNEST,
  
  "Enemy spawner" =
    SHAPE_SPAWNER,
  
  "Home base" =
    SHAPE_HOMEBASE
)

fill_values <- c(
  "Ground tower" =
    COL_TOWER,
  
  "Bird-nest tower" =
    COL_BIRDNEST,
  
  "Enemy spawner" =
    COL_SPAWNER,
  
  "Home base" =
    COL_HOMEBASE
)

colour_values <- c(
  "Ground tower" =
    "black",
  
  "Bird-nest tower" =
    "black",
  
  "Enemy spawner" =
    COL_SPAWNER,
  
  "Home base" =
    "black"
)


# --------------------------------------------------------------
# 27. COMPACT FINDING 3 FIGURE
# --------------------------------------------------------------

p_compact <- ggplot(
  tower_plot_data,
  aes(
    x = relative_x,
    y = relative_z
  )
) +
  
  geom_point(
    aes(
      shape = tower_type_label,
      fill = tower_type_label,
      colour = tower_type_label
    ),
    size = 3.2,
    stroke = 0.9,
    alpha = 0.90
  ) +
  
  scale_shape_manual(
    values = shape_values,
    drop = FALSE
  ) +
  
  scale_fill_manual(
    values = fill_values,
    drop = FALSE
  ) +
  
  scale_colour_manual(
    values = colour_values,
    drop = FALSE
  ) +
  
  facet_wrap(
    ~ panel,
    ncol = 6,
    scales = "fixed",
    drop = TRUE
  ) +
  
  coord_equal() +
  
  labs(
    x = "Relative X (m)",
    y = "Relative Z (m)",
    shape = "Placement type",
    fill = "Placement type",
    colour = "Placement type"
  ) +
  
  theme_minimal(
    base_size = 10
  ) +
  
  theme(
    panel.grid.minor =
      element_blank(),
    
    panel.grid.major =
      element_line(
        linewidth = 0.25
      ),
    
    strip.text =
      element_text(
        face = "bold",
        size = 9
      ),
    
    legend.position =
      "bottom",
    
    legend.box =
      "horizontal"
  )


ggsave(
  filename =
    file.path(
      DIR_FIG,
      "finding3_tower_placements_compact.png"
    ),
  
  plot =
    p_compact,
  
  width =
    14,
  
  height =
    12,
  
  dpi =
    300,
  
  bg =
    "white"
)


# --------------------------------------------------------------
# 28. SAVE FIGURE DATA
# --------------------------------------------------------------

write_csv(
  tower_plot_data,
  file.path(
    DIR_OUT,
    "finding3_tower_plot_data.csv"
  )
)


# --------------------------------------------------------------
# 29. SESSION INFO
# --------------------------------------------------------------

session_info_file <- file.path(
  DIR_OUT,
  "finding3_session_info.txt"
)

sink(
  session_info_file
)

cat(
  "Finding 3 tower-placement analysis\n\n"
)

cat(
  "Run:",
  as.character(
    Sys.time()
  ),
  "\n"
)

cat(
  "R:",
  R.version.string,
  "\n\n"
)

cat(
  "Scope:\n"
)

cat(
  "Tower-placement strategy only.\n"
)

cat(
  "No point-cloud or scanning metrics are analysed.\n\n"
)

cat(
  "Documented unavailable tower logs:\n"
)

cat(
  "P01-A: empty\n"
)

cat(
  "P02-A: empty\n"
)

cat(
  "P03-B: empty\n"
)

cat(
  "P14-A: not saved\n\n"
)

cat(
  "Valid tower sessions:",
  n_valid_sessions,
  "\n"
)

cat(
  "Unavailable tower sessions:",
  n_invalid_sessions,
  "\n"
)

cat(
  "Complete A/B tower pairs:",
  length(spatial_complete),
  "\n\n"
)

cat(
  "Predefined analysis family:\n"
)

cat(
  "1. Number of towers placed\n"
)

cat(
  "2. Number of spatial tower clusters\n"
)

cat(
  "3. DBSCAN noise fraction / isolated-placement fraction\n"
)

cat(
  "4. Median interval between tower placements\n"
)

cat(
  "5. Percentage of placements in first half of placement period\n\n"
)

cat(
  "DBSCAN:\n"
)

cat(
  "eps =",
  DBSCAN_EPS,
  "m\n"
)

cat(
  "minPts =",
  DBSCAN_MINP,
  "\n\n"
)

cat(
  "Timing:\n"
)

cat(
  "Placement period begins with first recorded tower timestamp.\n"
)

cat(
  "Placement period ends with last recorded tower timestamp.\n"
)

cat(
  "First-half percentage refers to the first half of this placement period.\n\n"
)

cat(
  "Statistics:\n"
)

cat(
  "Inferential unit = participant-level paired A/B observation\n"
)

cat(
  "Delta = B - A\n"
)

cat(
  "Wilcoxon exact = FALSE, correct = FALSE\n"
)

cat(
  "Bootstrap resamples =",
  BOOT_N,
  "\n"
)

cat(
  "Effect size = matched-pairs rank-biserial correlation\n"
)

cat(
  "BH/FDR family = exactly five predefined Finding 3 comparisons\n\n"
)

cat(
  "Figure:\n"
)

cat(
  "Coordinates translated relative to HomeBase where available.\n"
)

cat(
  "Translation does not change inter-tower distances.\n"
)

cat(
  "Square = ground tower\n"
)

cat(
  "Triangle = bird-nest tower (elevated)\n"
)

cat(
  "Cross = enemy spawner\n"
)

cat(
  "Circle = home base / round-start anchor\n"
)

sink()


# --------------------------------------------------------------
# 30. FINISH
# --------------------------------------------------------------

cat(
  "\n\n============================================\n"
)

cat(
  "FINDING 3 ANALYSIS COMPLETE\n"
)

cat(
  "============================================\n\n"
)

cat(
  "Expected complete pairs: 11\n"
)

cat(
  "Observed complete pairs:",
  length(spatial_complete),
  "\n\n"
)

cat(
  "Outputs:\n"
)

cat(
  "  finding3_final_summary.txt\n"
)

cat(
  "  finding3_final_analysis.csv\n"
)

cat(
  "  finding3_participant_deltas.csv\n"
)

cat(
  "  finding3_session_metrics.csv\n"
)

cat(
  "  finding3_paired_metrics.csv\n"
)

cat(
  "  table_finding3_final.tex\n"
)

cat(
  "  finding3_tower_plot_data.csv\n"
)

cat(
  "  finding3_session_info.txt\n"
)

cat(
  "  figures/finding3_tower_placements_compact.png\n\n"
)

cat(
  "Done.\n"
)