# ==============================================================
#
# Analyses:
#   1. Relative-height LiDAR return distribution
#   2. Device orientation variability
#   3. Temporal scanning continuity
#
# Unit of analysis:
#   Participant x configuration
#
# Inferential unit:
#   Participant-level paired A/B differences
#
# IMPORTANT:
# - No participant-level significance tests.
# - No tower-related measures.
# - Movement/orientation uses first 600 s from first valid timestamp.
# - Relative-height fractions use the COMPLETE valid point cloud.
# - No point-cloud subsampling is performed.
#
# Orientation:
# - Pitch SD = circular SD.
# - Yaw SD   = circular SD.
# - Roll SD  = circular SD.
# - Pitch and yaw changes used for fragmentation are calculated
#   using shortest circular angular differences.
#
# Statistical analysis:
# - Median A and B
# - Median paired difference B - A
# - Bootstrap 95% CI for median paired difference
# - Paired Wilcoxon signed-rank test
#       exact = FALSE
#       correct = FALSE
# - Matched-pairs rank-biserial correlation
# - Participant direction counts
# - BH/FDR correction across the same nine Finding 2 comparisons
#
# ==============================================================


# ==============================================================
# 0. PACKAGES
# ==============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(data.table)
})


# ==============================================================
# 1. SETTINGS
# ==============================================================

set.seed(42)

DIR_POINTCLOUDS <- "data/pointclouds"
DIR_MOVEMENT_BASE <- "scans_all"
DIR_OUT <- "out_finding2"

GAME_WINDOW_S <- 600
FRAG_THRESHOLD_DEG <- 15
BOOT_N <- 10000
WRAP_BOUNDARY_DEG <- 20

dir.create(
  DIR_OUT,
  showWarnings = FALSE,
  recursive = TRUE
)


# ==============================================================
# 2. MESSAGE HELPER
# ==============================================================

msg <- function(...) {
  cat(
    "[",
    format(Sys.time(), "%H:%M:%S"),
    "] ",
    ...,
    "\n",
    sep = ""
  )
}


# ==============================================================
# 3. INDEX POINT-CLOUD FILES
# ==============================================================

index_pc <- function(dir) {
  
  files <- list.files(
    dir,
    pattern = "\\.txt$",
    full.names = TRUE
  )
  
  out <- tibble(
    path = files
  )
  
  out <- out %>%
    mutate(
      file = basename(path)
    )
  
  out <- out %>%
    mutate(
      condition = case_when(
        str_detect(file, "Game-A") ~ "A",
        str_detect(file, "Game-B") ~ "B",
        TRUE ~ NA_character_
      )
    )
  
  out <- out %>%
    mutate(
      participant = str_extract(
        file,
        "P\\d+"
      )
    )
  
  out <- out %>%
    filter(
      !is.na(condition),
      !is.na(participant)
    )
  
  out <- out %>%
    mutate(
      participant = sprintf(
        "P%02d",
        as.integer(
          str_remove(
            participant,
            "P"
          )
        )
      )
    )
  
  return(out)
}


# ==============================================================
# 4. INDEX MOVEMENT FILES
# ==============================================================

index_mv <- function(base) {
  
  files <- list.files(
    base,
    pattern = "\\.txt$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  out <- tibble(
    path = files
  )
  
  out <- out %>%
    mutate(
      file = basename(path),
      game_dir = basename(dirname(path))
    )
  
  out <- out %>%
    mutate(
      condition = case_when(
        
        str_detect(
          game_dir,
          regex(
            "^Game-A$",
            ignore_case = TRUE
          )
        ) ~ "A",
        
        str_detect(
          game_dir,
          regex(
            "^Game-B$",
            ignore_case = TRUE
          )
        ) ~ "B",
        
        TRUE ~ NA_character_
      )
    )
  
  participant_match <- str_match(
    out$file,
    regex(
      "Participant\\s*(\\d+)",
      ignore_case = TRUE
    )
  )
  
  out$participant <- participant_match[, 2]
  
  out <- out %>%
    filter(
      !is.na(condition),
      !is.na(participant)
    )
  
  out <- out %>%
    mutate(
      participant = sprintf(
        "P%02d",
        as.integer(participant)
      )
    )
  
  return(out)
}


# ==============================================================
# 5. LOAD COMPLETE POINT CLOUD
#
# All valid XYZ returns are retained.
# No random/reservoir subsampling is performed.
# ==============================================================

load_pc <- function(path) {
  
  lines <- readLines(
    path,
    warn = FALSE
  )
  
  lines <- trimws(lines)
  
  lines <- lines[
    nchar(lines) > 0
  ]
  
  parts <- strsplit(
    lines,
    "\\s+"
  )
  
  part_lengths <- vapply(
    parts,
    length,
    integer(1)
  )
  
  keep <- which(
    part_lengths >= 3L
  )
  
  if (length(keep) == 0) {
    return(NULL)
  }
  
  numeric_rows <- lapply(
    parts[keep],
    function(v) {
      suppressWarnings(
        as.numeric(
          v[1:3]
        )
      )
    }
  )
  
  mat <- do.call(
    rbind,
    numeric_rows
  )
  
  colnames(mat) <- c(
    "x",
    "y",
    "z"
  )
  
  good <- apply(
    mat,
    1,
    function(r) {
      all(
        is.finite(r)
      )
    }
  )
  
  mat <- mat[
    good,
    ,
    drop = FALSE
  ]
  
  if (nrow(mat) == 0) {
    return(NULL)
  }
  
  return(
    as_tibble(mat)
  )
}


# ==============================================================
# 6. LOAD MOVEMENT / ORIENTATION
# ==============================================================

load_mv <- function(path) {
  
  df <- read.table(
    path,
    header = FALSE
  )
  
  if (ncol(df) < 10) {
    stop(
      paste(
        "Movement file has fewer than 10 columns:",
        path
      )
    )
  }
  
  out <- tibble(
    X = as.numeric(df[[1]]),
    Y = as.numeric(df[[2]]),
    Z = as.numeric(df[[3]]),
    pitch = as.numeric(df[[7]]),
    yaw = as.numeric(df[[8]]),
    roll = as.numeric(df[[9]]),
    t = as.numeric(df[[10]])
  )
  
  out <- out %>%
    filter(
      is.finite(X),
      is.finite(Y),
      is.finite(Z),
      is.finite(pitch),
      is.finite(yaw),
      is.finite(roll),
      is.finite(t)
    )
  
  out <- out %>%
    distinct(
      t,
      .keep_all = TRUE
    )
  
  out <- out %>%
    arrange(t)
  
  return(out)
}


# ==============================================================
# 7. RESTRICT TO FIRST 600 SECONDS
# ==============================================================

first_600s <- function(mv) {
  
  if (is.null(mv)) {
    return(NULL)
  }
  
  if (nrow(mv) == 0) {
    return(mv)
  }
  
  t0 <- min(
    mv$t,
    na.rm = TRUE
  )
  
  out <- mv %>%
    mutate(
      elapsed_s = t - t0
    )
  
  out <- out %>%
    filter(
      elapsed_s >= 0,
      elapsed_s <= GAME_WINDOW_S
    )
  
  return(out)
}


# ==============================================================
# 8. CIRCULAR STANDARD DEVIATION
# ==============================================================

circular_sd_deg <- function(angle_deg) {
  
  angle_deg <- angle_deg[
    is.finite(angle_deg)
  ]
  
  if (length(angle_deg) < 2) {
    return(NA_real_)
  }
  
  angle_rad <- (
    angle_deg *
      pi /
      180
  )
  
  vectors <- exp(
    1i * angle_rad
  )
  
  R_bar <- Mod(
    mean(vectors)
  )
  
  if (!is.finite(R_bar)) {
    return(NA_real_)
  }
  
  if (R_bar <= 0) {
    return(NA_real_)
  }
  
  R_bar <- min(
    R_bar,
    1
  )
  
  circ_sd_rad <- sqrt(
    -2 * log(R_bar)
  )
  
  circ_sd_deg_value <- (
    circ_sd_rad *
      180 /
      pi
  )
  
  return(
    circ_sd_deg_value
  )
}


# ==============================================================
# 9. SHORTEST CIRCULAR ANGULAR DIFFERENCE
# ==============================================================

circular_diff_deg <- function(angle_deg) {
  
  if (length(angle_deg) < 2) {
    return(numeric(0))
  }
  
  angle_rad <- (
    angle_deg *
      pi /
      180
  )
  
  raw_diff <- diff(
    angle_rad
  )
  
  shortest_diff <- atan2(
    sin(raw_diff),
    cos(raw_diff)
  )
  
  shortest_diff_deg <- abs(
    shortest_diff *
      180 /
      pi
  )
  
  return(
    shortest_diff_deg
  )
}


# ==============================================================
# 10. RELATIVE-HEIGHT RETURN DISTRIBUTION
#
# Internal VCI_* names are retained for compatibility.
#
# Session-relative height:
#
#   h = Y - Q05(Y)
#
# Bands:
#   Low-height:           0 <= h < 1 m
#   Intermediate-height:  1 <= h < 3 m
#   Elevated:             h >= 3 m
#
# These are relative-height return fractions rather than direct
# ecological classifications of vegetation strata.
# ==============================================================

compute_vci <- function(pts) {
  
  if (is.null(pts)) {
    return(
      tibble(
        VCI_Ground = NA_real_,
        VCI_Understory = NA_real_,
        VCI_MidCanopy = NA_real_
      )
    )
  }
  
  if (nrow(pts) == 0) {
    return(
      tibble(
        VCI_Ground = NA_real_,
        VCI_Understory = NA_real_,
        VCI_MidCanopy = NA_real_
      )
    )
  }
  
  baseline_y <- quantile(
    pts$y,
    0.05,
    na.rm = TRUE,
    names = FALSE
  )
  
  h <- (
    pts$y -
      baseline_y
  )
  
  band <- case_when(
    
    h >= 0 &
      h < 1 ~ "Ground",
    
    h >= 1 &
      h < 3 ~ "Understory",
    
    h >= 3 ~ "MidCanopy",
    
    TRUE ~ NA_character_
  )
  
  tab <- tibble(
    band = band
  )
  
  tab <- tab %>%
    filter(
      !is.na(band)
    )
  
  if (nrow(tab) == 0) {
    return(
      tibble(
        VCI_Ground = NA_real_,
        VCI_Understory = NA_real_,
        VCI_MidCanopy = NA_real_
      )
    )
  }
  
  tab <- tab %>%
    count(
      band,
      name = "n"
    )
  
  tab <- tab %>%
    mutate(
      frac = n / sum(n)
    )
  
  get_fraction <- function(name) {
    
    x <- tab$frac[
      tab$band == name
    ]
    
    if (length(x) == 0) {
      return(0)
    }
    
    return(
      x[1]
    )
  }
  
  result <- tibble(
    
    VCI_Ground =
      get_fraction(
        "Ground"
      ),
    
    VCI_Understory =
      get_fraction(
        "Understory"
      ),
    
    VCI_MidCanopy =
      get_fraction(
        "MidCanopy"
      )
  )
  
  return(result)
}


# ==============================================================
# 11. ORIENTATION DIAGNOSTICS
# ==============================================================

compute_orientation_diagnostics <- function(mv) {
  
  empty_result <- tibble(
    
    pitch_min = NA_real_,
    pitch_max = NA_real_,
    pitch_range = NA_real_,
    pitch_pct_lt20 = NA_real_,
    pitch_pct_gt340 = NA_real_,
    
    yaw_min = NA_real_,
    yaw_max = NA_real_,
    yaw_range = NA_real_,
    yaw_pct_lt20 = NA_real_,
    yaw_pct_gt340 = NA_real_,
    
    roll_min = NA_real_,
    roll_max = NA_real_,
    roll_range = NA_real_,
    roll_pct_lt20 = NA_real_,
    roll_pct_gt340 = NA_real_
  )
  
  if (is.null(mv)) {
    return(empty_result)
  }
  
  if (nrow(mv) == 0) {
    return(empty_result)
  }
  
  angle_diag <- function(x) {
    
    x <- x[
      is.finite(x)
    ]
    
    if (length(x) == 0) {
      return(
        c(
          min = NA_real_,
          max = NA_real_,
          range = NA_real_,
          pct_lt20 = NA_real_,
          pct_gt340 = NA_real_
        )
      )
    }
    
    x_min <- min(
      x,
      na.rm = TRUE
    )
    
    x_max <- max(
      x,
      na.rm = TRUE
    )
    
    x_range <- (
      x_max -
        x_min
    )
    
    pct_lt20 <- mean(
      x < WRAP_BOUNDARY_DEG,
      na.rm = TRUE
    ) * 100
    
    pct_gt340 <- mean(
      x > (
        360 -
          WRAP_BOUNDARY_DEG
      ),
      na.rm = TRUE
    ) * 100
    
    return(
      c(
        min = x_min,
        max = x_max,
        range = x_range,
        pct_lt20 = pct_lt20,
        pct_gt340 = pct_gt340
      )
    )
  }
  
  pitch_diag <- angle_diag(
    mv$pitch
  )
  
  yaw_diag <- angle_diag(
    mv$yaw
  )
  
  roll_diag <- angle_diag(
    mv$roll
  )
  
  result <- tibble(
    
    pitch_min =
      unname(
        pitch_diag["min"]
      ),
    
    pitch_max =
      unname(
        pitch_diag["max"]
      ),
    
    pitch_range =
      unname(
        pitch_diag["range"]
      ),
    
    pitch_pct_lt20 =
      unname(
        pitch_diag["pct_lt20"]
      ),
    
    pitch_pct_gt340 =
      unname(
        pitch_diag["pct_gt340"]
      ),
    
    yaw_min =
      unname(
        yaw_diag["min"]
      ),
    
    yaw_max =
      unname(
        yaw_diag["max"]
      ),
    
    yaw_range =
      unname(
        yaw_diag["range"]
      ),
    
    yaw_pct_lt20 =
      unname(
        yaw_diag["pct_lt20"]
      ),
    
    yaw_pct_gt340 =
      unname(
        yaw_diag["pct_gt340"]
      ),
    
    roll_min =
      unname(
        roll_diag["min"]
      ),
    
    roll_max =
      unname(
        roll_diag["max"]
      ),
    
    roll_range =
      unname(
        roll_diag["range"]
      ),
    
    roll_pct_lt20 =
      unname(
        roll_diag["pct_lt20"]
      ),
    
    roll_pct_gt340 =
      unname(
        roll_diag["pct_gt340"]
      )
  )
  
  return(result)
}


# ==============================================================
# 12. ORIENTATION VARIABILITY + TEMPORAL CONTINUITY
# ==============================================================

compute_scan <- function(mv) {
  
  if (is.null(mv)) {
    return(
      tibble(
        pitch_sd = NA_real_,
        yaw_sd = NA_real_,
        roll_sd = NA_real_,
        frag_rate = NA_real_,
        longest_run_frames = NA_integer_,
        longest_run_prop = NA_real_,
        longest_run_seconds = NA_real_,
        recorded_duration_s = NA_real_,
        n_frames = 0L
      )
    )
  }
  
  if (nrow(mv) < 2) {
    return(
      tibble(
        pitch_sd = NA_real_,
        yaw_sd = NA_real_,
        roll_sd = NA_real_,
        frag_rate = NA_real_,
        longest_run_frames = NA_integer_,
        longest_run_prop = NA_real_,
        longest_run_seconds = NA_real_,
        recorded_duration_s = NA_real_,
        n_frames = nrow(mv)
      )
    )
  }
  
  
  # ------------------------------------------------------------
  # ORIENTATION VARIABILITY
  # ------------------------------------------------------------
  
  pitch_sd_value <- circular_sd_deg(
    mv$pitch
  )
  
  yaw_sd_value <- circular_sd_deg(
    mv$yaw
  )
  
  roll_sd_value <- circular_sd_deg(
    mv$roll
  )
  
  
  # ------------------------------------------------------------
  # CONSECUTIVE ORIENTATION CHANGES
  # ------------------------------------------------------------
  
  dp <- circular_diff_deg(
    mv$pitch
  )
  
  dy <- circular_diff_deg(
    mv$yaw
  )
  
  
  # ------------------------------------------------------------
  # FRAGMENTATION
  #
  # A transition is an interruption when either pitch or yaw
  # changes by more than 15 degrees.
  # ------------------------------------------------------------
  
  jump <- (
    dp > FRAG_THRESHOLD_DEG
  ) | (
    dy > FRAG_THRESHOLD_DEG
  )
  
  valid_jump <- !is.na(
    jump
  )
  
  if (any(valid_jump)) {
    
    frag_rate_value <- mean(
      jump[
        valid_jump
      ]
    )
    
  } else {
    
    frag_rate_value <- NA_real_
  }
  
  
  # ------------------------------------------------------------
  # CONTINUOUS RUNS
  # ------------------------------------------------------------
  
  jump_for_runs <- jump
  
  jump_for_runs[
    is.na(jump_for_runs)
  ] <- TRUE
  
  run_id <- c(
    1L,
    1L + cumsum(
      jump_for_runs
    )
  )
  
  runs <- tibble(
    run_id = run_id,
    t = mv$t
  )
  
  runs <- runs %>%
    group_by(
      run_id
    ) %>%
    summarise(
      
      run_frames = n(),
      
      run_seconds = if (n() >= 2) {
        
        max(
          t,
          na.rm = TRUE
        ) -
          min(
            t,
            na.rm = TRUE
          )
        
      } else {
        
        0
      },
      
      .groups = "drop"
    )
  
  longest_run_frames_value <- max(
    runs$run_frames,
    na.rm = TRUE
  )
  
  longest_run_prop_value <- (
    longest_run_frames_value /
      nrow(mv)
  )
  
  longest_run_seconds_value <- max(
    runs$run_seconds,
    na.rm = TRUE
  )
  
  recorded_duration_value <- (
    max(
      mv$t,
      na.rm = TRUE
    ) -
      min(
        mv$t,
        na.rm = TRUE
      )
  )
  
  result <- tibble(
    
    pitch_sd =
      pitch_sd_value,
    
    yaw_sd =
      yaw_sd_value,
    
    roll_sd =
      roll_sd_value,
    
    frag_rate =
      frag_rate_value,
    
    longest_run_frames =
      longest_run_frames_value,
    
    longest_run_prop =
      longest_run_prop_value,
    
    longest_run_seconds =
      longest_run_seconds_value,
    
    recorded_duration_s =
      recorded_duration_value,
    
    n_frames =
      nrow(mv)
  )
  
  return(result)
}


# ==============================================================
# 13. INDEX DATA
# ==============================================================

msg(
  "Indexing files..."
)

pc_idx <- index_pc(
  DIR_POINTCLOUDS
)

mv_idx <- index_mv(
  DIR_MOVEMENT_BASE
)

pc_participants <- unique(
  pc_idx$participant
)

mv_participants <- unique(
  mv_idx$participant
)

participants <- sort(
  intersect(
    pc_participants,
    mv_participants
  )
)

msg(
  "Participants: ",
  paste(
    participants,
    collapse = ", "
  )
)


# ==============================================================
# 14. PROCESS PARTICIPANT x CONDITION
# ==============================================================

rows <- list()
diagnostic_rows <- list()

for (pid in participants) {
  
  for (cond in c("A", "B")) {
    
    msg(
      "Processing ",
      pid,
      " ",
      cond
    )
    
    
    # ----------------------------------------------------------
    # POINT CLOUD
    # ----------------------------------------------------------
    
    pc_match <- (
      pc_idx$participant == pid &
        pc_idx$condition == cond
    )
    
    pcp <- pc_idx$path[
      pc_match
    ]
    
    if (length(pcp) == 1) {
      
      pts <- tryCatch(
        
        load_pc(
          pcp
        ),
        
        error = function(e) {
          
          msg(
            "Point-cloud error: ",
            conditionMessage(e)
          )
          
          return(NULL)
        }
      )
      
      vci <- compute_vci(
        pts
      )
      
    } else {
      
      msg(
        "Point-cloud file count = ",
        length(pcp)
      )
      
      vci <- compute_vci(
        NULL
      )
    }
    
    
    # ----------------------------------------------------------
    # MOVEMENT / ORIENTATION
    # ----------------------------------------------------------
    
    mv_match <- (
      mv_idx$participant == pid &
        mv_idx$condition == cond
    )
    
    mvp <- mv_idx$path[
      mv_match
    ]
    
    if (length(mvp) == 1) {
      
      mv_raw <- tryCatch(
        
        load_mv(
          mvp
        ),
        
        error = function(e) {
          
          msg(
            "Movement error: ",
            conditionMessage(e)
          )
          
          return(NULL)
        }
      )
      
    } else {
      
      msg(
        "Movement file count = ",
        length(mvp)
      )
      
      mv_raw <- NULL
    }
    
    
    # ----------------------------------------------------------
    # RAW RECORDING DURATION
    # Diagnostic only.
    # ----------------------------------------------------------
    
    if (
      !is.null(mv_raw) &&
      nrow(mv_raw) >= 2
    ) {
      
      raw_duration <- (
        max(
          mv_raw$t,
          na.rm = TRUE
        ) -
          min(
            mv_raw$t,
            na.rm = TRUE
          )
      )
      
    } else {
      
      raw_duration <- NA_real_
    }
    
    
    # ----------------------------------------------------------
    # FIRST 600 SECONDS
    # ----------------------------------------------------------
    
    mv_game <- first_600s(
      mv_raw
    )
    
    
    # ----------------------------------------------------------
    # COMPUTE SCANNING METRICS
    # ----------------------------------------------------------
    
    scan <- compute_scan(
      mv_game
    )
    
    
    # ----------------------------------------------------------
    # ORIENTATION DIAGNOSTICS
    # ----------------------------------------------------------
    
    orientation_diag <- compute_orientation_diagnostics(
      mv_game
    )
    
    
    # ----------------------------------------------------------
    # STORE MAIN ROW
    # ----------------------------------------------------------
    
    main_row <- bind_cols(
      
      tibble(
        participant = pid,
        condition = cond,
        raw_recorded_duration_s = raw_duration,
        analysis_window_s = GAME_WINDOW_S
      ),
      
      vci,
      
      scan
    )
    
    rows[[length(rows) + 1]] <- main_row
    
    
    # ----------------------------------------------------------
    # STORE DIAGNOSTIC ROW
    # ----------------------------------------------------------
    
    if (is.null(mv_game)) {
      
      diagnostic_n_frames <- 0L
      
    } else {
      
      diagnostic_n_frames <- nrow(
        mv_game
      )
    }
    
    diagnostic_row <- bind_cols(
      
      tibble(
        participant = pid,
        condition = cond,
        n_frames = diagnostic_n_frames
      ),
      
      orientation_diag
    )
    
    diagnostic_rows[[length(diagnostic_rows) + 1]] <- diagnostic_row
  }
}


# ==============================================================
# 15. COMBINE RESULTS
# ==============================================================

sessions <- bind_rows(
  rows
)

orientation_diagnostics <- bind_rows(
  diagnostic_rows
)


# ==============================================================
# 16. SAVE PARTICIPANT-CONDITION METRICS
# ==============================================================

participant_metrics_path <- file.path(
  DIR_OUT,
  "finding2_participant_metrics.csv"
)

fwrite(
  sessions,
  participant_metrics_path
)


# ==============================================================
# 17. SAVE ORIENTATION DIAGNOSTICS
# ==============================================================

orientation_diagnostics_path <- file.path(
  DIR_OUT,
  "finding2_orientation_diagnostics.csv"
)

fwrite(
  orientation_diagnostics,
  orientation_diagnostics_path
)


# ==============================================================
# 18. PRIMARY FINDING 2 METRICS
# ==============================================================

metrics <- c(
  "VCI_Ground",
  "VCI_Understory",
  "VCI_MidCanopy",
  "pitch_sd",
  "yaw_sd",
  "roll_sd",
  "frag_rate",
  "longest_run_prop",
  "longest_run_seconds"
)


# ==============================================================
# 19. CREATE PAIRED WIDE DATA
# ==============================================================

wide <- sessions %>%
  
  select(
    participant,
    condition,
    all_of(metrics)
  ) %>%
  
  pivot_wider(
    names_from = condition,
    values_from = all_of(metrics),
    names_glue = "{.value}_{condition}"
  ) %>%
  
  arrange(
    participant
  )


# ==============================================================
# 20. BOOTSTRAP CI FOR MEDIAN PAIRED DIFFERENCE
# ==============================================================

boot_ci <- function(d) {
  
  d <- d[
    is.finite(d)
  ]
  
  if (length(d) < 2) {
    return(
      c(
        NA_real_,
        NA_real_
      )
    )
  }
  
  bootstrap_medians <- replicate(
    
    BOOT_N,
    
    {
      
      sampled_d <- sample(
        d,
        size = length(d),
        replace = TRUE
      )
      
      median(
        sampled_d,
        na.rm = TRUE
      )
    }
  )
  
  ci <- quantile(
    bootstrap_medians,
    probs = c(
      0.025,
      0.975
    ),
    na.rm = TRUE,
    names = FALSE
  )
  
  return(
    as.numeric(ci)
  )
}


# ==============================================================
# 21. MATCHED-PAIRS RANK-BISERIAL CORRELATION
# ==============================================================

rrb <- function(d) {
  
  d <- d[
    is.finite(d) &
      d != 0
  ]
  
  if (length(d) == 0) {
    return(NA_real_)
  }
  
  ranks <- rank(
    abs(d),
    ties.method = "average"
  )
  
  positive_sum <- sum(
    ranks[
      d > 0
    ]
  )
  
  negative_sum <- sum(
    ranks[
      d < 0
    ]
  )
  
  denominator <- (
    positive_sum +
      negative_sum
  )
  
  if (denominator == 0) {
    return(NA_real_)
  }
  
  result <- (
    positive_sum -
      negative_sum
  ) / denominator
  
  return(result)
}


# ==============================================================
# 22. GROUP-LEVEL PAIRED ANALYSIS
# ==============================================================

analyse <- function(metric) {
  
  A_column <- paste0(
    metric,
    "_A"
  )
  
  B_column <- paste0(
    metric,
    "_B"
  )
  
  A <- wide[[A_column]]
  B <- wide[[B_column]]
  
  complete <- (
    is.finite(A) &
      is.finite(B)
  )
  
  A <- A[
    complete
  ]
  
  B <- B[
    complete
  ]
  
  d <- B - A
  
  N <- length(
    d
  )
  
  if (N < 2) {
    
    result <- tibble(
      metric = metric,
      N = N,
      median_A = NA_real_,
      median_B = NA_real_,
      median_delta = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      wilcoxon_V = NA_real_,
      p = NA_real_,
      r_rb = NA_real_,
      B_higher = sum(d > 0),
      B_lower = sum(d < 0),
      ties = sum(d == 0)
    )
    
    return(result)
  }
  
  
  # ------------------------------------------------------------
  # BOOTSTRAP CI
  # ------------------------------------------------------------
  
  ci <- boot_ci(
    d
  )
  
  
  # ------------------------------------------------------------
  # PAIRED WILCOXON SIGNED-RANK TEST
  #
  # Same convention as Finding 1:
  # exact   = FALSE
  # correct = FALSE
  # ------------------------------------------------------------
  
  wilcox_result <- suppressWarnings(
    
    wilcox.test(
      B,
      A,
      paired = TRUE,
      alternative = "two.sided",
      exact = FALSE,
      correct = FALSE
    )
  )
  
  
  # ------------------------------------------------------------
  # EFFECT SIZE
  # ------------------------------------------------------------
  
  r_rb_value <- rrb(
    d
  )
  
  
  # ------------------------------------------------------------
  # RESULT
  # ------------------------------------------------------------
  
  result <- tibble(
    
    metric = metric,
    
    N = N,
    
    median_A = median(
      A,
      na.rm = TRUE
    ),
    
    median_B = median(
      B,
      na.rm = TRUE
    ),
    
    median_delta = median(
      d,
      na.rm = TRUE
    ),
    
    ci_low = ci[1],
    
    ci_high = ci[2],
    
    wilcoxon_V = unname(
      wilcox_result$statistic
    ),
    
    p = unname(
      wilcox_result$p.value
    ),
    
    r_rb = r_rb_value,
    
    B_higher = sum(
      d > 0
    ),
    
    B_lower = sum(
      d < 0
    ),
    
    ties = sum(
      d == 0
    )
  )
  
  return(result)
}


# ==============================================================
# 23. RUN ALL NINE FINDING 2 COMPARISONS
# ==============================================================

results_list <- lapply(
  metrics,
  analyse
)

results <- bind_rows(
  results_list
)


# ==============================================================
# 24. BH / FDR CORRECTION
#
# Fixed analytical family:
# all nine predefined Finding 2 comparisons.
# ==============================================================

results$p_fdr <- p.adjust(
  results$p,
  method = "BH"
)


# ==============================================================
# 25. SAVE PAIRED ANALYSIS
# ==============================================================

paired_analysis_path <- file.path(
  DIR_OUT,
  "finding2_paired_analysis.csv"
)

fwrite(
  results,
  paired_analysis_path
)


# ==============================================================
# 26. SAVE DIRECTION COUNTS
# ==============================================================

direction_counts <- results %>%
  
  select(
    metric,
    N,
    B_higher,
    B_lower,
    ties
  )

direction_counts_path <- file.path(
  DIR_OUT,
  "finding2_direction_counts.csv"
)

fwrite(
  direction_counts,
  direction_counts_path
)


# ==============================================================
# 27. SAVE PARTICIPANT-LEVEL PAIRED DIFFERENCES
# ==============================================================

participant_delta_rows <- list()

for (metric in metrics) {
  
  A_column <- paste0(
    metric,
    "_A"
  )
  
  B_column <- paste0(
    metric,
    "_B"
  )
  
  delta_table <- tibble(
    
    participant = wide$participant,
    
    metric = metric,
    
    A = wide[[A_column]],
    
    B = wide[[B_column]],
    
    delta = (
      wide[[B_column]] -
        wide[[A_column]]
    )
  )
  
  delta_table <- delta_table %>%
    filter(
      is.finite(A),
      is.finite(B)
    )
  
  delta_table <- delta_table %>%
    mutate(
      direction = case_when(
        delta > 0 ~ "B>A",
        delta < 0 ~ "B<A",
        TRUE ~ "Tie"
      )
    )
  
  participant_delta_rows[[length(participant_delta_rows) + 1]] <- delta_table
}

participant_deltas <- bind_rows(
  participant_delta_rows
)

participant_deltas_path <- file.path(
  DIR_OUT,
  "finding2_participant_deltas.csv"
)

fwrite(
  participant_deltas,
  participant_deltas_path
)


# ==============================================================
# 28. RECORDING DIAGNOSTICS
# ==============================================================

recording_diagnostics <- sessions %>%
  
  select(
    participant,
    condition,
    raw_recorded_duration_s,
    recorded_duration_s,
    n_frames
  ) %>%
  
  arrange(
    participant,
    condition
  )

recording_diagnostics_path <- file.path(
  DIR_OUT,
  "finding2_recording_diagnostics.csv"
)

fwrite(
  recording_diagnostics,
  recording_diagnostics_path
)


# ==============================================================
# 29. PAPER-FACING LABELS
#
# Internal VCI_* variable names are retained for compatibility,
# but the paper-facing terminology is deliberately conservative.
# ==============================================================

labels <- c(
  
  VCI_Ground =
    "Low-height returns (0--1 m)",
  
  VCI_Understory =
    "Intermediate-height returns (1--3 m)",
  
  VCI_MidCanopy =
    "Elevated returns ($\\geq$3 m)",
  
  pitch_sd =
    "Pitch circular SD ($^\\circ$)",
  
  yaw_sd =
    "Yaw circular SD ($^\\circ$)",
  
  roll_sd =
    "Roll circular SD ($^\\circ$)",
  
  frag_rate =
    "Fragmentation rate",
  
  longest_run_prop =
    "Longest-run proportion",
  
  longest_run_seconds =
    "Longest-run duration (s)"
)


# ==============================================================
# 30. FORMATTING HELPERS
# ==============================================================

fmt_p <- function(x) {
  
  if (is.na(x)) {
    return("--")
  }
  
  if (x < 0.001) {
    return("<.001")
  }
  
  sub(
    "^0",
    "",
    sprintf(
      "%.3f",
      x
    )
  )
}


fmt_fixed <- function(x, digits) {
  
  if (is.na(x)) {
    return("--")
  }
  
  sprintf(
    paste0(
      "%.",
      digits,
      "f"
    ),
    x
  )
}


fmt_signed <- function(x, digits) {
  
  if (is.na(x)) {
    return("--")
  }
  
  sprintf(
    paste0(
      "%+.",
      digits,
      "f"
    ),
    x
  )
}


metric_digits <- function(metric) {
  
  if (
    metric %in%
    c(
      "VCI_Ground",
      "VCI_Understory"
    )
  ) {
    return(3)
  }
  
  if (metric == "VCI_MidCanopy") {
    return(5)
  }
  
  if (
    metric %in%
    c(
      "pitch_sd",
      "yaw_sd",
      "roll_sd"
    )
  ) {
    return(2)
  }
  
  if (metric == "frag_rate") {
    return(5)
  }
  
  if (metric == "longest_run_prop") {
    return(3)
  }
  
  if (metric == "longest_run_seconds") {
    return(1)
  }
  
  return(3)
}


# ==============================================================
# 31. GENERATE LATEX TABLE
# ==============================================================

make_latex_row <- function(i) {
  
  metric_name <- results$metric[i]
  
  digits <- metric_digits(
    metric_name
  )
  
  label_text <- unname(
    labels[
      metric_name
    ]
  )
  
  delta_text <- fmt_signed(
    results$median_delta[i],
    digits
  )
  
  ci_text <- paste0(
    "[",
    fmt_fixed(
      results$ci_low[i],
      digits
    ),
    ", ",
    fmt_fixed(
      results$ci_high[i],
      digits
    ),
    "]"
  )
  
  rrb_text <- fmt_signed(
    results$r_rb[i],
    2
  )
  
  if (results$ties[i] > 0) {
    
    direction_text <- paste0(
      results$B_higher[i],
      " / ",
      results$B_lower[i],
      " / ",
      results$ties[i]
    )
    
  } else {
    
    direction_text <- paste0(
      results$B_higher[i],
      " / ",
      results$B_lower[i]
    )
  }
  
  row_text <- paste0(
    label_text,
    " & ",
    results$N[i],
    " & $",
    delta_text,
    "$ & $",
    ci_text,
    "$ & ",
    fmt_p(
      results$p[i]
    ),
    " & ",
    fmt_p(
      results$p_fdr[i]
    ),
    " & $",
    rrb_text,
    "$ & ",
    direction_text,
    " \\\\"
  )
  
  return(row_text)
}


latex_rows <- vapply(
  seq_len(
    nrow(results)
  ),
  make_latex_row,
  character(1)
)


tex <- c(
  
  "\\begin{table}[H]",
  
  "\\centering",
  
  paste0(
    "\\caption{Participant-level paired comparisons of relative-height ",
    "LiDAR return distribution, device-orientation variability, and ",
    "temporal scanning continuity between Configurations A and B. ",
    "$\\Delta$ denotes the median paired difference (B $-$ A). ",
    "Confidence intervals are bootstrap 95\\% CIs based on ",
    "10{,}000 resamples of participant-level paired differences. ",
    "$p$ values are from paired Wilcoxon signed-rank tests. ",
    "$p_{\\mathrm{FDR}}$ values are Benjamini--Hochberg adjusted ",
    "across the nine predefined Finding~2 comparisons. ",
    "$r_{rb}$ denotes matched-pairs rank-biserial correlation. ",
    "Direction counts indicate B$>$A / B$<$A, with ties reported ",
    "as a third value where present.}"
  ),
  
  "\\label{tab:finding2_scanning}",
  
  "\\small",
  
  "\\begin{tabular}{lrrrrrrr}",
  
  "\\toprule",
  
  paste0(
    "\\textbf{Metric} & ",
    "\\textbf{$N$} & ",
    "\\textbf{$\\Delta$} & ",
    "\\textbf{95\\% CI} & ",
    "\\textbf{$p$} & ",
    "\\textbf{$p_{\\mathrm{FDR}}$} & ",
    "\\textbf{$r_{rb}$} & ",
    "\\textbf{Direction} \\\\"
  ),
  
  "\\midrule",
  
  latex_rows[1:3],
  
  "\\addlinespace",
  
  latex_rows[4:6],
  
  "\\addlinespace",
  
  latex_rows[7:9],
  
  "\\bottomrule",
  
  "\\end{tabular}",
  
  "\\end{table}"
)


latex_table_path <- file.path(
  DIR_OUT,
  "finding2_table.tex"
)

writeLines(
  tex,
  latex_table_path
)


# ==============================================================
# 32. FINAL TEXT SUMMARY
# ==============================================================

summary_lines <- c(
  
  "Finding 2 / SQ2 final statistical analysis",
  
  "",
  
  "Analysis family: 9 predefined paired comparisons.",
  
  paste0(
    "Bootstrap resamples: ",
    BOOT_N
  ),
  
  "Wilcoxon convention: exact = FALSE, correct = FALSE.",
  
  paste0(
    "FDR correction: Benjamini-Hochberg across the ",
    "9 predefined Finding 2 comparisons only."
  ),
  
  paste0(
    "Point-cloud sampling: none; all valid XYZ returns ",
    "used for relative-height fractions."
  ),
  
  paste0(
    "Movement/orientation analysis window: first ",
    GAME_WINDOW_S,
    " seconds from first valid timestamp."
  ),
  
  ""
)


for (i in seq_len(nrow(results))) {
  
  metric_name <- results$metric[i]
  
  summary_line <- paste0(
    
    unname(
      labels[
        metric_name
      ]
    ),
    
    ": N=",
    results$N[i],
    
    ", median A=",
    round(
      results$median_A[i],
      5
    ),
    
    ", median B=",
    round(
      results$median_B[i],
      5
    ),
    
    ", median delta=",
    round(
      results$median_delta[i],
      5
    ),
    
    ", 95% CI [",
    round(
      results$ci_low[i],
      5
    ),
    
    ", ",
    round(
      results$ci_high[i],
      5
    ),
    
    "]",
    
    ", raw p=",
    round(
      results$p[i],
      4
    ),
    
    ", p_FDR=",
    round(
      results$p_fdr[i],
      4
    ),
    
    ", r_rb=",
    round(
      results$r_rb[i],
      3
    ),
    
    ", B>A=",
    results$B_higher[i],
    
    ", B<A=",
    results$B_lower[i],
    
    ", ties=",
    results$ties[i]
  )
  
  summary_lines <- c(
    summary_lines,
    summary_line
  )
}


summary_path <- file.path(
  DIR_OUT,
  "finding2_final_summary.txt"
)

writeLines(
  summary_lines,
  summary_path
)


# ==============================================================
# 33. SESSION INFORMATION
# ==============================================================

session_info <- c(
  
  paste0(
    "GAME_WINDOW_S = ",
    GAME_WINDOW_S
  ),
  
  paste0(
    "FRAG_THRESHOLD_DEG = ",
    FRAG_THRESHOLD_DEG
  ),
  
  paste0(
    "BOOT_N = ",
    BOOT_N
  ),
  
  paste0(
    "WRAP_BOUNDARY_DEG = ",
    WRAP_BOUNDARY_DEG
  ),
  
  paste0(
    "POINT_CLOUD_SAMPLING = none ",
    "(all valid XYZ returns used)"
  ),
  
  "",
  
  "Finding 2 scope:",
  
  "Relative-height LiDAR return distribution",
  
  "Orientation variability",
  
  "Temporal scanning continuity",
  
  "",
  
  "Relative-height treatment:",
  
  "Baseline = 5th percentile of point-cloud Y coordinates",
  
  "Low-height returns = 0 <= h < 1 m",
  
  "Intermediate-height returns = 1 <= h < 3 m",
  
  "Elevated returns = h >= 3 m",
  
  "All valid point-cloud returns are used.",
  
  paste0(
    "The point-cloud data have no timestamp, so relative-height ",
    "fractions are not restricted to the 600-s ",
    "movement/orientation window."
  ),
  
  "",
  
  "Orientation treatment:",
  
  "Pitch SD = circular standard deviation",
  
  "Yaw SD = circular standard deviation",
  
  "Roll SD = circular standard deviation",
  
  paste0(
    "Fragmentation pitch differences = ",
    "shortest circular angular difference"
  ),
  
  paste0(
    "Fragmentation yaw differences = ",
    "shortest circular angular difference"
  ),
  
  "",
  
  "No tower-related measures are analysed.",
  
  "No participant-level significance tests are performed.",
  
  paste0(
    "Movement/orientation metrics use the first ",
    GAME_WINDOW_S,
    " s from the first valid timestamp."
  ),
  
  "",
  
  "Statistical treatment:",
  
  paste0(
    "Inferential unit = participant-level paired ",
    "A/B observation"
  ),
  
  "Paired difference = B - A",
  
  paste0(
    "Wilcoxon signed-rank test: ",
    "exact = FALSE, correct = FALSE"
  ),
  
  paste0(
    "Bootstrap 95% CIs use ",
    BOOT_N,
    " participant-level resamples."
  ),
  
  paste0(
    "Effect size = matched-pairs ",
    "rank-biserial correlation"
  ),
  
  "Participant direction counts are reported.",
  
  "",
  
  "Multiple-comparison treatment:",
  
  paste0(
    "BH/FDR correction is applied across all nine ",
    "predefined Finding 2 comparisons."
  ),
  
  paste0(
    "The correction family was not changed based on ",
    "the observed statistical results."
  ),
  
  "",
  
  "Circular-angle treatment:",
  
  paste0(
    "Pitch, yaw, and roll variability use circular ",
    "standard deviations."
  ),
  
  paste0(
    "Pitch and yaw consecutive-frame changes use the ",
    "shortest circular angular difference."
  ),
  
  paste0(
    "This prevents transitions across the 0/360-degree ",
    "boundary from being treated as large angular changes."
  )
)


session_info_path <- file.path(
  DIR_OUT,
  "finding2_session_info.txt"
)

writeLines(
  session_info,
  session_info_path
)


# ==============================================================
# 34. FINAL CONSOLE OUTPUT
# ==============================================================

cat("\n")

cat(
  "============================================================\n"
)

cat(
  "FINAL FINDING 2 / SQ2 ANALYSIS\n"
)

cat(
  "============================================================\n"
)

cat("\n")

cat(
  "PAIRED ANALYSIS\n"
)

cat(
  "---------------\n"
)

print(
  results
)

cat("\n")

cat(
  "DIRECTION COUNTS\n"
)

cat(
  "----------------\n"
)

print(
  direction_counts
)

cat("\n")

cat(
  "FDR CHECK\n"
)

cat(
  "---------\n"
)

cat(
  "Number of predefined comparisons: ",
  nrow(results),
  "\n",
  sep = ""
)

cat(
  "Raw p < .05: ",
  sum(
    results$p < 0.05,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "FDR-adjusted p < .05: ",
  sum(
    results$p_fdr < 0.05,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


# ==============================================================
# 35. OUTPUT FILES
# ==============================================================

cat("\n")

cat(
  "============================================================\n"
)

cat(
  "OUTPUT FILES\n"
)

cat(
  "============================================================\n"
)

cat(
  "Directory: ",
  normalizePath(DIR_OUT),
  "\n\n",
  sep = ""
)

cat(
  "finding2_participant_metrics.csv\n"
)

cat(
  "finding2_paired_analysis.csv\n"
)

cat(
  "finding2_direction_counts.csv\n"
)

cat(
  "finding2_participant_deltas.csv\n"
)

cat(
  "finding2_recording_diagnostics.csv\n"
)

cat(
  "finding2_orientation_diagnostics.csv\n"
)

cat(
  "finding2_table.tex\n"
)

cat(
  "finding2_final_summary.txt\n"
)

cat(
  "finding2_session_info.txt\n"
)

cat(
  "\nFinding 2 final analysis complete.\n"
)