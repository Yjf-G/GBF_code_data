# =====================================================================
# GBF paper - main computation pipeline
# Algorithmic reference: the submitted Supplementary Information (A.1-A.5 methods, B.1-B.2 uncertainty and sensitivity analysis)
#
# Implemented formulas (one-to-one with the numbering of the Supplementary Information):
#   Formula 1   Normalization (0-100, with 100 always the better value; lower bound = 2.5 quantile, upper bound from the decision tree;
#               values better than best are recorded as 100, values worse than worst as 0; for indicators where lower is better the sign is flipped after normalization)
#   Formula 2   Target/goal value = equal-weight arithmetic mean of the indicators it contains (for n=1 it is the basic score at the indicator level)
#   Formula 3   MIS (Mean Index Score) = arithmetic mean of the values of the targets (or goals) (computed year by year)
#   Formula 4   EIS (Evenness Index Score) = area of the radar plot / area of the circle with the same perimeter x 100
#               S = Σ(1/n)π·r²,  L = 2(r_max−r_min) + (2π/n)Σr,  EIS = 400π·S/L²
#   Formula 5   DC  domestic change = |r_{t+1} - r_t| (r is the normalized value)
#   Formula 6   FC  foreign change = Euclidean distance of the year-on-year changes across countries sqrt(Sum_j Dr_j^2)
#   Formula 7   SyE  synergy effect S_t = (1/C) x (Sum_j DC_j - FC) (divided by the number of countries C)
#   Formula 8   SyD  synergy direction sgn(S_t) = sign of the mean signed change across countries compared with the expected direction;
#               +1 if they agree, -1 otherwise
#   Formula 9   CSC country synergy contribution = (DC_j - FC/C)/(S_t x sgn(S_t) x C) x 100%
#               the sum over countries = +100% or -100%; the country-level CSC = mean across indicators and years
#   Formula 10  PGR_sim = 100 x (Y_target year - v_2020) / v_2020
#   Formula 11  PGR_req: an explicit GBF threshold -> the required value at the target year is Y=100;
#               no threshold -> Y = v_2020 + |mean historical change| x n, with Y >= v_2020 (reversing a decline)
#   Formula 12  PGR gap = PGR_req - PGR_sim (positive = shortfall, negative = surplus); aggregated from the target level to the goal level
#   Formula 13  AR = mean(I),  I = 1{PGR_sim >= PGR_req and the bootstrap p of the annual rate of change < 0.05 (the rate of change is significantly != 0)}
#
# Trajectory projection (A.5):
#   high noise (SNR < 1) -> DLM; low noise (SNR >= 1) -> the parametric model family (linear / log /
#   quadratic / segmented), the best form is selected by AICc and extrapolated; the target year is 2030 (indicator/target)
#   or 2050 (goal).
# Bootstrap significance (A.5):
#   the model path after 2020 -> resampling -> linear fit of the annual rate of change -> 10,000 times
#   -> the 2.5%/97.5% quantiles form the 95% CI -> significant if the CI does not contain 0.
# Trajectory classification (A.5): rho = PGR_sim/PGR_req
#   Full achievement   rho >= 1   and a significantly positive rate of change
#   Favourable         0.5 <= rho < 1 and a significantly positive rate of change
#   Mixed              a significant change but 0 < rho < 0.5, or rho >= 1 but not significant
#   Unfavourable       rho <= 0, or not significant and rho < 1
#   Borderline case: the class is unchanged under a +/-0.1 perturbation of the rho thresholds -> confidence flag = stable
# =====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# Make sure the dlm package is installed before running (the DLM implementation of ref 10, Tittensor et al. 2014).
# It is installed automatically when missing; if the installation fails a warning is issued and the built-in equivalent Kalman filter is used (identical convention of the results).
if (!requireNamespace("dlm", quietly = TRUE)) {
  message("[setup] the dlm package is missing, installing it automatically ...")
  ok <- tryCatch({
    install.packages("dlm", repos = "https://cloud.r-project.org", quiet = TRUE)
    requireNamespace("dlm", quietly = TRUE)
  }, error = function(e) FALSE)
  if (isTRUE(ok)) {
    message("[setup] the dlm package was installed successfully")
  } else {
    warning("[setup] installation of the dlm package failed, the built-in equivalent DLM (local linear trend Kalman filter) will be used")
  }
}

options(scipen = 10, dplyr.summarise.inform = FALSE)

# ---------------------------------------------------------------------
# Section 0. Configuration
# ---------------------------------------------------------------------

# GBF targets with an explicit numerical threshold (Formula 11, branch 1: the required value at the target year is Y = 100).
# Based on the explicit quantitative requirements in the text of GBF decision 15/4:
#   T2  >=30% of ecosystems under effective restoration | T3  >=30% under effective protection | T6  introduction rate of invasive alien species down by >=50%
#   T7  nutrient and pesticide risk halved     | T16 food waste halved      | T18 harmful incentives reduced by >=$500bn per year
#   T19  mobilising >=$200bn per year
# The remaining targets follow branch 2 of Formula 11 (extrapolation of the mean historical change).
THRESHOLD_TARGETS <- c("T2", "T3", "T6", "T7", "T16", "T18", "T19")

TARGET_YEAR_TARGET <- 2030L   # target year of the indicators and the targets
TARGET_YEAR_GOAL   <- 2050L   # target year of the goals

# Regional indicators (27 RLI food/medicinal species; 28 RLI traded species; 29 RLI impact of invasive alien species)
# The raw data are continent/region series (meta.csv countries = "regions/global").
# They are loaded only for the global/regional analyses; they take no part in the country-level analyses (they do not enter the coverage denominator,
# the aggregation of target values, the projections, the AR, the synergy computations or any other country-level statistic).
REGIONAL_ONLY_INDICATORS <- c(27L, 28L, 29L)

# Unit conversion factors for the indicators: used when the unit of the raw data differs from the unit of the normalization bounds in meta.csv
#   (register a factor here if a similar case arises again; the values are converted on reading)
#   18 Extent of natural ecosystems: the raw data are a proportion (0-1) while the meta bounds are percentages
#   (13.57025/99.35725) - on 2026-10-06 the raw data were corrected directly by a factor of 100
#   (data/origin_data/indicator18_ExtNatEco.csv);
#   therefore no conversion factor has to be registered here.
#   so no conversion factor has to be registered here.
UNIT_SCALE <- c()

BASELINE_YEARS <- 2011:2020   # CBD reference baseline period (the four baseline parameters, Supplementary Table 3)
TRAJ_YEARS     <- 2011:2024   # annual series used for the trajectory analysis (SI B.2.7, 13 annual observations)
N_BOOT_RATE    <- 10000L      # A.5: 10,000 bootstrap resamples
# Noise criterion (Tittensor et al. 2014, the basis for the ref 10/11 split):
#   signal-to-noise ratio SNR = process variance / observation variance (MLE of a random walk plus noise model)
#   SNR >= 1 -> low noise -> PMA;  SNR < 1 -> high noise -> DLM
SNR_DLM_THRESHOLD <- 1

RHO_FULL <- 1.0               # A.5 rho threshold
RHO_FAV  <- 0.5
RHO_PERTURB <- 0.1            # perturbation used for the borderline sensitivity

# ---------------------------------------------------------------------
# Section 1. Metadata / inputs
# ---------------------------------------------------------------------

# Clean the bound strings: "0.6*", "100 +", "2.21764E+11" -> numeric
# Note: a trailing "+" means "at least" and is not part of the number; only the +/- inside an exponent is kept
parse_bound <- function(x) {
  x <- trimws(as.character(x))
  m <- regmatches(x, regexpr("[0-9]+(\\.[0-9]+)?([eE][+-]?[0-9]+)?", x))
  out <- rep(NA_real_, length(x))
  ok <- !is.na(m) & m != ""
  out[ok] <- suppressWarnings(as.numeric(m[ok]))
  out
}

build_ind_meta <- function(meta_path = "data/meta.csv",
                           i2t_path  = "data/indicator_to_target.csv") {
  if (!file.exists(meta_path)) stop("missing ", meta_path)
  meta <- read_csv(meta_path, show_col_types = FALSE) %>%
    rename(indicator = indicator_id) %>%
    mutate(
      indicator = as.integer(indicator),
      worst = parse_bound(worst),
      best  = parse_bound(best),
      # "+" means higher is better; "−"/"-" means lower is better (the sign is flipped after normalization, 100 is always the better value)
      direction = ifelse(grepl("^\\s*[+]", as.character(direction)), "higher", "lower"),
      used = used %in% c("Y", "y", "TRUE", "true", TRUE)
    )
  if (any((is.na(meta$worst) | is.na(meta$best)) & meta$used)) {
    warning("meta.csv contains worst/best bounds that cannot be parsed; the A.3 decision tree fallback will be used")
  }
  if (file.exists(i2t_path)) {
    i2t <- read_csv(i2t_path, show_col_types = FALSE) %>%
      mutate(indicator = as.integer(indicator_id),
             target = as.character(target))
    # The target column of the mapping table mixes two kinds of values: "A"–"D" are goal mappings, "1"–"23" are target mappings
    i2t <- i2t %>%
      mutate(layer = ifelse(target %in% c("A", "B", "C", "D"), "goal", "target"),
             layer_id = ifelse(layer == "goal", target, paste0("T", target)))
    meta <- meta %>%
      left_join(i2t %>% select(indicator, layer, layer_id), by = "indicator")
  }
  meta
}

# Locate the indicator file: data/origin_data/indicator<i>(_code).csv
find_indicator_file <- function(id, dir = "data/origin_data") {
  fs <- list.files(dir, pattern = paste0("^indicator", id, "(_|\\.).*\\.csv$"),
                   full.names = TRUE)
  if (length(fs) == 0) return(NA_character_)
  fs[1]
}

# Wide CSVs (year columns) and long CSVs (value/year columns) are read in the same way
# Returns the long format: country, year, indicator, value (NA is kept, before imputation)
read_indicator_raw <- function(path, id) {
  df <- suppressWarnings(read_csv(path, show_col_types = FALSE, progress = FALSE))
  names(df) <- trimws(names(df))
  if ("GeoAreaName" %in% names(df)) {
    df <- df %>% rename(country = GeoAreaName)
  } else {
    names(df)[1] <- "country"
  }
  if (all(c("year", "value") %in% names(df))) {
    out <- df %>%
      transmute(country = as.character(country),
                year = suppressWarnings(as.integer(year)),
                value = suppressWarnings(as.numeric(value)),
                indicator = as.integer(id))
  } else {
    year_cols <- grep("^[0-9]{4}$", names(df), value = TRUE)
    if (length(year_cols) == 0) stop("no year columns and no value/year columns found in ", path, "; the file must contain year columns or value/year columns")
    out <- df %>%
      select(country, all_of(year_cols)) %>%
      # The year columns may contain non-numeric text -> coerce to numeric before pivoting
      mutate(across(all_of(year_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
      pivot_longer(-country, names_to = "year", values_to = "value") %>%
      transmute(country = as.character(country),
                year = as.integer(year),
                value = suppressWarnings(as.numeric(value)),
                indicator = as.integer(id))
  }
  out %>%
    filter(!is.na(year), year >= 1900, year <= 2100, !is.na(country)) %>%
    # Several rows for the same (country, year) -> take the mean, so that the x values do not repeat during imputation
    group_by(country, year) %>%
    summarise(value = mean(value, na.rm = TRUE), indicator = first(indicator),
              .groups = "drop") %>%
    mutate(value = ifelse(is.nan(value), NA_real_, value))
}

load_all_raw <- function(meta = NULL, dir = "data/origin_data") {
  if (is.null(meta)) meta <- build_ind_meta()
  # The main analysis uses only the 44 indicators with used=Y in meta.csv (not all 53 candidates)
  # All historical years are kept: the imputation only fills 2011-2024, but the projection fit needs the full history to identify the pattern
  ids <- unique(meta$indicator[meta$used])
  map_dfr(ids, function(i) {
    f <- find_indicator_file(i, dir)
    if (is.na(f)) stop("no data file found for indicator ", i, " (", dir, ")")
    out <- read_indicator_raw(f, i)
    sc <- UNIT_SCALE[as.character(i)]
    if (length(sc) == 1 && !is.na(sc)) {
      out <- out %>% mutate(value = value * sc)
    }
    out
  })
}

# Country / group metadata
#   preferred: data/geoclass.xlsx (UN M49 country classification, 288 rows x 19 columns)
#     with Continent / Subcontinent / Region / country groups1-2 / ISO /
#     EU / Eurozone / developed-developing / LDC / LLDC / SIDS / transition countries
#   fallback: aggregation of the geographic group columns carried by the indicator files (when geoclass.xlsx is missing)
build_country_meta <- function(path = "data/geoclass.xlsx",
                               dir = "data/origin_data") {
  if (file.exists(path) && requireNamespace("readxl", quietly = TRUE)) {
    g <- readxl::read_excel(path)
    names(g) <- trimws(names(g))
    col_country <- "Geographical region"
    pick <- function(...) {
      nms <- c(...)
      hit <- intersect(nms, names(g))
      if (length(hit) == 0) return(rep(NA_character_, nrow(g)))
      as.character(g[[hit[1]]])
    }
    flag <- function(x, val) !is.na(x) & x == val
    out <- tibble(
      country     = as.character(g[[col_country]]),
      iso3        = pick("ISO CODES_3"),
      continent   = pick("Continent"),
      subcontinent = pick("Subcontinent"),
      region      = pick("Region"),
      group1      = pick("country groups1"),
      group2      = pick("country groups2"),
      is_nation   = flag(pick("United Nations"), "nation"),
      development = case_when(
        flag(pick("Developed and developing regions-developing"), "dl") ~ "developing",
        flag(pick("Developed and developing regions-developed"), "dd") ~ "developed",
        TRUE ~ NA_character_
      ),
      sids = flag(pick("Small island developing States"), "sid"),
      ldc  = flag(pick("Least developed countries"), "ldc"),
      llc  = flag(pick("Landlocked developing countries"), "llc"),
      oecd = flag(pick("country groups1"), "OECD"),
      eu   = flag(pick("European Union"), "EU")
    ) %>%
      filter(!is.na(country), country != "")
    return(out)
  }
  # Fallback: aggregate the geographic columns of the indicator files
  geo_cols <- c("Continent", "Subcontinent", "Region", "GeoAreaType",
                "Developing_regions", "Developed_regions", "Least_dpe",
                "Landlocked_dpi", "Small_island_dpi", "group1")
  fs <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
  parts <- map(fs, function(f) {
    df <- suppressWarnings(read_csv(f, show_col_types = FALSE, progress = FALSE,
                                    n_max = 1e5))
    names(df) <- trimws(names(df))
    if (!"GeoAreaName" %in% names(df)) return(NULL)
    keep <- intersect(geo_cols, names(df))
    if (length(keep) == 0) return(NULL)
    df %>% select(country = GeoAreaName, all_of(keep)) %>% distinct()
  })
  meta <- bind_rows(parts) %>%
    filter(!is.na(country), country != "") %>%
    group_by(country) %>%
    summarise(across(everything(), ~ {
      x <- .x[!.x %in% c("0", "", NA)]
      if (length(x) == 0) NA_character_ else x[1]
    }), .groups = "drop")
  meta %>%
    mutate(
      iso3 = NA_character_,
      continent = Continent,
      subcontinent = Subcontinent,
      region = Region,
      oecd = (group1 == "OECD"),          # oecd is computed first and then overwritten by group1
      group1 = NA_character_, group2 = NA_character_,
      is_nation = GeoAreaType == "nation",
      development = case_when(
        Developing_regions == "dl" ~ "developing",
        Developed_regions == "dd" ~ "developed",
        TRUE ~ NA_character_
      ),
      sids = Small_island_dpi == "sid",
      ldc  = Least_dpe == "ldc",
      llc  = Landlocked_dpi == "llc",
      eu = NA
    ) %>%
    select(country, iso3, continent, subcontinent, region, group1, group2,
           is_nation, development, sids, ldc, llc, oecd, eu)
}

# Country universe: only the 193 countries with United Nations == "nation" in geoclass.xlsx
#   non-country entities such as regions/continents/"Global average" enter no country-level statistic
#   (coverage denominator, C in Formula 13 (AR), C in Formula 7 (SyE), n_countries, rankings, etc.)
nation_universe <- function(country_meta) {
  country_meta %>% filter(is_nation %in% TRUE) %>% pull(country)
}

# Filter the data rows to the country universe; returns the filtered data plus a report
filter_to_nations <- function(df, country_meta, exclude_indicators = integer(0),
                             verbose = TRUE) {
  nat <- country_meta %>%
    filter(is_nation %in% TRUE) %>%
    transmute(key = tolower(country))
  keys <- tolower(unique(df$country))
  keep <- unique(nat$key)
  dropped <- setdiff(keys, keep)
  out <- df %>%
    mutate(.key = tolower(country)) %>%
    filter(.key %in% keep) %>%
    select(-.key)
  if (length(exclude_indicators) > 0) {
    out <- out %>% filter(!indicator %in% exclude_indicators)
  }
  if (verbose) {
    message(sprintf("[universe] %d country entities (geoclass United Nations=nation); %d non-country rows removed",
                    length(intersect(keys, keep)), length(dropped)))
    if (length(exclude_indicators) > 0) {
      message("[universe]   regional indicators take no part in the country-level analyses (loaded for the global/regional analyses only): ",
              paste(exclude_indicators, collapse = ", "))
    }
  }
  out
}

# ---------------------------------------------------------------------
# Section 2. Imputation of missing values (A.2)
#   0) First create the rows for the missing years: for wide tables without a year column (e.g. indicator 18 every five years, 39 every ten years,
#      26/30 every other year, etc.) there is no corresponding row in the panel at all, so the annual grid of the
#      impute_years window has to be completed first for the imputation to have anything to work on
#   1) Missing values at the start or end of a numeric series -> forward / backward fill
#   2) Missing values inside a numeric series -> linear interpolation
#   3) Categorical series -> most recent previous value (LOCF)
#   Anchor = all the observations the series already has (including the history outside the 2011-2024 window), not only
#   the points inside the period to be imputed (A.2; the imputation of missing years relies especially on the full history)
#   Retention rule for countries: indicator coverage >= 80%; countries with more than 20% of the indicators missing are excluded from MIS/EIS/rankings
# ---------------------------------------------------------------------

# Safe mean wrapper: returns NA instead of NaN when everything is NA (to keep NaN out of the downstream projections/tables)
mean_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else mean(x)
}

# Complete the rows for the missing years of the impute_years window for every (indicator, country)
# (for the wide-table indicators without year columns, e.g. 18/39/16/22/23/26/30/36/37).
# Historical years get no rows: outside the window the values are not imputed and stay NA, following A.2.
expand_impute_grid <- function(raw_df, impute_years = TRAJ_YEARS) {
  raw_df %>%
    group_by(indicator, country) %>%
    group_modify(~ {
      full <- sort(unique(c(.x$year, impute_years)))
      tibble(year = full) %>%
        left_join(.x %>% select(year, value), by = "year")
    }) %>%
    ungroup() %>%
    select(country, year, indicator, value)
}

# The imputation only fills the missing values inside the impute_years window (2011-2024);
# missing values outside the window stay NA (the historical observations are kept as they are, for the projection fit).
# The fill anchors refer to all the original observations of the series (including those outside the window), consistent with A.2.
impute_series <- function(year, value, categorical = FALSE, impute_years = NULL) {
  ord <- order(year)
  y <- year[ord]; v <- value[ord]
  was_na <- is.na(v)
  if (all(was_na)) return(v)
  if (categorical) {
    for (i in seq_along(v)) if (is.na(v[i]) && i > 1) v[i] <- v[i - 1]
    for (i in rev(seq_along(v))) if (is.na(v[i]) && i < length(v)) v[i] <- v[i + 1]
    if (!is.null(impute_years)) v[was_na & !(y %in% impute_years)] <- NA
    return(v)
  }
  obs <- which(!is.na(v))
  if (length(obs) == 0) return(v)
  # Endpoints: a missing tail is forward filled; a missing head is backward filled
  if (obs[1] > 1) v[1:(obs[1] - 1)] <- v[obs[1]]
  if (obs[length(obs)] < length(v)) v[(obs[length(obs)] + 1):length(v)] <- v[obs[length(obs)]]
  # Interior: linear interpolation (the anchors are all the observations, including those outside 2011-2024)
  nas <- which(is.na(v))
  if (length(nas) > 0) {
    v[nas] <- approx(y[!is.na(v)], v[!is.na(v)], xout = y[nas], rule = 1)$y
  }
  # Only the filled values inside the window are kept
  if (!is.null(impute_years)) v[was_na & !(y %in% impute_years)] <- NA
  v
}

# The imputation applies only to the 44 indicators with used=Y and only to the missing values inside the 2011-2024 window
# The inputs keep all historical years; only the missing values inside the 2011-2024 window are imputed
# Procedure: (1) create the rows for the missing years (expand_impute_grid) -> (2) impute using the full history as anchors
impute_all <- function(raw_df, categorical_ids = integer(0),
                       impute_years = TRAJ_YEARS, used_ids = NULL) {
  d <- raw_df
  if (!is.null(used_ids)) d <- d %>% filter(indicator %in% used_ids)
  d %>%
    expand_impute_grid(impute_years) %>%
    arrange(indicator, country, year) %>%
    group_by(indicator, country) %>%
    mutate(value = {
      first(indicator) %in% categorical_ids -> cat
      impute_series(year, value, categorical = cat, impute_years = impute_years)
    }) %>%
    ungroup()
}

# Cell-by-cell imputation detail (every imputation detail, for imputation_detail.csv / Supplementary Table 9)
#   status:
#     observed                   original observation, kept as it is
#     imputed                    a filled value imputed inside the window and retained
#     missing_no_anchor          no observation anchor inside the window, cannot be imputed
#     discarded_outside_window   a fill computed outside the window but not retained, following A.2 (kept for the record)
#   method: linear_interpolation / forward_fill / backward_fill / none
#   anchor_*: the basis of the imputation (all existing observations, including the history outside the window), the nearest previous / next anchor
#   value_raw the original value | value_imputed the imputed value | value_final the value that enters the panel
imputation_detail <- function(raw_df, imputed_df, impute_years = TRAJ_YEARS) {
  # Final value inside the window (including the rows created by expand)
  imp_win <- imputed_df %>%
    filter(year %in% impute_years) %>%
    select(indicator, country, year, value_final = value)
  raw_key <- raw_df %>%
    distinct(indicator, country, year, .keep_all = TRUE) %>%
    select(indicator, country, year, value_raw = value)
  det <- imp_win %>%
    left_join(raw_key, by = c("indicator", "country", "year")) %>%
    mutate(grid_added = is.na(value_raw) &
             !paste(indicator, country, year) %in%
             paste(raw_key$indicator, raw_key$country, raw_key$year))
  # Outside the window: the fills that were computed but discarded following A.2 are also recorded (raw missing, no created row, no final value)
  raw_out <- raw_key %>% filter(!year %in% impute_years, is.na(value_raw))
  if (nrow(raw_out) > 0) {
    det_out <- raw_out %>%
      select(indicator, country, year, value_raw) %>%
      mutate(value_final = NA_real_, grid_added = FALSE)
    det <- bind_rows(det, det_out)
  }
  # All existing observation anchors
  obs <- raw_df %>%
    filter(!is.na(value)) %>%
    group_by(indicator, country) %>%
    summarise(obs_year = list(year), obs_val = list(value), .groups = "drop")
  det <- det %>% left_join(obs, by = c("indicator", "country"))

  classify_one <- function(y, v_raw, v_fin, oy, ov) {
    n <- length(y)
    status <- character(n); method <- character(n)
    v_imp <- rep(NA_real_, n)
    a1y <- a1v <- a2y <- a2v <- rep(NA_real_, n)
    for (i in seq_len(n)) {
      if (!is.na(v_raw[i])) {
        status[i] <- "observed"; method[i] <- "none"
      } else {
        # the imputed value of the cell: constant extrapolation at the endpoints + linear interpolation inside = approx(rule = 2)
        imp_val <- if (length(oy) == 0) {
          NA_real_
        } else if (length(oy) == 1) {
          ov[1]
        } else {
          suppressWarnings(approx(oy, ov, xout = y[i], rule = 2)$y)
        }
        v_imp[i] <- imp_val
        if (is.na(imp_val)) {
          status[i] <- "missing_no_anchor"
          method[i] <- "none"
        } else {
          lo <- oy[oy < y[i]]; hi <- oy[oy > y[i]]
          if (length(lo) > 0 && length(hi) > 0) {
            method[i] <- "linear_interpolation"
            a1y[i] <- max(lo); a1v[i] <- ov[which.max(lo)]
            a2y[i] <- min(hi); a2v[i] <- ov[which.min(hi)]
          } else if (length(lo) > 0) {
            method[i] <- "forward_fill"
            a1y[i] <- max(lo); a1v[i] <- ov[which.max(lo)]
          } else {
            method[i] <- "backward_fill"
            a2y[i] <- min(hi); a2v[i] <- ov[which.min(hi)]
          }
          status[i] <- if (!is.na(v_fin[i])) {
            "imputed"
          } else {
            "discarded_outside_window"
          }
        }
      }
    }
    tibble(status = status, method = method, value_imputed = v_imp,
           anchor_year_before = a1y, anchor_value_before = a1v,
           anchor_year_after = a2y, anchor_value_after = a2v)
  }

  det %>%
    group_by(indicator, country) %>%
    group_modify(~ {
      oy <- if (length(.x$obs_year) == 0 || is.null(.x$obs_year[[1]])) {
        numeric(0)
      } else .x$obs_year[[1]]
      ov <- if (length(.x$obs_val) == 0 || is.null(.x$obs_val[[1]])) {
        numeric(0)
      } else .x$obs_val[[1]]
      bind_cols(
        .x %>% select(year, value_raw, value_final, grid_added),
        classify_one(.x$year, .x$value_raw, .x$value_final, oy, ov)
      )
    }) %>%
    ungroup() %>%
    select(indicator, country, year, status, method, value_raw, value_imputed,
           value_final, anchor_year_before, anchor_value_before,
           anchor_year_after, anchor_value_after, grid_added) %>%
    arrange(indicator, country, year)
}

# Country coverage: countries with >=80% indicator coverage are kept; more than 20% missing excludes them from MIS/EIS/rankings (A.2)
country_indicator_coverage <- function(imputed_df, meta) {
  # Denominator = the indicators actually available at the country level (the regional indicators 27-29 have no nation rows and do not enter the denominator)
  ids_used <- intersect(unique(meta$indicator[meta$used]),
                        unique(imputed_df$indicator))
  imputed_df %>%
    filter(indicator %in% ids_used) %>%
    group_by(country) %>%
    summarise(
      n_ind_observed = n_distinct(indicator[!is.na(value)]),
      n_ind_total    = length(ids_used),
      share_missing  = 1 - n_ind_observed / n_ind_total,
      .groups = "drop"
    ) %>%
    mutate(keep_index = share_missing <= 0.20)
}

# ---------------------------------------------------------------------
# Section 3. Normalization - Formula 1
#   x' = (x - worst)/(best - worst) x 100 (worst/best are the bounds of the 0 and 100 points;
#   for indicators where lower is better worst > best, so the sign is flipped naturally and after normalization 100 is always the better value)
#   better than best -> 100; worse than worst -> 0
#   Derivation of the bounds (A.3 decision tree, used as a fallback when meta is missing):
#     lower bound = the 2.5 quantile of the observed distribution;
#     upper bound = absolute threshold / LNOB universal value / scientific target / mean of the top five countries
# ---------------------------------------------------------------------

normalize_higher_better <- function(x, worst, best) {
  100 * (x - worst) / (best - worst)
}
normalize_lower_better <- function(x, worst, best) {
  100 * (worst - x) / (worst - best)
}
clamp100 <- function(s) pmin(100, pmax(0, s))

# A.3 decision tree fallback: returns c(worst, best)
derive_bounds <- function(values, direction = "higher", top_n = 5) {
  v <- values[!is.na(values)]
  worst <- as.numeric(quantile(v, 0.025, na.rm = TRUE))   # lower bound = the 2.5 quantile
  # Branch 4: mean of the top five countries (best observed values); branches 1-3 depend on the specific threshold and are provided by meta.csv
  srt <- sort(v, decreasing = (direction == "higher"))
  best <- mean(srt[seq_len(min(top_n, length(srt)))])
  c(worst = worst, best = best)
}

build_normalized <- function(imputed_df, meta = NULL, indicators = NULL) {
  if (is.null(meta)) meta <- build_ind_meta()
  # After the left join with the mapping table an indicator can have several rows -> the normalization bound table must be de-duplicated by indicator
  # The main analysis normalizes only the 44 indicators with used=Y; the indicators argument can override this as needed
  meta_u <- meta %>% distinct(indicator, .keep_all = TRUE)
  ids <- if (is.null(indicators)) meta_u$indicator[meta_u$used] else indicators
  m <- meta_u %>% filter(indicator %in% ids) %>%
    select(indicator, worst, best, direction)
  # When the meta bounds are missing, fall back to the A.3 decision tree
  need_fb <- setdiff(unique(imputed_df$indicator),
                     m$indicator[!is.na(m$worst) & !is.na(m$best)])
  fallback <- if (length(need_fb) > 0) {
    map_dfr(need_fb, function(id) {
      dir_i <- meta_u$direction[match(id, meta_u$indicator)]
      if (length(dir_i) == 0 || is.na(dir_i)) dir_i <- "higher"
      b <- derive_bounds(imputed_df$value[imputed_df$indicator == id], dir_i)
      tibble(indicator = id, worst = b[["worst"]], best = b[["best"]],
             direction = dir_i)
    })
  } else NULL
  m <- bind_rows(m %>% filter(!is.na(worst), !is.na(best)), fallback)

  imputed_df %>%
    inner_join(m, by = "indicator") %>%
    mutate(
      norm_raw = ifelse(direction == "higher",
                        normalize_higher_better(value, worst, best),
                        normalize_lower_better(value, worst, best)),
      norm = clamp100(norm_raw)          # the 0/100 truncation that comes with Formula 1
    ) %>%
    select(country, year, indicator, value, norm)
}

# ---------------------------------------------------------------------
# Section 4. Target/goal values - Formula 2: r = (1/n) Sum_j ind_j (equal weights)
#   An indicator mapped to several targets enters each of the targets it is mapped to
#   layer_id: targets "T1"–"T23"; goals "A"–"D"
# ---------------------------------------------------------------------

load_i2t <- function(i2t_path = "data/indicator_to_target.csv") {
  i2t <- read_csv(i2t_path, show_col_types = FALSE) %>%
    mutate(indicator = as.integer(indicator_id),
           target = as.character(target))
  i2t %>%
    mutate(layer = ifelse(target %in% c("A", "B", "C", "D"), "goal", "target"),
           layer_id = ifelse(layer == "goal", target, paste0("T", target))) %>%
    select(indicator, layer, layer_id) %>%
    distinct()                          # de-duplicate the mapping table, to avoid row inflation after the join
}

# Target <-> section/goal (target_section.csv)
load_target_section <- function(path = "data/target_section.csv") {
  read_csv(path, show_col_types = FALSE) %>%
    filter(!id %in% c("A", "B", "C", "D")) %>%
    transmute(layer_id = paste0("T", id),
              section = section,
              goal = goal,
              section_label = theme)
}

aggregate_goal_target <- function(norm_df, meta = NULL,
                                  i2t_path = "data/indicator_to_target.csv") {
  if (is.null(meta)) meta <- build_ind_meta()
  i2t <- load_i2t(i2t_path)

  norm_df %>%
    inner_join(i2t, by = "indicator", relationship = "many-to-many") %>%
    group_by(country, year, layer, layer_id) %>%
    summarise(index = mean_na(norm),                # Formula 2 (all NA -> NA)
              n_ind = n_distinct(indicator),
              .groups = "drop")
}

# ---------------------------------------------------------------------
# Section 5. MIS - Formula 3: arithmetic mean of the values of the targets (or goals) (computed year by year)
#   Countries with more than 20% of the indicators missing do not enter MIS (nor EIS/rankings)
# ---------------------------------------------------------------------

compute_mis <- function(layer_df, coverage_df, layer = "goal") {
  layer_df %>%
    filter(layer == .env$layer) %>%
    inner_join(coverage_df %>% filter(keep_index) %>% select(country),
               by = "country") %>%
    group_by(country, year) %>%
    summarise(MIS = mean_na(index),                 # Formula 3 (all NA -> NA)
              n_layers = sum(!is.na(index)), .groups = "drop")
}

# ---------------------------------------------------------------------
# Section 6. EIS - Formula 4 (Liu et al. 2021; Qi et al. 2024)
  # S = Sum_j (1/n).pi.r_j^2                 area of the radar plot (sum of the sector areas)
  # L = 2(r_max - r_min) + (2*pi/n).Sum_j r_j  perimeter estimate (range term + mean radius term)
  # EIS = S / (pi.(L/2*pi)^2) x 100 = 400*pi*S / L^2   ratio to the area of the circle with the same perimeter
  # r_j = the j-th target (or goal) value of the country; EIS = 100 when all values are equal
# ---------------------------------------------------------------------

eis_from_values <- function(r) {
  r <- r[!is.na(r)]
  n <- length(r)
  if (n < 2) return(NA_real_)
  S <- sum((1 / n) * pi * r^2)                            # area of the radar plot
  L <- 2 * (max(r) - min(r)) + (2 * pi / n) * sum(r)       # perimeter estimate
  if (L <= 0) return(NA_real_)
  400 * pi * S / L^2                                      # ratio to the area of the circle with the same perimeter x 100
}

compute_eis <- function(layer_df, coverage_df, layer = "goal") {
  layer_df %>%
    filter(layer == .env$layer) %>%
    inner_join(coverage_df %>% filter(keep_index) %>% select(country),
               by = "country") %>%
    group_by(country, year) %>%
    summarise(EIS = eis_from_values(index),
              n_layers = sum(!is.na(index)), .groups = "drop")
}

# Annual MIS/EIS series (computed at both the goal and the target resolution)
compute_mis_eis_annual <- function(layer_df, coverage_df) {
  # The annual MIS/EIS series is limited to 2011-2024 (SI B.2.7, 13 annual observations)
  layer_df <- layer_df %>% filter(year %in% TRAJ_YEARS)
  bind_rows(
    map_dfr(c("goal", "target"), function(ly) {
      mis <- compute_mis(layer_df, coverage_df, ly) %>%
        rename(!!paste0("MIS_", ly) := MIS)
      eis <- compute_eis(layer_df, coverage_df, ly) %>%
        rename(!!paste0("EIS_", ly) := EIS)
      mis %>% left_join(eis %>% select(country, year, starts_with("EIS")),
                        by = c("country", "year"))
    })
  ) %>%
    group_by(country, year) %>%
    # The goal rows and the target rows are joined with an offset: take the first non-NA value of each group (the original .x[1] would turn
    # all of MIS_target/EIS_target into NA)
    summarise(across(starts_with(c("MIS", "EIS")), ~ {
      v <- .x[!is.na(.x)]
      if (length(v) > 0) v[1] else NA_real_
    }), .groups = "drop") %>%
    select(country, year, any_of(c("MIS_goal", "MIS_target", "EIS_goal",
                                   "EIS_target"))) %>%
    # The main resolution = the goal level (the worked example of Formulas 3/4 uses the four goal values)
    mutate(MIS = MIS_goal, EIS = EIS_goal) %>%
    select(country, year, MIS_goal, MIS_target, EIS_goal, EIS_target, MIS, EIS)
}

# ---------------------------------------------------------------------
# Section 6b. Regional and global aggregation of the target/goal time series (the country level is written out as well)
#   country level = the (country, year, layer_id) series of Formula 2;
#   regional level = the equal-weight mean of the non-missing index values within the group (no population weighting, consistent
#            with the equal-weight convention of the country level); global level = the equal-weight mean over all countries.
#   The grouping convention follows Supplementary Table 7 / SI B.2.7: continent / subcontinent / region /
#   development / SIDS / LDC / LLDC / OECD / EU, plus Global.
# ---------------------------------------------------------------------

geo_group_membership <- function(country_meta) {
  bind_rows(
    country_meta %>% filter(!is.na(continent)) %>%
      transmute(country, group_type = "continent", group = continent),
    country_meta %>% filter(!is.na(subcontinent)) %>%
      transmute(country, group_type = "subcontinent", group = subcontinent),
    country_meta %>% filter(!is.na(region)) %>%
      transmute(country, group_type = "region", group = region),
    country_meta %>% filter(!is.na(development)) %>%
      transmute(country, group_type = "development", group = development),
    country_meta %>% filter(sids %in% TRUE) %>%
      transmute(country, group_type = "sids", group = "SIDS"),
    country_meta %>% filter(ldc %in% TRUE) %>%
      transmute(country, group_type = "ldc", group = "LDC"),
    country_meta %>% filter(llc %in% TRUE) %>%
      transmute(country, group_type = "lldc", group = "LLDC"),
    country_meta %>% filter(oecd %in% TRUE) %>%
      transmute(country, group_type = "oecd", group = "OECD"),
    if ("eu" %in% names(country_meta)) {
      country_meta %>% filter(eu %in% TRUE) %>%
        transmute(country, group_type = "eu", group = "EU")
    } else NULL
  ) %>% distinct()
}

# Target/goal values (Formula 2) as regional and global time series
#   With with_projection=TRUE the year-by-year projection is also computed, giving the complete time series
#   up to the target year (2030/2050) (region_ts / global_ts)
aggregate_layer_geography <- function(layer_df, country_meta,
                                      with_projection = TRUE) {
  gm <- geo_group_membership(country_meta)
  region <- layer_df %>%
    inner_join(gm, by = "country") %>%
    group_by(group_type, group, year, layer, layer_id) %>%
    summarise(index = mean_na(index),
              n_countries = n_distinct(country[!is.na(index)]),
              .groups = "drop")
  global <- layer_df %>%
    group_by(year, layer, layer_id) %>%
    summarise(index = mean_na(index),
              n_countries = n_distinct(country[!is.na(index)]),
              .groups = "drop") %>%
    mutate(group_type = "global", group = "Global") %>%
    select(group_type, group, year, layer, layer_id, index, n_countries)
  region_ts <- NULL; global_ts <- NULL
  if (with_projection) {
    rpaths <- project_geography_paths(region, c("group_type", "group"))
    gpaths <- project_geography_paths(global, c("group_type", "group"))
    region_ts <- goal_target_ts_full(region, rpaths,
                                     id_cols = c("group_type", "group"))
    global_ts <- goal_target_ts_full(global, gpaths,
                                     id_cols = c("group_type", "group"))
  }
  list(region = region, global = global,
       region_ts = region_ts, global_ts = global_ts)
}

# Regional and global aggregation of the annual MIS/EIS series
#   Convention: the equal-weight mean of the country-level MIS/EIS within the group (consistent with Supplementary Table 7 / SI B.2.7);
#   the country universe = the countries that pass the coverage filter and enter MIS/EIS (the same convention as the country level).
#   Segments are aggregated separately when the value_type column (historical/projected) is present.
mis_eis_geography <- function(mis_eis_country, country_meta) {
  gm <- geo_group_membership(country_meta)
  cols <- c("MIS_goal", "MIS_target", "EIS_goal", "EIS_target", "MIS", "EIS")
  cols <- intersect(cols, names(mis_eis_country))
  d <- mis_eis_country
  if (!"value_type" %in% names(d)) d <- d %>% mutate(value_type = "historical")
  region <- d %>%
    inner_join(gm, by = "country") %>%
    group_by(group_type, group, year, value_type) %>%
    summarise(across(all_of(cols), mean_na),
              n_countries = n_distinct(country),
              .groups = "drop")
  global <- d %>%
    group_by(year, value_type) %>%
    summarise(across(all_of(cols), mean_na),
              n_countries = n_distinct(country),
              .groups = "drop") %>%
    mutate(group_type = "global", group = "Global") %>%
    select(group_type, group, year, value_type, all_of(cols), n_countries)
  list(region = region, global = global)
}

# Projected segment of MIS/EIS: computed from the projected target/goal values with Formulas 3/4
#   the goal-level series runs to 2050 and the target-level series to 2030, so the two resolutions of MIS/EIS extend
#   to their own target year; only the countries passing the coverage rule are counted (the same convention as the historical segment).
mis_eis_projected <- function(goal_target_ts_country, coverage_df) {
  ts <- goal_target_ts_country %>%
    filter(value_type == "projected", !is.na(index)) %>%
    semi_join(coverage_df %>% filter(keep_index) %>% select(country),
              by = "country")
  goal_part <- ts %>%
    filter(layer == "goal") %>%
    group_by(country, year) %>%
    summarise(MIS_goal = mean_na(index),
              EIS_goal = eis_from_values(index), .groups = "drop")
  target_part <- ts %>%
    filter(layer == "target") %>%
    group_by(country, year) %>%
    summarise(MIS_target = mean_na(index),
              EIS_target = eis_from_values(index), .groups = "drop")
  full_join(goal_part, target_part, by = c("country", "year")) %>%
    mutate(MIS = MIS_goal, EIS = EIS_goal, value_type = "projected")
}

# Complete MIS/EIS time series (country level): the 2011-2024 history plus the projected segment up to the target year
mis_eis_ts_country_full <- function(mis_eis_hist, goal_target_ts_country,
                                    coverage_df) {
  proj <- mis_eis_projected(goal_target_ts_country, coverage_df) %>%
    anti_join(mis_eis_hist %>% select(country, year), by = c("country", "year"))
  bind_rows(
    mis_eis_hist %>% mutate(value_type = "historical"),
    proj
  ) %>%
    select(country, year, value_type, MIS_goal, MIS_target, EIS_goal,
           EIS_target, MIS, EIS) %>%
    arrange(country, year)
}

# ---------------------------------------------------------------------
# Section 7. Baseline (2011-2020) - the four parameters of Supplementary Table 3: mean, SD, trend slope, range
#   with the R2/p of the linear regression (needed for Fig 1 / Supplementary Table 3)
# ---------------------------------------------------------------------

# Safe linear trend statistics: y ~ year
#   Degeneracy guard: returns NA when the sample has < 3 points or year has no variation (e.g. only one valid year
#   inside the baseline period), so that summary(fit)$coefficients[2,4] cannot be indexed out of bounds
lm_trend_stats <- function(y, year) {
  ok <- is.finite(y) & is.finite(year)
  y <- y[ok]; year <- year[ok]
  if (length(y) < 3 || n_distinct(year) < 2) {
    return(list(slope = NA_real_, r2 = NA_real_, p_value = NA_real_))
  }
  fit <- tryCatch(lm(y ~ year), error = function(e) NULL)
  if (is.null(fit)) return(list(slope = NA_real_, r2 = NA_real_, p_value = NA_real_))
  cm <- tryCatch(coef(summary(fit)), error = function(e) NULL)
  if (is.null(cm) || nrow(cm) < 2) {
    return(list(slope = NA_real_, r2 = NA_real_, p_value = NA_real_))
  }
  list(slope = unname(cm[2, 1]),
       r2 = tryCatch(summary(fit)$r.squared, error = function(e) NA_real_),
       p_value = unname(cm[2, 4]))
}

baseline_params <- function(data_df, value_col = "index",
                            years = BASELINE_YEARS) {
  data_df %>%
    filter(year %in% years) %>%
    group_by(layer, layer_id) %>%
    group_modify(~ {
      d <- .x %>% filter(!is.na(.data[[value_col]]))
      st <- lm_trend_stats(d[[value_col]], d$year)
      tibble(
        n_countries = n_distinct(d$country),
        mean_v  = mean(d[[value_col]], na.rm = TRUE),
        sd_v    = sd(d[[value_col]], na.rm = TRUE),
        range_v = diff(range(d[[value_col]], na.rm = TRUE)),
        slope   = st$slope, r2 = st$r2, p_value = st$p_value
      )
    }) %>%
    ungroup()
}

# Baseline at the indicator level (the base layer of Supplementary Tables 1/3)
baseline_params_indicator <- function(norm_df, years = BASELINE_YEARS) {
  norm_df %>%
    filter(year %in% years) %>%
    group_by(indicator) %>%
    group_modify(~ {
      d <- .x %>% filter(!is.na(norm))
      st <- lm_trend_stats(d$norm, d$year)
      tibble(
        n_countries = n_distinct(d$country),
        mean_v  = mean(d$norm, na.rm = TRUE),
        sd_v    = sd(d$norm, na.rm = TRUE),
        range_v = diff(range(d$norm, na.rm = TRUE)),
        slope   = st$slope, r2 = st$r2, p_value = st$p_value
      )
    }) %>%
    ungroup()
}

# The four baseline parameters per country (used for the between-continent tests; Supplementary Table 3 / Supplementary Fig. 1)
per_country_baseline_params <- function(data_df, value_col = "index",
                                        years = BASELINE_YEARS) {
  data_df %>%
    filter(year %in% years) %>%
    group_by(layer, layer_id, country) %>%
    group_modify(~ {
      d <- .x %>% filter(!is.na(.data[[value_col]]))
      st <- lm_trend_stats(d[[value_col]], d$year)
      tibble(
        mean_v  = mean(d[[value_col]], na.rm = TRUE),
        sd_v    = sd(d[[value_col]], na.rm = TRUE),
        range_v = diff(range(d[[value_col]], na.rm = TRUE)),
        slope   = st$slope
      )
    }) %>%
    ungroup()
}

# Kruskal-Wallis between continents + Bonferroni-corrected Mann-Whitney U post-hoc (Supplementary Table 3)
#   the post-hoc is only run when the omnibus p < 0.05 (omnibus-gated post-hoc of Supplementary Table 3)
baseline_kw_posthoc <- function(per_country_df, country_meta,
                                params = c("mean_v", "sd_v", "slope", "range_v")) {
  d <- per_country_df %>%
    inner_join(country_meta %>% select(country, continent), by = "country") %>%
    filter(!is.na(continent))
  kw <- d %>%
    pivot_longer(all_of(params), names_to = "param", values_to = "v") %>%
    group_by(layer, layer_id, param) %>%
    summarise(
      kw_H = {
        ok <- !is.na(v)
        tryCatch({
          if (sum(ok) >= 5 && n_distinct(continent[ok]) >= 2)
            unname(kruskal.test(v[ok] ~ factor(continent[ok]))$statistic) else NA_real_
        }, error = function(e) NA_real_)
      },
      kw_p = {
        ok <- !is.na(v)
        tryCatch({
          if (sum(ok) >= 5 && n_distinct(continent[ok]) >= 2)
            kruskal.test(v[ok] ~ factor(continent[ok]))$p.value else NA_real_
        }, error = function(e) NA_real_)
      },
      n = sum(!is.na(v)),
      .groups = "drop"
    )
  posthoc <- d %>%
    pivot_longer(all_of(params), names_to = "param", values_to = "v") %>%
    group_by(layer, layer_id, param) %>%
    group_modify(~ {
      # the post-hoc is only run when the omnibus test is significant
      ok <- !is.na(.x$v)
      if (sum(ok) < 5 || n_distinct(.x$continent[ok]) < 2)
        return(tibble(cont1 = character(0), cont2 = character(0),
                      mean_diff = numeric(0), u_p_raw = numeric(0),
                      u_p_bonf = numeric(0)))
      p_omnibus <- tryCatch(
        kruskal.test(.x$v[ok] ~ factor(.x$continent[ok]))$p.value,
        error = function(e) NA_real_)
      # the post-hoc is only run when the omnibus test is significant; it is skipped when the omnibus p is not available (NA)
      if (is.na(p_omnibus) || p_omnibus >= 0.05)
        return(tibble(cont1 = character(0), cont2 = character(0),
                      mean_diff = numeric(0), u_p_raw = numeric(0),
                      u_p_bonf = numeric(0)))
      lv <- sort(unique(.x$continent[ok]))
      prs <- combn(lv, 2, simplify = FALSE)
      res <- map_dfr(prs, function(pr) {
        s <- .x$continent %in% pr & ok
        p <- tryCatch(wilcox.test(.x$v[s] ~ factor(.x$continent[s]))$p.value,
                      error = function(e) NA_real_)
        md <- mean(.x$v[s & .x$continent == pr[1]], na.rm = TRUE) -
              mean(.x$v[s & .x$continent == pr[2]], na.rm = TRUE)
        tibble(cont1 = pr[1], cont2 = pr[2], mean_diff = md, u_p_raw = p)
      })
      res %>% mutate(u_p_bonf = pmin(1, u_p_raw * n()))
    }) %>%
    ungroup()
  list(kw = kw, posthoc_bonferroni = posthoc)
}

# ---------------------------------------------------------------------
# Section 8. Trajectory projection - A.5
#   high noise (SNR < 1) -> DLM; low noise (SNR >= 1) -> PMA
#   PMA = a genuine multi-model averaging (Tittensor et al. 2014, Science
#   346:241-244, doi:10.1126/science.1257484; Burnham & Anderson 2002):
#     candidate model family = 18 parametric models (identical to the original analysis code) + quad + segmented
#     multi-start nls fit -> AICc weights -> the confidence set with dAICc<2 -> weighted average prediction
#     unconditional variance = Sum w_i(sigma_i^2 + (pred_i - pred_avg)^2)
#   Target year: 2030 (indicator/target); 2050 (goal)
#   The fit uses the centred year x = year - 2000 (a pure shift, numerically stable; the prediction uses the same transformation)
# ---------------------------------------------------------------------

PMA_YEAR_CENTER <- 2000

# Model family of Tittensor et al. 2014 (18 parametric models, identical to the original analysis code)
# Signature f(x, a, b[, c[, d]]); poly1/quad/segmented use lm, the others use nls
pma_model_funs <- list(
  poly1            = function(x, a, b) a + b * x,
  hyperbolic       = function(x, a, b) a / (b + x),
  michaelis_menten = function(x, a, b) (a * x) / (b + x),
  holling3         = function(x, a, b) (a * x^2) / (b + x^2),
  holling4         = function(x, a, b, c) (a * x^2) / (b + c * x + x^2),
  rational         = function(x, a, b, c) (a + b * x) / (1 + c * x),
  exponential      = function(x, a, b) a * exp(b * x),
  monomolecular    = function(x, a, b) a - a * exp(b * x),
  ricker           = function(x, a, b) a * x * exp(b * x),
  gompertz         = function(x, a, b, c) exp(a * exp(b * (x - c))),
  chapman_richards = function(x, a, b, c) a * (1 - exp(b * x))^c,
  unnamed_exp      = function(x, a, b, c) a * exp(b / (x + c)),
  weibull          = function(x, a, b, c, d) a * (1 - exp(b * (x - c)))^d,
  power_law        = function(x, a, b) a * x^b,
  unnamed_power    = function(x, a, b) a * (1 + x)^b,
  asymptotic       = function(x, a, b, c) a - b * c^x,
  shepherd         = function(x, a, b, c) (a * x) / (b + x^c),
  hassell          = function(x, a, b, c) (a * x) / ((b + x)^c),
  # The quadratic form of the SI 20-form candidate family (A.5) is included among the candidates as well
  quad             = function(x, a, b, c) a + b * x + c * x^2
)

# Environment for the nls formulas: lets the formulas resolve the model functions by name (fixes the problem of the functions not being found inside a list in the old code)
pma_env <- new.env(parent = globalenv())
for (nm in names(pma_model_funs)) assign(nm, pma_model_funs[[nm]], envir = pma_env)

# Two-parameter subset (only the two-parameter models are used when n < 10, to avoid overfitting; identical to the original analysis code)
PMA_TWO_PARAM <- c("poly1", "hyperbolic", "michaelis_menten", "holling3",
                   "exponential", "monomolecular", "ricker", "power_law",
                   "unnamed_power")

# Multi-start nls fit of a single model; returns the best fit or NULL
# Performance constraint: maxiter is strictly limited and only genuinely converged fits are accepted;
# two sets of starting points are used (unit and order-of-magnitude starts), so that the high-dimensional models do not hang
fit_nls_multi <- function(model_name, x, value, k) {
  par_names <- letters[seq_len(k)]
  fml_txt <- paste("value ~", model_name, "(x,", paste(par_names, collapse = ", "), ")")
  fml <- as.formula(fml_txt)
  environment(fml) <- pma_env
  df <- data.frame(value = value, x = x)
  # Set of starting points: unit / order-of-magnitude / inverse order of magnitude / small parameters
  starts <- list(
    setNames(rep(1, k), par_names),
    setNames(c(mean(value), rep(0.1, k))[seq_len(k)], par_names),
    setNames(c(max(value), rep(-0.1, k))[seq_len(k)], par_names),
    setNames(rep(0.01, k), par_names)
  )
  best <- NULL
  for (st in starts) {
    fit <- tryCatch(
      nls(fml, data = df, start = as.list(st),
          control = nls.control(maxiter = 40, warnOnly = FALSE, minFactor = 1 / 1024)),
      error = function(e) NULL)
    if (is.null(fit)) next
    # Only genuinely converged fits are accepted (in warnOnly mode a fit can be falsely convergent)
    if (!isTRUE(fit$convInfo$isConv)) next
    a <- tryCatch(AIC(fit), error = function(e) NA_real_)
    if (!is.finite(a)) next
    if (is.null(best) || a < best$aic) best <- list(fit = fit, aic = a)
    # When the unit starting point has already converged to a good solution, the remaining starting points are skipped (speed)
    if (!is.null(best) && st[[1]] == 1 && best$aic < a + 2) break
  }
  best
}

# AICc: small-sample correction (Burnham & Anderson 2002); k = number of parameters + residual variance
aicc_of <- function(fit, k, n) {
  a <- tryCatch(AIC(fit), error = function(e) NA_real_)
  if (!is.finite(a)) return(NA_real_)
  if (n - k - 1 > 0) a + 2 * k * (k + 1) / (n - k - 1) else Inf
}

# PMA fit: all the candidates -> AICc weights -> the confidence set with dAICc<2 (if it is empty, the best single model is taken)
fit_pma <- function(year, value) {
  d <- tibble(year = year, value = value) %>%
    filter(is.finite(value), is.finite(year)) %>%
    arrange(year) %>% distinct(year, .keep_all = TRUE)
  n <- nrow(d)
  if (n < 5) return(NULL)
  x <- d$year - PMA_YEAR_CENTER          # centring
  y <- d$value

  cands <- list()   # name -> list(kind, fit, k, breakpoint)
  # --- lm fits: poly1 / quad / segmented (these always converge and serve as a fallback) ---
  m1 <- tryCatch(lm(y ~ x), error = function(e) NULL)
  if (!is.null(m1)) cands$poly1 <- list(kind = "lm", fit = m1, k = 3, name = "poly1")
  m2 <- tryCatch(lm(y ~ poly(x, 2, raw = TRUE)), error = function(e) NULL)
  if (!is.null(m2)) cands$quad <- list(kind = "lm", fit = m2, k = 4, name = "quad")
  if (n >= 8) {
    # segmented: continuous piecewise linear, the breakpoint is searched over the interior year grid
    best_seg <- NULL
    for (bp in x[-c(1, n)]) {
      m <- tryCatch(lm(y ~ x + pmax(0, x - bp)), error = function(e) NULL)
      if (is.null(m)) next
      a <- aicc_of(m, 4, n)
      if (is.finite(a) && (is.null(best_seg) || a < best_seg$aic)) {
        best_seg <- list(fit = m, aic = a, bp = bp + PMA_YEAR_CENTER)
      }
    }
    if (!is.null(best_seg)) {
      cands$segmented <- list(kind = "segmented", fit = best_seg$fit, name = "segmented",
                              k = 4, breakpoint = as.integer(best_seg$bp))
    }
  }
  # --- nls fits: the family of 18 parametric models ---
  nls_names <- setdiff(names(pma_model_funs), c("poly1", "quad"))
  if (n < 10) nls_names <- intersect(nls_names, PMA_TWO_PARAM)
  for (nm in nls_names) {
    k <- length(formals(pma_model_funs[[nm]])) - 1
    if (k + 1 >= n) next
    r <- fit_nls_multi(nm, x, y, k)
    if (!is.null(r)) {
      cands[[nm]] <- list(kind = "nls", fit = r$fit, k = k + 1, name = nm)
    }
  }
  if (length(cands) == 0) return(NULL)

  aics <- map_dbl(cands, function(cc) aicc_of(cc$fit, cc$k, n))
  if (!any(is.finite(aics))) return(NULL)
  delta <- aics - min(aics[is.finite(aics)])
  w <- exp(-0.5 * pmin(delta, 700))
  w[!is.finite(delta)] <- 0
  if (sum(w) <= 0) return(NULL)
  w <- w / sum(w)
  # Confidence set: dAICc < 2 (Burnham & Anderson); if it is empty the best model is taken
  in_set <- is.finite(delta) & delta < 2
  if (!any(in_set)) in_set[which.min(aics)] <- TRUE
  w_set <- w[in_set]; w_set <- w_set / sum(w_set)

  list(n = n, candidates = cands[in_set], weights = w_set,
       all_aicc = aics, all_weights = w,
       top = names(which.max(w)),
       top_weight = max(w),
       n_avg = sum(in_set),
       breakpoint = {
         seg <- cands$segmented
         if (!is.null(seg) && "segmented" %in% names(cands[in_set])) seg$breakpoint
         else NA_integer_
       })
}

# Prediction variance for nls: delta method
#   var = sigma^2 + g' V g,  g = numerical gradient of the model with respect to the parameters, V = parameter covariance
#   (includes the parameter uncertainty, so the extrapolation no longer underestimates; the sigma_i^2 of the unconditional variance of Burnham & Anderson)
nls_pred_var <- function(fit, x_new, model_name) {
  s2 <- tryCatch(summary(fit)$sigma^2, error = function(e) NA_real_)
  V <- tryCatch(vcov(fit), error = function(e) NULL)
  if (is.null(V) || any(!is.finite(V))) return(rep(s2, length(x_new)))
  theta <- coef(fit)
  mfun <- pma_model_funs[[model_name]]
  if (is.null(mfun)) return(rep(s2, length(x_new)))
  # central-difference numerical gradient
  grad <- matrix(NA_real_, nrow = length(x_new), ncol = length(theta))
  for (j in seq_along(theta)) {
    h <- max(abs(theta[j]) * 1e-5, 1e-8)
    tp <- theta; tm <- theta
    tp[j] <- tp[j] + h; tm[j] <- tm[j] - h
    fp <- tryCatch(do.call(mfun, c(list(x_new), as.list(tp))),
                   error = function(e) NULL)
    fm <- tryCatch(do.call(mfun, c(list(x_new), as.list(tm))),
                   error = function(e) NULL)
    if (is.null(fp) || is.null(fm)) return(rep(s2, length(x_new)))
    grad[, j] <- (fp - fm) / (2 * h)
  }
  v <- s2 + rowSums((grad %*% V) * grad)
  v
}

# Point prediction and prediction variance of a single candidate at a new x
predict_candidate <- function(cc, x_new) {
  nd <- data.frame(x = x_new)
  pred <- tryCatch(predict(cc$fit, newdata = nd), error = function(e) NULL)
  if (is.null(pred) || any(!is.finite(pred))) return(NULL)
  pred <- as.numeric(pred)
  # Within-model variance: lm uses the prediction variance (SE^2 + sigma^2); nls uses the delta method (sigma^2 + g'Vg)
  if (cc$kind == "nls") {
    v <- nls_pred_var(cc$fit, x_new, cc$name)
  } else {
    pr <- tryCatch(predict(cc$fit, newdata = nd, se.fit = TRUE),
                   error = function(e) NULL)
    if (is.null(pr)) {
      s2 <- tryCatch(summary(cc$fit)$sigma^2, error = function(e) NA_real_)
      v <- rep(s2, length(pred))
    } else {
      v <- pr$se.fit^2 + summary(cc$fit)$sigma^2
    }
  }
  v <- as.numeric(v)
  v[!is.finite(v) | v < 0] <- NA_real_
  list(pred = pred, var = v)
}

# Model-averaged prediction: pred_avg = Sum w_i pred_i
# Unconditional variance (Burnham & Anderson 2002): Sum w_i(sigma_i^2 + (pred_i - pred_avg)^2)
predict_pma <- function(pma, year_new) {
  x_new <- year_new - PMA_YEAR_CENTER
  preds <- matrix(NA_real_, nrow = length(x_new), ncol = length(pma$candidates))
  vars <- matrix(NA_real_, nrow = length(x_new), ncol = length(pma$candidates))
  valid <- logical(length(pma$candidates))
  for (j in seq_along(pma$candidates)) {
    r <- predict_candidate(pma$candidates[[j]], x_new)
    if (!is.null(r)) {
      preds[, j] <- r$pred
      vars[, j] <- r$var
      valid[j] <- TRUE
    }
  }
  if (!any(valid)) return(NULL)
  w <- pma$weights[valid]; w <- w / sum(w)
  preds <- preds[, valid, drop = FALSE]
  vars <- vars[, valid, drop = FALSE]
  avg <- as.numeric(preds %*% w)
  uncond <- map_dbl(seq_along(avg), function(i) {
    sum(w * (vars[i, ] + (preds[i, ] - avg[i])^2), na.rm = TRUE)
  })
  list(mean = avg, var = uncond)
}

# Signal-to-noise ratio SNR (the noise split criterion of Tittensor et al. 2014, defined as in the original analysis code)
#   random walk plus noise model (dlmModPoly(order=1)) MLE: SNR = process variance w / observation variance v
#   fallback when the estimation fails: SNR = var(y) / (var(diff(y))/2)
#   SNR >= SNR_DLM_THRESHOLD -> low noise -> PMA; otherwise high noise -> DLM
calculate_snr <- function(ts_data) {
  y <- as.numeric(ts_data)
  y <- y[is.finite(y)]
  if (length(y) < 3) return(Inf)
  snr <- tryCatch({
    if (requireNamespace("dlm", quietly = TRUE)) {
      build_rw <- function(par) {
        m <- dlm::dlmModPoly(order = 1)
        m$V[] <- exp(par[1]); m$W[] <- exp(par[2])
        m
      }
      fit <- dlm::dlmMLE(y, parm = c(0, 0), build = build_rw)
      if (is.null(fit) || !is.finite(fit$value)) stop("snr mle failed")
      v <- exp(fit$par[1]); w <- exp(fit$par[2])
      if (!is.finite(v) || !is.finite(w) || v <= 0 || w <= 0) stop("bad var")
      w / v
    } else {
      # built-in local level model MLE (equivalent to dlmModPoly(order=1))
      kalman_ll <- function(logpar, y) {
        V <- exp(logpar[1]); W <- exp(logpar[2])
        m <- y[1]; C <- V * 10; ll <- 0
        for (t in seq_along(y)) {
          a <- m; R <- C + W
          f <- a; Q <- R + V; e <- y[t] - f
          ll <- ll - 0.5 * (log(2 * pi * Q) + e^2 / Q)
          if (!is.finite(ll)) return(1e10)
          K <- R / Q; m <- a + K * e; C <- R - K * R
        }
        -ll
      }
      mle <- optim(c(0, 0), kalman_ll, y = y,
                   control = list(maxit = 200, reltol = 1e-8))
      v <- exp(mle$par[1]); w <- exp(mle$par[2])
      if (!is.finite(v) || !is.finite(w) || v <= 0 || w <= 0) stop("bad var")
      w / v
    }
  }, error = function(e) {
    # simple fallback (the convention of the original analysis code)
    signal_var <- var(y)
    noise_var <- var(diff(y)) / 2
    if (is.na(signal_var) || is.na(noise_var) || noise_var <= 0) Inf
    else signal_var / noise_var
  })
  if (length(snr) != 1 || is.na(snr)) snr <- Inf
  snr
}

# DLM (local linear trend, Tittensor et al. 2014 ref 10 = dlmModPoly(order=2))
#   the dlm package is used when available; otherwise the built-in equivalent Kalman filter (the same state-space model):
#     mu_t  = mu_{t-1} + nu_{t-1} + w1_t   (level)
#     nu_t  = nu_{t-1} + w2_t              (slope)
#     y_t   = mu_t + v_t                   (observation)
#   the variances (V, W1, W2) are estimated by maximum likelihood
fit_dlm <- function(value) {
  y <- as.numeric(value)
  y <- y[is.finite(y)]
  if (length(y) < 4) return(NULL)
  if (requireNamespace("dlm", quietly = TRUE)) {
    build <- function(p) {
      m <- dlm::dlmModPoly(order = 2)
      m$V[] <- exp(p[1]); m$W[] <- exp(p[2])
      m
    }
    mle <- tryCatch(dlm::dlmMLE(y, parm = c(0, 0), build = build),
                    error = function(e) NULL)
    if (is.null(mle)) return(NULL)
    return(tryCatch(build(mle$par), error = function(e) NULL))
  }
  # --- built-in local linear trend DLM: Kalman filter + MLE ---
  kalman_ll <- function(logpar, y) {
    V <- exp(logpar[1]); W1 <- exp(logpar[2]); W2 <- exp(logpar[3])
    n <- length(y)
    m <- c(y[1], 0)                      # state: (level, slope)
    C <- diag(c(V * 10, V * 10))
    G <- matrix(c(1, 0, 1, 1), 2, 2)     # state transition
    W <- diag(c(W1, W2))
    ll <- 0
    for (t in seq_len(n)) {
      a <- G %*% m
      R <- G %*% C %*% t(G) + W
      f <- a[1]
      Q <- R[1, 1] + V
      e <- y[t] - f
      ll <- ll - 0.5 * (log(2 * pi * Q) + e^2 / Q)
      if (!is.finite(ll)) return(1e10)
      K <- R[, 1] / Q
      m <- a + K * e
      C <- R - K %*% t(R[, 1, drop = FALSE])
      C <- (C + t(C)) / 2                # keep it symmetric
    }
    -ll
  }
  mle <- tryCatch(
    optim(c(0, 0, 0), kalman_ll, y = y,
          control = list(maxit = 200, reltol = 1e-8)),
    error = function(e) NULL)
  if (is.null(mle) || !is.finite(mle$value)) return(NULL)
  # returns a lightweight object compatible with the dlm package
  structure(list(type = "builtin_local_linear_trend",
                 V = exp(mle$par[1]), W = exp(mle$par[2:3]), y = y),
            class = "builtin_dlm")
}

# DLM filter + h-step prediction (built-in object; objects of the dlm package follow the original path of project_entity)
#   returns the endpoint (f, R) and the step-by-step path for steps 1..h (f_path, R_path), for the complete time series
builtin_dlm_forecast <- function(mod, h) {
  y <- mod$y; V <- mod$V; W <- mod$W
  n <- length(y)
  m <- c(y[1], 0); C <- diag(c(V * 10, V * 10))
  G <- matrix(c(1, 0, 1, 1), 2, 2); Wm <- diag(W)
  for (t in seq_len(n)) {
    a <- G %*% m; R <- G %*% C %*% t(G) + Wm
    f <- a[1]; Q <- R[1, 1] + V; e <- y[t] - f
    K <- R[, 1] / Q
    m <- a + K * e; C <- R - K %*% t(R[, 1, drop = FALSE]); C <- (C + t(C)) / 2
  }
  # h-step prediction: accumulated step by step
  f_path <- rep(NA_real_, h); R_path <- rep(NA_real_, h)
  a <- G %*% m; R <- G %*% C %*% t(G) + Wm
  f_path[1] <- as.numeric(a[1]); R_path[1] <- R[1, 1] + V
  if (h > 1) {
    for (i in 2:h) {
      a <- G %*% a; R <- G %*% R %*% t(G) + Wm
      f_path[i] <- as.numeric(a[1]); R_path[i] <- R[1, 1] + V
    }
  }
  list(f = f_path[h], R = R_path[h], f_path = f_path, R_path = R_path)
}

# Project a single series -> value at the target year + 95% CI + model averaging information
# With return_fit=TRUE the PMA fit object is attached as attr(, "fit"), so that the path bootstrap can reuse it
# attr(, "path"): the year-by-year prediction path (from the year after the end of the series to the target year), for the complete time series
#   written out: year, value, ci_lo, ci_hi
project_entity <- function(year, value, target_year, return_fit = FALSE) {
  d <- tibble(year = year, value = value) %>%
    filter(!is.na(value)) %>% arrange(year)
  path_empty <- tibble(year = integer(0), value = numeric(0),
                       ci_lo = numeric(0), ci_hi = numeric(0))
  out_na <- tibble(projected = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                   model_type = NA_character_, model_form = NA_character_,
                   n_models = NA_integer_, top_weight = NA_real_,
                   aicc = NA_real_, breakpoint = NA_integer_,
                   cv = NA_real_, snr = NA_real_)
  attr(out_na, "path") <- path_empty
  if (nrow(d) < 5) return(out_na)
  cv <- sd(d$value) / abs(mean(d$value))
  snr <- calculate_snr(d$value)          # the noise criterion of Tittensor et al.
  h <- target_year - max(d$year)
  # Years of the year-by-year path: after the end of the series -> the target year (empty when the series already reaches the target year)
  path_years <- if (h >= 1) seq.int(max(d$year) + 1, target_year) else integer(0)

  # High noise (SNR < 1) -> DLM (Tittensor ref 10; the dlm package or the built-in equivalent implementation)
  if (snr < SNR_DLM_THRESHOLD) {
    mod <- fit_dlm(d$value)
    if (!is.null(mod)) {
      res <- tryCatch({
        if (inherits(mod, "builtin_dlm")) {
          fc <- builtin_dlm_forecast(mod, max(h, 1))
          f_path <- fc$f_path; R_path <- fc$R_path
          list(f_last = fc$f, se = sqrt(fc$R), f_path = f_path,
               R_path = R_path)
        } else {
          filt <- dlm::dlmFilter(d$value, mod)
          fc <- dlm::dlmForecast(filt, nAhead = max(h, 1))
          f_last <- as.numeric(tail(fc$f, 1)[[1]])[1]
          # variance of the observation level = Z R Z' + V; for a local linear trend Z = (1, 0) -> R[1,1] + V
          R_last <- tail(fc$R, 1)[[1]]
          obs_var <- as.numeric(R_last[1, 1]) + as.numeric(mod$V)[1]
          f_path <- vapply(fc$f, function(x) as.numeric(x)[1], numeric(1))
          R_path <- vapply(fc$R, function(x)
            as.numeric(x[1, 1]) + as.numeric(mod$V)[1], numeric(1))
          list(f_last = f_last[1], se = sqrt(max(obs_var, 0)),
               f_path = f_path, R_path = R_path)
        }
      }, error = function(e) NULL)
      if (!is.null(res) && is.finite(res$f_last) && length(res$f_last) == 1) {
        out <- tibble(projected = as.numeric(res$f_last)[1],
                      ci_lo = as.numeric(res$f_last)[1] - 1.96 * as.numeric(res$se)[1],
                      ci_hi = as.numeric(res$f_last)[1] + 1.96 * as.numeric(res$se)[1],
                      model_type = "DLM", model_form = "DLM",
                      n_models = 1L, top_weight = 1,
                      aicc = NA_real_, breakpoint = NA_integer_,
                      cv = cv, snr = snr)
        if (length(path_years) > 0) {
          k <- seq_along(path_years)
          se_p <- sqrt(pmax(res$R_path[k], 0))
          attr(out, "path") <- tibble(
            year = path_years, value = res$f_path[k],
            ci_lo = res$f_path[k] - 1.96 * se_p,
            ci_hi = res$f_path[k] + 1.96 * se_p)
        } else {
          attr(out, "path") <- path_empty
        }
        return(out)
      }
    }
  }
  # Low noise -> PMA: multi-model averaging
  pma <- fit_pma(d$year, d$value)
  if (is.null(pma)) return(out_na %>% mutate(cv = cv))
  pr <- predict_pma(pma, target_year)
  if (is.null(pr)) return(out_na %>% mutate(cv = cv))
  m <- pr$mean[length(pr$mean)]
  v <- pr$var[length(pr$var)]
  se <- sqrt(max(v, 0))
  out <- tibble(projected = m,
                ci_lo = m - 1.96 * se, ci_hi = m + 1.96 * se,
                model_type = "PMA",
                model_form = pma$top,
                n_models = pma$n_avg,
                top_weight = round(pma$top_weight, 3),
                aicc = unname(pma$all_aicc[pma$top]),
                breakpoint = pma$breakpoint, cv = cv, snr = snr)
  if (length(path_years) > 0) {
    pr_path <- predict_pma(pma, path_years)
    if (!is.null(pr_path)) {
      se_p <- sqrt(pmax(pr_path$var, 0))
      attr(out, "path") <- tibble(
        year = path_years, value = pr_path$mean,
        ci_lo = pr_path$mean - 1.96 * se_p,
        ci_hi = pr_path$mean + 1.96 * se_p)
    }
  } else {
    attr(out, "path") <- path_empty
  }
  if (return_fit) attr(out, "fit") <- pma
  out
}

# Model-averaged path: the year-by-year model values after 2020 (2021..target year) -> bootstrap the annual rate of change
model_path_after_2020 <- function(pma, from_year = 2021L, target_year) {
  ys <- from_year:target_year
  pr <- predict_pma(pma, ys)
  if (is.null(pr)) return(tibble(year = ys, value = rep(NA_real_, length(ys))))
  tibble(year = ys, value = pr$mean)
}

# Projection + bootstrap of the annual rate of change of the post-2020 model path (A.5, called as one unit)
# The PMA fit is done only once and shared by the projection and the path (the same model average)
# Returns list(proj = the projection and rate of change in one row, path = the year-by-year prediction path, for the complete time series)
project_with_rate <- function(year, value, target_year, n_boot = N_BOOT_RATE) {
  d <- tibble(year = year, value = value) %>% filter(!is.na(value)) %>% arrange(year)
  pr <- project_entity(d$year, d$value, target_year, return_fit = TRUE)
  ts_path <- attr(pr, "path")
  if (is.null(ts_path)) {
    ts_path <- tibble(year = integer(0), value = numeric(0),
                      ci_lo = numeric(0), ci_hi = numeric(0))
  }
  if (is.na(pr$projected)) {
    proj <- bind_cols(pr, tibble(rate = NA_real_, rate_ci_lo = NA_real_,
                                 rate_ci_hi = NA_real_, p_boot = NA_real_,
                                 significant = NA))
    return(list(proj = proj, path = ts_path))
  }
  pma <- attr(pr, "fit")
  if (pr$model_type == "DLM" || is.null(pma)) {
    ok <- !is.na(d$value)
    slope <- coef(lm(d$value[ok] ~ d$year[ok]))[2]
    path <- tibble(year = (2021L):target_year,
                   value = pr$projected - slope * rev(seq_along(2021L:target_year)))
  } else {
    path <- model_path_after_2020(pma, 2021L, target_year)
  }
  br <- bootstrap_rate_of_change(path$year, path$value, n_boot = n_boot)
  out <- bind_cols(pr, br %>% rename(rate = rate,
                                     rate_ci_lo = ci_lo, rate_ci_hi = ci_hi,
                                     significant = significant))
  attr(out, "fit") <- NULL
  list(proj = out, path = ts_path)
}

# Projection without bootstrap: for writing out the time series of the regional/global series
project_path_only <- function(year, value, target_year) {
  pr <- project_entity(year, value, target_year)
  list(proj = pr, path = attr(pr, "path"))
}

# Year-by-year prediction path of the regional/global series (the same model family as A.5; no bootstrap)
#   id_cols: the columns identifying the entity, e.g. c("group_type", "group")
project_geography_paths <- function(layer_geo_df, id_cols) {
  keys <- layer_geo_df %>%
    distinct(across(all_of(c(id_cols, "layer", "layer_id"))))
  map_dfr(seq_len(nrow(keys)), function(k) {
    key <- keys[k, , drop = FALSE]
    d <- layer_geo_df %>%
      semi_join(key, by = c(id_cols, "layer", "layer_id")) %>%
      filter(!is.na(index)) %>%
      arrange(year)
    if (nrow(d) < 5) return(NULL)
    ty <- if (key$layer[1] == "goal") TARGET_YEAR_GOAL else TARGET_YEAR_TARGET
    pr <- project_path_only(d$year, d$index, ty)
    if (is.null(pr$path) || nrow(pr$path) == 0) return(NULL)
    out <- pr$path %>%
      mutate(layer = key$layer[1], layer_id = key$layer_id[1],
             model_type = pr$proj$model_type[1],
             model_form = pr$proj$model_form[1])
    for (cc in id_cols) out[[cc]] <- key[[cc]][1]
    out
  })
}

# Complete target/goal time series: the existing data (observed + imputed) + the year-by-year projected segment, up to the target year
#   the hist segment = the rows of layer_df that have a value (from the first data year);
#   the projected segment = the model path fills the years without data (up to 2030/2050).
goal_target_ts_full <- function(layer_df, proj_paths, id_cols = "country") {
  extra <- intersect(c("n_ind", "n_countries"), names(layer_df))
  hist <- layer_df %>%
    filter(!is.na(index)) %>%
    select(all_of(id_cols), year, layer, layer_id, index, any_of(extra)) %>%
    mutate(value_type = "historical", ci_lo = NA_real_, ci_hi = NA_real_,
           model_type = NA_character_, model_form = NA_character_)
  if (is.null(proj_paths) || nrow(proj_paths) == 0) {
    return(hist %>% arrange(!!!rlang::syms(id_cols), layer, layer_id, year))
  }
  proj <- proj_paths %>%
    anti_join(hist %>% select(all_of(id_cols), year, layer, layer_id),
              by = c(id_cols, "year", "layer", "layer_id")) %>%
    select(all_of(id_cols), year, layer, layer_id,
           index = value, ci_lo, ci_hi, model_type, model_form) %>%
    mutate(value_type = "projected",
           model_type = as.character(model_type),
           model_form = as.character(model_form))
  bind_rows(hist, proj) %>%
    arrange(!!!rlang::syms(id_cols), layer, layer_id, year)
}

# ---------------------------------------------------------------------
# Section 9. Bootstrap of the annual rate of change (A.5)
#   the model path after 2020 -> resampling -> linear fit of the annual rate of change -> 10,000 times
#   -> the 2.5%/97.5% quantiles form the 95% CI -> significant if it does not contain 0
# ---------------------------------------------------------------------

bootstrap_rate_of_change <- function(path_years, path_values,
                                     n_boot = N_BOOT_RATE, seed = 1) {
  ok <- is.finite(path_values)
  y <- path_years[ok]; v <- path_values[ok]
  n <- length(v)
  if (n < 3) {
    return(tibble(rate = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                  p_boot = NA_real_, significant = NA))
  }
  set.seed(seed)
  # vectorised bootstrap: each resample draws n points and the slope is computed from the sufficient statistics
  idx <- matrix(sample.int(n, n_boot * n, replace = TRUE), nrow = n_boot)
  Y <- matrix(v[idx], nrow = n_boot)
  X <- matrix(y[idx], nrow = n_boot)
  mx <- rowMeans(X); my <- rowMeans(Y)
  sxx <- rowMeans(X * X) - mx^2
  sxy <- rowMeans(X * Y) - mx * my
  rates <- ifelse(sxx > 0, sxy / sxx, NA_real_)
  ci <- quantile(rates, c(0.025, 0.975), na.rm = TRUE)
  rate <- coef(lm(v ~ y))[2]
  # two-sided bootstrap p value
  p_boot <- 2 * min(mean(rates <= 0, na.rm = TRUE), mean(rates >= 0, na.rm = TRUE))
  p_boot <- min(1, p_boot)
  tibble(rate = rate,
         ci_lo = unname(ci[1]), ci_hi = unname(ci[2]),
         p_boot = p_boot,
         significant = !(ci[1] <= 0 && ci[2] >= 0))
}

# ---------------------------------------------------------------------
# Section 10. Formulas 10-13 - the PGR family and the achievement rate AR
#   Formula 10 PGR_sim = 100×(Ŷ − v_2020)/v_2020
#   Formula 11 PGR_req: threshold target Y=100; otherwise Y = v_2020 + |m|.n (Y >= v_2020)
#   Formula 12 PGR gap = PGR_req - PGR_sim (positive = shortfall; negative = surplus)
#   Formula 13 I = 1{PGR_sim >= PGR_req and a rate of change significantly != 0 (bootstrap p < 0.05)}; AR = mean(I)
# ---------------------------------------------------------------------

required_target_value <- function(v2020, hist_annual_change, years_to_target,
                                  has_threshold) {
  ifelse(has_threshold, 100,
         # extrapolate |mean historical change| and make sure Y >= the historical value (reversing a downward trend)
         pmax(v2020 + abs(hist_annual_change) * years_to_target, v2020))
}

compute_pgr <- function(proj_df) {
  # proj_df must contain: country, layer, layer_id, baseline_value(=v_2020),
  #               projected, projected_year, hist_annual_change
  proj_df %>%
    mutate(
      has_threshold = layer_id %in% THRESHOLD_TARGETS,
      required_value = required_target_value(baseline_value, hist_annual_change,
                                             projected_year - 2020, has_threshold),
      pgr_sim = ifelse(baseline_value > 0,
                       100 * (projected - baseline_value) / baseline_value, NA_real_),
      pgr_req = ifelse(baseline_value > 0,
                       100 * (required_value - baseline_value) / baseline_value, NA_real_),
      pgr_gap = pgr_req - pgr_sim           # Formula 12: positive = shortfall, negative = surplus
    )
}

# Formula 13: the indicator function I and the achievement rate AR
#   I = 1 if and only if PGR_sim >= PGR_req and the bootstrap p of the annual rate of change < 0.05 (significant != 0)
compute_ar <- function(pgr_df) {
  pgr_df %>%
    mutate(
      I = as.integer(!is.na(pgr_sim) & !is.na(pgr_req) &
                       pgr_sim >= pgr_req &
                       !is.na(p_boot) & p_boot < 0.05)
    )
}

# AR aggregated by layer (the mean over the target/goal/indicator levels gives the overall achievement rate)
summarise_ar <- function(ar_df, by = c("layer", "layer_id")) {
  ar_df %>%
    group_by(across(all_of(by))) %>%
    summarise(AR = mean(I, na.rm = TRUE), n_countries = sum(!is.na(I)),
              .groups = "drop")
}

# Bootstrap interval of the AR over the resampling of countries (SI A.5: the 50% and 95% quantile intervals)
ar_bootstrap_ci <- function(ar_df, n_boot = 2000, seed = 1) {
  set.seed(seed)
  ar_df %>%
    group_by(layer, layer_id) %>%
    group_modify(~ {
      x <- .x$I[!is.na(.x$I)]
      n <- length(x)
      if (n == 0) return(tibble(AR = NA_real_, ci50_lo = NA_real_, ci50_hi = NA_real_,
                                ci95_lo = NA_real_, ci95_hi = NA_real_, n = 0L))
      boots <- replicate(n_boot, mean(sample(x, n, replace = TRUE)))
      q50 <- quantile(boots, c(0.25, 0.75), na.rm = TRUE)
      q95 <- quantile(boots, c(0.025, 0.975), na.rm = TRUE)
      tibble(AR = mean(x),
             ci50_lo = unname(q50[1]), ci50_hi = unname(q50[2]),
             ci95_lo = unname(q95[1]), ci95_hi = unname(q95[2]), n = n)
    }) %>%
    ungroup()
}

# ---------------------------------------------------------------------
# Section 11. A.5 trajectory classification - rho = PGR_sim / PGR_req
#   Full achievement   rho >= 1 and a significantly positive rate of change
#   Favourable         0.5 <= rho < 1 and a significantly positive rate of change
#   Mixed              a significant change but 0 < rho < 0.5, or rho >= 1 but not significant
#   Unfavourable       rho <= 0, or not significant and rho < 1
# ---------------------------------------------------------------------

classify_rho <- function(rho, significant_positive) {
  dplyr::case_when(
    is.na(rho)                          ~ NA_character_,
    rho >= RHO_FULL & significant_positive ~ "Full achievement",
    rho >= RHO_FAV  & significant_positive ~ "Favourable",
    (significant_positive & rho > 0 & rho < RHO_FAV) |
      (rho >= RHO_FULL & !significant_positive)        ~ "Mixed",
    rho <= 0 | !significant_positive                   ~ "Unfavourable",
    TRUE                                               ~ "Mixed"
  )
}

# Trajectory classification by layer: rho = mean(PGR_sim)/mean(PGR_req); a significantly positive change is decided by the bootstrap
classify_trajectory <- function(pgr_df) {
  pgr_df %>%
    filter(!is.na(pgr_sim), !is.na(pgr_req)) %>%
    group_by(layer, layer_id) %>%
    summarise(
      pgr_sim_mean = mean(pgr_sim), pgr_sim_sd = sd(pgr_sim),
      pgr_req_mean = mean(pgr_req), pgr_req_sd = sd(pgr_req),
      rate_mean = mean(rate, na.rm = TRUE),
      sig_share = mean(significant %in% TRUE & rate > 0, na.rm = TRUE),
      p_boot = mean(p_boot, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      rho = ifelse(pgr_req_mean != 0, pgr_sim_mean / pgr_req_mean,
                   ifelse(pgr_sim_mean > 0, Inf, -Inf)),
      significant_positive = sig_share > 0.5 & rate_mean > 0 &
        !is.na(p_boot) & p_boot < 0.05,
      # A.5 special case: PGR_req = 0 (zero requirement, rho undefined) and non-negative growth -> recorded as achieved
      class = ifelse(pgr_req_mean == 0 & pgr_sim_mean >= 0,
                     "Full achievement",
                     classify_rho(rho, significant_positive))
    )
}

# Confidence flag for the borderline cases: the class is unchanged under a +/-0.1 perturbation of the rho thresholds (0.5, 1) -> stable
classify_rho_perturbed <- function(rho, significant_positive,
                                   rho_full = RHO_FULL, rho_fav = RHO_FAV) {
  dplyr::case_when(
    is.na(rho) ~ NA_character_,
    rho >= rho_full & significant_positive ~ "Full achievement",
    rho >= rho_fav  & significant_positive ~ "Favourable",
    (significant_positive & rho > 0 & rho < rho_fav) |
      (rho >= rho_full & !significant_positive) ~ "Mixed",
    rho <= 0 | !significant_positive           ~ "Unfavourable",
    TRUE ~ "Mixed"
  )
}

trajectory_class_with_flag <- function(class_df) {
  class_df %>%
    mutate(
      class_lo = classify_rho_perturbed(rho, significant_positive,
                                        RHO_FULL - RHO_PERTURB, RHO_FAV - RHO_PERTURB),
      class_hi = classify_rho_perturbed(rho, significant_positive,
                                        RHO_FULL + RHO_PERTURB, RHO_FAV + RHO_PERTURB),
      class_confidence = ifelse(class == class_lo & class == class_hi,
                                "stable", "boundary")
    )
}

# ---------------------------------------------------------------------
# Section 12. Formulas 5-9 - trend alignment and country contributions (Liu et al. 2024)
#   all the quantities are computed on the normalized values r in [0,100]
#   eq. (9)  d_t^c = sqrt((r_{t+1}^c - r_t^c)^2) = |r_{t+1}^c - r_t^c|  domestic change
#   eq. (10) d_t = sqrt(Sum_c (r_{t+1}^c - r_t^c)^2)                      foreign change
#   eq. (11) S_t = (1/C) x (Sum_c d_t^c - d_t)                          synergy effect
#   eq. (12) g_t = (1/C) Sum_c (r_{t+1}^c - r_t^c);
#          sgn(S_t) = +1 if sgn(g_t) agrees with the expected direction of the indicator, otherwise -1   synergy direction
#   eq. (13) CSC_t^c = (d_t^c - d_t/C)/(S_t x sgn(S_t) x C) x 100%        country contribution
#          the sum over countries = +100% (the gross synergy is positive) or -100% (it is negative);
#          the country-level CSC = mean across indicators and years (Supplementary Table 8)
#   Anchors of the worked example of the docx: d=(10,-5,10) -> d_t=15, S_t=(25-15)/3=3.33,
#          g_t=+5 -> sgn(S_t)=+1, CSC=(50%,0%,50%); for the reversed segment sgn(S_t)=-1,
#          CSC=(−50%,0%,−50%)
# ---------------------------------------------------------------------

synergy_decompose <- function(norm_df, meta = NULL, coverage_df = NULL) {
  if (is.null(meta)) meta <- build_ind_meta()
  # A.2 retention rule: countries with more than 20% of the indicators missing do not enter the rankings (the country contributions of the CSC / the synergy are likewise
  # excluded, with the same convention as MIS/EIS/rankings). It applies when coverage_df provides keep_index.
  if (!is.null(coverage_df)) {
    norm_df <- norm_df %>%
      semi_join(coverage_df %>% filter(keep_index) %>% select(country),
                by = "country")
  }
  dir_map <- meta %>% filter(used) %>%
    distinct(indicator, .keep_all = TRUE) %>%
    transmute(indicator, expected_dir = ifelse(direction == "higher", 1, -1))

  # year-on-year pairs (t, t+1), limited to the 2011-2024 analysis window
  pairs <- norm_df %>%
    filter(year %in% TRAJ_YEARS) %>%
    arrange(indicator, country, year) %>%
    group_by(indicator, country) %>%
    mutate(r_next = lead(norm), year_next = lead(year)) %>%
    filter(!is.na(r_next), year_next == year + 1) %>%
    ungroup() %>%
    mutate(delta = r_next - norm,        # signed change
           d_abs = abs(delta))           # eq. (9) d_t^c = |Dr|

  by_ind_year <- pairs %>%
    group_by(indicator, year) %>%
    summarise(
      C = n_distinct(country),
      d_foreign = sqrt(sum(delta^2)),    # eq. (10) Euclidean distance
      sum_d_dom = sum(d_abs),            # Σ_c d_t^c
      g_t = mean(delta),                 # eq. (12) mean signed change
      .groups = "drop"
    ) %>%
    mutate(S_t = (sum_d_dom - d_foreign) / C) %>%   # eq. (11) divided by the number of countries C
    inner_join(dir_map, by = "indicator") %>%
    mutate(
      sgn_g = sign(g_t),
      # eq. (12): sgn(g_t) agrees with the expected direction -> +1, otherwise -1
      sgn_S = ifelse(sgn_g == expected_dir, 1L, -1L)
    )

  # eq. (13) country synergy contribution CSC
  csc <- pairs %>%
    select(indicator, year, country, DC = d_abs, delta) %>%
    inner_join(by_ind_year %>% select(indicator, year, C, d_foreign, S_t, sgn_S),
               by = c("indicator", "year")) %>%
    mutate(
      CSC_pct = ifelse(S_t != 0,
                       (DC - d_foreign / C) / (S_t * sgn_S * C) * 100,
                       NA_real_)
    )

  # Country-level CSC: mean across indicators and years (Supplementary Table 8)
  csc_country <- csc %>%
    group_by(country) %>%
    summarise(
      CSC_mean = mean(CSC_pct, na.rm = TRUE),
      CSC_sd   = sd(CSC_pct, na.rm = TRUE),
      DC_mean  = mean(DC, na.rm = TRUE),
      foreign_ref_mean = mean(d_foreign / C, na.rm = TRUE),
      n_obs = n(),
      .groups = "drop"
    ) %>%
    mutate(rank = rank(-CSC_mean, ties.method = "min"),
           percentile = 100 * (1 - (rank - 1) / n())) %>%
    arrange(rank)

  list(dc = pairs %>% select(indicator, year, country, DC = d_abs),
       fc = by_ind_year %>% select(indicator, year, d_foreign),
       sye = by_ind_year %>% select(indicator, year, S_t, g_t),
       syd = by_ind_year %>% select(indicator, year, S_t, sgn_S, g_t, expected_dir),
       csc = csc,
       csc_country = csc_country)
}

# ---------------------------------------------------------------------
# Section 13. Projection at the indicator level (needed for Supplementary Table 10 / the achievement rate at the indicator level of Formula 13)
# ---------------------------------------------------------------------

project_indicator_level <- function(norm_df, meta = NULL, n_boot = N_BOOT_RATE) {
  if (is.null(meta)) meta <- build_ind_meta()
  ids <- intersect(unique(meta$indicator[meta$used]),
                   unique(norm_df$indicator))
  map_dfr(ids, function(id) {
    d <- norm_df %>% filter(indicator == id)
    countries <- unique(d$country)
    map_dfr(countries, function(cc) {
      dd <- d %>% filter(country == cc, !is.na(norm))   # the fit uses the complete historical series
      if (sum(!is.na(dd$norm)) < 5) return(NULL)
      pr <- project_with_rate(dd$year, dd$norm, TARGET_YEAR_TARGET, n_boot = n_boot)
      if (is.na(pr$proj$projected)) return(NULL)
      base <- dd %>% filter(year == 2020)
      v2020 <- if (nrow(base) > 0) base$norm[1] else NA_real_
      hist_ch <- dd %>% filter(year %in% BASELINE_YEARS) %>%
        arrange(year) %>%
        summarise(m = mean(diff(norm), na.rm = TRUE)) %>% pull(m)
      bind_cols(tibble(country = cc, indicator = id,
                       baseline_value = v2020,
                       hist_annual_change = hist_ch),
                pr$proj %>% mutate(projected_year = TARGET_YEAR_TARGET))
    })
  })
}

# ---------------------------------------------------------------------
# Section 14. Main entry point: the whole main chain runs with one call
# Checkpoint version: increment it when the logic changes, so that the old chk_*.rds files are invalidated automatically
PIPELINE_CHECKPOINT_VERSION <- "v8-goal-target-ts-full"

run_pipeline <- function(n_boot = N_BOOT_RATE, verbose = TRUE,
                         with_indicator_projection = TRUE,
                         checkpoint_dir = "data/intermediate",
                         resume = TRUE) {
  log_msg <- if (verbose) function(...) message("[pipeline] ", ...) else function(...) invisible(NULL)

  # Every stage is written to disk immediately as an RDS (resume after an interruption):
  #   resume = TRUE and the checkpoint exists and the version matches -> it is loaded and the stage is skipped
  #   otherwise the stage is computed and saved immediately, so that an error does not lose the completed stages
  if (!is.null(checkpoint_dir) && !dir.exists(checkpoint_dir)) {
    dir.create(checkpoint_dir, recursive = TRUE)
  }
  # Checkpoint version check: if the version does not match, everything is recomputed
  ver_f <- file.path(checkpoint_dir, "chk_00_version.rds")
  ver_write_later <- FALSE
  if (resume && file.exists(ver_f)) {
    if (!identical(readRDS(ver_f), PIPELINE_CHECKPOINT_VERSION)) {
      log_msg("    [checkpoint] version changed (", PIPELINE_CHECKPOINT_VERSION,
              "), the old checkpoints are invalid, everything is recomputed")
      resume <- FALSE
      # the version number is only updated after the full recomputation has succeeded (if an error occurs the old version number is kept,
      # so the next run is forced to recompute and old and new checkpoints are never mixed)
      ver_write_later <- TRUE
    }
  } else if (!is.null(checkpoint_dir)) {
    saveRDS(PIPELINE_CHECKPOINT_VERSION, ver_f)
  }
  stage <- function(name, compute) {
    f <- if (is.null(checkpoint_dir)) NULL else
      file.path(checkpoint_dir, paste0("chk_", name, ".rds"))
    if (resume && !is.null(f) && file.exists(f)) {
      log_msg("    [checkpoint] ", name, " already exists, it is loaded and the computation is skipped")
      return(readRDS(f))
    }
    val <- compute()
    if (!is.null(f)) {
      saveRDS(val, f)
      log_msg("    [checkpoint] ", name, " → ", f)
    }
    val
  }

  log_msg("1/9 metadata and raw data (the 44 indicators with used=Y, all historical years) ...")
  meta <- stage("01_meta", function() build_ind_meta())
  raw_all <- stage("02_raw", function() load_all_raw(meta))   # 44 indicators, full history
  country_meta <- stage("03_country_meta", function() build_country_meta())

  log_msg("2/9 imputation (all 44 indicators imputed with one common rule; rows created for the missing years + anchors from the full history; only 2011-2024 filled) ...")
  # Impute first: all 44 indicators, all entities, handled uniformly by the same rule
  imputed_all <- stage("04_imputed_all", function()
    impute_all(raw_all, impute_years = TRAJ_YEARS,
               used_ids = unique(meta$indicator[meta$used])))
  # Cell-by-cell imputation detail (the basis of imputation_detail.csv / Supplementary Table 9)
  impute_detail <- stage("04b_impute_detail", function()
    imputation_detail(raw_all %>% filter(indicator %in% unique(meta$indicator[meta$used])),
                      imputed_all, impute_years = TRAJ_YEARS))

  log_msg("2b/9 filtering after the imputation (country universe + routing of the regional indicators) ...")
  # Then filter: the country level = the 193 countries with geoclass United Nations == "nation",
  #          and without the regional indicators 27-29
  imputed <- stage("05_imputed_country", function()
    filter_to_nations(imputed_all, country_meta,
                      exclude_indicators = REGIONAL_ONLY_INDICATORS,
                      verbose = verbose))
  # Routing of the regional indicators 27-29: they are used only for the global/regional analyses
  imputed_regional <- stage("05b_imputed_regional", function()
    imputed_all %>% filter(indicator %in% REGIONAL_ONLY_INDICATORS))
  raw <- raw_all

  log_msg("3/9 country coverage filter (>=80% of the indicators) ...")
  coverage <- stage("05_coverage", function() country_indicator_coverage(imputed, meta))

  log_msg("4/9 normalization Formula 1 ...")
  norm <- stage("06_norm", function() build_normalized(imputed, meta))

  # Normalization of the regional indicators (27-29, for the global/regional analyses)
  norm_regional <- stage("06b_norm_regional", function() {
    if (nrow(imputed_regional) == 0) return(imputed_regional)
    build_normalized(imputed_regional, meta, indicators = REGIONAL_ONLY_INDICATORS)
  })

  log_msg("5/9 target/goal values Formula 2 ...")
  layer <- stage("07_layer", function() aggregate_goal_target(norm, meta))

  # Regional and global aggregation of the target/goal time series (the country level = the layer itself)
  layer_geo <- stage("07b_layer_geography", function()
    aggregate_layer_geography(layer, country_meta))

  log_msg("6/9 MIS Formula 3 / EIS Formula 4 (annual, 2011-2024) ...")
  mis_eis <- stage("08_mis_eis", function() compute_mis_eis_annual(layer, coverage))

  # Regional and global aggregation of the annual MIS/EIS series
  mis_eis_geo <- stage("08b_mis_eis_geography", function()
    mis_eis_geography(mis_eis, country_meta))

  log_msg("7/9 baseline four parameters over 2011-2020 (Supplementary Table 3) ...")
  baseline <- stage("09_baseline", function() {
    list(indicator = baseline_params_indicator(norm),
         target_goal = baseline_params(layer),
         per_country = per_country_baseline_params(layer))
  })

  log_msg("8/9 projection + bootstrap + Formulas 10-13 ...")
  proj <- stage("10_projections", function() {
    # mean historical change (branch 2 of Formula 11): the mean of the year-on-year differences of the baseline period
    layer_hist <- layer %>%
      arrange(layer, layer_id, country, year) %>%
      group_by(layer, layer_id, country) %>%
      filter(year %in% BASELINE_YEARS) %>%
      summarise(hist_annual_change = mean(diff(index), na.rm = TRUE),
                baseline_value = index[year == 2020][1],
                .groups = "drop")

    keys <- layer_hist %>% distinct(layer, layer_id, country)
    proj_list <- vector("list", nrow(keys))
    path_list <- vector("list", nrow(keys))
    for (k in seq_len(nrow(keys))) {
      ly <- keys$layer[k]; lid <- keys$layer_id[k]; cc <- keys$country[k]
      # the projection fit uses the complete historical series (to identify the pattern); the target year is extrapolated
      d <- layer %>% filter(layer == ly, layer_id == lid,
                            country == cc, !is.na(index))
      ty <- if (ly == "goal") TARGET_YEAR_GOAL else TARGET_YEAR_TARGET
      pr <- project_with_rate(d$year, d$index, ty, n_boot = n_boot)
      proj_list[[k]] <- pr$proj %>%
        mutate(layer = ly, layer_id = lid, country = cc, projected_year = ty)
      if (nrow(pr$path) > 0) {
        path_list[[k]] <- pr$path %>%
          mutate(layer = ly, layer_id = lid, country = cc,
                 model_type = pr$proj$model_type[1],
                 model_form = pr$proj$model_form[1])
      }
    }
    proj_rows <- bind_rows(proj_list) %>%
      left_join(layer_hist %>% select(layer, layer_id, country,
                                      baseline_value, hist_annual_change),
                by = c("layer", "layer_id", "country"))
    proj_paths <- bind_rows(path_list)
    # Complete country-level target/goal time series: the existing data -> the end of the projection (2030/2050)
    ts_country <- goal_target_ts_full(layer, proj_paths, id_cols = "country")
    list(proj_rows = proj_rows, layer_hist = layer_hist,
         proj_paths = proj_paths, goal_target_ts = ts_country)
  })
  proj_rows <- proj$proj_rows

  # Complete MIS/EIS time series: the 2011-2024 history + the projected segment up to the target year (2030/2050)
  mis_eis_ts <- stage("10b_mis_eis_ts", function() {
    ts_c <- mis_eis_ts_country_full(mis_eis, proj$goal_target_ts, coverage)
    geo_ts <- mis_eis_geography(ts_c, country_meta)
    list(country = ts_c, region = geo_ts$region, global = geo_ts$global)
  })

  pgr_ar <- stage("11_pgr_ar", function() {
    pgr <- compute_pgr(proj_rows)
    ar  <- compute_ar(pgr)
    list(pgr = pgr, ar = ar,
         ar_summary = summarise_ar(ar),
         ar_ci = ar_bootstrap_ci(ar),
         trajectory_class = trajectory_class_with_flag(classify_trajectory(pgr)))
  })

  log_msg("9/9 Formulas 5-9 synergy/contributions ...")
  # The country contributions of the CSC are likewise removed according to the coverage rule (more than 20% missing excludes them from the rankings, A.2)
  syn <- stage("12_synergy", function() synergy_decompose(norm, meta, coverage))

  ind_proj <- NULL; ind_ar <- NULL
  if (with_indicator_projection) {
    log_msg("    projection at the indicator level (Supplementary Table 10 / the AR at the indicator level) ...")
    ind <- stage("13_ind_projection", function() {
      ip <- project_indicator_level(norm, meta, n_boot = n_boot)
      if (is.null(ip) || nrow(ip) == 0) return(list(proj = NULL, ar = NULL))
      # Indicator threshold: an explicit GBF threshold for any of the mapped targets makes it a threshold indicator (Formula 11)
      thr_ind <- load_i2t() %>%
        filter(layer == "target", layer_id %in% THRESHOLD_TARGETS) %>%
        distinct(indicator) %>% pull(indicator)
      ind_pgr <- ip %>%
        mutate(layer = "indicator",
               layer_id = paste0("ind", indicator)) %>%
        compute_pgr() %>%
        mutate(has_threshold = indicator %in% thr_ind,
               required_value = required_target_value(
                 baseline_value, hist_annual_change, projected_year - 2020,
                 has_threshold),
               pgr_req = ifelse(baseline_value > 0,
                                100 * (required_value - baseline_value) /
                                  baseline_value, NA_real_),
               pgr_gap = pgr_req - pgr_sim)
      list(proj = ind_pgr, ar = compute_ar(ind_pgr))
    })
    ind_proj <- ind$proj; ind_ar <- ind$ar
  }

  # The checkpoint version number is written after all the stages have completed
  if (ver_write_later && !is.null(checkpoint_dir)) {
    saveRDS(PIPELINE_CHECKPOINT_VERSION, ver_f)
  }

  list(
    meta = meta,
    country_meta = country_meta,
    raw = raw_all,
    imputed_all = imputed_all,
    imputation_detail = impute_detail,
    imputed = imputed,
    imputed_regional = imputed_regional,
    norm_regional = norm_regional,
    coverage = coverage,
    norm_indicator = norm,
    norm_layer = layer,
    layer_geography = layer_geo,
    goal_target_ts = proj$goal_target_ts,
    projection_paths = proj$proj_paths,
    mis_eis_annual = mis_eis,
    mis_eis_geography = mis_eis_geo,
    mis_eis_ts = mis_eis_ts,
    baseline_indicator = baseline$indicator,
    baseline_target_goal = baseline$target_goal,
    baseline_per_country = baseline$per_country,
    projections = proj_rows,
    pgr = pgr_ar$pgr,
    ar = pgr_ar$ar,
    ar_summary = pgr_ar$ar_summary,
    ar_ci = pgr_ar$ar_ci,
    trajectory_class = pgr_ar$trajectory_class,
    indicator_projection = ind_proj,
    indicator_ar = ind_ar,
    synergy = syn
  )
}

# ---------------------------------------------------------------------
# Section 15. Writing the results to disk
# ---------------------------------------------------------------------

save_intermediates <- function(res, dir = "data/intermediate") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  save_one <- function(obj, name) {
    if (is.null(obj)) return(invisible(NULL))
    saveRDS(obj, file.path(dir, paste0(name, ".rds")))
    if (is.data.frame(obj)) write_csv(obj, file.path(dir, paste0(name, ".csv")))
  }
  save_one(res$country_meta, "country_meta")
  save_one(res$imputed_all, "imputed_all")
  save_one(res$imputation_detail, "imputation_detail")
  save_one(res$imputed_regional, "imputed_regional")
  save_one(res$norm_regional, "norm_regional")
  save_one(res$norm_indicator, "norm_indicator")
  save_one(res$norm_layer, "norm_layer")
  # Target/goal time series: country / regional / global (existing data -> end of the projection 2030/2050)
  if (!is.null(res$goal_target_ts)) {
    save_one(res$goal_target_ts, "goal_target_ts_country")
  } else if (!is.null(res$norm_layer)) {
    save_one(res$norm_layer, "goal_target_ts_country")
  }
  if (!is.null(res$layer_geography)) {
    save_one(res$layer_geography$region_ts, "goal_target_ts_region")
    save_one(res$layer_geography$global_ts, "goal_target_ts_global")
  }
  # Annual MIS/EIS series: country / regional / global (2011-2024, the B.2.7 window)
  save_one(res$mis_eis_annual, "mis_eis_annual")
  if (!is.null(res$mis_eis_geography)) {
    save_one(res$mis_eis_geography$region, "mis_eis_annual_region")
    save_one(res$mis_eis_geography$global, "mis_eis_annual_global")
  }
  # Complete MIS/EIS time series: country / regional / global (history + projected segment up to the target year)
  if (!is.null(res$mis_eis_ts)) {
    save_one(res$mis_eis_ts$country, "mis_eis_ts_country")
    save_one(res$mis_eis_ts$region, "mis_eis_ts_region")
    save_one(res$mis_eis_ts$global, "mis_eis_ts_global")
  }
  save_one(res$projections, "projections")
  save_one(res$pgr, "pgr")
  save_one(res$ar, "ar")
  save_one(res$ar_summary, "ar_summary")
  save_one(res$ar_ci, "ar_ci")
  save_one(res$trajectory_class, "trajectory_class")
  save_one(res$indicator_projection, "indicator_projection")
  save_one(res$synergy$csc_country, "csc_country")
  save_one(res$baseline_target_goal, "baseline_target_goal")
  save_one(res$baseline_per_country, "baseline_per_country")
  save_one(res$coverage, "country_coverage")
  message("[pipeline] intermediates saved to ", dir)
}

# Restore the res list of run_pipeline() from the checkpoints (so that run_remaining() can be called without arguments)
load_res_from_checkpoints <- function(checkpoint_dir = "data/intermediate") {
  get1 <- function(name) {
    f <- file.path(checkpoint_dir, paste0("chk_", name, ".rds"))
    if (!file.exists(f)) return(NULL)
    readRDS(f)
  }
  if (is.null(get1("01_meta"))) {
    stop("no checkpoint found (", checkpoint_dir, "/chk_*.rds); please run run_pipeline() first")
  }
  baseline <- get1("09_baseline")
  proj     <- get1("10_projections")
  pgr_ar   <- get1("11_pgr_ar")
  ind      <- get1("13_ind_projection")
  list(
    meta                  = get1("01_meta"),
    raw                   = get1("02_raw"),
    country_meta          = get1("03_country_meta"),
    imputed_all           = get1("04_imputed_all"),
    imputation_detail     = get1("04b_impute_detail"),
    imputed               = get1("05_imputed_country"),
    imputed_regional      = get1("05b_imputed_regional"),
    coverage              = get1("05_coverage"),
    norm_indicator        = get1("06_norm"),
    norm_regional         = get1("06b_norm_regional"),
    norm_layer            = get1("07_layer"),
    layer_geography       = get1("07b_layer_geography"),
    mis_eis_annual        = get1("08_mis_eis"),
    mis_eis_geography     = get1("08b_mis_eis_geography"),
    mis_eis_ts            = get1("10b_mis_eis_ts"),
    baseline_indicator    = if (!is.null(baseline)) baseline$indicator else NULL,
    baseline_target_goal  = if (!is.null(baseline)) baseline$target_goal else NULL,
    baseline_per_country  = if (!is.null(baseline)) baseline$per_country else NULL,
    projections           = if (!is.null(proj)) proj$proj_rows else NULL,
    projection_paths      = if (!is.null(proj)) proj$proj_paths else NULL,
    goal_target_ts        = if (!is.null(proj)) proj$goal_target_ts else NULL,
    pgr                   = if (!is.null(pgr_ar)) pgr_ar$pgr else NULL,
    ar                    = if (!is.null(pgr_ar)) pgr_ar$ar else NULL,
    ar_summary            = if (!is.null(pgr_ar)) pgr_ar$ar_summary else NULL,
    ar_ci                 = if (!is.null(pgr_ar)) pgr_ar$ar_ci else NULL,
    trajectory_class      = if (!is.null(pgr_ar)) pgr_ar$trajectory_class else NULL,
    indicator_projection  = if (!is.null(ind)) ind$proj else NULL,
    indicator_ar          = if (!is.null(ind)) ind$ar else NULL,
    synergy               = get1("12_synergy")
  )
}

# Usage: sourcing this file runs the whole pipeline and writes the results to disk (default n_boot = 10000).
#   To load the functions without running them: set GBF_NO_AUTORUN <- TRUE before sourcing,
#   or call res <- run_pipeline(); save_intermediates(res) explicitly
if (!exists("GBF_NO_AUTORUN") || !isTRUE(GBF_NO_AUTORUN)) {
  res <- run_pipeline()
  save_intermediates(res)
}
