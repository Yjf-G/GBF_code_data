# =====================================================================
# GBF paper - derived analyses module
# Algorithmic reference: the Supplementary Information of the submission (2026-09-30 version)
# Companion: gbf_compute_pipeline.R (run res <- run_pipeline() first)
#
# Products (one-to-one with the numbering of the Supplementary Information):
#   Table S1 -> Supplementary Table 10    projection model per indicator (series type CV / selected form / AICc / breakpoint year)
#   Table S2 -> Supplementary Table 3     baseline regression + geographic heterogeneity (slope/R2/p + KW H + MWU-Bonferroni
#               + global BH-FDR, the note of Supplementary Table 3)
#   Table S3 -> Supplementary Table 4     post-2020 trajectories, achievement and classification (mean(SD) of PGR_sim/req / rho /
#               bootstrap p / class / confidence flag, A.5)
#   Table S4 -> Supplementary Table 7     annual MIS/EIS series by region and development group (2011-2024, B.2.7)
#   Table S5 -> Supplementary Table 8     the country synergy contributions CSC (country rows, domestic/foreign decomposition)
#   Table S6 -> Supplementary Table 9     imputation records per indicator (A.2)
#   b6_*        section/goal summaries (trend slope / R2 / AR + bootstrap intervals + KW; Fig. 2c)
#   b7_*        trajectory analyses (SI B.2.7): share of monotonicity / acceleration / endpoint bias / COVID-2020 /
#               leave-one-year-out / regional annual paths
#   c1_*        50%/95% bootstrap intervals of the achievement rate (A.5; overall/target/goal/section)
#   uncertainty__a1..a5  uncertainty, Algorithms 1-5 (SI B.1/B.2.1)
#   uncertainty__c3  weight sensitivity (SI B.2.2; Dirichlet 1,000 / one PCA component / equal weights; dropping indicator 46)
#   (dropped-candidate sensitivity) - not implemented (the 9 candidates were dropped because of the years/coverage criteria)
#   c5_*        endpoint vs trajectory robustness (SI B.2.3/B.2.4; extrapolation test / leave-one-year-out)
#   Fig 1/3/4   figure source data (slope x R2 quadrant / MIS-EIS KW + year-on-year paired / SyE binomial test)
#
# Run:
#   source("gbf_compute_pipeline.R")
#   source("gbf_compute_remaining.R")
#   res  <- run_pipeline()
#   rest <- run_remaining(res)
# =====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# Make sure the main pipeline is loaded (this module depends on load_res_from_checkpoints / fit_pma /
# calculate_snr and the like). When this module is sourced on its own, gbf_compute_pipeline.R is loaded automatically.
if (!all(vapply(c("load_res_from_checkpoints", "fit_pma", "calculate_snr"),
                exists, logical(1), mode = "function"))) {
  pf <- "gbf_compute_pipeline.R"
  if (!file.exists(pf)) pf <- file.path(getwd(), "gbf_compute_pipeline.R")
  if (!file.exists(pf)) {
    stop("gbf_compute_pipeline.R not found; please source the main pipeline before sourcing this module")
  }
  source(pf)   # note: the end of the main pipeline contains the automatic run_pipeline() entry point
}

options(scipen = 10, dplyr.summarise.inform = FALSE)

OUTPUT_DIR <- "data/output"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Safe wrapper for the tests: returns NA instead of aborting when the data are degenerate
safe_stat <- function(expr) tryCatch(expr, error = function(e) NA_real_)
safe_kwp <- function(x, g) safe_stat(kruskal.test(x ~ factor(g))$p.value)
safe_kwh <- function(x, g) safe_stat(unname(kruskal.test(x ~ factor(g))$statistic))
safe_wilcox <- function(...) safe_stat(wilcox.test(...)$p.value)
safe_ks <- function(x, y) safe_stat(ks.test(x, y)$p.value)

# Quantiles of the B.2.7 acceleration flag (the place where the SI threshold is still to be computed: the distribution quantile across indicators)
ACCEL_QUANTILES <- c(0.75, 0.90)

# ---------------------------------------------------------------------
# Section 1. Supplementary Table 10 (internal product table_S1) - the projection model per indicator
#   series type (CV) / selected form / AICc / breakpoint year; the piecewise description of OHI and the like
#   the AICc form selection is done on the global mean series of the indicator (2011-2024)
# ---------------------------------------------------------------------

table_s1_projection_models <- function(norm_df, norm_regional = NULL) {
  # The projection model is fitted on the complete historical series (identical to the projection of run_pipeline);
  # the regional indicators 27-29 enter the table through their regional series
  # PMA is a multi-model averaging: selected_form records the model with the highest weight,
  # notes records the number of models in the confidence set and the highest weight
  if (!is.null(norm_regional) && nrow(norm_regional) > 0) {
    norm_df <- bind_rows(norm_df, norm_regional %>% select(any_of(names(norm_df))))
  }
  norm_df %>%
    group_by(indicator, year) %>%
    summarise(v = mean(norm, na.rm = TRUE), .groups = "drop") %>%
    group_by(indicator) %>%
    group_modify(~ {
      d <- .x %>% arrange(year) %>% filter(is.finite(v))
      if (nrow(d) < 5) {
        return(tibble(series_type = NA_character_, cv = NA_real_,
                      snr = NA_real_,
                      selected_form = NA_character_, aicc = NA_real_,
                      n_models_averaged = NA_integer_, top_weight = NA_real_,
                      breakpoint_year = NA_integer_, notes = "series too short"))
      }
      cv <- sd(d$v) / abs(mean(d$v))
      snr <- calculate_snr(d$v)      # the noise criterion of Tittensor et al.
      series_type <- if (snr < SNR_DLM_THRESHOLD) "high-noise (DLM)" else "low-noise (PMA)"
      if (snr < SNR_DLM_THRESHOLD) {
        return(tibble(series_type = series_type, cv = round(cv, 3),
                      snr = round(snr, 3),
                      selected_form = "DLM", aicc = NA_real_,
                      n_models_averaged = 1L, top_weight = 1,
                      breakpoint_year = NA_integer_, notes = "dynamic linear model"))
      }
      pma <- fit_pma(d$year, d$v)
      if (is.null(pma)) {
        return(tibble(series_type = series_type, cv = round(cv, 3),
                      snr = round(snr, 3),
                      selected_form = NA_character_, aicc = NA_real_,
                      n_models_averaged = NA_integer_, top_weight = NA_real_,
                      breakpoint_year = NA_integer_, notes = "fit failed"))
      }
      notes <- paste0("PMA averaged over ", pma$n_avg, " models (delta AICc<2)")
      if (!is.na(pma$breakpoint)) {
        notes <- paste0(notes, "; segmented breakpoint ", pma$breakpoint)
      }
      tibble(series_type = series_type, cv = round(cv, 3),
             snr = round(snr, 3),
             selected_form = pma$top,
             aicc = round(unname(pma$all_aicc[pma$top]), 2),
             n_models_averaged = pma$n_avg,
             top_weight = round(pma$top_weight, 3),
             breakpoint_year = pma$breakpoint, notes = notes)
    }) %>%
    ungroup()
}

# ---------------------------------------------------------------------
# Section 1b. Global/regional summary of the regional indicators (27-29)
#   Used only in the global/regional analyses; they take no part in any country-level statistic
# ---------------------------------------------------------------------

regional_indicator_summary <- function(norm_regional) {
  if (is.null(norm_regional) || nrow(norm_regional) == 0) return(NULL)
  # The entity rows are the regions (GeoAreaName is a continent/region name) -> aggregated by entity x year
  by_entity <- norm_regional %>%
    group_by(indicator, entity = country, year) %>%
    summarise(norm = mean(norm, na.rm = TRUE), .groups = "drop")
  global <- norm_regional %>%
    group_by(indicator, year) %>%
    summarise(entity = "Global", norm = mean(norm, na.rm = TRUE),
              n_regions = n_distinct(country), .groups = "drop")
  bind_rows(by_entity %>% mutate(n_regions = NA_integer_), global) %>%
    arrange(indicator, entity, year)
}

# ---------------------------------------------------------------------
# Section 2. Supplementary Table 3 (internal product table_S2) - baseline regression + geographic heterogeneity
#   slope / R² / p + KW H + post-hoc (mean diff, Bonferroni p) + BH-FDR
#   FDR: Benjamini-Hochberg global correction over all the reported pairwise tests (the note of Supplementary Table 3)
# ---------------------------------------------------------------------

bh_fdr <- function(p) {
  ok <- !is.na(p)
  q <- rep(NA_real_, length(p))
  if (sum(ok) > 0) q[ok] <- p.adjust(p[ok], method = "BH")
  q
}

table_s2_baseline_geo <- function(res) {
  geo <- baseline_kw_posthoc(res$baseline_per_country, res$country_meta)
  posthoc <- geo$posthoc_bonferroni %>%
    mutate(u_p_fdr = bh_fdr(u_p_raw))   # global BH-FDR check (the note of Supplementary Table 3)
  list(
    baseline_regression = res$baseline_target_goal %>%
      select(layer, layer_id, n_countries, slope, r2, p_value, mean_v, sd_v, range_v),
    kw = geo$kw,
    posthoc = posthoc
  )
}

# ---------------------------------------------------------------------
# Section 3. Supplementary Table 4 (internal product table_S3) - post-2020 trajectories, achievement and classification (A.5)
#   target/goal; section; PGR_sim mean (SD); PGR_req mean (SD); ρ;
#   bootstrap p; class; class confidence flag
# ---------------------------------------------------------------------

table_s3_trajectories <- function(res) {
  ts <- load_target_section()
  sec_map <- bind_rows(
    ts %>% select(layer_id, section),
    tibble(layer_id = c("A", "B", "C", "D"),
           section = "Goals")
  )
  res$trajectory_class %>%
    left_join(sec_map, by = "layer_id") %>%
    transmute(layer, layer_id, section,
              pgr_sim_mean = round(pgr_sim_mean, 1),
              pgr_sim_sd = round(pgr_sim_sd, 1),
              pgr_req_mean = round(pgr_req_mean, 1),
              pgr_req_sd = round(pgr_req_sd, 1),
              rho = round(rho, 2),
              bootstrap_p = signif(p_boot, 3),
              class, class_confidence, n)
}

# ---------------------------------------------------------------------
# Section 4. Supplementary Table 7 (internal product table_S4) - MIS/EIS by region and development group, annual series 2011-2024
#   group; n; MIS(2011,2014,2017,2020,2024); the same for EIS; the summary columns have one set of MIS/EIS each:
#   delta annualized; delta headline (2011->2024); KS test p (the 2011 vs 2024 country distribution within the group);
#   share of monotonic trajectories (the country series over 2011-2024 without a reversal of sign)
# ---------------------------------------------------------------------

# a constant sign of the year-on-year differences = monotonic (a plateau is allowed)
is_monotonic <- function(v) {
  d <- diff(v[!is.na(v)])
  d <- d[d != 0]
  if (length(d) == 0) return(TRUE)
  all(d > 0) || all(d < 0)
}

group_membership <- function(country_meta) {
  # The grouping convention comes from data/geoclass.xlsx (UN M49):
  #   continent / subcontinent / region / developed-developing / SIDS / LDC / LLDC / OECD / EU
  bind_rows(
    country_meta %>% filter(!is.na(continent)) %>%
      transmute(country, group = paste0("continent: ", continent)),
    country_meta %>% filter(!is.na(subcontinent)) %>%
      transmute(country, group = paste0("subcontinent: ", subcontinent)),
    country_meta %>% filter(!is.na(region)) %>%
      transmute(country, group = paste0("region: ", region)),
    country_meta %>% filter(!is.na(development)) %>%
      transmute(country, group = paste0("development: ", development)),
    country_meta %>% filter(sids %in% TRUE) %>%
      transmute(country, group = "SIDS"),
    country_meta %>% filter(ldc %in% TRUE) %>%
      transmute(country, group = "LDC"),
    country_meta %>% filter(llc %in% TRUE) %>%
      transmute(country, group = "LLDC"),
    country_meta %>% filter(oecd %in% TRUE) %>%
      transmute(country, group = "OECD"),
    if ("eu" %in% names(country_meta)) {
      country_meta %>% filter(eu %in% TRUE) %>% transmute(country, group = "EU")
    } else NULL
  ) %>% distinct()
}

table_s4_mis_eis_groups <- function(res) {
  gm <- group_membership(res$country_meta)
  mis <- res$mis_eis_annual %>% select(country, year, MIS, EIS)
  yrs <- sort(unique(mis$year))   # the complete 2011-2024 series
  mis %>%
    inner_join(gm, by = "country") %>%
    group_by(group) %>%
    group_modify(~ {
      d <- .x
      n_cty <- n_distinct(d$country)
      wide_mis <- d %>% select(country, year, MIS) %>%
        pivot_wider(names_from = year, values_from = MIS, names_prefix = "MIS_")
      wide_eis <- d %>% select(country, year, EIS) %>%
        pivot_wider(names_from = year, values_from = EIS, names_prefix = "EIS_")
      pick <- function(w, prefix, yr) {
        col <- paste0(prefix, yr)
        if (col %in% names(w)) mean(w[[col]], na.rm = TRUE) else NA_real_
      }
      # delta annualized: the mean over the countries of (2024-2011)/13; delta headline: the difference of the endpoints of the annual means
      ann <- d %>% select(country, year, MIS) %>%
        pivot_wider(names_from = year, values_from = MIS) %>%
        mutate(delta = .data[["2024"]] - .data[["2011"]])
      d_ann <- mean(ann$delta, na.rm = TRUE) / 13
      ann_e <- d %>% select(country, year, EIS) %>%
        pivot_wider(names_from = year, values_from = EIS) %>%
        mutate(delta = .data[["2024"]] - .data[["2011"]])
      d_ann_eis <- mean(ann_e$delta, na.rm = TRUE) / 13
      mis2011 <- mean(d$MIS[d$year == 2011], na.rm = TRUE)
      mis2024 <- mean(d$MIS[d$year == 2024], na.rm = TRUE)
      eis2011 <- mean(d$EIS[d$year == 2011], na.rm = TRUE)
      eis2024 <- mean(d$EIS[d$year == 2024], na.rm = TRUE)
      # KS: the 2011 vs 2024 country distribution within the group
      ks_p <- tryCatch(
        safe_ks(d$MIS[d$year == 2011], d$MIS[d$year == 2024]),
        error = function(e) NA_real_)
      ks_p_eis <- tryCatch(
        safe_ks(d$EIS[d$year == 2011], d$EIS[d$year == 2024]),
        error = function(e) NA_real_)
      mono_share <- d %>% select(country, year, MIS) %>%
        pivot_wider(names_from = year, values_from = MIS) %>%
        rowwise() %>%
        mutate(mono = is_monotonic(c_across(starts_with("20")))) %>%
        ungroup() %>%
        summarise(m = mean(mono)) %>% pull(m)
      mono_share_eis <- d %>% select(country, year, EIS) %>%
        pivot_wider(names_from = year, values_from = EIS) %>%
        rowwise() %>%
        mutate(mono = is_monotonic(c_across(starts_with("20")))) %>%
        ungroup() %>%
        summarise(m = mean(mono)) %>% pull(m)
      out_mis <- set_names(lapply(yrs, function(y) pick(wide_mis, "MIS_", y)),
                          paste0("MIS_", yrs))
      out_eis <- set_names(lapply(yrs, function(y) pick(wide_eis, "EIS_", y)),
                          paste0("EIS_", yrs))
      bind_cols(tibble(n = n_cty), as_tibble(out_mis), as_tibble(out_eis),
                tibble(
                  delta_annualized = d_ann,
                  delta_headline = mis2024 - mis2011,
                  ks_p = signif(ks_p, 3),
                  monotonic_share = mono_share,
                  delta_annualized_eis = d_ann_eis,
                  delta_headline_eis = eis2024 - eis2011,
                  ks_p_eis = signif(ks_p_eis, 3),
                  monotonic_share_eis = mono_share_eis
                ))
    }) %>%
    ungroup()
}

# ---------------------------------------------------------------------
# Section 5. Supplementary Table 8 (internal product table_S5) - the country synergy contributions CSC (Formula 9)
#   country; CSC mean (SD); rank; percentile; domestic/foreign decomposition; the China row is enforced
# ---------------------------------------------------------------------

table_s5_csc <- function(res) {
  csc_cty <- res$synergy$csc_country
  china <- csc_cty %>% filter(grepl("^china$", tolower(country)))
  out <- csc_cty %>%
    transmute(country,
              CSC_mean = round(CSC_mean, 2),
              CSC_sd = round(CSC_sd, 2),
              rank, percentile = round(percentile, 1),
              DC_mean = round(DC_mean, 2),
              foreign_ref_mean = round(foreign_ref_mean, 2),
              n_obs)
  # The China row has to be present
  if (nrow(china) == 0) {
    out <- bind_rows(out, tibble(country = "China", CSC_mean = NA_real_,
                                 CSC_sd = NA_real_, rank = NA_integer_,
                                 percentile = NA_real_, DC_mean = NA_real_,
                                 foreign_ref_mean = NA_real_, n_obs = NA_integer_))
  }
  out
}

# ---------------------------------------------------------------------
# Section 6. Supplementary Table 9 (internal product table_S6) - imputation records per indicator (A.2)
#   aggregated from the cell-by-cell detail of imputation_detail:
#   observed cells / imputed cells (linear interpolation / forward fill / backward fill) / cells that cannot be imputed /
#   rows created for the missing years / fills computed outside the window but not retained following A.2.
#   The indicators without year columns (18 every five years, 39 every ten years, 26/30 every other year, etc.) are
#   correctly reflected here as rows created for the missing years + imputation from the full-history anchors, and are no longer misreported as "imputed = FALSE".
# ---------------------------------------------------------------------

table_s6_imputation <- function(imputation_detail, raw_df = NULL) {
  det <- imputation_detail
  # Aggregation convention = all the cells inside the 2011-2024 window (observed / imputed / not imputable)
  # n_regions: the raw data contain country rows as well as region-scale rows (continents/subregions/global), so the number of regions is counted rather than the number of countries
  win <- det %>% filter(year %in% TRAJ_YEARS)
  per <- win %>%
    group_by(indicator) %>%
    summarise(
      n_regions = n_distinct(country),
      window_years = n_distinct(year),
      n_observed_cells = sum(status == "observed"),
      n_imputed_cells = sum(status == "imputed"),
      n_imputed_linear = sum(method == "linear_interpolation"),
      n_imputed_forward = sum(method == "forward_fill"),
      n_imputed_backward = sum(method == "backward_fill"),
      n_missing_cells = sum(status == "missing_no_anchor"),
      n_grid_added_rows = sum(grid_added),
      .groups = "drop"
    )
  disc <- det %>%
    filter(status == "discarded_outside_window") %>%
    count(indicator, name = "n_discarded_outside_window")
  per <- per %>%
    left_join(disc, by = "indicator") %>%
    mutate(n_discarded_outside_window = replace_na(n_discarded_outside_window, 0L))
  if (!is.null(raw_df)) {
    per <- per %>%
      left_join(raw_df %>% group_by(indicator) %>%
                  summarise(period = paste0(min(year[!is.na(value)]),
                                            "-", max(year[!is.na(value)])),
                            .groups = "drop"),
                by = "indicator")
  } else {
    per <- per %>% mutate(period = as.character(NA))
  }
  per %>%
    mutate(
      data_type = "numeric",
      imputed = n_imputed_cells > 0,
      method = case_when(
        n_imputed_cells == 0 ~ "none",
        TRUE ~ paste0("endpoints: forward/backward fill (LOCF); ",
                      "interior: linear interpolation; ",
                      "anchors = all available observations (incl. outside 2011-2024)")
      )
    ) %>%
    select(indicator, period, n_regions, window_years,
           n_observed_cells, n_imputed_cells, n_imputed_linear,
           n_imputed_forward, n_imputed_backward, n_missing_cells,
           n_grid_added_rows, n_discarded_outside_window,
           data_type, imputed, method)
}

# ---------------------------------------------------------------------
# Section 7. section/goal summaries (b6_*; Fig. 2c)
#   aggregate trend slope, R2, achievement rate + bootstrap intervals
#   + Kruskal-Wallis test of the differences between sections
# ---------------------------------------------------------------------

section_goal_summary <- function(res, n_boot = 2000, seed = 1) {
  ts <- load_target_section()
  set.seed(seed)
  # Section level: slope/R2/AR of the target level
  ar_t <- res$ar_summary %>% filter(layer == "target") %>%
    select(layer_id, AR)
  joined <- res$trajectory_class %>%
    filter(layer == "target") %>%
    left_join(ts %>% select(layer_id, section, goal), by = "layer_id") %>%
    left_join(res$baseline_target_goal %>% filter(layer == "target") %>%
                select(layer_id, slope, r2, p_value), by = "layer_id") %>%
    left_join(ar_t, by = "layer_id")
  by_section <- joined %>%
    group_by(section) %>%
    summarise(
      n_targets = n(),
      slope_mean = mean(slope, na.rm = TRUE),
      r2_mean = mean(r2, na.rm = TRUE),
      AR_mean = mean(AR, na.rm = TRUE),
      rho_mean = mean(rho, na.rm = TRUE),
      .groups = "drop"
    )
  # Bootstrap interval of the AR at the section level (the targets are resampled)
  ar_t <- res$ar_summary %>% filter(layer == "target") %>%
    left_join(ts %>% select(layer_id, section), by = "layer_id")
  ar_sec_ci <- ar_t %>%
    group_by(section) %>%
    group_modify(~ {
      x <- .x$AR[!is.na(.x$AR)]
      if (length(x) == 0) return(tibble(AR = NA_real_, ci95_lo = NA_real_,
                                        ci95_hi = NA_real_))
      boots <- replicate(n_boot, mean(sample(x, length(x), replace = TRUE)))
      q <- quantile(boots, c(0.025, 0.975), na.rm = TRUE)
      tibble(AR = mean(x), ci95_lo = unname(q[1]), ci95_hi = unname(q[2]))
    }) %>%
    ungroup()
  by_section <- by_section %>% left_join(ar_sec_ci, by = "section")

  # Goal level
  by_goal <- joined %>%
    group_by(goal) %>%
    summarise(
      n_targets = n(),
      slope_mean = mean(slope, na.rm = TRUE),
      r2_mean = mean(r2, na.rm = TRUE),
      rho_mean = mean(rho, na.rm = TRUE),
      .groups = "drop"
    )
  ar_g <- res$ar_summary %>% filter(layer == "target") %>%
    left_join(ts %>% select(layer_id, goal), by = "layer_id")
  ar_goal_ci <- ar_g %>%
    group_by(goal) %>%
    group_modify(~ {
      x <- .x$AR[!is.na(.x$AR)]
      if (length(x) == 0) return(tibble(AR = NA_real_, ci95_lo = NA_real_,
                                        ci95_hi = NA_real_))
      boots <- replicate(n_boot, mean(sample(x, length(x), replace = TRUE)))
      q <- quantile(boots, c(0.025, 0.975), na.rm = TRUE)
      tibble(AR = mean(x), ci95_lo = unname(q[1]), ci95_hi = unname(q[2]))
    }) %>%
    ungroup()
  by_goal <- by_goal %>% left_join(ar_goal_ci, by = "goal")

  # Kruskal-Wallis test of the differences between sections
  kw_section <- joined %>%
    summarise(
      kw_slope_H = {
        ok <- !is.na(slope)
        if (sum(ok) >= 5 && n_distinct(section[ok]) >= 2)
          safe_kwh(slope[ok], section[ok]) else NA_real_
      },
      kw_slope_p = {
        ok <- !is.na(slope)
        if (sum(ok) >= 5 && n_distinct(section[ok]) >= 2)
          safe_kwp(slope[ok], section[ok]) else NA_real_
      },
      kw_rho_H = {
        ok <- !is.na(rho)
        if (sum(ok) >= 5 && n_distinct(section[ok]) >= 2)
          safe_kwh(rho[ok], section[ok]) else NA_real_
      },
      kw_rho_p = {
        ok <- !is.na(rho)
        if (sum(ok) >= 5 && n_distinct(section[ok]) >= 2)
          safe_kwp(rho[ok], section[ok]) else NA_real_
      }
    )
  list(by_section = by_section, by_goal = by_goal, kw_section = kw_section)
}

# ---------------------------------------------------------------------
# Section 8. trajectory analyses (b7_*; SI B.2.7, annual series 2011-2024)
# ---------------------------------------------------------------------

# 8.1 Share of monotonicity: statistics of the sign reversals of the derivative, per indicator and country
monotonicity_analysis <- function(norm_df) {
  per_ic <- norm_df %>%
    filter(year %in% TRAJ_YEARS) %>%
    select(country, year, indicator, norm) %>%
    pivot_wider(names_from = year, values_from = norm) %>%
    rowwise() %>%
    mutate(
      mono = is_monotonic(c_across(starts_with("20"))),
      n_reversals = {
        d <- diff(c_across(starts_with("20")))
        d <- d[!is.na(d) & d != 0]
        if (length(d) <= 1) 0L else sum(sign(d[-1]) != sign(d[-length(d)]))
      }
    ) %>%
    ungroup()
  by_indicator <- per_ic %>%
    group_by(indicator) %>%
    summarise(mono_share = mean(mono, na.rm = TRUE),
              one_reversal_share = mean(n_reversals == 1, na.rm = TRUE),
              multi_reversal_share = mean(n_reversals > 1, na.rm = TRUE),
              .groups = "drop")
  overall <- tibble(
    mono_share = mean(per_ic$mono, na.rm = TRUE),
    one_reversal_share = mean(per_ic$n_reversals == 1, na.rm = TRUE),
    multi_reversal_share = mean(per_ic$n_reversals > 1, na.rm = TRUE),
    n_ind = n_distinct(per_ic$indicator)
  )
  list(overall = overall, by_indicator = by_indicator)
}

# 8.2 Acceleration vs linear: the quadratic slope at 2024 minus the linear slope
acceleration_analysis <- function(norm_df, quantiles = ACCEL_QUANTILES) {
  acc <- norm_df %>%
    filter(year %in% TRAJ_YEARS) %>%
    group_by(indicator, year) %>%
    summarise(v = mean(norm, na.rm = TRUE), .groups = "drop") %>%
    group_by(indicator) %>%
    group_modify(~ {
      d <- .x %>% filter(is.finite(v))
      if (nrow(d) < 5) return(tibble(slope_linear = NA_real_,
                                     slope_quad_2024 = NA_real_,
                                     acceleration = NA_real_))
      m1 <- lm(v ~ year, data = d)
      m2 <- lm(v ~ poly(year, 2, raw = TRUE), data = d)
      s1 <- coef(m1)[2]
      # the instantaneous slope of the quadratic model at 2024
      s2 <- coef(m2)[2] + 2 * coef(m2)[3] * 2024
      tibble(slope_linear = unname(s1),
             slope_quad_2024 = unname(s2),
             acceleration = unname(s2 - s1))
    }) %>%
    ungroup()
  thr <- quantile(acc$acceleration, quantiles, na.rm = TRUE)
  list(
    by_indicator = acc %>%
      mutate(flag_q75 = abs(acceleration) > thr[1],
             flag_q90 = abs(acceleration) > thr[2]),
    thresholds = tibble(quantile = names(thr), threshold = as.numeric(thr))
  )
}

# 8.3 Endpoint vs trajectory fit: the endpoint difference vs the linear extrapolation of 2011-2023 to 2024
endpoint_vs_trajectory <- function(mis_eis_annual) {
  wide <- mis_eis_annual %>%
    select(country, year, MIS) %>%
    filter(year %in% TRAJ_YEARS) %>%
    pivot_wider(names_from = year, values_from = MIS)
  yr_cols <- intersect(as.character(TRAJ_YEARS), names(wide))
  fit_cols <- intersect(as.character(2011:2023), names(wide))
  res <- map_dfr(seq_len(nrow(wide)), function(i) {
    vals <- as.numeric(wide[i, fit_cols])
    ok <- !is.na(vals)
    traj_2024 <- if (sum(ok) >= 3) {
      fit <- lm(vals[ok] ~ as.numeric(fit_cols)[ok])
      unname(coef(fit)[1] + coef(fit)[2] * 2024)
    } else NA_real_
    obs11 <- as.numeric(wide[i, "2011"])
    obs24 <- as.numeric(wide[i, "2024"])
    tibble(endpoint_diff = obs24 - obs11,
           traj_diff = traj_2024 - obs11,
           endpoint_bias = (obs24 - obs11) - (traj_2024 - obs11),
           mae = abs(obs24 - traj_2024))
  })
  tibble(endpoint_bias_mean = mean(res$endpoint_bias, na.rm = TRUE),
         endpoint_bias_sd = sd(res$endpoint_bias, na.rm = TRUE),
         mae = mean(res$mae, na.rm = TRUE),
         n = sum(!is.na(res$endpoint_bias)))
}

# 8.4 COVID-2020 perturbation: recompute the MIS(2011,2024) headline without 2020
covid_perturbation <- function(mis_eis_annual) {
  full <- mis_eis_annual %>%
    filter(year %in% c(2011, 2024)) %>%
    group_by(year) %>%
    summarise(m = mean(MIS, na.rm = TRUE), .groups = "drop")
  headline_full <- full$m[full$year == 2024] - full$m[full$year == 2011]
  # Removing 2020: the endpoints stay the same, the comparison uses the linear trend of 2011-2024 (without 2020) extrapolated
  traj <- mis_eis_annual %>%
    select(country, year, MIS) %>%
    filter(year %in% TRAJ_YEARS) %>%
    group_by(country) %>%
    group_modify(~ {
      d_full <- .x %>% filter(!is.na(MIS))
      d_ex <- d_full %>% filter(year != 2020)
      fit_full <- if (nrow(d_full) >= 3) lm(MIS ~ year, data = d_full) else NULL
      fit_ex <- if (nrow(d_ex) >= 3) lm(MIS ~ year, data = d_ex) else NULL
      tibble(
        trend_full = if (!is.null(fit_full)) coef(fit_full)[2] else NA_real_,
        trend_ex2020 = if (!is.null(fit_ex)) coef(fit_ex)[2] else NA_real_
      )
    }) %>%
    ungroup()
  tibble(
    headline_change_full = headline_full,
    trend_mean_full = mean(traj$trend_full, na.rm = TRUE),
    trend_mean_ex2020 = mean(traj$trend_ex2020, na.rm = TRUE),
    delta_trend = mean(traj$trend_ex2020 - traj$trend_full, na.rm = TRUE),
    n_countries = nrow(traj)
  )
}

# 8.5 Regional annual paths (inset of Fig 3)
regional_annual_paths <- function(res) {
  gm <- group_membership(res$country_meta)
  res$mis_eis_annual %>%
    inner_join(gm, by = "country") %>%
    group_by(group, year) %>%
    summarise(MIS = mean(MIS, na.rm = TRUE),
              EIS = mean(EIS, na.rm = TRUE),
              n = n(), .groups = "drop")
}

# ---------------------------------------------------------------------
# Section 9. bootstrap intervals of the achievement rate (c1_*; SI A.5)
# ---------------------------------------------------------------------

achievement_ci_report <- function(res, n_boot = 2000) {
  ts <- load_target_section()
  per_target <- res$ar_ci %>% filter(layer == "target")
  per_goal <- res$ar_ci %>% filter(layer == "goal")
  overall <- res$ar %>%
    group_by(layer) %>%
    group_modify(~ {
      x <- .x$I[!is.na(.x$I)]
      n <- length(x)
      boots <- replicate(n_boot, mean(sample(x, n, replace = TRUE)))
      q50 <- quantile(boots, c(0.25, 0.75), na.rm = TRUE)
      q95 <- quantile(boots, c(0.025, 0.975), na.rm = TRUE)
      tibble(AR = mean(x),
             ci50_lo = unname(q50[1]), ci50_hi = unname(q50[2]),
             ci95_lo = unname(q95[1]), ci95_hi = unname(q95[2]), n = n)
    }) %>%
    ungroup()
  by_section <- res$ar_summary %>%
    filter(layer == "target") %>%
    left_join(ts %>% select(layer_id, section), by = "layer_id") %>%
    group_by(section) %>%
    summarise(AR = mean(AR, na.rm = TRUE), n_targets = n(), .groups = "drop")
  list(overall = overall, per_target = per_target, per_goal = per_goal,
       by_section = by_section)
}

# ---------------------------------------------------------------------
# Section 10. uncertainty, Algorithms 1-5 (SI B.1/B.2.1)
# ---------------------------------------------------------------------

# Algorithm 1 - exhaustive enumeration of the indicator subsets: for the targets/goals with >=2 indicators, one composite value per k-subset (the pseudocode score_c of the SI)
#   score_c = f(X_c): first the units are combined with Formula 2 (the mean of the available indicators inside the subset), then the mean across the units
#   (identical to the global convention of the paper, mean_na). norm_df passes the country level plus the regional level stacked together; the regional indicators
#   (27-29, regional rows only) enter the combination through their regional values. For k = n there is only one combination -> a single horizontal line, not a box plot.
algo1_indicator_subset <- function(norm_df, i2t, panel_year = 2014,
                                   max_combos_per_k = 200) {
  set.seed(42)
  # The targets T1-T23 and the goals A-D are included together (those with >=2 indicators; goal C has only 1 indicator and is excluded by n>=2)
  map_i2t <- i2t %>% distinct(indicator, layer_id)
  multi <- map_i2t %>% count(layer_id) %>% filter(n >= 2)
  map_dfr(multi$layer_id, function(tgt) {
    inds <- map_i2t$indicator[map_i2t$layer_id == tgt]
    map_dfr(seq_along(inds), function(k) {
      combos <- combn(inds, k, simplify = FALSE)
      if (length(combos) > max_combos_per_k) {
        combos <- combos[sort(sample.int(length(combos), max_combos_per_k))]
      }
      map_dfr(combos, function(cmb) {
        v <- norm_df %>%
          filter(year == panel_year, indicator %in% cmb) %>%
          group_by(country) %>%
          summarise(x = mean(norm, na.rm = TRUE), .groups = "drop") %>%
          pull(x)
        tibble(layer_id = tgt, k = k, n_combos = length(combos),
               n_units = sum(is.finite(v)),
               value = mean(v, na.rm = TRUE))
      })
    })
  })
}

# Algorithm 2 — bootstrap-based imputation uncertainty (SI Algorithm 2)
#   Stage 1  LOO residuals: after removing each observation, the value is predicted again with the imputation rule, eps = x - x_hat;
#            the residuals are stratified by the gap context g - g = the span between the two anchors of the target point (right bound - left bound - 1),
#            i.e. the length of the run of consecutive missing values in which the target sits; LOO uses the same convention as the missing cells
#   Stage 2  For every imputed cell inside a window: x* = f_hat + resampled residual (B times), the CI [L,U] = the quantiles,
#            width = U − L; stability = P(sign(x*) = sign(f̂))
#            the residual pools are stratified by (indicator, g); when the pool < min_pool the fallback goes stepwise to (indicator, all g) and then to the global normalised pool
#   Stage 3  Calibration coverage = P(x in [L,U]) (observed points, the intervals are centred on the LOO prediction);
#            the mean width and the entropy are reported by group g (the differential entropy of the bootstrap distribution under a normal approximation 0.5*log(2*pi*e*sigma^2))
#   Sensitivity  One factor at a time around the baseline (B=1000, alpha=0.05, Linear, LOO):
#            B in {200,500,2000} (the first B samples of the B_max resamples are taken, which has the same distribution as an independent rerun);
#            alpha in {0.10,0.05,0.01}, imputer in {Linear,FFill,BFill,Spline},
#            rolling-window CV (the anchors are restricted to the target year +/- roll_window years, which affects only the residual estimation)
#   Input raw_df: the long table of the raw observations (indicator, country, year, value), value=NA means missing;
#   returns list(per_indicator, gap_length, sensitivity) - put() writes the 3 CSVs
algo2_imputation <- function(raw_df, impute_years = TRAJ_YEARS,
                             B = 1000, alpha = 0.05, B_max = 2000,
                             min_pool = 30, roll_window = 5, seed = 42) {
  set.seed(seed)
  B_grid <- c(200, 500, 1000, 2000); B_grid <- B_grid[B_grid <= B_max]
  a_grid <- c(0.10, 0.05, 0.01)
  methods <- c("Linear", "FFill", "BFill", "Spline")

  raw_df <- raw_df %>% distinct(indicator, country, year, .keep_all = TRUE)

  # Deterministic imputation f: given the observations (yo, vo) (yo is sorted), predict the target year yt
  #   Linear = the main method of A.2 (linear interpolation inside, forward/backward fill at the endpoints); the others are alternative imputers
  f_predict <- function(yo, vo, yt, method) {
    n <- length(yo)
    if (n == 0) return(rep(NA_real_, length(yt)))
    vapply(yt, function(y0) {
      if (method == "FFill") {
        left <- which(yo <= y0)
        if (length(left)) vo[max(left)] else vo[1]
      } else if (method == "BFill") {
        right <- which(yo >= y0)
        if (length(right)) vo[min(right)] else vo[n]
      } else if (method == "Spline") {
        if (n == 1 || y0 <= yo[1]) vo[1]
        else if (y0 >= yo[n]) vo[n]
        else spline(yo, vo, xout = y0)$y
      } else {  # Linear
        if (n == 1) vo[1]
        else {
          left <- which(yo < y0); right <- which(yo > y0)
          if (length(left) && length(right)) {
            i1 <- max(left); i2 <- min(right)
            vo[i1] + (vo[i2] - vo[i1]) * (y0 - yo[i1]) / (yo[i2] - yo[i1])
          } else if (length(left)) vo[max(left)]
          else vo[min(right)]
        }
      }
    }, numeric(1))
  }

  # Gap length g = (right bound - left bound - 1): the length of the run of consecutive missing values in which the target sits
  #   yo = the reference observation years; y0 = the target year; tl_min/tl_max = the first and last year of the timeline of the series (the bound on the side without observations)
  gap_len <- function(yo, y0, tl_min, tl_max) {
    lo <- max(c(yo[yo < y0], tl_min - 1))
    hi <- min(c(yo[yo > y0], tl_max + 1))
    hi - lo - 1
  }

  # ---- Series list: one row per (indicator, country) ----
  ser_all <- raw_df %>%
    group_by(indicator, country) %>%
    summarise(
      yo = list(sort(year[!is.na(value)])),
      vo = list(value[!is.na(value)][order(year[!is.na(value)])]),
      tl_min = min(c(year, impute_years)), tl_max = max(c(year, impute_years)),
      miss_win = list(setdiff(impute_years, year[!is.na(value)])),
      .groups = "drop")
  ser  <- ser_all %>% filter(lengths(yo) >= 2)   # LOO needs >=2 observations
  ser1 <- ser_all %>% filter(lengths(yo) == 1)   # a single observation: Stage 2 only (constant fill)

  # Scale of each indicator (IQR): makes the widths comparable across indicators
  scale_i <- raw_df %>% filter(!is.na(value)) %>%
    group_by(indicator) %>%
    summarise(scale = {
      q <- quantile(value, c(0.25, 0.75), na.rm = TRUE, names = FALSE)
      s <- q[2] - q[1]
      if (!is.finite(s) || s <= 0) max(abs(value), 1) else s
    }, .groups = "drop")

  # ---- Stage 1: LOO residuals (each observation x method; the rolling anchor set applies to Linear only) ----
  loo_all <- map_dfr(seq_len(nrow(ser)), function(k) {
    yo <- ser$yo[[k]]; vo <- ser$vo[[k]]; n <- length(yo)
    g <- vapply(seq_len(n), function(j)
      gap_len(yo[-j], yo[j], ser$tl_min[k], ser$tl_max[k]), numeric(1))
    map_dfr(methods, function(m) {
      eps_full <- vapply(seq_len(n), function(j)
        vo[j] - f_predict(yo[-j], vo[-j], yo[j], m), numeric(1))
      eps_roll <- if (m == "Linear") {
        vapply(seq_len(n), function(j) {
          keep <- abs(yo[-j] - yo[j]) <= roll_window
          if (!any(keep)) NA_real_
          else vo[j] - f_predict(yo[-j][keep], vo[-j][keep], yo[j], m)
        }, numeric(1))
      } else rep(NA_real_, n)
      tibble(indicator = ser$indicator[k], country = ser$country[k],
             year = yo, g = g, method = m, eps = eps_full, eps_roll = eps_roll)
    })
  })

  # ---- Residual pools (a single split; the key is method|indicator|g) ----
  fin <- is.finite(loo_all$eps);  finr <- is.finite(loo_all$eps_roll)
  kg <- paste(loo_all$method, loo_all$indicator, loo_all$g, sep = "|")
  ki <- paste(loo_all$method, loo_all$indicator, sep = "|")
  pool_g   <- split(loo_all$eps[fin],       kg[fin])
  pool_i   <- split(loo_all$eps[fin],       ki[fin])
  pool_g_r <- split(loo_all$eps_roll[finr], kg[finr])
  pool_i_r <- split(loo_all$eps_roll[finr], ki[finr])
  gn <- loo_all %>% filter(fin) %>%
    left_join(scale_i, by = "indicator") %>% mutate(e = eps / scale)
  glob   <- split(gn$e, gn$method)
  gnr <- loo_all %>% filter(finr) %>%
    left_join(scale_i, by = "indicator") %>% mutate(e = eps_roll / scale)
  glob_r <- split(gnr$e, gnr$method)

  pool_lookup <- function(m, i_id, g0, use_roll) {
    k <- paste(m, i_id, g0, sep = "|")
    p <- if (use_roll) pool_g_r[[k]] else pool_g[[k]]
    if (!is.null(p) && length(p) >= min_pool) return(p)
    k2 <- paste(m, i_id, sep = "|")
    p <- if (use_roll) pool_i_r[[k2]] else pool_i[[k2]]
    if (!is.null(p) && length(p) >= min_pool) return(p)
    s <- scale_i$scale[scale_i$indicator == i_id]
    if (!length(s)) s <- 1
    gn0 <- if (use_roll) glob_r[[m]] else glob[[m]]
    gn0 * s[1]
  }

  # Resampling cache: one B_max resampling per (m, i, g, roll); the B grid takes a prefix of it (the same distribution)
  rs_cache <- new.env(parent = emptyenv())
  get_draws <- function(m, i_id, g0, use_roll, bb) {
    ck <- paste(m, i_id, g0, as.integer(use_roll), sep = "|")
    if (!exists(ck, rs_cache, inherits = FALSE))
      assign(ck, sample(pool_lookup(m, i_id, g0, use_roll), B_max, replace = TRUE),
             rs_cache)
    get(ck, rs_cache, inherits = FALSE)[seq_len(bb)]
  }
  # Cache of the sorted draws (for the stability, vectorised with findInterval in O(log B))
  srt_cache <- new.env(parent = emptyenv())
  get_sorted <- function(m, i_id, g0, use_roll, bb) {
    ck <- paste(m, i_id, g0, as.integer(use_roll), bb, sep = "|")
    if (!exists(ck, srt_cache, inherits = FALSE))
      assign(ck, sort(get_draws(m, i_id, g0, use_roll, bb)), srt_cache)
    get(ck, srt_cache, inherits = FALSE)
  }
  # Cache of the quantiles/variances: the width and the entropy are shared by all the cells within the same pool (m,i,g,roll,B,alpha)
  q_cache <- new.env(parent = emptyenv())
  get_q <- function(m, i_id, g0, use_roll, bb, aa) {
    ck <- paste(m, i_id, g0, as.integer(use_roll), bb, aa, sep = "|")
    if (!exists(ck, q_cache, inherits = FALSE)) {
      dr <- get_draws(m, i_id, g0, use_roll, bb)
      assign(ck, c(quantile(dr, c(aa / 2, 1 - aa / 2), names = FALSE),
                   var = var(dr)), q_cache)
    }
    get(ck, q_cache, inherits = FALSE)
  }

  # ---- Cell list (input of Stage 2, including the single-observation series), with the point estimates f_hat of the 4 methods precomputed ----
  ent_all <- bind_rows(lapply(list(ser, ser1), function(sx) {
    map_dfr(seq_len(nrow(sx)), function(k) {
      miss <- sx$miss_win[[k]]
      if (!length(miss)) return(NULL)
      g0 <- vapply(miss, function(y0)
        gap_len(sx$yo[[k]], y0, sx$tl_min[k], sx$tl_max[k]), numeric(1))
      fmat <- vapply(methods, function(m)
        f_predict(sx$yo[[k]], sx$vo[[k]], miss, m), numeric(length(miss)))
      # fmat: a length(miss) x 4 matrix; with a single missing year vapply degenerates into a named vector of length 4
      if (is.null(dim(fmat)))
        fmat <- matrix(fmat, nrow = 1, dimnames = list(NULL, methods))
      tibble(indicator = sx$indicator[k], country = sx$country[k],
             year = miss, g = g0, fhat = asplit(fmat, 1))
    })
  }))
  # fhat list column (each element a named numeric: Linear/FFill/BFill/Spline) -> take the column of the single method
  fhat_of <- function(m) vapply(ent_all$fhat, `[[`, numeric(1), m)

  # ---- Stage 2: bootstrap for every imputed cell (given the method/anchor set/B/alpha) ----
  # stability = P(sign(f_hat+eps*) = sign(f_hat)): when f_hat>0, eps* > -f_hat; when f_hat<0, eps* < -f_hat
  #   vectorised with the sorted draws + findInterval (O(log B) per cell)
  stab_of <- function(fhat_v, m, i_id, g0, use_roll, bb) {
    sdr <- get_sorted(m, i_id, g0, use_roll, bb)
    Bb <- length(sdr)
    ifelse(fhat_v > 0,
           (Bb - findInterval(-fhat_v, sdr)) / Bb,
           ifelse(fhat_v < 0,
                  findInterval(-fhat_v, sdr, left.open = TRUE) / Bb, 0))
  }
  entry_one <- function(m, use_roll, bb, aa) {
    fhat <- fhat_of(m)
    n <- nrow(ent_all)
    width <- stab <- entv <- numeric(n)
    key <- paste(m, ent_all$indicator, ent_all$g, as.integer(use_roll), sep = "|")
    for (k in unique(key)) {
      jj <- which(key == k)
      i_id <- ent_all$indicator[jj[1]]; g0 <- ent_all$g[jj[1]]
      qv <- get_q(m, i_id, g0, use_roll, bb, aa)
      width[jj] <- qv[2] - qv[1]
      stab[jj]  <- stab_of(fhat[jj], m, i_id, g0, use_roll, bb)
      entv[jj]  <- 0.5 * log(2 * pi * exp(1) * max(qv["var"], .Machine$double.eps))
    }
    ent_all %>% mutate(fhat = fhat, width = width, stab = stab, ent = entv,
                       method = m, use_roll = use_roll, B = bb, alpha = aa)
  }

  # Coverage (the observed points, per indicator): eps in [q_alpha/2, q_1-alpha/2] (the interval = f_hat_LOO + the pool quantiles)
  cov_one <- function(m, use_roll = FALSE, bb = B, aa = alpha) {
    e_col <- if (use_roll) "eps_roll" else "eps"
    d0 <- loo_all %>% filter(method == m, is.finite(.data[[e_col]]))
    e0 <- d0[[e_col]]
    key <- paste(m, d0$indicator, d0$g, as.integer(use_roll), sep = "|")
    covered <- numeric(nrow(d0))
    for (k in unique(key)) {
      jj <- which(key == k)
      qv <- get_q(m, d0$indicator[jj[1]], d0$g[jj[1]], use_roll, bb, aa)
      covered[jj] <- as.numeric(e0[jj] >= qv[1] & e0[jj] <= qv[2])
    }
    tibble(indicator = d0$indicator, covered = covered) %>%
      group_by(indicator) %>%
      summarise(coverage = mean(covered), n_obs = n(), .groups = "drop") %>%
      mutate(method = m, use_roll = use_roll)
  }

  message("[algo2] Stage 1: LOO residuals (", nrow(ser), " series) ...")
  message("[algo2] Stage 2: baseline Linear/LOO on B x alpha grid ...")
  base <- bind_rows(
    lapply(B_grid, function(bb)
      bind_rows(lapply(a_grid, function(aa) entry_one("Linear", FALSE, bb, aa))))) %>%
    left_join(scale_i, by = "indicator") %>%
    mutate(width_rel = width / scale)
  base_row <- base %>% filter(B == .env$B, alpha == .env$alpha)

  message("[algo2] sensitivity: alternative imputers ...")
  alt <- bind_rows(lapply(setdiff(methods, "Linear"),
                          function(m) entry_one(m, FALSE, B, alpha))) %>%
    left_join(scale_i, by = "indicator") %>%
    mutate(width_rel = width / scale)
  message("[algo2] sensitivity: rolling-window CV ...")
  cv_roll <- entry_one("Linear", TRUE, B, alpha) %>%
    left_join(scale_i, by = "indicator") %>%
    mutate(width_rel = width / scale)

  message("[algo2] Stage 3: calibration coverage ...")
  cov_all <- bind_rows(
    cov_one("Linear"),
    map_dfr(setdiff(methods, "Linear"), cov_one),
    cov_one("Linear", use_roll = TRUE))
  cov_B <- map_dfr(setdiff(B_grid, B), function(bb) {
    cov_one("Linear", FALSE, bb, alpha) %>%
      summarise(coverage = mean(coverage)) %>% mutate(B = bb)
  })
  cov_A <- map_dfr(setdiff(a_grid, alpha), function(aa) {
    cov_one("Linear", FALSE, B, aa) %>%
      summarise(coverage = mean(coverage)) %>% mutate(alpha = aa)
  })
  cov_mean <- function(m, use_roll = FALSE)
    mean((cov_all %>% filter(method == m, .data$use_roll == .env$use_roll))$coverage)

  # ---- Output 1: per indicator (the baseline Linear/LOO, B, alpha) ----
  per_indicator <- base_row %>%
    group_by(indicator) %>%
    summarise(n_imputed = n(),
              mean_gap = mean(g),
              mean_ci_width = mean(width),
              mean_ci_width_rel = mean(width_rel),
              median_ci_width_rel = median(width_rel),
              mean_stability = mean(stab),
              entropy = mean(ent),
              .groups = "drop") %>%
    left_join(cov_all %>% filter(method == "Linear", !use_roll) %>%
                select(indicator, coverage, n_obs), by = "indicator") %>%
    left_join(scale_i, by = "indicator")

  # ---- Output 2: grouped by gap length g (the baseline; across indicators the relative width is used) ----
  gap_length <- base_row %>%
    group_by(g) %>%
    summarise(n_entries = n(),
              mean_ci_width_rel = mean(width_rel),
              median_ci_width_rel = median(width_rel),
              mean_stability = mean(stab),
              entropy = mean(ent),
              .groups = "drop") %>%
    arrange(g)

  # ---- Output 3: sensitivity (one factor at a time around the baseline) ----
  agg <- function(d, factor, level, cov) {
    d %>%
      summarise(mean_ci_width_rel = mean(width_rel),
                mean_stability = mean(stab),
                entropy = mean(ent),
                n_entries = n(), .groups = "drop") %>%
      mutate(factor = factor, level = as.character(level),
             coverage = cov, .before = 1)
  }
  bl <- sprintf("B=%d, alpha=%.2f, Linear, LOO", B, alpha)
  sensitivity <- bind_rows(
    agg(base_row, "baseline", bl, cov_mean("Linear")),
    bind_rows(lapply(setdiff(B_grid, B), function(bb)
      agg(base %>% filter(method == "Linear", B == bb, alpha == .env$alpha),
          "B", bb, cov_B$coverage[cov_B$B == bb]))),
    bind_rows(lapply(setdiff(a_grid, alpha), function(aa)
      agg(base %>% filter(method == "Linear", B == .env$B, alpha == aa),
          "alpha", aa, cov_A$coverage[cov_A$alpha == aa]))),
    bind_rows(lapply(setdiff(methods, "Linear"), function(m)
      agg(alt %>% filter(method == m), "imputer", m, cov_mean(m)))),
    agg(cv_roll, "cv", sprintf("rolling-window (+/- %d yr)", roll_window),
        cov_mean("Linear", use_roll = TRUE))
  )

  message("[algo2] done: ", nrow(per_indicator), " indicators, ",
          nrow(gap_length), " gap groups, ", nrow(sensitivity), " sensitivity rows")
  list(per_indicator = per_indicator,
       gap_length = gap_length,
       sensitivity = sensitivity)
}

# Algorithm 3 - time-series length: block-bootstrap resampling of the window + one-step prediction dispersion of three models
#   The pseudocode of the SI is followed line by line:
#     L = {L_min, ..., floor(T/2)}; M = {ETS, auto-ARIMA, Random Walk}; B draws with replacement of the blocks
#       // the upper bound floor(T/2): a window close to the full length leaves only a handful of blocks, so sigma collapses spuriously (indicator 27, L=30: 3 blocks)
#     sigma_{i,L} = std_{b,m}(y_hat)  // the dimensions b and m are pooled; the per-model sigma is reported separately
#     CI_{i,L} = 1.96·σ_{i,L}
#     L* = Kneedle(sigma(L)) (the maximum deviation from the chord to the endpoints after normalisation; the curve is smoothed with a 3-point moving average first)
#       // the knee is only taken for convex (diminishing returns) curves; flat/linear/concave curves have no knee -> fall back to L^der
#     L^der = the first point with |Delta sigma| <= eps; eps = eps_der*sigma(L_min) (the uncertainty enters the plateau)
#   Two-step speed-up (numerically equivalent to fitting the blocks one by one):
#     1) Cache: for every (i, L) the three model predictions are computed in advance for all T-L+1 possible contiguous blocks, and the B draws are looked up
#     2) Only the history up to year_max is used (the raw data of indicator 2 contain scenario values for 2025-2050, which do not enter the window analysis)
#   The regional indicators 27-29 are merged in through the norm_regional rows (the Global summary rows are removed; the annual mean convention is the same as at the country level)
#   returns list(curves, optimal, cross_model) -> uncertainty__a3_window_length__*.csv
algo3_window_length <- function(norm_df, B = 1000, L_min = 5, year_max = 2024,
                                eps_der = 0.05, seed = 42) {
  set.seed(seed)
  has_forecast <- requireNamespace("forecast", quietly = TRUE)
  MODS <- c("ETS", "ARIMA", "RW")

  # Global mean series: the annual mean after removing the Global summary rows, keeping only the history up to year_max
  series_all <- norm_df %>%
    filter(as.character(country) != "Global", year <= year_max) %>%
    group_by(indicator, year) %>%
    summarise(v = mean(norm, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(v)) %>%
    arrange(indicator, year) %>%
    group_split(indicator)

  # One-step prediction of a single block and a single model: RW = the last value; without the forecast package both ETS and ARIMA fall back to the linear trend
  fit1 <- function(train, model) {
    if (model == "RW") return(train[length(train)])
    if (!has_forecast) {
      tt <- seq_along(train)
      cf <- coef(lm(train ~ tt))
      return(unname(cf[1] + cf[2] * (length(train) + 1)))
    }
    tryCatch({
      if (model == "ETS") {
        as.numeric(forecast::forecast(forecast::ets(train), h = 1)$mean)
      } else {
        as.numeric(forecast::forecast(
          forecast::auto.arima(train, stepwise = TRUE, approximation = TRUE),
          h = 1)$mean)
      }
    }, error = function(e) mean(train))
  }

  # Simplified implementation of Kneedle (Satopaa et al. 2011): the maximum deviation from the chord to the endpoints after normalisation
  #   Well-posedness guard: a flat curve (range/median < flat_tol) or a non-convex one (linear/concave) has no knee of diminishing returns
  #   -> NA is returned and the caller falls back to L^der
  knee_L <- function(L, sig, flat_tol = 0.2, knee_tol = 0.05) {
    n <- length(L)
    if (n < 3) return(NA_real_)
    rng <- max(sig) - min(sig)
    if (!is.finite(rng) || rng < 1e-12) return(NA_real_)
    if (rng / median(sig) < flat_tol) return(NA_real_)
    x <- (L - L[1]) / (L[n] - L[1])
    y <- (sig - min(sig)) / rng
    if (n >= 5) {  # 3-point moving average to reduce the noise
      y2 <- y
      for (i in 2:(n - 1)) y2[i] <- mean(y[(i - 1):(i + 1)])
      y <- y2
    }
    d <- y - (y[1] + (y[n] - y[1]) * x)
    if (mean(d) >= 0 || max(abs(d)) < knee_tol) return(NA_real_)
    L[which.min(d)]
  }
  # Derivative-based alternative: the first point with |sigma(L)-sigma(L-1)| <= eps
  der_L <- function(L, sig, eps) {
    if (length(L) < 2) return(NA_real_)
    dd <- abs(diff(sig))
    idx <- which(dd <= eps)
    if (length(idx) == 0) return(L[length(L)])
    L[idx[1] + 1]
  }

  blocks <- map(series_all, function(srs) {
    id <- srs$indicator[1]
    y <- srs$v; n <- length(y)
    L_grid <- seq(L_min, max(L_min, floor(n / 2)))
    if (length(L_grid) < 1) return(NULL)
    sigma <- sigma_m <- matrix(NA_real_, length(L_grid), length(MODS),
                               dimnames = list(NULL, MODS))
    n_blk <- integer(length(L_grid))
    for (li in seq_along(L_grid)) {
      L <- L_grid[li]
      starts <- seq_len(n - L + 1)
      n_blk[li] <- length(starts)
      # cache the three model predictions for all the possible blocks
      Fmat <- vapply(MODS, function(m)
        vapply(starts, function(st) fit1(y[st:(st + L - 1)], m), numeric(1)),
        numeric(length(starts)))
      # B draws with replacement of the blocks (pseudocode lines 3-7)
      draws <- Fmat[sample(starts, B, replace = TRUE), , drop = FALSE]
      sigma_m[li, ] <- apply(draws, 2, sd)      # the per-model std_b
      sigma[li, ]   <- sd(as.vector(draws))     # the pooled std_{b,m} (line 8)
    }
    sigma_pooled <- sigma[, 1]
    eps <- eps_der * sigma_pooled[1]
    L_knee <- knee_L(L_grid, sigma_pooled)
    L_der <- der_L(L_grid, sigma_pooled, eps)
    L_opt <- if (is.na(L_knee)) L_der else L_knee
    L_opt_m <- vapply(MODS, function(m) {
      lm_ <- knee_L(L_grid, sigma_m[, m])
      if (is.na(lm_)) der_L(L_grid, sigma_m[, m], eps_der * sigma_m[1, m]) else lm_
    }, numeric(1))
    opt_at <- function(Lv, sig) if (is.na(Lv)) NA_real_ else sig[L_grid == Lv][1]
    list(
      curves = tibble(indicator = id, L = L_grid, n_blocks = n_blk,
                      sigma = sigma_pooled, ci95_halfwidth = 1.96 * sigma_pooled,
                      sigma_ets = sigma_m[, "ETS"], sigma_arima = sigma_m[, "ARIMA"],
                      sigma_rw = sigma_m[, "RW"]),
      optimal = tibble(indicator = id, n_years = n, L_opt = L_opt, L_der = L_der,
                       knee_found = !is.na(L_knee),
                       sigma_at_L_opt = opt_at(L_opt, sigma_pooled),
                       L_opt_ets = unname(L_opt_m["ETS"]),
                       L_opt_arima = unname(L_opt_m["ARIMA"]),
                       L_opt_rw = unname(L_opt_m["RW"]),
                       model_spread = max(L_opt_m) - min(L_opt_m)),
      cross_model = tibble(indicator = id, model = MODS,
                           L_opt_kneedle = unname(L_opt_m),
                           L_opt_der = vapply(MODS, function(m)
                             der_L(L_grid, sigma_m[, m], eps_der * sigma_m[1, m]),
                             numeric(1)),
                           sigma_at_L_opt = vapply(MODS, function(m)
                             opt_at(L_opt_m[m], sigma_m[, m]), numeric(1)))
    )
  })
  blocks <- compact(blocks)
  curves <- map_dfr(blocks, "curves")
  optimal <- map_dfr(blocks, "optimal") %>% arrange(indicator)
  cross_model <- map_dfr(blocks, "cross_model") %>% arrange(indicator, model)

  message("[algo3] done: ", nrow(optimal), " indicators, L_opt mean ",
          sprintf("%.2f", mean(optimal$L_opt, na.rm = TRUE)),
          " ± ", sprintf("%.2f", sd(optimal$L_opt, na.rm = TRUE)),
          " (median ", median(optimal$L_opt, na.rm = TRUE),
          ", range ", min(optimal$L_opt, na.rm = TRUE), "-",
          max(optimal$L_opt, na.rm = TRUE),
          "); convex knee ", sum(optimal$knee_found), "/", nrow(optimal),
          "; model spread mean ",
          sprintf("%.2f", mean(optimal$model_spread, na.rm = TRUE)))
  list(curves = curves, optimal = optimal, cross_model = cross_model)
}

# Algorithm 4 - country coverage: the scaling law sigma_i(n) per indicator + the stability of the full-size bootstrap
#   + the country bootstrap of the headline (it feeds directly into the conclusions; the keep_index 130-country convention)
#   Phase 1: sigma_i(n) = a_i/sqrt(n) + b_i (the fit is only a summary of the 1/sqrt(n) decay)
#   Phase 2: the share of years stable within +/-5% of the indicator mean in the full-panel bootstrap (the n=30 reference quantifies the inadequacy of a small panel)
#   Phase 3: the country bootstrap distribution of the headline at the goal level
algo4_country_coverage <- function(norm_df, i2t, coverage_df, layer_df,
                                   eps = 0.05,
                                   n_grid = c(5, 10, 20, 40, 80, 120, 160),
                                   B = 50, B_stab = 200, B_head = 500) {
  set.seed(42)
  keep_cty <- coverage_df %>% filter(keep_index) %>% pull(country) %>% unique()
  inds <- sort(unique(norm_df$indicator))

  # -- Phase 1: the subsampling scaling law sigma_i(n) per indicator --
  sigma_by_n <- map_dfr(inds, function(i) {
    d <- norm_df %>% filter(indicator == .env$i, year %in% TRAJ_YEARS)
    cs_all <- unique(d$country)
    years <- sort(unique(d$year))
    map_dfr(n_grid, function(n) {
      if (length(cs_all) < n) return(NULL)
      ses <- map_dbl(years, function(yr) {
        vals <- replicate(B, {
          cs <- sample(cs_all, n)
          mean(d$norm[d$year == yr & d$country %in% cs], na.rm = TRUE)
        })
        sd(vals, na.rm = TRUE)
      })
      tibble(indicator = i, n = n, sigma = mean(ses, na.rm = TRUE))
    })
  })
  sigma_fits <- sigma_by_n %>% group_by(indicator) %>% group_modify(~ {
    if (nrow(.x) < 3) return(tibble(a = NA_real_, b = NA_real_, r2 = NA_real_))
    fit <- tryCatch(nls(sigma ~ a / sqrt(n) + b, data = .x,
                        start = list(a = 10, b = 0)),
                    error = function(e) NULL)
    if (is.null(fit)) return(tibble(a = NA_real_, b = NA_real_, r2 = NA_real_))
    pred <- predict(fit)
    r2 <- 1 - sum((.x$sigma - pred)^2) / sum((.x$sigma - mean(.x$sigma))^2)
    tibble(a = coef(fit)[1], b = coef(fit)[2], r2 = r2)
  }) %>% ungroup()

  # -- Phase 2: the stability of the full-size bootstrap (n_fixed=30 as the small-panel reference) --
  stab_share <- function(n_fixed = NULL) {
    map_dbl(inds, function(i) {
      d <- norm_df %>% filter(indicator == .env$i, year %in% TRAJ_YEARS)
      years <- sort(unique(d$year))
      mean(map_lgl(years, function(yr) {
        v <- d$norm[d$year == yr]; v <- v[!is.na(v)]
        n <- length(v)
        if (n < 2) return(NA)
        nn <- if (is.null(n_fixed)) n else min(n_fixed, n)
        mu <- mean(v)
        hit <- mean(replicate(B_stab, {
          xs <- v[sample.int(n, nn, replace = TRUE)]
          abs(mean(xs) - mu) <= eps * abs(mu)
        }))
        hit >= 0.5
      }), na.rm = TRUE)
    })
  }
  prop_stable <- tibble(indicator = inds,
                        prop_stable = stab_share(NULL),
                        prop_stable_n30 = stab_share(30))

  # -- Phase 3: the country bootstrap of the headline (goal level, Formula 3) --
  mis <- compute_mis(layer_df, coverage_df, "goal")
  w <- mis %>% filter(year %in% c(2011, 2024)) %>%
    select(country, year, MIS) %>%
    tidyr::pivot_wider(names_from = year, values_from = MIS) %>%
    filter(!is.na(`2011`), !is.na(`2024`))
  n_c <- nrow(w)
  full_headline <- mean(w$`2024`) - mean(w$`2011`)
  headline_boot <- map_dfr(seq_len(B_head), function(b) {
    idx <- sample.int(n_c, n_c, replace = TRUE)
    tibble(draw = b, net_change = mean(w$`2024`[idx]) - mean(w$`2011`[idx]))
  })
  headline_summary <- tibble(full = full_headline,
                             boot_mean = mean(headline_boot$net_change),
                             boot_sd = sd(headline_boot$net_change),
                             n_countries = n_c, B = B_head)

  list(sigma_by_n = sigma_by_n, sigma_fits = sigma_fits,
       prop_stable = prop_stable,
       headline_boot = headline_boot, headline_summary = headline_summary)
}

# Algorithm 5 - perturbation of a single indicator by +/-10%: |Sx| = (DeltaC/C)/0.1
#   Statistic: for each indicator j, all of its norm values are multiplied by (1 +/- delta) and the composite value C of each target is recomputed
#   (Formula 2); Sx = the mean over all the (country, year, target) units with C>0 of
#   (DeltaC/C)/delta; |Sx| is the mean of the two directions. The target units that do not contain j have DeltaC=0 and still enter the mean
#   (a system-level effect: the relative change of the average composite value).
#   One row is output per indicator (a target with several mappings does not produce duplicate rows).
algo5_perturbation <- function(norm_df, i2t) {
  map_i2t <- i2t %>% filter(layer == "target") %>% distinct(indicator, layer_id)
  inds <- unique(map_i2t$indicator)
  base <- norm_df %>%
    inner_join(map_i2t, by = "indicator", relationship = "many-to-many") %>%
    group_by(country, year, layer_id) %>%
    summarise(C = mean(norm, na.rm = TRUE), .groups = "drop")
  map_dfr(inds, function(id) {
    sx <- map_dfr(c(1.1, 0.9), function(f) {
      pert <- norm_df %>% filter(indicator == .env$id) %>%
        mutate(norm = norm * f)
      bind_rows(norm_df %>% filter(indicator != .env$id), pert) %>%
        inner_join(map_i2t, by = "indicator", relationship = "many-to-many") %>%
        group_by(country, year, layer_id) %>%
        summarise(C_p = mean(norm, na.rm = TRUE), .groups = "drop") %>%
        inner_join(base, by = c("country", "year", "layer_id")) %>%
        filter(C > 0) %>%
        summarise(Sx = mean(((C_p - C) / C) / (f - 1), na.rm = TRUE),
                  .groups = "drop")
    })
    tibble(indicator = id, Sx_mean_abs = mean(abs(sx$Sx), na.rm = TRUE))
  })
}

# ---------------------------------------------------------------------
# Section 11. weight sensitivity (SI B.2.2; keep_index countries + the target/goal two levels)
#   Dirichlet(1,000 draws) / one PCA component / equal weights; includes dropping indicator 46 (BCH)
#   Composite value = target value (F2) -> MIS (F3); the headline of the goal level = the convention of the paper (3.55)
#   Implementation: matrix algebra (M %*% W computes all the draws in one go), not a loop over the draws
# ---------------------------------------------------------------------

weight_sensitivity <- function(norm_df, i2t, coverage_df, n_draws = 1000, seed = 42) {
  set.seed(seed)
  # B1 convention: only the keep_index countries (the same set of countries as compute_mis / the official MIS)
  keep_cty <- coverage_df %>% filter(keep_index) %>% pull(country) %>% unique()
  keys <- norm_df %>%
    filter(year %in% TRAJ_YEARS, country %in% keep_cty) %>%
    distinct(country, year) %>% mutate(uid = row_number())
  n_u <- nrow(keys)

  # For every layer (target/goal): M (n_u x k, NA = not observed), O (observation mask 0/1)
  build_layer_data <- function(map_i2t) {
    lids <- sort(unique(map_i2t$layer_id))
    ld <- map(lids, function(lid) {
      inds <- sort(unique(map_i2t$indicator[map_i2t$layer_id == lid]))
      k <- length(inds)
      M <- matrix(NA_real_, n_u, k)
      sub <- norm_df %>%
        filter(indicator %in% inds, year %in% TRAJ_YEARS, country %in% keep_cty) %>%
        select(country, year, indicator, norm) %>%
        inner_join(keys, by = c("country", "year"))
      if (nrow(sub) > 0) {
        M[cbind(sub$uid, match(sub$indicator, inds))] <- sub$norm
      }
      list(inds = inds, M = M, O = ifelse(is.na(M), 0, 1))
    })
    names(ld) <- lids
    ld
  }

  # Weight list (a k_l x D matrix per layer) -> MIS matrix (n_u x D)
  #   within a layer: Sum w.x(observed)/Sum w(observed); across layers: the equal-weight arithmetic mean (F3)
  mis_from_weights <- function(W_list, layer_data) {
    lids <- names(layer_data)
    lids <- lids[vapply(lids, function(lid) {
      !is.null(W_list[[lid]]) && nrow(W_list[[lid]]) > 0
    }, logical(1))]
    if (length(lids) == 0) return(NULL)
    D <- ncol(W_list[[lids[1]]])
    acc <- matrix(0, n_u, D); cnt <- matrix(0, n_u, D)
    for (lid in lids) {
      W <- W_list[[lid]]
      ld <- layer_data[[lid]]
      M0 <- ld$M; M0[is.na(M0)] <- 0
      num <- M0 %*% W
      den <- ld$O %*% W
      comp <- ifelse(den > 0, num / den, NA_real_)
      acc <- acc + ifelse(is.na(comp), 0, comp)
      cnt <- cnt + ifelse(is.na(comp), 0, 1)
    }
    ifelse(cnt > 0, acc / cnt, NA_real_)
  }
  headline_from_mis <- function(mis) {
    yr <- keys$year
    colMeans(mis[yr == 2024, , drop = FALSE], na.rm = TRUE) -
      colMeans(mis[yr == 2011, , drop = FALSE], na.rm = TRUE)
  }

  # The four weighting schemes -> the change of the headline of the layer (mean MIS 2024 - mean MIS 2011)
  scheme_headlines <- function(layer_data) {
    # --- equal weights (k x 1 per layer) ---
    eq_W <- map(layer_data, function(ld) {
      k <- length(ld$inds)
      if (k == 0) return(NULL)
      matrix(1 / k, nrow = k, ncol = 1)
    })
    headline_eq <- headline_from_mis(mis_from_weights(eq_W, layer_data))[1]

    # --- Dirichlet random weights (k x n_draws per layer, columns normalised) ---
    dir_W <- map(layer_data, function(ld) {
      k <- length(ld$inds)
      if (k == 0) return(NULL)
      W <- matrix(rgamma(k * n_draws, 1), nrow = k)
      sweep(W, 2, colSums(W), `/`)
    })
    dir_head <- headline_from_mis(mis_from_weights(dir_W, layer_data))

    # --- weights from one PCA component (k x 1 per layer) ---
    pca_W <- map(layer_data, function(ld) {
      k <- length(ld$inds)
      if (k == 0) return(NULL)
      if (k == 1) return(matrix(1, 1, 1))
      X <- ld$M
      ok <- rowSums(!is.na(X)) == k
      if (sum(ok) < 5) return(matrix(1 / k, k, 1))
      pc <- tryCatch(prcomp(X[ok, , drop = FALSE], center = TRUE, scale. = TRUE),
                     error = function(e) NULL)
      if (is.null(pc)) return(matrix(1 / k, k, 1))
      w <- abs(pc$rotation[, 1]); w <- w / sum(w)
      matrix(w, ncol = 1)
    })
    headline_pca <- headline_from_mis(mis_from_weights(pca_W, layer_data))[1]

    # --- equal weights after dropping indicator 46 (BCH, constructed by the authors) ---
    layer_data_no46 <- map(layer_data, function(ld) {
      keep <- ld$inds != 46
      if (!any(keep)) return(NULL)          # the layer contains only 46 -> the whole layer is dropped
      list(inds = ld$inds[keep],
           M = ld$M[, keep, drop = FALSE],
           O = ld$O[, keep, drop = FALSE])
    })
    names(layer_data_no46) <- names(layer_data)
    W_no46 <- map(layer_data_no46, function(ld) {
      if (is.null(ld) || length(ld$inds) == 0) return(NULL)
      matrix(1 / length(ld$inds), nrow = length(ld$inds), ncol = 1)
    })
    headline_no46 <- headline_from_mis(mis_from_weights(W_no46, layer_data_no46))[1]

    c(equal = headline_eq,
      dirichlet_mean = mean(dir_head, na.rm = TRUE),
      dirichlet_sd = sd(dir_head, na.rm = TRUE),
      PCA_first_component = headline_pca,
      equal_excluding_BCH46 = headline_no46)
  }

  # The target level and the goal level are each computed once (the goal level = the net change convention of the paper, 3.55 in the main text)
  map_dfr(c("target", "goal"), function(ly) {
    map_i2t <- i2t %>% filter(layer == .env$ly) %>% distinct(indicator, layer_id)
    h <- scheme_headlines(build_layer_data(map_i2t))
    tibble(scheme = names(h), level = ly, net_change = unname(h))
  })
}

# ---------------------------------------------------------------------
# Section 12. sensitivity to adding the dropped candidates back: not implemented
#   The 9 dropped candidates (Supplementary Table 1, used = N) were all excluded because they did not meet the
#   selection criteria (A.1: too few years, limited country coverage, etc.), so the question of whether the composite
#   results are robust after adding them back does not arise, and this module computes nothing for it.

# ---------------------------------------------------------------------
# Section 13. endpoint vs trajectory robustness (c5_*; SI B.2.3/B.2.4)
#   (1) Extrapolation test: the linear trend of 2011-2023 is extrapolated to 2024 and compared with the measured value -> endpoint bias
#   (2) Leave-one-year-out: the years {2012..2024} are removed one by one and the MIS/EIS headline is recomputed -> sensitivity to the year
# ---------------------------------------------------------------------

endpoint_robustness <- function(layer_df, coverage_df) {
  # (1) Extrapolation test (per country, per indicator)
  extrap <- layer_df %>%
    filter(layer == "target") %>%
    group_by(layer_id, country) %>%
    group_modify(~ {
      d <- .x %>% filter(year %in% TRAJ_YEARS, !is.na(index)) %>% arrange(year)
      if (nrow(d) < 5) return(tibble(mad = NA_real_))
      d23 <- d %>% filter(year <= 2023)
      fit <- lm(index ~ year, data = d23)
      pred24 <- unname(coef(fit)[1] + coef(fit)[2] * 2024)
      obs24 <- d$index[d$year == 2024][1]
      tibble(mad = abs(obs24 - pred24))
    }) %>%
    ungroup()
  extrapolation <- tibble(
    mean_absolute_deviation = mean(extrap$mad, na.rm = TRUE),
    sd = sd(extrap$mad, na.rm = TRUE),
    n = sum(!is.na(extrap$mad))
  )

  # (2) Leave-one-year-out (MIS/EIS headline)
  headline <- function(yrs) {
    mis <- compute_mis(layer_df %>% filter(year %in% yrs), coverage_df, "goal")
    eis <- compute_eis(layer_df %>% filter(year %in% yrs), coverage_df, "goal")
    tibble(
      MIS_change = mean(mis$MIS[mis$year == max(yrs)], na.rm = TRUE) -
        mean(mis$MIS[mis$year == min(yrs)], na.rm = TRUE),
      EIS_change = mean(eis$EIS[eis$year == max(yrs)], na.rm = TRUE) -
        mean(eis$EIS[eis$year == min(yrs)], na.rm = TRUE)
    )
  }
  h_full <- headline(TRAJ_YEARS)
  loo <- map_dfr(setdiff(TRAJ_YEARS, 2011), function(y) {
    h <- headline(setdiff(TRAJ_YEARS, y))
    tibble(omitted_year = y,
           MIS_change = h$MIS_change,
           EIS_change = h$EIS_change,
           dMIS = h$MIS_change - h_full$MIS_change,
           dEIS = h$EIS_change - h_full$EIS_change)
  })
  list(
    extrapolation_test = extrapolation,
    extrapolation_series = extrap %>% filter(!is.na(mad)) %>%
      rename(abs_deviation = mad),
    leave_one_year_out = loo,
    loo_max = tibble(
      max_abs_dMIS = max(abs(loo$dMIS), na.rm = TRUE),
      max_abs_dEIS = max(abs(loo$dEIS), na.rm = TRUE),
      influential_year_MIS = loo$omitted_year[which.max(abs(loo$dMIS))],
      influential_year_EIS = loo$omitted_year[which.max(abs(loo$dEIS))]
    )
  )
}

# ---------------------------------------------------------------------
# Section 14. Figure source data for Figs 1/3/4
# ---------------------------------------------------------------------

# Fig 1a/b: the year-on-year regression slope over 2011-2020 + 95% CI, aggregated by section
fig1_slopes <- function(res) {
  ts <- load_target_section()
  res$baseline_target_goal %>%
    select(layer, layer_id, n_countries, slope, r2, p_value) %>%
    mutate(ci_half = 1.96 * abs(slope) / pmax(sqrt(pmax(n_countries, 2)), 1),
           slope_ci_lo = slope - ci_half,
           slope_ci_hi = slope + ci_half) %>%
    left_join(bind_rows(ts %>% select(layer_id, section),
                        tibble(layer_id = c("A", "B", "C", "D"), section = "Goals")),
              by = "layer_id")
}

# Fig 1c: the slope x R2 quadrant of the 23 targets + the sections
fig1c_quadrant <- function(res) {
  ts <- load_target_section()
  res$baseline_target_goal %>%
    filter(layer == "target") %>%
    select(layer_id, slope, r2) %>%
    left_join(ts %>% select(layer_id, section), by = "layer_id") %>%
    mutate(quadrant = case_when(
      slope >= 0 & r2 >= median(r2, na.rm = TRUE) ~ "steep positive",
      slope >= 0 ~ "shallow positive",
      slope < 0 & r2 >= median(r2, na.rm = TRUE) ~ "steep negative",
      TRUE ~ "shallow negative"
    ))
}

# Fig 3: KW test of MIS/EIS by continent/development group (2011, 2024) + the year-on-year paired tests
fig3_tests <- function(res) {
  gm <- group_membership(res$country_meta)
  d <- res$mis_eis_annual %>%
    inner_join(gm, by = "country") %>%
    filter(year %in% c(2011, 2024))
  kw <- d %>%
    group_by(year) %>%
    summarise(
      kw_MIS_H = {
        ok <- !is.na(MIS)
        if (sum(ok) >= 5 && n_distinct(group[ok]) >= 2)
          safe_kwh(MIS[ok], group[ok]) else NA_real_
      },
      kw_MIS_p = {
        ok <- !is.na(MIS)
        if (sum(ok) >= 5 && n_distinct(group[ok]) >= 2)
          safe_kwp(MIS[ok], group[ok]) else NA_real_
      },
      kw_EIS_p = {
        ok <- !is.na(EIS)
        if (sum(ok) >= 5 && n_distinct(group[ok]) >= 2)
          safe_kwp(EIS[ok], group[ok]) else NA_real_
      },
      .groups = "drop"
    )
  paired <- map_dfr(c("MIS_target", "MIS_goal", "EIS_target", "EIS_goal"),
                    function(m) {
    wide <- res$mis_eis_annual %>%
      select(country, year, all_of(m)) %>%
      filter(year %in% TRAJ_YEARS) %>%
      pivot_wider(names_from = year, values_from = all_of(m))
    map_dfr(2012:2024, function(y) {
      a <- as.numeric(wide[[as.character(y - 1)]])
      b <- as.numeric(wide[[as.character(y)]])
      ok <- !is.na(a) & !is.na(b)
      p <- if (sum(ok) >= 5)
        safe_wilcox(a[ok], b[ok], paired = TRUE) else NA_real_
      tibble(metric = m, year_pair = paste0(y - 1, "-", y),
             mean_diff = mean(b[ok] - a[ok], na.rm = TRUE), p_paired = p)
    })
  })
  list(kw = kw, paired_year_on_year = paired)
}

# Fig 4: binomial test of SyE (the share of (indicator, year) cells with a positive synergy vs 50%)
fig4_synergy_test <- function(res) {
  syd <- res$synergy$syd %>%
    mutate(synergy_positive_aligned = S_t * sgn_S > 0)   # the sign of the gross synergy
  n <- sum(!is.na(syd$synergy_positive_aligned))
  k <- sum(syd$synergy_positive_aligned, na.rm = TRUE)
  bt <- binom.test(k, n, p = 0.5)
  tibble(n_periods = n, n_positive = k,
         share_positive = k / n,
         p_binom = bt$p.value,
         ci_lo = bt$conf.int[1], ci_hi = bt$conf.int[2])
}

# ---------------------------------------------------------------------
# Section 15. Main entry point
# ---------------------------------------------------------------------

# When res = NULL it is restored automatically from the data/intermediate/chk_*.rds checkpoints
run_remaining <- function(res = NULL, save_outputs = TRUE, output_dir = OUTPUT_DIR,
                          with_uncertainty = TRUE, n_boot = 2000, verbose = TRUE) {
  log_msg <- if (verbose) function(...) message("[remaining] ", ...) else function(...) invisible(NULL)
  if (is.null(res)) {
    log_msg("res was not provided, restoring it from the checkpoints (load_res_from_checkpoints) ...")
    res <- load_res_from_checkpoints()
  }
  i2t <- load_i2t()
  if (save_outputs && !dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Every block is written to disk as soon as it is computed, so that an error does not lose the completed tables/figure data
  out <- list()
  put <- function(name, value) {
    out[[name]] <<- value
    if (save_outputs) {
      f <- file.path(output_dir, paste0(name, ".csv"))
      if (is.data.frame(value)) {
        invisible(write_csv(value, f))
        log_msg("    [saved] ", f)
      } else if (is.list(value)) {
        iwalk(value, function(sub, sub_nm) {
          if (is.data.frame(sub)) {
            invisible(write_csv(sub, file.path(output_dir, paste0(name, "__", sub_nm, ".csv"))))
          }
        })
        log_msg("    [saved] ", output_dir, "/", name, "__*.csv")
      }
    }
    invisible(value)
  }

  log_msg("Supplementary Table 10: the projection model per indicator ...")
  put("table_S1_projection_models",
      table_s1_projection_models(res$norm_indicator, res$norm_regional))

  if (!is.null(res$norm_regional) && nrow(res$norm_regional) > 0) {
    log_msg("global/regional summary of the regional indicators (27-29) ...")
    put("regional_indicator_summary", regional_indicator_summary(res$norm_regional))
  }

  log_msg("Supplementary Table 3: baseline regression + geographic heterogeneity + FDR ...")
  s2 <- table_s2_baseline_geo(res)
  put("table_S2_baseline_regression", s2$baseline_regression)
  put("table_S2_kw", s2$kw)
  put("table_S2_posthoc", s2$posthoc)

  log_msg("Supplementary Table 4: trajectories/achievement/classification ...")
  put("table_S3_trajectories", table_s3_trajectories(res))

  log_msg("Supplementary Table 7: annual MIS/EIS by group ...")
  put("table_S4_mis_eis_groups", table_s4_mis_eis_groups(res))

  log_msg("Supplementary Table 8: the country synergy contributions CSC ...")
  put("table_S5_csc", table_s5_csc(res))

  log_msg("Supplementary Table 9: imputation records (cell-by-cell detail imputation_detail.csv) ...")
  # The imputation records cover all 44 indicators (including the regional ones 27-29), aggregated from the cell-by-cell detail
  imp_det <- res$imputation_detail
  if (is.null(imp_det)) {
    imp_det <- imputation_detail(
      if (!is.null(res$raw)) res$raw else res$imputed,
      if (!is.null(res$imputed_all)) res$imputed_all else res$imputed)
  }
  put("imputation_detail", imp_det)
  put("table_S6_imputation", table_s6_imputation(imp_det, res$raw))

  log_msg("target/goal time series (country/region/global, existing data -> end of the projection) ...")
  put("goal_target_ts_country",
      if (!is.null(res$goal_target_ts)) res$goal_target_ts else res$norm_layer)
  if (!is.null(res$layer_geography)) {
    put("goal_target_ts_region", res$layer_geography$region_ts)
    put("goal_target_ts_global", res$layer_geography$global_ts)
  }

  log_msg("annual MIS/EIS series (country/region/global) ...")
  put("mis_eis_annual", res$mis_eis_annual)
  if (!is.null(res$mis_eis_geography)) {
    put("mis_eis_annual_region", res$mis_eis_geography$region)
    put("mis_eis_annual_global", res$mis_eis_geography$global)
  }
  if (!is.null(res$mis_eis_ts)) {
    put("mis_eis_ts_country", res$mis_eis_ts$country)
    put("mis_eis_ts_region", res$mis_eis_ts$region)
    put("mis_eis_ts_global", res$mis_eis_ts$global)
  }

  log_msg("section/goal summaries ...")
  b6 <- section_goal_summary(res, n_boot = n_boot)
  put("b6_section_summary", b6$by_section)
  put("b6_goal_summary", b6$by_goal)
  put("b6_kw_section", b6$kw_section)

  log_msg("trajectory analyses (monotonicity/acceleration/endpoint/COVID/regions) ...")
  mono <- monotonicity_analysis(res$norm_indicator)
  put("b7_monotonicity_overall", mono$overall)
  put("b7_monotonicity_by_indicator", mono$by_indicator)
  accel <- acceleration_analysis(res$norm_indicator)
  put("b7_acceleration", accel$by_indicator)
  put("b7_acceleration_thresholds", accel$thresholds)
  put("b7_endpoint_vs_trajectory", endpoint_vs_trajectory(res$mis_eis_annual))
  put("b7_covid_perturbation", covid_perturbation(res$mis_eis_annual))
  put("b7_regional_paths", regional_annual_paths(res))

  log_msg("intervals of the achievement rate ...")
  c1 <- achievement_ci_report(res, n_boot = n_boot)
  put("c1_achievement_ci", c1)

  log_msg("endpoint robustness ...")
  c5 <- endpoint_robustness(res$norm_layer, res$coverage)
  put("c5_extrapolation_test", c5$extrapolation_test)
  put("c5_extrapolation_series", c5$extrapolation_series)
  put("c5_leave_one_year_out", c5$leave_one_year_out)
  put("c5_loo_max", c5$loo_max)

  log_msg("figure source data for Figs 1/3/4 ...")
  put("fig1_slopes", fig1_slopes(res))
  put("fig1c_quadrant", fig1c_quadrant(res))
  f3 <- fig3_tests(res)
  put("fig3_kw", f3$kw)
  put("fig3_paired", f3$paired_year_on_year)
  put("fig4_synergy_binomial", fig4_synergy_test(res))

  if (with_uncertainty) {
    log_msg("Algorithm 1: exhaustive enumeration of the indicator subsets ...")
    put("uncertainty__a1_indicator_subset",
        algo1_indicator_subset(bind_rows(res$norm_indicator, res$norm_regional), i2t))
    log_msg("Algorithm 2: imputation comparison ...")
    put("uncertainty__a2_imputation",
        algo2_imputation(if (!is.null(res$raw)) res$raw else res$imputed_all))
    log_msg("Algorithm 3: window length (block-bootstrap + three models) ...")
    put("uncertainty__a3_window_length",
        algo3_window_length(bind_rows(res$norm_indicator, res$norm_regional)))
    log_msg("Algorithm 4: country coverage ...")
    put("uncertainty__a4_country_coverage",
        algo4_country_coverage(res$norm_indicator, i2t, res$coverage, res$norm_layer))
    log_msg("Algorithm 5: single-indicator perturbation ...")
    # 27-29 are regional indicators (norm_regional only) and are merged in with the same convention as a1/a3 before the perturbation
    put("uncertainty__a5_perturbation",
        algo5_perturbation(bind_rows(res$norm_indicator, res$norm_regional), i2t))
    log_msg("weight sensitivity (Dirichlet/PCA/equal weights/dropping BCH) ...")
    put("uncertainty__c3_weight_sensitivity",
        weight_sensitivity(res$norm_indicator, i2t, res$coverage))
  }

  if (save_outputs) {
    saveRDS(out, file.path(output_dir, "results_remaining.rds"))
    log_msg("all blocks saved to ", output_dir)
  }
  out
}

# Running the remaining computations: sourcing this file executes them automatically;
#   to load the functions without running them, set GBF_NO_AUTORUN <- TRUE before sourcing
if (!exists("GBF_NO_AUTORUN") || !isTRUE(GBF_NO_AUTORUN)) {
  run_remaining()
}

# ---------------------------------------------------------------------
# Section 16. Smoke test: numerical verification of the SI worked examples (no data needed)
#   1) Formula 4 EIS: (50,50)=100.0, (70,30)=73.7, (80,60,40,20)=62.8
#   2) Formula 1 normalization: OHI 64.3->50.0, EcoFootprint 15->21.3
#   3) Formulas 10-13: PGR_sim=37.5%, PGR_req=150%, AR(3 countries)=33.3%
#   4) Formulas 5-9: FC=15, SyE=10, CSC=(50,0,50) and =100%
#   5) Formula 3 MIS=60; the imputation example (43,45,50)
# ---------------------------------------------------------------------

smoke_test <- function() {
  message("[smoke] Formula 4 EIS ...")
  stopifnot(abs(eis_from_values(c(50, 50)) - 100.0) < 0.05)
  stopifnot(abs(eis_from_values(c(70, 30)) - 73.7) < 0.05)
  stopifnot(abs(eis_from_values(c(30, 70)) - 73.7) < 0.05)
  stopifnot(abs(eis_from_values(c(50, 50, 50, 50)) - 100.0) < 0.05)
  stopifnot(abs(eis_from_values(c(80, 60, 40, 20)) - 62.8) < 0.05)
  message("  EIS(50,50)=100.0 / (70,30)=73.7 / (80,60,40,20)=62.8 ✓")

  message("[smoke] Formula 1 normalization ...")
  stopifnot(abs(100 * (64.3 - 28.6) / (100 - 28.6) - 50.0) < 0.05)
  stopifnot(abs(100 * (18.89 - 15) / (18.89 - 0.61) - 21.3) < 0.05)
  stopifnot(clamp100(100 * (110 - 28.6) / (100 - 28.6)) == 100)
  stopifnot(clamp100(100 * (25 - 28.6) / (100 - 28.6)) == 0)
  message("  OHI 64.3→50.0 / EcoFootprint 15→21.3 ✓")

  message("[smoke] Formula 2/3 ...")
  stopifnot(abs(mean(c(80, 60)) - 70) < 1e-9)
  stopifnot(abs(mean(c(70, 30, 50, 90)) - 60) < 1e-9)
  message("  F2 (80,60)→70 / F3 MIS(70,30,50,90)=60 ✓")

  message("[smoke] Formulas 10-13 (the worked example of Target 3 with three countries) ...")
  pgr_sim <- c(A = 100 * (110 - 50) / 50, B = 100 * (85 - 40) / 40, C = 100 * (60 - 30) / 30)
  pgr_req <- c(A = 100 * (100 - 50) / 50, B = 100 * (100 - 40) / 40, C = 100 * (100 - 30) / 30)
  stopifnot(abs(pgr_sim["A"] - 120) < 1e-9, abs(pgr_req["A"] - 100) < 1e-9)
  stopifnot(abs(pgr_sim["B"] - 112.5) < 1e-9, abs(pgr_req["B"] - 150) < 1e-9)
  stopifnot(abs(pgr_sim["C"] - 100) < 1e-9, abs(pgr_req["C"] - 233.3) < 0.05)
  pgr_gap <- pgr_req - pgr_sim
  stopifnot(abs(pgr_gap["A"] - (-20)) < 1e-9)   # surplus
  stopifnot(abs(pgr_gap["B"] - 37.5) < 1e-9)    # shortfall
  # Branch without a threshold: Y = 40 + 1.0*10 = 50 -> PGR_req = 25%; PGR = 25 - 37.5 = -12.5
  stopifnot(abs(100 * (50 - 40) / 40 - 25) < 1e-9)
  stopifnot(abs((25 - 37.5) - (-12.5)) < 1e-9)
  # I = 1{PGR_sim >= PGR_req and the bootstrap p<0.05}:
  #   A=1 (p=0.02, 120>=100); B=0 (112.5<150); C=0 (p=0.07, not significant)
  p_boot <- c(A = 0.02, B = 0.03, C = 0.07)
  I <- as.integer(pgr_sim >= pgr_req & p_boot < 0.05)
  stopifnot(identical(unname(I), c(1L, 0L, 0L)))
  stopifnot(abs(mean(I) * 100 - 33.3) < 0.1)    # AR = 33.3%
  message("  PGR_sim/req/gap + AR(Target 3)=33.3% ✓")

  message("[smoke] the synergy worked example of Formulas 5-9 ...")
  # eq. (9)  d_t^c = |Dr|            eq. (10) d_t = sqrt(Sum Dr^2)
  # eq. (11) S_t = (1/C)(Sum d_t^c - d_t)   eq. (13) CSC = (d_t^c - d_t/C)/(S_t*sgnS*C)x100%
  syn <- function(d, expected_dir = 1) {
    C <- length(d)
    d_dom <- abs(d)
    d_for <- sqrt(sum(d^2))
    S_t <- (sum(d_dom) - d_for) / C
    g_t <- mean(d)
    sgn_S <- ifelse(sign(g_t) == expected_dir, 1, -1)
    csc <- (d_dom - d_for / C) / (S_t * sgn_S * C) * 100
    list(d_for = d_for, S_t = S_t, g_t = g_t, sgn_S = sgn_S, csc = csc)
  }
  p1 <- syn(c(10, -5, 10))
  stopifnot(abs(p1$d_for - 15) < 1e-9)               # docx: d_t = 15
  stopifnot(abs(p1$S_t - 10 / 3) < 1e-9)             # docx: S_t = (25−15)/3 = 3.33
  stopifnot(abs(p1$g_t - 5) < 1e-9, p1$sgn_S == 1)   # docx: g_t = +5, sgn(S_t) = +1
  stopifnot(all(abs(p1$csc - c(50, 0, 50)) < 1e-9))  # docx: CSC = 50%, 0%, 50%
  stopifnot(abs(sum(p1$csc) - 100) < 1e-6)
  p2 <- syn(c(-10, -5, -10))
  stopifnot(abs(p2$S_t - 10 / 3) < 1e-9)             # docx: S_t is also = 3.33
  stopifnot(abs(p2$g_t + 25 / 3) < 1e-9, p2$sgn_S == -1)  # docx: g_t = −8.33, sgnS = −1
  stopifnot(all(abs(p2$csc - c(-50, 0, -50)) < 1e-9))
  stopifnot(abs(sum(p2$csc) + 100) < 1e-6)           # docx: Σ = −100%
  message("  d_t=15 / S_t=3.33 / CSC=(50,0,50)sum=100% and (-50,0,-50)sum=-100% OK")

  message("[smoke] the A.2 imputation example (41,NA,NA,47,50,NA) ...")
  v <- impute_series(2015:2020, c(41, NA, NA, 47, 50, NA))
  stopifnot(abs(v[2] - 43.0) < 1e-9, abs(v[3] - 45.0) < 1e-9, abs(v[6] - 50) < 1e-9)
  message("  43.0 / 45.0 / 50 ✓")

  message("[smoke] A.2 creation of the rows for the missing years + imputation from the full-history anchors (sparse grids of the 18/39 type) ...")
  # The observations are only in 2010/2015/2020 (every five years): the missing years inside the window must be created and then imputed from the full-history anchors
  raw <- tibble(country = "X", year = c(2010, 2015, 2020),
                value = c(0.9, 0.8, 0.7), indicator = 18L)
  imp <- impute_all(raw, impute_years = 2011:2024, used_ids = 18L)
  getv <- function(y) imp$value[imp$year == y]
  stopifnot(nrow(imp) == 15)                    # the 2010 observation + the rows created for the 2011-2024 window
  stopifnot(abs(getv(2011) - 0.88) < 1e-9,      # 2010->2015 linear
            abs(getv(2013) - 0.84) < 1e-9,
            abs(getv(2014) - 0.82) < 1e-9,
            abs(getv(2016) - 0.78) < 1e-9,      # 2015->2020 linear
            abs(getv(2019) - 0.72) < 1e-9,
            abs(getv(2021) - 0.70) < 1e-9,      # endpoint forward fill
            abs(getv(2024) - 0.70) < 1e-9)
  det <- imputation_detail(raw, imp)
  d11 <- det %>% filter(year == 2011)
  d21 <- det %>% filter(year == 2021)
  stopifnot(d11$status == "imputed", d11$method == "linear_interpolation",
            d11$anchor_year_before == 2010, d11$anchor_year_after == 2015,
            d21$status == "imputed", d21$method == "forward_fill",
            d21$anchor_year_before == 2020)
  # Full-history anchors: observations outside the window can also support the imputation inside it
  raw2 <- tibble(country = "X", year = 1990, value = 5, indicator = 39L)
  imp2 <- impute_all(raw2, impute_years = 2011:2024, used_ids = 39L)
  stopifnot(all(abs(imp2$value[imp2$year %in% 2011:2024] - 5) < 1e-9))
  message("  rows created for the missing years + linear/endpoint fills + the anchors recorded OK")

  message("[smoke] the non-missing regression of the target-level MIS/EIS of Formulas 3/4 ...")
  layer_t <- tibble(
    country = rep("X", 16), year = rep(2011:2012, each = 8),
    layer = rep(rep(c("goal", "target"), each = 4), 2),
    layer_id = rep(c("A", "B", "C", "D", "T1", "T2", "T3", "T4"), 2),
    index = c(70, 30, 50, 90, 60, 40, 80, 20, 70, 30, 50, 90, 60, 40, 80, 20),
    n_ind = 1L)
  cov_t <- tibble(country = "X", keep_index = TRUE)
  me <- compute_mis_eis_annual(layer_t, cov_t)
  stopifnot(!is.na(me$MIS_target[1]), !is.na(me$EIS_target[1]),
            !is.na(me$MIS_goal[1]), !is.na(me$EIS_goal[1]))
  stopifnot(abs(me$MIS_goal[1] - 60) < 1e-9,    # mean(70,30,50,90)
            abs(me$MIS_target[1] - 50) < 1e-9)  # mean(60,40,80,20)
  message("  MIS_goal=60 / MIS_target=50 are both non-NA OK")

  message("[smoke] the complete target/goal time series (existing data -> end of the projection) ...")
  layer_x <- tibble(country = rep("X", 10), year = 2011:2020,
                    layer = "target", layer_id = "T1",
                    index = 40 + 2 * (0:9), n_ind = 1L)
  paths_x <- tibble(country = "X", year = 2021:2030,
                    value = 58 + 2 * (0:9), ci_lo = 56 + 2 * (0:9),
                    ci_hi = 60 + 2 * (0:9),
                    layer = "target", layer_id = "T1",
                    model_type = "PMA", model_form = "linear")
  ts_x <- goal_target_ts_full(layer_x, paths_x, id_cols = "country")
  stopifnot(nrow(ts_x) == 20)   # the 2011-2020 history + the 2021-2030 projection
  stopifnot(all(ts_x$value_type[ts_x$year <= 2020] == "historical"),
            all(ts_x$value_type[ts_x$year > 2020] == "projected"),
            max(ts_x$year) == 2030,
            abs(ts_x$index[ts_x$year == 2030] - 76) < 1e-9)
  message("  the 2011-2020 history + the 2021-2030 projection are continuous up to the target year OK")

  message("[smoke] the A.5 rho classification ...")
  stopifnot(classify_rho(1.2, TRUE) == "Full achievement")
  stopifnot(classify_rho(0.7, TRUE) == "Favourable")
  stopifnot(classify_rho(0.3, TRUE) == "Mixed")
  stopifnot(classify_rho(1.2, FALSE) == "Mixed")
  stopifnot(classify_rho(0.3, FALSE) == "Unfavourable")
  stopifnot(classify_rho(-0.1, FALSE) == "Unfavourable")
  message("  the four rho classes OK")

  message("SMOKE TEST PASSED")
  invisible(TRUE)
}

# Automatic smoke test (uncomment to enable):
# smoke_test()
