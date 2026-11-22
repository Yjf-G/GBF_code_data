# =====================================================================
# endfill_sensitivity.R - sensitivity to endpoint forward filling (the LOCF tail)  v2
#   (reviewer comment: forward filling to 2024 creates an artificial zero change at the end of the outcome series,
#     which affects the year-pair Wilcoxon tests and the trajectory fits; report the true last observation year + rerun truncated)
#   1) endfill_last_observed.csv - the true last observation year of each indicator (global / median and minimum across countries)
#      plus the number of endpoint-filled cells and the tail length (used for the new columns of Table 9)
#   2) Truncation rule: drop the imputed cells with method=forward_fill (verified to lie entirely after the last
#      observation of each series); internal linear interpolations and the backward fill after the first year are kept
#   3) The rerun convention = **the difference on the common indicator set** (to avoid compositional drift):
#      for every country and pair (y0,y1) take only the indicators with values in both years to obtain Delta norm -> target-level Delta -> MIS Delta;
#      for EIS, rebuild the target values of the two years on the same common indicator set and take the difference.
#      On the full panel this convention must reproduce table_6 cell by cell (checked before the truncated run).
#   4) Analysis: a. year-pair Wilcoxon, 52 tests (full vs truncated)
#            b. endpoint net change 2011->2024 (common-set convention)
#            c. indicator-level slopes, full vs truncated + DPSIR class summaries + action/outcome comparison
#   Written to: data/output/endfill_*.csv
# Run from the repository root:  Rscript fig_code/endfill_sensitivity.R
# =====================================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
})

put <- function(df, name) {
  write_csv(df, file.path("data/output", paste0(name, ".csv")))
  message("[out] data/output/", name, ".csv  (", nrow(df), " rows)")
}

# ---- Formula 2/3/4 (identical formula by formula to gbf_compute_pipeline.R) ----
mean_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else mean(x)
}
eis_from_values <- function(r) {
  r <- r[!is.na(r)]
  n <- length(r)
  if (n < 2) return(NA_real_)
  S <- sum((1 / n) * pi * r^2)
  L <- 2 * (max(r) - min(r)) + (2 * pi / n) * sum(r)
  if (L <= 0) return(NA_real_)
  400 * pi * S / L^2
}

# ---- 0. Data ----
norm_raw <- read_csv("data/intermediate/norm_indicator.csv", show_col_types = FALSE)
imp <- read_csv("data/intermediate/imputation_detail.csv", show_col_types = FALSE)
cov <- read_csv("data/intermediate/country_coverage.csv", show_col_types = FALSE)
i2t <- read_csv("data/indicator_to_target.csv", show_col_types = FALSE)
dps <- read_csv("data/indicator_dpsir.csv", show_col_types = FALSE)
keep <- cov$country[cov$keep_index]

map <- i2t %>%
  transmute(indicator = as.integer(indicator_id),
            layer = ifelse(grepl("^[0-9]+$", target), "target", "goal"),
            layer_id = ifelse(layer == "target", paste0("T", target), target))

norm_w <- norm_raw %>%
  filter(year %in% 2011:2024, country %in% keep) %>%
  left_join(imp %>% select(indicator, country, year, status, method),
            by = c("indicator", "country", "year")) %>%
  mutate(end_fill = !is.na(method) & method == "forward_fill" & status == "imputed")

# ---- 1. True last observation year per indicator + the tail ----
obs_c <- imp %>% filter(status == "observed") %>%
  group_by(indicator, country) %>%
  summarise(last_obs_c = max(year), .groups = "drop")
last_obs <- obs_c %>% group_by(indicator) %>%
  summarise(global_last_obs = max(last_obs_c),
            median_last_obs = median(last_obs_c),
            min_last_obs = min(last_obs_c),
            n_countries = n_distinct(country), .groups = "drop")
tail_by_c <- norm_w %>% filter(end_fill) %>% count(indicator, country, name = "tail_n")
tails <- norm_w %>% group_by(indicator) %>%
  summarise(n_cells_window = n(),
            n_end_filled = sum(end_fill),
            share_end_filled = mean(end_fill),
            share_2021_2024_end_filled = mean(end_fill[year >= 2021]),
            .groups = "drop") %>%
  left_join(tail_by_c %>% group_by(indicator) %>%
              summarise(max_tail_years = max(tail_n), .groups = "drop"),
            by = "indicator") %>%
  mutate(max_tail_years = ifelse(is.na(max_tail_years), 0L, max_tail_years))
tbl_last <- last_obs %>% left_join(tails, by = "indicator") %>%
  left_join(dps %>% transmute(indicator = indicator_id, indicator_name,
                              dpsir_class = as.character(dpsir_class), used),
            by = "indicator") %>%
  relocate(indicator_name, dpsir_class, used, .after = indicator) %>%
  arrange(global_last_obs, indicator)
put(tbl_last, "endfill_last_observed")

# ---- 2. The two difference conventions ----
# common: the common indicator set at the indicator level (indicators with values in both years), fixed composition, no compositional drift
# literal: the convention of the paper - the year-specific sets are combined first and the difference is taken afterwards
pair_common <- function(nw, ly, y0, y1) {
  a <- nw %>% filter(year == y0) %>% select(country, indicator, v0 = norm)
  b <- nw %>% filter(year == y1) %>% select(country, indicator, v1 = norm)
  inner_join(a, b, by = c("country", "indicator")) %>%
    filter(!is.na(v0), !is.na(v1)) %>%
    inner_join(map %>% filter(layer == ly), by = "indicator",
               relationship = "many-to-many")
}
mis_delta_common <- function(nw, ly, y0, y1) {
  pair_common(nw, ly, y0, y1) %>%
    group_by(country, layer_id) %>%
    summarise(d = mean(v1 - v0), .groups = "drop") %>%
    group_by(country) %>%
    summarise(d = mean(d), .groups = "drop")
}
eis_delta_common <- function(nw, ly, y0, y1) {
  pair_common(nw, ly, y0, y1) %>%
    group_by(country, layer_id) %>%
    summarise(v0 = mean_na(v0), v1 = mean_na(v1), .groups = "drop") %>%
    group_by(country) %>%
    summarise(d = eis_from_values(v1) - eis_from_values(v0), .groups = "drop")
}
comp_series <- function(nw, ly, met) {
  ld <- nw %>%
    inner_join(map %>% filter(layer == ly), by = "indicator",
               relationship = "many-to-many") %>%
    group_by(country, year, layer_id) %>%
    summarise(index = mean_na(norm), .groups = "drop")
  if (met == "MIS") {
    ld %>% group_by(country, year) %>%
      summarise(value = mean_na(index), .groups = "drop")
  } else {
    ld %>% group_by(country, year) %>%
      summarise(value = eis_from_values(index), .groups = "drop")
  }
}
delta_literal <- function(nw, s, y0, y1) {
  ly <- if (grepl("target", s)) "target" else "goal"
  met <- sub(paste0("_", ly), "", s)
  comp_series(nw, ly, met) %>%
    filter(year %in% c(y0, y1)) %>%
    pivot_wider(names_from = year, values_from = value) %>% drop_na() %>%
    transmute(d = .data[[as.character(y1)]] - .data[[as.character(y0)]]) %>%
    pull(d)
}

# ---- 3. Year-pair Wilcoxon (same convention as make_tables.R) ----
wtest <- function(d, series, y0, variant) {
  d <- d[is.finite(d) & abs(d) > 1e-12]
  if (length(d) < 3) return(NULL)
  wt <- wilcox.test(d, exact = FALSE, correct = TRUE)
  wtc <- suppressWarnings(tryCatch(
    wilcox.test(d, exact = FALSE, correct = TRUE, conf.int = TRUE),
    error = function(e) NULL))
  n <- length(d)
  tibble(series = series, y0 = y0, variant = variant, n = n,
         mean_delta = mean(d), W = unname(wt$statistic), p = wt$p.value,
         r_rank_biserial = 2 * unname(wt$statistic) / (n * (n + 1) / 2) - 1,
         hl_est = if (is.null(wtc)) NA_real_ else unname(wtc$estimate),
         hl_lo = if (is.null(wtc)) NA_real_ else unname(wtc$conf.int[1]),
         hl_hi = if (is.null(wtc)) NA_real_ else unname(wtc$conf.int[2]))
}

delta_of <- function(nw, s, y0, y1, mode = "common") {
  if (mode == "literal") return(delta_literal(nw, s, y0, y1))
  ly <- if (grepl("target", s)) "target" else "goal"
  met <- sub(paste0("_", ly), "", s)
  if (met == "MIS") mis_delta_common(nw, ly, y0, y1)$d else eis_delta_common(nw, ly, y0, y1)$d
}

nw_trunc <- norm_w %>% filter(!end_fill)
SER <- c("MIS_target", "MIS_goal", "EIS_target", "EIS_goal")
GRID <- list(c("full_literal", "full_pairwise"),
             c("truncated_literal", "truncated_pairwise"))
run_grid <- function(nw, mode, variant) {
  bind_rows(lapply(SER, function(s) bind_rows(lapply(2011:2023, function(y)
    wtest(delta_of(nw, s, y, y + 1, mode), s, y, variant)))))
}
yos <- bind_rows(run_grid(norm_w, "literal", "full_literal"),
                 run_grid(norm_w, "common", "full_common"),
                 run_grid(nw_trunc, "literal", "truncated_literal"),
                 run_grid(nw_trunc, "common", "truncated_common"))

# Check 1: full_literal must reproduce table_6 cell by cell (the convention of the paper)
# Check 2: the difference between full_common and table_6 (arising only from the small compositional differences of the missing cells)
t6 <- read_csv("data/tables/table_6_yoy_paired_tests.csv", show_col_types = FALSE) %>%
  select(series, y0, W_ref = W, p_ref = p, mean_ref = mean_delta)
chk <- yos %>% filter(variant == "full_literal") %>%
  inner_join(t6, by = c("series", "y0"))
max_dW <- max(abs(chk$W - chk$W_ref)); max_dp <- max(abs(chk$p - chk$p_ref))
message(sprintf("[check] full_literal reproduces table_6: max|dW|=%.3g  max|dp|=%.3g",
                max_dW, max_dp))
stopifnot(max_dW == 0, max_dp < 1e-12)
chk2 <- yos %>% filter(variant == "full_common") %>%
  inner_join(t6, by = c("series", "y0"))
message(sprintf("[note] full_common vs table_6: max|d mean|=%.3g  max|dW|=%.3g  max|dp|=%.3g",
                max(abs(chk2$mean_delta - chk2$mean_ref)),
                max(abs(chk2$W - chk2$W_ref)), max(abs(chk2$p - chk2$p_ref))))

yos <- yos %>%
  mutate(sig = p < 0.05, direction = ifelse(mean_delta > 0, "up", "down")) %>%
  left_join(t6 %>% select(series, y0, p_ref, W_ref), by = c("series", "y0"))
put(yos, "endfill_yoy_wilcox")

yoy_cmp <- yos %>%
  select(series, y0, variant, n, mean_delta, p, r_rank_biserial, hl_est,
         sig, direction) %>%
  pivot_wider(names_from = variant,
              values_from = c(n, mean_delta, p, r_rank_biserial, hl_est,
                              sig, direction)) %>%
  mutate(sig_change_common = sig_full_common != sig_truncated_common |
           direction_full_common != direction_truncated_common,
         sig_change_literal = sig_full_literal != sig_truncated_literal |
           direction_full_literal != direction_truncated_literal)
put(yoy_cmp, "endfill_yoy_wilcox_diff")

# ---- 4. Endpoint net change 2011->2024 (both conventions) ----
net <- bind_rows(
  bind_rows(lapply(SER, function(s)
    wtest(delta_of(norm_w, s, 2011, 2024, "literal"), s, 2011, "full_literal"))),
  bind_rows(lapply(SER, function(s)
    wtest(delta_of(norm_w, s, 2011, 2024, "common"), s, 2011, "full_common"))),
  bind_rows(lapply(SER, function(s)
    wtest(delta_of(nw_trunc, s, 2011, 2024, "common"), s, 2011, "truncated_common")))) %>%
  mutate(year_pair = "2011-2024") %>%
  select(series, year_pair, variant, n, mean_delta, W, p, r_rank_biserial,
         hl_est, hl_lo, hl_hi)
put(net, "endfill_net_change")

# ---- 5. Indicator-level slopes / net change, full vs truncated + class summaries + action/outcome ----
ind_stats <- function(nw, variant) {
  nw %>% group_by(indicator, year) %>%
    summarise(norm = mean_na(norm), .groups = "drop") %>%
    group_by(indicator) %>%
    summarise(n_years = sum(!is.na(norm)),
              v_first = norm[which.min(year)],
              v_last = norm[which.max(year)],
              first_year = min(year), last_year = max(year),
              net_own = v_last - v_first,
              slope = unname(coef(lm(norm ~ year))[2]),
              r2 = summary(lm(norm ~ year))$r.squared, .groups = "drop") %>%
    mutate(variant = variant)
}
sl <- bind_rows(ind_stats(norm_w, "full"), ind_stats(nw_trunc, "truncated")) %>%
  left_join(dps %>% transmute(indicator = indicator_id,
                              dpsir_class = as.character(dpsir_class),
                              side = ifelse(dpsir_class %in% c("Process", "Response"),
                                            "action", "outcome")),
            by = "indicator") %>%
  left_join(tbl_last %>% select(indicator, global_last_obs, n_end_filled,
                                max_tail_years), by = "indicator")
put(sl, "endfill_slopes")

sl_cmp <- sl %>%
  select(indicator, dpsir_class, side, global_last_obs, max_tail_years,
         variant, slope, r2, net_own) %>%
  pivot_wider(names_from = variant, values_from = c(slope, r2, net_own)) %>%
  mutate(slope_diff = slope_truncated - slope_full)
put(sl_cmp, "endfill_slopes_diff")

# Class summaries (median of the truncated slopes and of the net change over each indicator's own window)
cls_sum <- sl_cmp %>% group_by(dpsir_class, side) %>%
  summarise(n = n(),
            med_slope_full = median(slope_full),
            med_slope_trunc = median(slope_truncated),
            med_net_own_trunc = median(net_own_truncated),
            .groups = "drop")
put(cls_sum, "endfill_class_summary")

# action vs outcome (MWU rerun on the truncated slopes)
mwu_pair <- function(x1, x2) {
  wt <- wilcox.test(x1, x2, exact = FALSE)
  c(median1 = median(x1), median2 = median(x2), p = wt$p.value,
    rrb = 2 * unname(wt$statistic) / (length(x1) * length(x2)) - 1)
}
ao_full <- mwu_pair(sl_cmp$slope_full[sl_cmp$side == "action"],
                    sl_cmp$slope_full[sl_cmp$side == "outcome"])
ao_tr <- mwu_pair(sl_cmp$slope_truncated[sl_cmp$side == "action"],
                  sl_cmp$slope_truncated[sl_cmp$side == "outcome"])
ao_row <- function(v, label) tibble(
  analysis = label, action_median = unname(v["median1"]),
  outcome_median = unname(v["median2"]), mwu_p = unname(v["p"]),
  rank_biserial = unname(v["rrb"]))
ao_summary <- bind_rows(ao_row(ao_full, "action_vs_outcome_slope_full"),
                        ao_row(ao_tr, "action_vs_outcome_slope_truncated"))
put(ao_summary, "endfill_action_outcome")

message("ENDFILL SENSITIVITY DONE")
