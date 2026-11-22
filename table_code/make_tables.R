# =====================================================================
# table_code/make_tables.R — Supplementary Tables 3-10 (numbering of the submitted Supplementary Information)
#   Generated here: 3 (baseline regression) / 4 (trajectory classes) / 5 (Global MIS/EIS 2011-2024) /
#                   6 (year-on-year paired tests) / 7 (MIS/EIS by group) / 8 (CSC) /
#                   9 (imputation records) / 10 (projection models)
#   Compiled manually and therefore not generated here: 1 (list of the 53 indicators) / 2 (element alignment)
#   Committed and not rewritten here: 11 (countries excluded for insufficient data coverage,
#                   table_11_excluded_countries.csv, assembled from the country-coverage products of the pipeline)
#   Source data: data/output/table_S*.csv + data/intermediate/*.csv
#   Output: data/tables/table_N_*.csv
#   Note: Tables 1 and 2 are the indicator list/element alignment (compiled manually, see the Supplementary
#         Information) and are not part of the computational pipeline. The internal product names
#         table_S1-S6 (data/output/) predate the final table numbering; the mapping is given below.
#   Run from the repository root:  Rscript table_code/make_tables.R
# =====================================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
})
OUT <- "data/tables"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- Table 3 — baseline regression + KW(slope) + leading post-hoc pair ----
# Definition (consistent with docx Supplementary Table 3 / figure legend, 2026-11-05):
#   slope/R²/p = OLS of the country-mean baseline series 2011–2020 (n=10 years); the pipeline
#     r2/p in table_S2_baseline_regression are from the pooled country regression (same slope, different R²/p) and are overwritten here.
#   leading pair = the largest |mean diff| among the 10 continent pairs (largest within-layer contrast), not the smallest p;
#     layers with a non-significant omnibus test have no post-hoc row, corresponding to an empty value.
b  <- read_csv("data/output/table_S2_baseline_regression.csv", show_col_types = FALSE)
gts <- read_csv("data/intermediate/goal_target_ts_global.csv", show_col_types = FALSE) %>%
  filter(value_type == "historical", group_type == "global",
         year >= 2011, year <= 2020)
fit_gl <- gts %>% group_by(layer, layer_id) %>% group_modify(~ {
  m <- lm(index ~ year, data = .x)
  tibble(slope = unname(coef(m)[2]),
         r2 = summary(m)$r.squared,
         p_value = coef(summary(m))[2, 4])
}) %>% ungroup()
b <- b %>% select(-slope, -r2, -p_value) %>%
  left_join(fit_gl, by = c("layer", "layer_id")) %>%
  relocate(slope, r2, p_value, .after = n_countries)
kw <- read_csv("data/output/table_S2_kw.csv", show_col_types = FALSE) %>%
  filter(param == "slope") %>% select(layer, layer_id, kw_H_slope = kw_H, kw_p_slope = kw_p)
ph <- read_csv("data/output/table_S2_posthoc.csv", show_col_types = FALSE) %>%
  filter(param == "slope")
lead <- ph %>% group_by(layer, layer_id) %>%
  slice_max(abs(mean_diff), n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(layer, layer_id,
            leading_pair = paste(cont1, cont2, sep = " vs "),
            leading_mean_diff = mean_diff, leading_p_bonf = u_p_bonf)
t3 <- b %>% left_join(kw, by = c("layer", "layer_id")) %>%
  left_join(lead, by = c("layer", "layer_id"))
write_csv(t3, file.path(OUT, "table_3_baseline_regression.csv"))

# ---- Table 4 — post-2020 trajectories, achievement and classes ----
# The Class column uses the global-level determination of the manuscript: when the required growth
#   rate is zero (T6), a non-negative simulated growth is recorded as "Full achievement" (SI A.5).
# The Achievement classes (%) column gives the country-level class shares (Full/Fav./Mixed/Unfav.,
#   two decimals) over the country-layer pairs with defined growth rates, computed with the same
#   convention as Extended Data Fig. 1 (classify_rho of gbf_compute_pipeline.R + the PGR_req = 0
#   special case); the four displayed shares are adjusted by 0.01 where needed so that they sum to
#   100.00, as in the manuscript table. n_full/n_fav/n_mixed/n_unfav give the underlying counts.
GBF_NO_AUTORUN <- TRUE
source("gbf_compute_pipeline.R")
pgr <- read_csv("data/intermediate/pgr.csv", show_col_types = FALSE)
CLS <- c("Full achievement", "Favourable", "Mixed", "Unfavourable")
cnt <- pgr %>%
  filter(!is.na(pgr_sim), !is.na(pgr_req)) %>%
  mutate(
    sig_pos = significant %in% TRUE & rate > 0 & !is.na(p_boot) & p_boot < 0.05,
    rho = ifelse(pgr_req != 0, pgr_sim / pgr_req, ifelse(pgr_sim >= 0, Inf, -Inf)),
    class = ifelse(pgr_req == 0 & pgr_sim >= 0, "Full achievement",
                   classify_rho(rho, sig_pos))) %>%
  count(layer, layer_id, class) %>%
  complete(nesting(layer, layer_id), class = CLS, fill = list(n = 0)) %>%
  mutate(class = factor(class, levels = CLS)) %>%
  arrange(layer, layer_id, class) %>%
  group_by(layer, layer_id) %>%
  summarise(n_full  = n[class == CLS[1]], n_fav   = n[class == CLS[2]],
            n_mixed = n[class == CLS[3]], n_unfav = n[class == CLS[4]],
            .groups = "drop")
share_str <- function(f, v, m, u) {
  raw <- 100 * c(f, v, m, u) / (f + v + m + u)
  d <- round(raw, 2)
  if (abs(sum(d) - 100) > 1e-9) {
    i <- which.max(round(abs(d - raw), 9))   # ties resolved to the first class
    d[i] <- d[i] + ifelse(sum(d) < 100, 0.01, -0.01)
  }
  paste(sprintf("%.2f", d), collapse = "/")
}
t4 <- read_csv("data/output/table_S3_trajectories.csv", show_col_types = FALSE) %>%
  left_join(cnt, by = c("layer", "layer_id")) %>%
  mutate(
    class = ifelse(is.na(class) & pgr_req_mean == 0 & pgr_sim_mean >= 0,
                   "Full achievement", class),
    class_shares_pct = mapply(share_str, n_full, n_fav, n_mixed, n_unfav)) %>%
  relocate(n_full, n_fav, n_mixed, n_unfav, class_shares_pct, .before = n)
write_csv(t4, file.path(OUT, "table_4_trajectories_classes.csv"))

# ---- Table 5 — Global MIS and EIS during 2011–2024 ----
g <- read_csv("data/intermediate/mis_eis_annual_global.csv", show_col_types = FALSE) %>%
  arrange(year)
vals <- g %>% transmute(Year = year,
                        MIS_goal = round(MIS_goal, 1), MIS_target = round(MIS_target, 1),
                        EIS_goal = round(EIS_goal, 1), EIS_target = round(EIS_target, 1),
                        n = as.character(n_countries))
net_cell <- function(disp_net, disp_base) {
  # Consistent with the manuscript Table 5: the net change is the difference of the displayed
  #   values, and the percentage is computed from the same displayed values (integer per cent
  #   when at least 1%, one decimal below that): +3.5 (+11%), +6.2 (+15%), +0.3 (+0.5%), +1.2 (+2%)
  pct <- 100 * disp_net / disp_base
  if (abs(pct) >= 1) sprintf("%+.1f (+%d%%)", disp_net, round(pct))
  else sprintf("%+.1f (+%.1f%%)", disp_net, round(pct, 1))
}
COLS <- c("MIS_goal", "MIS_target", "EIS_goal", "EIS_target")
v24 <- vals %>% filter(Year == 2024) %>% select(all_of(COLS)) %>% as.data.frame()
v11 <- vals %>% filter(Year == 2011) %>% select(all_of(COLS)) %>% as.data.frame()
disp_net <- as.list(v24 - v11)
net_row <- tibble(
  Year = "Net change",
  MIS_goal   = net_cell(disp_net$MIS_goal,   v11$MIS_goal),
  MIS_target = net_cell(disp_net$MIS_target, v11$MIS_target),
  EIS_goal   = net_cell(disp_net$EIS_goal,   v11$EIS_goal),
  EIS_target = net_cell(disp_net$EIS_target, v11$EIS_target),
  n = "/")
t5 <- bind_rows(vals %>% mutate(Year = as.character(Year),
                                # one decimal in every cell, as in the manuscript table ("41.0", not "41")
                                across(all_of(COLS), ~ sprintf("%.1f", .x))), net_row)
write_csv(t5, file.path(OUT, "table_5_global_mis_eis.csv"))

# ---- Table 6 — year-on-year paired tests (SI B.2.7) ----
# 2026-11-10: effect sizes added in response to reviewer comments (reported alongside the Wilcoxon results):
#   r_rank_biserial = rank-biserial correlation for matched pairs (2W/n(n+1) − 1), −1..1, sign = direction of the median change;
#   hl_est/lo/hi   = Hodges–Lehmann pseudo-median (median of the Walsh averages of the paired differences) + 95% CI,
#                    computed on the same basis as the p value (normal approximation + continuity correction).
me <- read_csv("data/intermediate/mis_eis_annual.csv", show_col_types = FALSE)
MET <- c("MIS_target", "MIS_goal", "EIS_target", "EIS_goal")
wtest <- function(s, y0) {
  w <- me %>% filter(year %in% c(y0, y0 + 1)) %>%
    select(country, year, value = all_of(s)) %>%
    pivot_wider(names_from = year, values_from = value) %>% drop_na()
  d <- w[[as.character(y0 + 1)]] - w[[as.character(y0)]]
  d <- d[abs(d) > 1e-12]
  wt <- wilcox.test(d, exact = FALSE, correct = TRUE)
  wtc <- suppressWarnings(tryCatch(
    wilcox.test(d, exact = FALSE, correct = TRUE, conf.int = TRUE),
    error = function(e) NULL))
  n <- length(d)
  tibble(series = s, y0 = y0, n = n,
         mean_delta = mean(d), W = unname(wt$statistic), p = wt$p.value,
         r_rank_biserial = 2 * unname(wt$statistic) / (n * (n + 1) / 2) - 1,
         hl_est = if (is.null(wtc)) NA_real_ else unname(wtc$estimate),
         hl_lo = if (is.null(wtc)) NA_real_ else unname(wtc$conf.int[1]),
         hl_hi = if (is.null(wtc)) NA_real_ else unname(wtc$conf.int[2]))
}
t6 <- bind_rows(lapply(MET, function(s) bind_rows(lapply(2011:2023, function(y) wtest(s, y)))))
write_csv(t6, file.path(OUT, "table_6_yoy_paired_tests.csv"))

# ---- Table 7 — MIS/EIS by region and development group ----
cm_raw <- read_csv("data/intermediate/country_meta.csv", show_col_types = FALSE)
cm <- cm_raw %>% filter(is_nation %in% TRUE)
mem <- bind_rows(
  cm %>% filter(!is.na(continent))    %>% transmute(country, group = paste0("continent: ", continent)),
  cm %>% filter(!is.na(subcontinent)) %>% transmute(country, group = paste0("subcontinent: ", subcontinent)),
  cm %>% filter(!is.na(region))       %>% transmute(country, group = paste0("region: ", region)),
  cm %>% filter(!is.na(development))  %>% transmute(country, group = paste0("development: ", development)),
  cm %>% filter(sids %in% TRUE)  %>% transmute(country, group = "SIDS"),
  cm %>% filter(ldc %in% TRUE)   %>% transmute(country, group = "LDC"),
  cm %>% filter(llc %in% TRUE)   %>% transmute(country, group = "LLDC"),
  cm %>% filter(oecd %in% TRUE)  %>% transmute(country, group = "OECD"),
  cm %>% filter(eu %in% TRUE)    %>% transmute(country, group = "EU")) %>% distinct()
is_monotonic <- function(v) {
  d <- diff(v[!is.na(v)]); d <- d[d != 0]
  if (length(d) == 0) return(TRUE)
  all(d > 0) || all(d < 0)
}
YRS <- c(2011, 2014, 2017, 2020, 2024)
t7 <- me %>% inner_join(mem, by = "country") %>% group_by(group) %>%
  group_modify(~ {
    d <- .x
    pick <- function(metric, y) mean(d[[metric]][d$year == y], na.rm = TRUE)
    block <- function(metric) {
      wide <- d %>% select(country, year, value = all_of(metric)) %>%
        pivot_wider(names_from = year, values_from = value)
      ann <- mean(wide[["2024"]] - wide[["2011"]], na.rm = TRUE) / 13
      net <- pick(metric, 2024) - pick(metric, 2011)
      ks <- suppressWarnings(tryCatch(
        ks.test(d[[metric]][d$year == 2011], d[[metric]][d$year == 2024])$p.value,
        error = function(e) NA_real_))
      mono <- wide %>% rowwise() %>%
        mutate(m = is_monotonic(c_across(starts_with("20")))) %>%
        ungroup() %>% summarise(m = mean(m)) %>% pull(m)
      c(setNames(vapply(YRS, function(y) pick(metric, y), numeric(1)), paste0(metric, "_", YRS)),
        setNames(c(ann, net, ks, mono), paste0(metric, c("_delta_ann", "_delta_net", "_ks_p", "_mono"))))
    }
    bind_cols(tibble(n = n_distinct(d$country)),
              as_tibble(as.list(c(block("MIS"), block("EIS")))))
  }) %>% ungroup()
write_csv(t7, file.path(OUT, "table_7_mis_eis_by_group.csv"))

# ---- Table 8 / 9 / 10 — direct renumbered copies of pipeline products ----
file.copy("data/output/table_S5_csc.csv", file.path(OUT, "table_8_csc.csv"), overwrite = TRUE)
file.copy("data/output/table_S1_projection_models.csv", file.path(OUT, "table_10_projection_models.csv"), overwrite = TRUE)
# Table 9: imputation records + "true last observed year" and endpoint-fill columns added in response to reviewer comments (2026-11-10)
#   (product of endfill_sensitivity.R; this script must be run first)
t9 <- read_csv("data/output/table_S6_imputation.csv", show_col_types = FALSE)
lo <- read_csv("data/output/endfill_last_observed.csv", show_col_types = FALSE) %>%
  select(indicator, global_last_obs, median_last_obs, min_last_obs,
         n_end_filled, max_tail_years, share_2021_2024_end_filled)
t9 <- t9 %>% left_join(lo, by = "indicator")
write_csv(t9, file.path(OUT, "table_9_imputation_records.csv"))

message("tables 3-10 written to ", OUT)
