# =====================================================================
# inherited_trend_test.R - test of inherited-trend continuation (reviewer comment: the 2021-2024 data can be read as
#   business-as-usual extrapolations of the Aichi-era trends; the projections should therefore be framed as counterfactual extrapolations)
#   Test: for each indicator / each target-goal layer / the global MIS and EIS series, compare
#     the segment slopes of 2011-2019 (tail of the Aichi period, excluding COVID year 2020) and 2021-2024 (post-GBF window) -
#   if the two slopes are highly correlated and the median ratio is close to 1, then continuation as usual holds.
#   Additionally: the historical length actually used by the projection fits (answering the "13 points on the limit" remark).
# Written to: data/output/inherited_trend_{indicator,layer,summary,fit_window}.csv
# Run from the repository root:  Rscript fig_code/inherited_trend_test.R
# =====================================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
})

put <- function(df, name) {
  write_csv(df, file.path("data/output", paste0(name, ".csv")))
  message("[out] data/output/", name, ".csv  (", nrow(df), " rows)")
}
fit_seg <- function(df, y0, y1, min_pts = 3) {
  s <- df %>% filter(year >= y0, year <= y1)
  if (nrow(s) < min_pts) return(tibble(slope = NA_real_, r2 = NA_real_, n_yr = nrow(s)))
  m <- lm(norm ~ year, data = s)
  tibble(slope = unname(coef(m)[2]), r2 = summary(m)$r.squared, n_yr = nrow(s))
}

# ---- 0. Data ----
norm <- read_csv("data/intermediate/norm_indicator.csv", show_col_types = FALSE)
imp  <- read_csv("data/intermediate/imputation_detail.csv", show_col_types = FALSE)
cov  <- read_csv("data/intermediate/country_coverage.csv", show_col_types = FALSE)
dps  <- read_csv("data/indicator_dpsir.csv", show_col_types = FALSE)
reg  <- read_csv("data/output/regional_indicator_summary.csv", show_col_types = FALSE)
gts  <- read_csv("data/intermediate/goal_target_ts_global.csv", show_col_types = FALSE)
meg  <- read_csv("data/intermediate/mis_eis_annual_global.csv", show_col_types = FALSE)
keep <- cov$country[cov$keep_index]

# ---- 1. Indicator level: global country-mean series (27-29 use the Global regional series) ----
#    both the full panel and the truncated version (LOCF tail removed) are reported, to guard against artefacts of the padding
norm_w <- norm %>% filter(country %in% keep, year %in% 2011:2024) %>%
  left_join(imp %>% select(indicator, country, year, status, method),
            by = c("indicator", "country", "year")) %>%
  mutate(end_fill = !is.na(method) & method == "forward_fill" & status == "imputed")
mean_series <- function(nw) {
  nw %>% group_by(indicator, year) %>%
    summarise(norm = mean(norm, na.rm = TRUE), .groups = "drop")
}
g_reg <- reg %>% filter(indicator %in% c(27, 28, 29), entity == "Global") %>%
  group_by(indicator, year) %>%
  summarise(norm = mean(norm, na.rm = TRUE), .groups = "drop")
g <- bind_rows(mean_series(norm_w), g_reg)
g_t <- bind_rows(mean_series(norm_w %>% filter(!end_fill)), g_reg)
g_hist <- bind_rows(
  norm %>% filter(country %in% keep) %>%
    group_by(indicator, year) %>%
    summarise(norm = mean(norm, na.rm = TRUE), .groups = "drop"),
  g_reg)
hist_stats <- g_hist %>% group_by(indicator) %>%
  summarise(hist_first = min(year), hist_last = max(year),
            hist_years = n_distinct(year), .groups = "drop")

ind <- bind_rows(lapply(sort(unique(g$indicator)), function(i) {
  s <- g %>% filter(indicator == i)
  a <- fit_seg(s, 2011, 2019) %>% rename_with(~ paste0("pre_", .x))
  b <- fit_seg(s, 2021, 2024) %>% rename_with(~ paste0("post_", .x))
  cbind(tibble(indicator = i), a, b, tibble(n_yr_window = nrow(s)))
})) %>%
  left_join(hist_stats, by = "indicator") %>%
  left_join(dps %>% transmute(indicator = indicator_id, indicator_name,
                              dpsir_class = as.character(dpsir_class),
                              side = ifelse(dpsir_class %in% c("Process", "Response"),
                                            "action", "outcome")),
            by = "indicator")
put(ind, "inherited_trend_indicator")

# ---- 2. Target/goal level (paper convention: the historical segment of goal_target_ts_global) ----
gl <- gts %>% filter(value_type == "historical", group_type == "global") %>%
  select(layer, layer_id, year, norm = index)
lay <- bind_rows(lapply(sort(unique(gl$layer_id)), function(t) {
  s <- gl %>% filter(layer_id == t)
  a <- fit_seg(s, 2011, 2019) %>% rename_with(~ paste0("pre_", .x))
  b <- fit_seg(s, 2021, 2024) %>% rename_with(~ paste0("post_", .x))
  cbind(tibble(layer_id = t, layer = s$layer[1]), a, b)
}))
put(lay, "inherited_trend_layer")

# ---- 3. Segment slopes of the global MIS/EIS series ----
me <- meg %>% select(year, MIS_target, MIS_goal, EIS_target, EIS_goal) %>%
  pivot_longer(-year, names_to = "series", values_to = "norm")
me_seg <- bind_rows(lapply(unique(me$series), function(s) {
  d <- me %>% filter(series == s)
  a <- fit_seg(d, 2011, 2019) %>% rename_with(~ paste0("pre_", .x))
  b <- fit_seg(d, 2021, 2024) %>% rename_with(~ paste0("post_", .x))
  cbind(tibble(series = s), a, b)
}))
put(me_seg, "inherited_trend_mis_eis")

# ---- 4. Summary: continuity statistics ----
summarise_pair <- function(df, name) {
  d <- df %>% filter(is.finite(pre_slope), is.finite(post_slope))
  ct <- cor.test(d$pre_slope, d$post_slope, method = "spearman", exact = FALSE)
  wt <- wilcox.test(d$post_slope, d$pre_slope, paired = TRUE, exact = FALSE)
  tibble(
    unit = name, n = nrow(d),
    spearman_rho = unname(ct$estimate), spearman_p = ct$p.value,
    sign_agree = mean(sign(d$pre_slope) == sign(d$post_slope)),
    med_pre_slope = median(d$pre_slope), med_post_slope = median(d$post_slope),
    med_ratio = median(d$post_slope / d$pre_slope, na.rm = TRUE),
    paired_wilcox_p = wt$p.value,
    paired_rank_biserial = 2 * unname(wt$statistic) / (nrow(d) * (nrow(d) + 1) / 2) - 1)
}
sm <- bind_rows(
  summarise_pair(ind, "indicator_all"),
  summarise_pair(ind %>% filter(side == "action"), "indicator_action"),
  summarise_pair(ind %>% filter(side == "outcome"), "indicator_outcome"),
  summarise_pair(lay %>% filter(layer == "target"), "target"),
  summarise_pair(lay %>% filter(layer == "goal"), "goal"))
put(sm, "inherited_trend_summary")

# ---- 4b. Extrapolation deviation test: observed 2021-2024 values vs the trend line extrapolated from 2011-2019 ----
#   dev_z = mean(observed - extrapolated) / SD of the Aichi-period residuals; |dev_z| <= 1 means consistency with the inherited trend
dev_one <- function(df) {
  pre <- df %>% filter(year >= 2011, year <= 2019)
  post <- df %>% filter(year >= 2021, year <= 2024)
  if (nrow(pre) < 4 || nrow(post) < 2) return(NULL)
  m <- lm(norm ~ year, data = pre)
  sd_res <- sd(residuals(m))
  pred <- predict(m, newdata = post)
  dev <- post$norm - pred
  tibble(mean_dev = mean(dev), mean_abs_dev = mean(abs(dev)),
         dev_z = if (sd_res > 0) mean(dev) / sd_res else NA_real_,
         sd_res_pre = sd_res)
}
dev_panel <- function(gg, variant) {
  bind_rows(lapply(sort(unique(gg$indicator)), function(i) {
    d <- dev_one(gg %>% filter(indicator == i))
    if (is.null(d)) return(NULL)
    cbind(tibble(indicator = i, variant = variant), d)
  }))
}
dev_ind <- bind_rows(dev_panel(g, "full"), dev_panel(g_t, "truncated")) %>%
  left_join(dps %>% transmute(indicator = indicator_id, dpsir_class = as.character(dpsir_class),
                              side = ifelse(dpsir_class %in% c("Process", "Response"),
                                            "action", "outcome")), by = "indicator")
put(dev_ind, "inherited_trend_deviation")

dev_lay <- bind_rows(lapply(sort(unique(gl$layer_id)), function(t) {
  d <- dev_one(gl %>% filter(layer_id == t) %>% rename(indicator = layer_id))
  if (is.null(d)) return(NULL)
  cbind(tibble(layer_id = t), d)
}))
put(dev_lay, "inherited_trend_deviation_layer")

dev_sum <- bind_rows(lapply(c("full", "truncated"), function(v) {
  di <- dev_ind %>% filter(variant == v)
  tibble(
    unit = c("indicator_all", "indicator_action", "indicator_outcome"),
    variant = v,
    n = c(nrow(di), sum(di$side == "action"), sum(di$side == "outcome")),
    med_dev_z = c(median(di$dev_z, na.rm = TRUE),
                  median(di$dev_z[di$side == "action"], na.rm = TRUE),
                  median(di$dev_z[di$side == "outcome"], na.rm = TRUE)),
    share_within_1sd = c(mean(abs(di$dev_z) <= 1, na.rm = TRUE),
                         mean(abs(di$dev_z[di$side == "action"]) <= 1, na.rm = TRUE),
                         mean(abs(di$dev_z[di$side == "outcome"]) <= 1, na.rm = TRUE)),
    share_below = c(mean(di$dev_z < 0, na.rm = TRUE),
                    mean(di$dev_z[di$side == "action"] < 0, na.rm = TRUE),
                    mean(di$dev_z[di$side == "outcome"] < 0, na.rm = TRUE)))
}))
dev_sum <- bind_rows(dev_sum, tibble(
  unit = "layer_all", variant = "full", n = nrow(dev_lay),
  med_dev_z = median(dev_lay$dev_z, na.rm = TRUE),
  share_within_1sd = mean(abs(dev_lay$dev_z) <= 1, na.rm = TRUE),
  share_below = mean(dev_lay$dev_z < 0, na.rm = TRUE)))
put(dev_sum, "inherited_trend_deviation_summary")

# ---- 5. Historical length of the projection fits (answering "13 points on the limit") ----
fw <- ind %>% transmute(indicator, indicator_name, dpsir_class,
                        hist_first, hist_last, hist_years,
                        fit_window = hist_years) %>%
  arrange(hist_years)
fw_sum <- tibble(
  stat = c("median_hist_years", "q25_hist_years", "q75_hist_years",
           "min_hist_years", "n_series_lt15y", "n_series_lt10y", "n_total"),
  value = c(median(fw$hist_years), quantile(fw$hist_years, 0.25),
            quantile(fw$hist_years, 0.75), min(fw$hist_years),
            sum(fw$hist_years < 15), sum(fw$hist_years < 10), nrow(fw)))
put(fw, "inherited_trend_fit_window")
put(fw_sum, "inherited_trend_fit_window_summary")

print(as.data.frame(sm))
print(as.data.frame(dev_sum))
print(as.data.frame(fw_sum))
message("INHERITED TREND TEST DONE")
