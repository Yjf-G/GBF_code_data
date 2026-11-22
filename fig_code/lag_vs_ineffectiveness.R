# =====================================================================
# lag_vs_ineffectiveness.R - evidence for discriminating a lag from ineffective action (paper parks)
#   (reviewer comment: the Discussion cites 5-15 year lags, which implies that part of the divergence is expected,
#     yet it also argues that outputs mask failure - one has to explain how a lag is told apart from ineffective action)
#   Two of the three possible criteria can be provided from the existing data:
#   1) Temporal criterion (legacy-stock test): the age structure of the gains of action plus the stock of legacy institutions
#      (coverage and management plans accumulated in 2000-2010, 14-24 years old by 2024) - has ecological response already
#      been observed? Do the segment slopes of long State series before/after 2010 turn upwards?
#   2) Within-target adoption-vs-response criterion: the action side rises within a target while the state/pressure side
#      shows no response (T5/T6/T10 etc.) = the paper-park signature at the indicator level (rearranged from the indicator-type diagnostics)
#   3) Effect-size criterion (effectiveness indicators): missing from this framework - it is discussed only in the text
#      as a limitation and a recommendation, and is not computed
# Written to: data/output/lag_test_*.csv
# Run from the repository root:  Rscript fig_code/lag_vs_ineffectiveness.R
# =====================================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
})

put <- function(df, name) {
  write_csv(df, file.path("data/output", paste0(name, ".csv")))
  message("[out] data/output/", name, ".csv  (", nrow(df), " rows)")
}
fit_win <- function(df, y0, y1) {
  s <- df %>% filter(year >= y0, year <= y1)
  if (nrow(s) < 3) return(tibble(slope = NA_real_, r2 = NA_real_, n_yr = nrow(s)))
  m <- lm(norm ~ year, data = s)
  tibble(slope = unname(coef(m)[2]), r2 = summary(m)$r.squared, n_yr = nrow(s))
}

# ---- 0. Data (same conventions as the indicator-type and endfill analyses) ----
#   g_full = full panel (fixed composition; used for the gain staging and the stock age)
#   g_t    = truncated panel (LOCF tail removed; used for the State segment-upturn test, to avoid spurious upturns)
norm <- read_csv("data/intermediate/norm_indicator.csv", show_col_types = FALSE)
imp  <- read_csv("data/intermediate/imputation_detail.csv", show_col_types = FALSE)
cov  <- read_csv("data/intermediate/country_coverage.csv", show_col_types = FALSE)
dps  <- read_csv("data/indicator_dpsir.csv", show_col_types = FALSE)
reg  <- read_csv("data/output/regional_indicator_summary.csv", show_col_types = FALSE)
keep <- cov$country[cov$keep_index]

norm_w <- norm %>% filter(country %in% keep) %>%
  left_join(imp %>% select(indicator, country, year, status, method),
            by = c("indicator", "country", "year")) %>%
  mutate(end_fill = !is.na(method) & method == "forward_fill" & status == "imputed")
mean_series <- function(nw) {
  nw %>% group_by(indicator, year) %>%
    summarise(norm = mean(norm, na.rm = TRUE), .groups = "drop")
}
reg_full <- reg %>% filter(indicator %in% c(27, 28, 29), entity == "Global") %>%
  group_by(indicator, year) %>%
  summarise(norm = mean(norm, na.rm = TRUE), .groups = "drop")
g_full <- bind_rows(mean_series(norm_w), reg_full)
g <- bind_rows(mean_series(norm_w %>% filter(!end_fill)),
               reg_full %>% filter(year <= 2020))   # truncation (27-29 cut at the last observation, 2020)

# ---- 1. Segment slopes before/after 2010 for long State series (has the legacy stock already produced a response?) ----
STATE <- dps %>% filter(used == "Y", dpsir_class == "State") %>% pull(indicator_id)
seg <- bind_rows(lapply(STATE, function(i) {
  s <- g %>% filter(indicator == i, year >= 1990)
  if (nrow(s) < 8) return(NULL)
  a <- fit_win(s, min(s$year), 2010) %>% rename_with(~ paste0("pre_", .x))
  b <- fit_win(s, 2011, max(s$year)) %>% rename_with(~ paste0("post_", .x))
  cbind(tibble(indicator = i, first = min(s$year), last = max(s$year)), a, b)
})) %>%
  left_join(dps %>% transmute(indicator = indicator_id, indicator_name,
                              dpsir_class = as.character(dpsir_class)), by = "indicator") %>%
  mutate(upturn = post_slope > pre_slope)   # a steeper later segment = first signs of a response
put(seg, "lag_test_state_segments")

# ---- 2. Age structure of the action gains (full panel: fixed composition; LOCF dilutes late gains,
#         so the true share of "young" gains can only be higher - this favours the lag concession, and the text is written accordingly) ----
# 2a. Gains of the Response/Process type series over 2011-2024, decomposed into periods
tact <- g_full %>%
  left_join(dps %>% transmute(indicator = indicator_id,
                              dpsir_class = as.character(dpsir_class)), by = "indicator") %>%
  filter(dpsir_class %in% c("Response", "Process")) %>%
  group_by(dpsir_class, year) %>%
  summarise(norm = mean(norm, na.rm = TRUE), .groups = "drop")
gain_age <- bind_rows(lapply(split(tact, tact$dpsir_class), function(s) {
  v <- function(y) s$norm[s$year == y][1]
  net <- v(2024) - v(2011)
  tibble(dpsir_class = s$dpsir_class[1], net_2011_2024 = net,
         gain_by_2015 = v(2015) - v(2011),
         gain_by_2019 = v(2019) - v(2011),
         share_by_2015 = (v(2015) - v(2011)) / net,
         share_by_2019 = (v(2019) - v(2011)) / net)
}))
put(gain_age, "lag_test_gain_age")

# 2b. The institutional-stock indicators themselves: level in 2010 / level in the final year (age structure of the stock; full panel)
STOCK <- c(4, 19, 20, 21, 22, 23, 24, 36)   # coverage / protection / management plans / certification
stock <- bind_rows(lapply(STOCK, function(i) {
  s <- g_full %>% filter(indicator == i)
  v <- function(y) { val <- s$norm[s$year == y]; if (length(val)) val[1] else NA_real_ }
  vend <- v(max(s$year))
  tibble(indicator = i, end_year = max(s$year), v2000 = v(2000), v2010 = v(2010),
         v_end = vend,
         share_stock_by_2010 = v(2010) / vend,
         gain_2000_2010 = v(2010) - v(2000),
         gain_2010_end = vend - v(2010))
})) %>%
  left_join(dps %>% transmute(indicator = indicator_id, indicator_name), by = "indicator")
put(stock, "lag_test_stock_age")

# ---- 3. Within-target adoption-vs-response (computed from the indicator typing of data/indicator_dpsir.csv) ----
#   per target x type: the mean normalized global series of the target's indicators of that type (2011-2024),
#   and the slope of a linear fit over the same window (same convention as the indicator-type analyses)
i2t <- read_csv("data/indicator_to_target.csv", show_col_types = FALSE)
g_nat <- norm %>%
  group_by(indicator, year) %>%
  summarise(norm = mean(norm, na.rm = TRUE), .groups = "drop")
g_reg <- reg %>%
  filter(indicator %in% c(27, 28, 29), entity == "Global") %>%
  group_by(indicator, year) %>%
  summarise(norm = mean(norm, na.rm = TRUE), .groups = "drop")
g_all <- bind_rows(g_nat, g_reg) %>%
  filter(year >= 2011, year <= 2024)
comp <- i2t %>%
  filter(grepl("^[0-9]+$", target)) %>%
  mutate(target = as.integer(target)) %>%
  inner_join(dps %>% transmute(indicator_id,
                               dpsir_class = factor(dpsir_class,
                                                    levels = c("Process", "Response",
                                                               "Pressure", "State"))),
             by = "indicator_id")
cells <- comp %>%
  inner_join(g_all, by = c("indicator_id" = "indicator")) %>%
  group_by(target, dpsir_class, year) %>%
  summarise(norm = mean(norm, na.rm = TRUE), n_ind = n_distinct(indicator_id),
            .groups = "drop")
within <- bind_rows(lapply(split(cells, list(cells$target, cells$dpsir_class)),
                           function(s) {
  if (nrow(s) == 0) return(NULL)
  m <- lm(norm ~ year, data = s)
  tibble(target = s$target[1], layer_id = paste0("T", s$target[1]),
         dpsir_class = s$dpsir_class[1], n_ind = s$n_ind[1],
         slope = unname(coef(m)[2]), r2 = summary(m)$r.squared)
})) %>%
  filter(dpsir_class %in% c("Response", "Process", "State", "Pressure")) %>%
  select(target = layer_id, dpsir_class, n_ind, slope) %>%
  pivot_wider(names_from = dpsir_class, values_from = c(slope, n_ind)) %>%
  filter((!is.na(slope_Response) | !is.na(slope_Process)) &
           (!is.na(slope_State) | !is.na(slope_Pressure)))
put(within, "lag_test_within_target")

# ---- 4. Summary ----
sm <- tibble(
  stat = c("state_n", "state_upturn_n",
           "gain_share_by_2015_Response", "gain_share_by_2019_Response",
           "gain_share_by_2015_Process", "gain_share_by_2019_Process",
           "stock_share_by_2010_median",
           "within_target_pairs"),
  value = c(nrow(seg), sum(seg$upturn, na.rm = TRUE),
            gain_age$share_by_2015[gain_age$dpsir_class == "Response"],
            gain_age$share_by_2019[gain_age$dpsir_class == "Response"],
            gain_age$share_by_2015[gain_age$dpsir_class == "Process"],
            gain_age$share_by_2019[gain_age$dpsir_class == "Process"],
            median(stock$share_stock_by_2010, na.rm = TRUE),
            nrow(within)))
put(sm, "lag_test_summary")
print(as.data.frame(sm))
print(as.data.frame(seg %>% select(indicator, indicator_name, pre_slope, post_slope, upturn)))
message("LAG VS INEFFECTIVENESS DONE")
