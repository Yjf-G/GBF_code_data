# =====================================================================
# make_fig_data.R - prepares the figure source data for all figures (main Figs. 1-4, Extended Data Fig. 1, Supplementary Figs. 1-13)
# Outputs: data/fig_data/<figure name>/ one folder per figure, with CSVs and a README
# Depends on: the results of run_pipeline() and run_remaining() (data/intermediate, data/output)
# Run from the repository root:  Rscript fig_code/make_fig_data.R
# =====================================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(tibble)
})
options(scipen = 10, dplyr.summarise.inform = FALSE)

# Main pipeline functions (loaded only, not run)
GBF_NO_AUTORUN <- TRUE
source("gbf_compute_pipeline.R")

FIG_DATA <- "data/fig_data"
dir.create(FIG_DATA, recursive = TRUE, showWarnings = FALSE)
new_folder <- function(name) {
  d <- file.path(FIG_DATA, name)
  if (dir.exists(d)) unlink(d, recursive = TRUE)
  dir.create(d, recursive = TRUE)
  d
}
save_csv <- function(df, folder, file) {
  write_csv(df, file.path(folder, file))
  message("  [data] ", file.path(folder, file), " (", nrow(df), " rows)")
}
save_readme <- function(folder, lines) {
  writeLines(lines, file.path(folder, "README.md"))
}

# ---------------- Common inputs ----------------
norm_layer <- read_csv("data/intermediate/norm_layer.csv", show_col_types = FALSE)
pgr        <- read_csv("data/intermediate/pgr.csv", show_col_types = FALSE)
cov        <- read_csv("data/intermediate/country_coverage.csv", show_col_types = FALSE)
traj       <- read_csv("data/intermediate/trajectory_class.csv", show_col_types = FALSE)
csc_cty    <- read_csv("data/intermediate/csc_country.csv", show_col_types = FALSE)
cmeta      <- read_csv("data/intermediate/country_meta.csv", show_col_types = FALSE)
me_glob    <- read_csv("data/intermediate/mis_eis_annual_global.csv", show_col_types = FALSE)
bpc        <- read_csv("data/intermediate/baseline_per_country.csv", show_col_types = FALSE)
syn        <- readRDS("data/intermediate/chk_12_synergy.rds")

c1o <- read_csv("data/output/c1_achievement_ci__overall.csv", show_col_types = FALSE)
c1g <- read_csv("data/output/c1_achievement_ci__per_goal.csv", show_col_types = FALSE)
c1s <- read_csv("data/output/c1_achievement_ci__by_section.csv", show_col_types = FALSE)
b6s <- read_csv("data/output/b6_section_summary.csv", show_col_types = FALSE)
f3p <- read_csv("data/output/fig3_paired.csv", show_col_types = FALSE)
f3k <- read_csv("data/output/fig3_kw.csv", show_col_types = FALSE)
rpath <- read_csv("data/output/b7_regional_paths.csv", show_col_types = FALSE)
a1 <- read_csv("data/output/uncertainty__a1_indicator_subset.csv", show_col_types = FALSE)
a2 <- read_csv("data/output/uncertainty__a2_imputation__gap_length.csv", show_col_types = FALSE)
a2pi <- read_csv("data/output/uncertainty__a2_imputation__per_indicator.csv", show_col_types = FALSE)
a2s <- read_csv("data/output/uncertainty__a2_imputation__sensitivity.csv", show_col_types = FALSE)
a3o <- read_csv("data/output/uncertainty__a3_window_length__optimal.csv", show_col_types = FALSE)
a3c <- read_csv("data/output/uncertainty__a3_window_length__curves.csv", show_col_types = FALSE)
a3x <- read_csv("data/output/uncertainty__a3_window_length__cross_model.csv", show_col_types = FALSE)
a4s <- read_csv("data/output/uncertainty__a4_country_coverage__sigma_by_n.csv", show_col_types = FALSE)
a4f <- read_csv("data/output/uncertainty__a4_country_coverage__sigma_fits.csv", show_col_types = FALSE)
a4p <- read_csv("data/output/uncertainty__a4_country_coverage__prop_stable.csv", show_col_types = FALSE)
a4h <- read_csv("data/output/uncertainty__a4_country_coverage__headline_boot.csv", show_col_types = FALSE)
a4hs <- read_csv("data/output/uncertainty__a4_country_coverage__headline_summary.csv", show_col_types = FALSE)
a5 <- read_csv("data/output/uncertainty__a5_perturbation.csv", show_col_types = FALSE)
c3 <- read_csv("data/output/uncertainty__c3_weight_sensitivity.csv", show_col_types = FALSE)
s2kw <- read_csv("data/output/table_S2_kw.csv", show_col_types = FALSE)
s2p  <- read_csv("data/output/table_S2_posthoc.csv", show_col_types = FALSE)
s1   <- read_csv("data/output/table_S1_projection_models.csv", show_col_types = FALSE)

SEC_LEVELS <- c("reducing_threats", "meeting_needs", "tools_solutions")
SEC_LABELS <- c(reducing_threats = "Reducing threats",
                meeting_needs = "Meeting people's needs",
                tools_solutions = "Tools and solutions")
SEC_COLORS <- c(reducing_threats = "#0072B2", meeting_needs = "#009E73",
                tools_solutions = "#E69F00")

sec_map <- read_csv("data/target_section.csv", show_col_types = FALSE) %>%
  transmute(layer_id = paste0("T", id), section) %>%
  bind_rows(tibble(layer_id = c("A", "B", "C", "D"), section = "Goals"))

# Baseline fit of the global mean series over 2011-2020 (identical convention to the main text / Table 3)
global_fit <- function(df) {
  df %>% group_by(layer, layer_id, year) %>%
    summarise(v = mean(index, na.rm = TRUE), .groups = "drop") %>%
    filter(year %in% 2011:2020, is.finite(v)) %>%
    group_by(layer, layer_id) %>%
    group_modify(~ {
      fit <- lm(v ~ year, data = .x)
      s <- summary(fit)$coefficients
      tibble(slope = s[2, 1], se = s[2, 2],
             ci_lo = s[2, 1] - 1.96 * s[2, 2], ci_hi = s[2, 1] + 1.96 * s[2, 2],
             r2 = summary(fit)$r.squared, p_value = s[2, 4])
    }) %>% ungroup()
}

# =====================================================================
# Fig 1 - baseline trends: (a) slopes of the 23 targets with CI (point size R2, colour section, filled if significant);
#         (b) global baseline trend lines for A-D;  (c) slope x R2 quadrant
# =====================================================================
message("Fig 1 ...")
d1 <- new_folder("fig1")
fits <- global_fit(norm_layer)
f1a <- fits %>% filter(layer == "target") %>%
  left_join(sec_map, by = "layer_id") %>%
  mutate(label = layer_id,
         significant = slope > 0 & p_value < 0.05 | slope < 0 & p_value < 0.05,
         direction = ifelse(slope >= 0, "up", "down")) %>%
  select(layer_id, label, section, slope, se, ci_lo, ci_hi, r2, p_value,
         significant, direction)
save_csv(f1a, d1, "fig1a_target_slopes.csv")

f1b <- norm_layer %>% filter(layer == "goal", year %in% 2011:2020, !is.na(index)) %>%
  group_by(layer_id, year) %>%
  summarise(index = mean(index, na.rm = TRUE), .groups = "drop") %>%
  left_join(fits %>% filter(layer == "goal") %>% select(layer_id, slope, p_value),
            by = "layer_id") %>%
  mutate(slope_sign = ifelse(slope > 0, "+", "−"))
save_csv(f1b, d1, "fig1b_goal_trends.csv")

f1c <- bind_rows(
  f1a %>% select(layer_id, section, slope, r2, p_value, significant),
  fits %>% filter(layer == "goal") %>%
    left_join(sec_map, by = "layer_id") %>%
    mutate(significant = p_value < 0.05) %>%
    select(layer_id, section, slope, r2, p_value, significant)
)
save_csv(f1c, d1, "fig1c_quadrant.csv")
save_readme(d1, c(
  "# Fig 1 figure source data",
  "",
  "- fig1a_target_slopes.csv: linear fit of the global mean series over 2011-2020 for the 23 targets (slope with 95% CI, R2, p); significant = the slope is significant (p<0.05); section = the GBF section",
  "- fig1b_goal_trends.csv: annual index of the global means for goals A-D (2011-2020) + the sign of the slope (+/-)",
  "- fig1c_quadrant.csv: slope x R2 quadrant data (23 targets + goals A-D, section=Goals)",
  "",
  "Definition: equal-weight mean of the country indices (Formula 2) -> annual global mean -> lm(index ~ year)."))

# =====================================================================
# Fig 2 - trajectories and achievement: (a,b) trend x achievement classification; (c) achievement rate (target/section);
#         (d,e) PGR_req per country (2030/2050); (f) distribution of PGR_req
# =====================================================================
message("Fig 2 ...")
d2 <- new_folder("fig2")
f2c <- bind_rows(
  c1o %>% mutate(group_type = "overall", group = toupper(layer)),
  c1g %>% mutate(group_type = "goal", group = layer_id),
  c1s %>% mutate(group_type = "section", group = section) %>%
    left_join(b6s %>% select(section, ci95_lo, ci95_hi), by = "section")
) %>% select(group_type, group, AR, ci50_lo, ci50_hi, ci95_lo, ci95_hi, n)
save_csv(f2c, d2, "fig2c_achievement.csv")

iso <- cmeta %>% select(country, iso3) %>% distinct()
f2de <- pgr %>%
  filter(!is.na(pgr_sim), !is.na(pgr_req)) %>%
  mutate(pgr = pgr_req - pgr_sim) %>%   # Formula 12: PGR gap = PGR_req - PGR_sim
  group_by(layer, country) %>%
  summarise(pgr_mean = mean(pgr), .groups = "drop") %>%
  # A.2 coverage filter: countries with >20% of the indicators missing are excluded (keep_index, same convention as MIS/EIS/ranking)
  semi_join(cov %>% filter(keep_index) %>% select(country), by = "country") %>%
  pivot_wider(names_from = layer, values_from = pgr_mean,
              names_prefix = "pgr_") %>%
  left_join(iso, by = "country")
save_csv(f2de, d2, "fig2de_pgr_by_country.csv")

mk_hist <- function(x, layer) {
  x <- x[is.finite(x)]
  t <- sign(x) * log10(1 + abs(x))       # signed pseudo-log: accommodates negative values (surplus)
  tr <- range(t)
  if (diff(tr) == 0) tr <- tr + c(-0.3, 0.3)
  breaks <- seq(tr[1], tr[2], length.out = 11)   # 10 equal-width bins (in pseudo-log space)
  h <- hist(t, breaks = breaks, plot = FALSE)
  tibble(layer = layer,
         bin_lo = h$breaks[-length(h$breaks)],
         bin_hi = h$breaks[-1],
         share = h$counts / sum(h$counts),
         mean_val = mean(x, na.rm = TRUE),
         median_val = median(x, na.rm = TRUE))
}
f2f <- bind_rows(
  mk_hist(f2de$pgr_target[!is.na(f2de$pgr_target)], "target"),
  mk_hist(f2de$pgr_goal[!is.na(f2de$pgr_goal)], "goal"))
save_csv(f2f, d2, "fig2f_pgr_distribution.csv")
save_readme(d2, c(
  "# Fig 2 figure source data",
  "",
  "- fig2c_achievement.csv: achievement rate AR with 50%/95% bootstrap intervals (overall, goals A-D, the three sections)",
  "- fig2de_pgr_by_country.csv: mean growth-rate shortfall PGR = PGR_req - PGR_sim (Formula 12) at the target level (2030) and the goal level (2050) for the countries passing the coverage rule (>=80% of the indicators, A.2), with iso3 for the map join",
  "- fig2f_pgr_distribution.csv: distribution of the mean PGR across countries (shares of 10 equal-width bins on a signed pseudo-log scale, mean, median)"))

# =====================================================================
# Fig 3 - MIS/EIS by year: (a) global dual-axis series; (b) regional trajectories; (c) year-on-year paired tests
# =====================================================================
message("Fig 3 ...")
d3 <- new_folder("fig3")
me_h <- read_csv("data/intermediate/mis_eis_annual.csv", show_col_types = FALSE)
METRICS <- c("MIS_target", "MIS_goal", "EIS_target", "EIS_goal")

# ---- Region / development-group membership (convention of the caption) ----
REG_LV <- c("Africa", "Asia", "Europe", "Oceania", "LAC", "North A.")
GRP_LV <- c("DD", "DL", "LDC", "LLDC", "OECD", "SIDS")
mem <- cmeta %>%
  mutate(
    region6 = case_when(
      subcontinent == "Latin America and the Caribbean" ~ "LAC",
      subcontinent == "Northern America" ~ "North A.",
      continent %in% c("Africa", "Asia", "Europe", "Oceania") ~ continent,
      TRUE ~ NA_character_),
    g_DD = development == "developed",
    g_DL = development == "developing",
    g_LDC = ldc %in% TRUE, g_LLDC = llc %in% TRUE,
    g_OECD = oecd %in% TRUE, g_SIDS = sids %in% TRUE) %>%
  select(country, region6, starts_with("g_"))

val_long <- me_h %>% filter(year %in% c(2011, 2024)) %>%
  select(country, year, all_of(METRICS)) %>%
  pivot_longer(all_of(METRICS), names_to = "metric", values_to = "value") %>%
  left_join(mem, by = "country")

# ---- Dunn pairwise comparisons (rank-sum method, Bonferroni correction) ----
dunn_pairs <- function(x, g, levels) {
  out <- list(); k <- 0
  for (i in seq_len(length(levels) - 1)) {
    for (j in (i + 1):length(levels)) {
      k <- k + 1
      xi <- x[g == levels[i]]; xj <- x[g == levels[j]]
      ok <- !is.na(xi) & !is.na(xj)
      ni <- sum(!is.na(xi)); nj <- sum(!is.na(xj))
      if (ni < 2 || nj < 2) {
        p <- NA_real_
      } else {
        xi <- xi[!is.na(xi)]; xj <- xj[!is.na(xj)]
        r <- rank(c(xi, xj))
        z <- (mean(r[seq_along(xi)]) - mean(r[-seq_along(xi)])) /
          sqrt((length(r) * (length(r) + 1) / 12) * (1 / ni + 1 / nj))
        p <- 2 * pnorm(-abs(z))
      }
      out[[k]] <- tibble(g1 = levels[i], g2 = levels[j], p_raw = p)
    }
  }
  bind_rows(out) %>% mutate(p_bonf = pmin(1, p_raw * n()))
}

mk_triangle <- function(metric, colname, levels) {
  bind_rows(map_dfr(c(2011, 2024), function(yr) {
    d <- val_long %>% filter(metric == !!metric, year == yr,
                             .data[[colname]] %in% levels)
    x <- d$value; g <- factor(d[[colname]], levels = levels)
    Dunn <- dunn_pairs(x, g, levels)
    means <- tapply(x, g, mean, na.rm = TRUE)
    crossing(g1 = levels, g2 = levels) %>%
      filter(g1 != g2) %>%
      mutate(year = yr,
             diff = as.numeric(means[g2]) - as.numeric(means[g1])) %>%
      left_join(Dunn %>% select(g1, g2, p_bonf), by = c("g1", "g2")) %>%
      mutate(metric = metric)
  }))
}

tri_region <- bind_rows(map_dfr(METRICS, ~ mk_triangle(.x, "region6", REG_LV)))
save_csv(tri_region, d3, "fig3_tri_regions.csv")

tri_group <- bind_rows(map_dfr(METRICS, function(m) {
  # Development groups: each group has its own membership set (the memberships overlap)
  d <- val_long %>% filter(metric == m, year %in% c(2011, 2024))
  gdf <- bind_rows(map_dfr(GRP_LV, function(gl) {
    d %>% filter(.data[[paste0("g_", gl)]] %in% TRUE) %>%
      mutate(grp = gl) %>% select(country, year, value, grp, metric)
  }))
  bind_rows(map_dfr(c(2011, 2024), function(yr) {
    dd <- gdf %>% filter(year == yr)
    x <- dd$value; g <- factor(dd$grp, levels = GRP_LV)
    Dunn <- dunn_pairs(x, g, GRP_LV)
    means <- tapply(x, g, mean, na.rm = TRUE)
    crossing(g1 = GRP_LV, g2 = GRP_LV) %>% filter(g1 != g2) %>%
      mutate(year = yr, diff = as.numeric(means[g2]) - as.numeric(means[g1])) %>%
      left_join(Dunn %>% select(g1, g2, p_bonf), by = c("g1", "g2")) %>%
      mutate(metric = m)
  }))
}))
save_csv(tri_group, d3, "fig3_tri_groups.csv")

mk_ks_region <- function(metric, colname, levels) {
  map_dfr(levels, function(gl) {
    d <- val_long %>% filter(metric == !!metric, .data[[colname]] == gl,
                             year %in% c(2011, 2024))
    v11 <- d$value[d$year == 2011]; v24 <- d$value[d$year == 2024]
    v11 <- v11[!is.na(v11)]; v24 <- v24[!is.na(v24)]
    tibble(metric = metric, group = gl, group_type = "region",
           delta_mean = mean(v24) - mean(v11),
           ks_p = if (length(v11) >= 3 && length(v24) >= 3)
             suppressWarnings(ks.test(v11, v24)$p.value) else NA_real_,
           n = length(v24))
  })
}
ks_region <- map_dfr(METRICS, ~ mk_ks_region(.x, "region6", REG_LV))
ks_grp <- bind_rows(map_dfr(METRICS, function(m) {
  d <- val_long %>% filter(metric == m, year %in% c(2011, 2024))
  map_dfr(GRP_LV, function(gl) {
    dd <- d %>% filter(.data[[paste0("g_", gl)]] %in% TRUE)
    v11 <- dd$value[dd$year == 2011]; v24 <- dd$value[dd$year == 2024]
    v11 <- v11[!is.na(v11)]; v24 <- v24[!is.na(v24)]
    tibble(metric = m, group = gl, group_type = "group",
           delta_mean = mean(v24) - mean(v11),
           ks_p = if (length(v11) >= 3 && length(v24) >= 3)
             suppressWarnings(ks.test(v11, v24)$p.value) else NA_real_,
           n = length(v24))
  })
}))
ks_all <- bind_rows(ks_region, ks_grp)
save_csv(ks_all, d3, "fig3_ks.csv")

map_delta <- val_long %>%
  filter(year %in% c(2011, 2024)) %>%
  select(country, year, metric, value) %>%
  pivot_wider(names_from = year, values_from = value, names_prefix = "v") %>%
  mutate(delta = v2024 - v2011) %>%
  left_join(cmeta %>% select(country, iso3) %>% distinct(), by = "country") %>%
  select(metric, country, iso3, delta)
save_csv(map_delta, d3, "fig3_map_delta.csv")

save_csv(f3p, d3, "fig3_paired_year_on_year.csv")
save_readme(d3, c(
  "# Fig 3 figure source data (Spatiotemporal dynamics of MIS and EIS, 2011-2024)",
  "",
  "- fig3_tri_regions.csv: triangular matrix of the pairwise differences between regions (6 regions: Africa/Asia/Europe/Oceania/LAC/North A.);",
  "  diff = column-group mean - row-group mean (X-Y, identical to the values in the figure; 2011 in blue in the upper triangle, 2024 in red in the lower triangle); p_bonf = Dunn post-hoc test (Bonferroni, coloured by significance in the figure)",
  "- fig3_tri_groups.csv: triangular matrix of the same format for the development-status groups (DD/DL/LDC/LLDC/OECD/SIDS)",
  "- fig3_ks.csv: KS test per group (2011 vs 2024 country distributions) and delta_mean = 2024 mean - 2011 mean (the \"ks\" column; the fill colour of the figure gives the significance of the KS p)",
  "- fig3_map_delta.csv: change of each country in 2024 relative to 2011 (blue->red = negative->positive)",
  "- fig3_paired_year_on_year.csv: year-on-year paired tests (paired Wilcoxon for 4 metrics x 13 year pairs: metric/year pair/mean_diff/p_paired)",
  "",
  "Definition: country-level MIS/EIS (both the goal and the target resolution); group means = equal-weight mean over the countries in the group."))

# =====================================================================
# Fig 4 - collective action: (a) GA ranges per indicator; (b) share of years with a positive GA; (c) CAC map per country
# =====================================================================
message("Fig 4 ...")
d4 <- new_folder("fig4")
syd <- syn$syd %>% mutate(GA = S_t * sgn_S)
f4a <- syd %>% group_by(indicator) %>%
  summarise(ga_min = min(GA), ga_max = max(GA), ga_mean = mean(GA),
            .groups = "drop")
save_csv(f4a, d4, "fig4a_ga_ranges.csv")

f4b <- syd %>% group_by(indicator) %>%
  summarise(share_positive = mean(GA > 0), n_years = n(),
            .groups = "drop") %>%
  arrange(share_positive)
save_csv(f4b, d4, "fig4b_share_positive_years.csv")

f4c <- csc_cty %>%
  left_join(cmeta %>% select(country, iso3) %>% distinct(), by = "country") %>%
  mutate(key_country = country %in% c("France", "Brazil", "Cuba", "China"))
save_csv(f4c, d4, "fig4c_csc_by_country.csv")
save_readme(d4, c(
  "# Fig 4 figure source data",
  "",
  "- fig4a_ga_ranges.csv: annual range of the aligned synergy effect SyE*SyD = S_t*sgn(S_t) for each indicator (min/max/mean, Formulas 7-8)",
  "- fig4b_share_positive_years.csv: share of the years with SyE*SyD>0 for each indicator (with the 50% threshold line)",
  "- fig4c_csc_by_country.csv: mean country synergy contribution CSC (Formula 9, 130 countries passing the coverage rule), with iso3 and a key_country flag column (country names are not annotated in the figure)"))

# =====================================================================
# Supplementary Fig. 1 - geographic heterogeneity of the baseline parameters (section/subcontinent/development status/SIDS/LDC)
# =====================================================================
message("Supplementary Fig. 1 ...")
de1 <- new_folder("sup_fig1")
geo <- cmeta %>% select(country, continent, subcontinent, region,
                        development, sids, ldc, llc, oecd)
e1 <- bpc %>%
  left_join(geo, by = "country") %>%
  pivot_longer(c(mean_v, sd_v, range_v, slope),
               names_to = "param", values_to = "value")
save_csv(e1, de1, "supfig1_baseline_params_by_country.csv")
# KW test tables: continent (s2kw) + the 6 status groups (DD/DL/LDC/LLDC/SIDS/OECD, a country can belong to several groups), omnibus
PARAMS_E1 <- c("mean_v", "sd_v", "slope", "range_v")
STATUS_LV <- c("DD", "DL", "LDC", "LLDC", "SIDS", "OECD")
e1s <- e1 %>%
  filter(!is.na(value)) %>%
  transmute(layer, layer_id, country, param, value,
            DD   = development == "developed",
            DL   = development == "developing",
            LDC  = ldc, LLDC = llc, SIDS = sids, OECD = oecd) %>%
  pivot_longer(all_of(STATUS_LV), names_to = "grp", values_to = "member") %>%
  filter(member) %>% select(-member)
safe_kw_e1 <- function(v, g) {
  ok <- !is.na(v)
  if (sum(ok) < 5 || n_distinct(g[ok]) < 2) return(list(H = NA_real_, p = NA_real_))
  kt <- tryCatch(kruskal.test(v[ok] ~ factor(g[ok])), error = function(e) NULL)
  if (is.null(kt)) list(H = NA_real_, p = NA_real_)
  else list(H = unname(kt$statistic), p = kt$p.value)
}
e1s_kw <- e1s %>%
  group_by(layer, layer_id, param) %>%
  group_modify(~ {
    kw <- safe_kw_e1(.x$value, .x$grp)
    tibble(kw_H = kw$H, kw_p = kw$p,
           n = n_distinct(.x$country[!is.na(.x$value)]))
  }) %>%
  ungroup() %>% mutate(grouping = "status")
e1kw <- bind_rows(
  s2kw %>% filter(param %in% PARAMS_E1) %>% mutate(grouping = "continent") %>%
    select(layer, layer_id, grouping, param, kw_H, kw_p, n),
  e1s_kw %>% select(layer, layer_id, grouping, param, kw_H, kw_p, n)
) %>% arrange(layer, layer_id, grouping, param)
save_csv(e1kw, de1, "supfig1_kw_tests.csv")
# Pairwise differences between continents at the goal level (post-hoc, averaged over A-D; kept as source data)
e1p <- s2p %>% filter(layer == "goal", param %in% PARAMS_E1) %>%
  group_by(param, cont1, cont2) %>%
  summarise(mean_diff = mean(mean_diff, na.rm = TRUE),
            p_bonf_median = median(u_p_bonf, na.rm = TRUE),
            n_sig_of_4 = sum(u_p_bonf < 0.05, na.rm = TRUE),
            .groups = "drop")
save_csv(e1p, de1, "supfig1_posthoc_goals_continent.csv")
# Pairwise tests between continents (significance matrices of panels c/d; B.2 post-hoc, values only for omnibus-significant cases)
CONT_LV_E1 <- c("Africa", "Americas", "Asia", "Europe", "Oceania")
ids_e1 <- bind_rows(tibble(layer = "target", layer_id = paste0("T", 1:23)),
                    tibble(layer = "goal", layer_id = c("A", "B", "C", "D")))
pairs_c <- as.data.frame(t(combn(CONT_LV_E1, 2)), stringsAsFactors = FALSE)
names(pairs_c) <- c("grp1", "grp2")
pw_c <- expand_grid(ii = seq_len(nrow(ids_e1)),
                    param = PARAMS_E1,
                    pi = seq_len(nrow(pairs_c))) %>%
  mutate(layer = ids_e1$layer[ii], layer_id = ids_e1$layer_id[ii],
         grp1 = pairs_c$grp1[pi], grp2 = pairs_c$grp2[pi]) %>%
  select(layer, layer_id, param, grp1, grp2) %>%
  left_join(s2p %>%
              select(layer, layer_id, param, cont1, cont2, u_p_raw, u_p_bonf) %>%
              rename(grp1 = cont1, grp2 = cont2,
                     p_raw = u_p_raw, p_bonf = u_p_bonf),
            by = c("layer", "layer_id", "param", "grp1", "grp2")) %>%
  left_join(e1kw %>% filter(grouping == "continent") %>%
              select(layer, layer_id, param, omnibus_p = kw_p),
            by = c("layer", "layer_id", "param")) %>%
  mutate(sig = !is.na(p_bonf) & p_bonf < 0.05) %>%
  arrange(layer, layer_id, param, grp1, grp2)
save_csv(pw_c, de1, "supfig1_continent_pairwise_tests.csv")
# Pairwise standardised gaps (10 continent pairs + 15 status pairs; heatmaps of caption panels c/d)
#   for each (layer_id, param): standardise the country values z -> group means -> pairwise |dz| -> average over the targets/goals;
#   n_sig: for continent pairs = the number of significant Bonferroni x 10 post-hoc results (B.2, only for omnibus-significant cases);
#           for status pairs = the number of significant MW Bonferroni x 15 results (again only for omnibus-significant cases)
e1z <- e1 %>%
  filter(!is.na(value), !is.na(continent)) %>%
  group_by(layer, layer_id, param) %>%
  mutate(z = as.numeric(scale(value))) %>%
  ungroup()
gaps_c <- e1z %>%
  group_by(layer, layer_id, param, continent) %>%
  summarise(mz = mean(z), .groups = "drop_last") %>%
  group_modify(~ {
    lv <- sort(unique(.x$continent))
    prs <- combn(lv, 2, simplify = FALSE)
    map_dfr(prs, function(pr) {
      tibble(grp1 = pr[1], grp2 = pr[2],
             gap_z = abs(.x$mz[.x$continent == pr[1]] -
                           .x$mz[.x$continent == pr[2]]))
    })
  }) %>%
  group_by(layer, param, grp1, grp2) %>%
  summarise(mean_gap_z = mean(gap_z), n_units = n(), .groups = "drop")
nsig_c <- s2p %>%
  filter(param %in% PARAMS_E1) %>%
  group_by(layer, param, cont1, cont2) %>%
  summarise(n_sig = sum(u_p_bonf < 0.05, na.rm = TRUE), .groups = "drop") %>%
  rename(grp1 = cont1, grp2 = cont2)
gaps_c <- gaps_c %>%
  left_join(nsig_c, by = c("layer", "param", "grp1", "grp2")) %>%
  mutate(n_sig = ifelse(is.na(n_sig), 0L, as.integer(n_sig)),
         pair_type = "continent")
e1s_z <- e1z %>%
  transmute(layer, layer_id, param, country, z,
            DD   = development == "developed",
            DL   = development == "developing",
            LDC  = ldc, LLDC = llc, SIDS = sids, OECD = oecd) %>%
  pivot_longer(all_of(STATUS_LV), names_to = "grp", values_to = "member") %>%
  filter(member) %>% select(-member)
mzs <- e1s_z %>%
  group_by(layer, layer_id, param, grp) %>%
  summarise(mz = mean(z), .groups = "drop")
# Pairwise tests between status groups (significance matrices of panels e/f; the n_sig of gaps_s also comes from here)
#   MW pairwise with Bonferroni x 15; only omnibus-significant cases count as tested (B.2)
pw <- map_dfr(combn(STATUS_LV, 2, simplify = FALSE), function(pr) {
  e1s_z %>% filter(grp %in% pr) %>%
    group_by(layer, layer_id, param) %>%
    group_modify(~ {
      v <- .x$z; g <- factor(.x$grp, levels = pr)
      if (n_distinct(g) < 2) return(tibble(p_raw = NA_real_))
      tibble(p_raw = tryCatch(suppressWarnings(wilcox.test(v ~ g)$p.value),
                              error = function(e) NA_real_))
    }) %>%
    ungroup() %>%
    mutate(grp1 = pr[1], grp2 = pr[2])
}) %>%
  left_join(e1s_kw %>% select(layer, layer_id, param, omnibus_p = kw_p),
            by = c("layer", "layer_id", "param")) %>%
  mutate(p_bonf = pmin(1, p_raw * 15),
         sig = !is.na(omnibus_p) & omnibus_p < 0.05 &
               !is.na(p_bonf) & p_bonf < 0.05) %>%
  select(layer, layer_id, param, grp1, grp2, p_raw, p_bonf, omnibus_p, sig) %>%
  arrange(layer, layer_id, param, grp1, grp2)
save_csv(pw, de1, "supfig1_status_pairwise_tests.csv")
gaps_s <- map_dfr(combn(STATUS_LV, 2, simplify = FALSE), function(pr) {
  mzs %>%
    filter(grp %in% pr) %>%
    mutate(side = ifelse(grp == pr[1], "g1", "g2")) %>%
    select(layer, layer_id, param, side, mz) %>%
    pivot_wider(names_from = side, values_from = mz) %>%
    filter(!is.na(g1), !is.na(g2)) %>%
    mutate(gap_z = abs(g1 - g2), grp1 = pr[1], grp2 = pr[2]) %>%
    left_join(pw %>% select(layer, layer_id, param, grp1, grp2, sig),
              by = c("layer", "layer_id", "param", "grp1", "grp2")) %>%
    group_by(layer, param) %>%
    summarise(mean_gap_z = mean(gap_z), n_units = n(),
              n_sig = sum(sig, na.rm = TRUE), .groups = "drop") %>%
    mutate(pair_type = "status", grp1 = pr[1], grp2 = pr[2])
})
e1g <- bind_rows(gaps_c %>% select(layer, pair_type, param, grp1, grp2,
                                   mean_gap_z, n_units, n_sig),
                 gaps_s %>% select(layer, pair_type, param, grp1, grp2,
                                   mean_gap_z, n_units, n_sig)) %>%
  arrange(layer, pair_type, param, grp1, grp2)
save_csv(e1g, de1, "supfig1_pair_gaps.csv")
save_readme(de1, c(
  "# Supplementary Fig. 1 figure source data",
  "",
  "- supfig1_baseline_params_by_country.csv: the four baseline parameters per country over 2011-2020 (mean/SD/slope/range, B.3) x target/goal, with the continent/subcontinent/region/development/SIDS/LDC/LLDC/OECD group assignments",
  "- supfig1_kw_tests.csv: grouped Kruskal-Wallis tests (kw_H, kw_p; grouping = continent (5 continents) / status (the six groups DD/DL/LDC/LLDC/SIDS/OECD, a country can belong to several groups); significance heatmaps of caption panels a/b and e/f)",
  "- supfig1_posthoc_goals_continent.csv: pairwise differences between continents at the goal level (mean_diff averaged over A-D, median Bonferroni p, number of significant differences among the 4 goals)",
  "- supfig1_pair_gaps.csv: pairwise standardised gaps |dz| (10 continent pairs + 15 status pairs) averaged over the targets/goals + the significant counts (auxiliary source data, not used directly in the current figure)",
  "- supfig1_continent_pairwise_tests.csv: MW tests for the 10 continent pairs, per (target/goal x parameter) (p_raw, p_bonf = u_p_bonf, omnibus_p, sig; significance matrices of caption panels c/d)",
  "- supfig1_status_pairwise_tests.csv: MW tests for the 15 status pairs, per (target/goal x parameter) (p_raw, p_bonf = p_raw*15, omnibus_p, sig; significance matrices of caption panels e/f)"))

# =====================================================================
# Supplementary Figs. 8, 10-11 - content unique to the former Supplementary Fig. 12 + endpoint robustness (SI B.1.4 / B.2.2 / B.2.3)
# =====================================================================
message("Supplementary Figs. 8, 10-11 ...")
ds6 <- new_folder("sup_fig8")
save_csv(a4s, ds6, "supfig8_sigma_curves.csv")
save_csv(a4p, ds6, "supfig8_prop_stable.csv")
save_csv(a4h, ds6, "supfig8_net_change_bootstrap.csv")
save_csv(a4hs, ds6, "supfig8_net_change_summary.csv")
save_readme(ds6, c(
  "# Supplementary Fig. 8 figure source data",
  "",
  "- supfig8_sigma_curves.csv: per-indicator subsampling SD curve sigma_i(n) (Algorithm 4, phase 1)",
  "- supfig8_prop_stable.csv: prop_stable, the share of years stable within +/-5% in the full-size bootstrap, and the 30-country comparison prop_stable_n30",
  "- supfig8_net_change_bootstrap.csv: the country bootstrap draws of the net MIS change (goal level; column net_change)",
  "- supfig8_net_change_summary.csv: bootstrap summary (full / boot_mean / boot_sd / n_countries / B)"))
ds7 <- new_folder("sup_fig10")
save_csv(c3, ds7, "supfig10_weight_sensitivity.csv")
save_readme(ds7, c(
  "# Supplementary Fig. 10 figure source data",
  "",
  "- supfig10_weight_sensitivity.csv: change in MIS between 2011 and 2024 under each weighting scheme (SI B.2.2) (level = the target/goal two",
  "  levels keep_index with the 130-country convention; column net_change = the MIS change of the paper); dirichlet_sd is the standard deviation of 1,000 Dirichlet draws"))
ds8 <- new_folder("sup_fig11")
save_csv(read_csv("data/output/c5_leave_one_year_out.csv", show_col_types = FALSE),
         ds8, "supfig11_leave_one_year_out.csv")
save_csv(read_csv("data/output/c5_extrapolation_test.csv", show_col_types = FALSE),
         ds8, "supfig11_extrapolation_test.csv")
save_csv(read_csv("data/output/c5_extrapolation_series.csv", show_col_types = FALSE),
         ds8, "supfig11_extrapolation_series.csv")
save_readme(ds8, c(
  "# Supplementary Fig. 11 figure source data",
  "",
  "- supfig11_leave_one_year_out.csv: change in the net MIS/EIS after leaving out one year (SI B.2.3)",
  "- supfig11_extrapolation_test.csv: mean absolute deviation of the linear trend over 2011-2023 extrapolated to 2024",
  "- supfig11_extrapolation_series.csv: absolute extrapolation deviation per country-target series (used for the distribution plot)"))

# =====================================================================
# Supplementary Fig. 12 - representative model fits (one example per model form)
# =====================================================================
message("Supplementary Fig. 12 ...")
de4 <- new_folder("sup_fig12")
fit_series <- function(id, name, engine = c("pma", "dlm")) {
  engine <- match.arg(engine)
  d <- read_csv("data/intermediate/norm_indicator.csv", show_col_types = FALSE) %>%
    filter(indicator == id, !is.na(norm)) %>%
    group_by(year) %>%
    summarise(v = mean(norm, na.rm = TRUE), .groups = "drop") %>%
    arrange(year)
  snr <- calculate_snr(d$v)
  out <- NULL
  if (engine == "dlm") {
    mod <- fit_dlm(d$v)
    if (is.null(mod)) return(NULL)
    h <- max(TARGET_YEAR_TARGET - max(d$year), 1)
    fut <- seq.int(max(d$year) + 1, TARGET_YEAR_TARGET)
    if (inherits(mod, "builtin_dlm")) {
      y <- mod$y; V <- mod$V; W <- mod$W; n <- length(y)
      m <- c(y[1], 0); C <- diag(c(V * 10, V * 10))
      G <- matrix(c(1, 0, 1, 1), 2, 2); Wm <- diag(W)
      f_in <- rep(NA_real_, n); R_in <- rep(NA_real_, n)
      for (t in seq_len(n)) {
        a <- G %*% m; R <- G %*% C %*% t(G) + Wm
        f0 <- as.numeric(a[1]); Q <- R[1, 1] + V
        e <- y[t] - f0; K <- R[, 1] / Q
        m <- as.numeric(a + K * e)
        C <- R - K %*% t(R[, 1, drop = FALSE]); C <- (C + t(C)) / 2
        f_in[t] <- m[1]; R_in[t] <- max(C[1, 1], 0)   # filtered level + state variance
      }
      fc <- builtin_dlm_forecast(mod, h)
      f_all <- c(f_in, fc$f_path); R_all <- c(R_in, fc$R_path)
    } else {
      filt <- dlm::dlmFilter(d$v, mod)
      f_in <- as.numeric(filt$m[-1, 1])          # filtered level
      C_list <- dlm::dlmSvd2var(filt$U.C, filt$D.C)
      R_in <- vapply(C_list[-1], function(C) max(C[1, 1], 0), numeric(1))
      fc <- dlm::dlmForecast(filt, nAhead = h)
      f_path <- vapply(fc$f, function(z) as.numeric(z)[1], numeric(1))
      R_path <- vapply(fc$R, function(R) R[1, 1] + mod$V, numeric(1))
      f_all <- c(f_in, f_path); R_all <- c(R_in, R_path)
    }
    yr <- c(d$year, fut)
    out <- tibble(indicator = id, indicator_name = name, year = yr,
                  fitted = f_all,
                  ci_lo = f_all - 1.96 * sqrt(pmax(R_all, 0)),
                  ci_hi = f_all + 1.96 * sqrt(pmax(R_all, 0)),
                  model_form = "DLM", snr = snr)
  } else {
    pma <- fit_pma(d$year, d$v)
    if (is.null(pma)) return(NULL)
    yrs <- seq(min(d$year), TARGET_YEAR_TARGET)
    pr <- predict_pma(pma, yrs)
    out <- tibble(indicator = id, indicator_name = name, year = yrs,
                  fitted = pr$mean,
                  ci_lo = pr$mean - 1.96 * sqrt(pmax(pr$var, 0)),
                  ci_hi = pr$mean + 1.96 * sqrt(pmax(pr$var, 0)),
                  model_form = pma$top, snr = snr)
  }
  obs <- d %>% transmute(indicator = id, indicator_name = name, year, observed = v)
  out %>% left_join(obs, by = c("indicator", "indicator_name", "year")) %>%
    mutate(segment = case_when(
      !is.na(observed) ~ "observed",
      year <= 2020 ~ "fitted",
      TRUE ~ "extrapolation"))
}
e4 <- {
  # One example per model form (Table 10 selected_form): poly1/power_law/quad/asymptotic/
  # exponential/segmented + DLM for the high-noise tracks
  meta_names <- read_csv("data/meta.csv", show_col_types = FALSE) %>%
    distinct(indicator_id, indicator_name)
  pick <- function(id) meta_names$indicator_name[match(id, meta_names$indicator_id)]
  message("  Sup Fig. 12: 7 representative forms - 11 linear / 13 power law / 39 quad / 21 asymptotic /",
          " 49 exponential / 3 segmented(OHI) / 14 DLM")
  bind_rows(
    fit_series(11L, pick(11L)),
    fit_series(13L, pick(13L)),
    fit_series(39L, pick(39L)),
    fit_series(21L, pick(21L)),
    fit_series(49L, pick(49L)),
    fit_series(3L,  pick(3L)),
    fit_series(14L, pick(14L), engine = "dlm"))
}
save_csv(e4, de4, "supfig12_model_fits.csv")
save_readme(de4, c(
  "# Supplementary Fig. 12 figure source data",
  "",
  "- supfig12_model_fits.csv: fit and extrapolation of the global mean series (fitted with 95% CI)",
  "- one example per model form: poly1(11) / power_law(13) / quad(39) / asymptotic(21) /",
  "  exponential(49) / segmented(3, OHI, breakpoint) / DLM(14, high-noise track)",
  "- segment distinguishes observed/fitted/extrapolation; 2020 is the boundary of the GBF era; the target year is 2030"))

# =====================================================================
# Extended Data Fig. 1 - counts of the trend/achievement status per target/goal x country and the transitions between them
#   country-level convention as in A.5: rho_c = PGR_sim/PGR_req; the trend follows the bootstrap p of each country's rate
# =====================================================================
message("Extended Data Fig. 1 ...")
de5 <- new_folder("ed_fig1")
pc <- pgr %>%
  filter(!is.na(pgr_sim), !is.na(pgr_req)) %>%
  mutate(
    sig_pos = significant %in% TRUE & rate > 0 & !is.na(p_boot) & p_boot < 0.05,
    rho = ifelse(pgr_req != 0, pgr_sim / pgr_req,
                 ifelse(pgr_sim >= 0, Inf, -Inf)),
    # A.5 zero-requirement special case (same-layer convention): PGR_req=0 and non-negative growth -> recorded as achieved
    class = ifelse(pgr_req == 0 & pgr_sim >= 0, "Full achievement",
                   classify_rho(rho, sig_pos)),
    achievement = case_when(
      class == "Full achievement" ~ "Ac",
      class == "Unfavourable"     ~ "Unac",
      TRUE                        ~ "Unc"),
    trend = case_when(
      rate > 0 & !is.na(p_boot) & p_boot < 0.05 ~ "Up",
      rate < 0 & !is.na(p_boot) & p_boot < 0.05 ~ "Down",
      TRUE ~ "NS"))
e5a <- pc %>% count(layer, layer_id, class) %>%
  complete(nesting(layer, layer_id),
           class = c("Full achievement", "Favourable", "Mixed", "Unfavourable"),
           fill = list(n = 0))
save_csv(e5a, de5, "edfig1a_achievement_counts.csv")
e5b <- pc %>% count(layer, layer_id, trend) %>%
  pivot_wider(names_from = trend, values_from = n, values_fill = 0) %>%
  mutate(n_total = Up + NS + Down) %>%
  select(layer, layer_id, Up, NS, Down, n_total)
save_csv(e5b, de5, "edfig1b_trend_counts.csv")
e5c <- pc %>% count(layer, layer_id, trend, class, name = "n_countries") %>%
  complete(nesting(layer, layer_id), trend = c("Up", "NS", "Down"),
           class = c("Full achievement", "Favourable", "Mixed", "Unfavourable"),
           fill = list(n_countries = 0))
save_csv(e5c, de5, "edfig1c_trend_x_achievement.csv")
e5d <- pc %>%
  mutate(situation = case_when(
    class == "Full achievement" & trend == "Up"   ~ "Ac_Up_sustained",
    class == "Full achievement" & trend == "Down" ~ "Ac_Down_at_risk",
    class == "Favourable" & trend == "Up"   ~ "Fav_Up_toward_Ac",
    class == "Favourable" & trend == "Down" ~ "Fav_Down_at_risk",
    class == "Mixed" & trend == "Up"   ~ "Mixed_Up_toward_Ac",
    class == "Mixed" & trend == "Down" ~ "Mixed_Down_deepening",
    class == "Unfavourable" & trend == "Up"   ~ "Unfav_Up_toward_Ac",
    class == "Unfavourable" & trend == "Down" ~ "Unfav_Down_deepening",
    TRUE                                    ~ "NS_stable")) %>%
  count(layer, layer_id, situation, name = "n_countries") %>%
  complete(nesting(layer, layer_id),
           situation = c("Ac_Up_sustained", "Ac_Down_at_risk",
                         "Fav_Up_toward_Ac", "Fav_Down_at_risk",
                         "Mixed_Up_toward_Ac", "Mixed_Down_deepening",
                         "Unfav_Up_toward_Ac", "Unfav_Down_deepening",
                         "NS_stable"),
           fill = list(n_countries = 0))
save_csv(e5d, de5, "edfig1d_conversion.csv")
save_readme(de5, c(
  "# Extended Data Fig. 1 figure source data",
  "",
  "- edfig1a_achievement_counts.csv: number of countries per target/goal x country in the four achievement classes (A.5: Full achievement/Favourable/Mixed/Unfavourable) (classification per country by rho; the PGR_req=0 special case as above)",
  "- edfig1b_trend_counts.csv: number of countries per target/goal x country by trend status (Up/NS/Down) (bootstrap p<0.05 of each country rate)",
  "- edfig1c_trend_x_achievement.csv: cross-tabulated country counts of the two above, trend x four classes (the full transition table)",
  "- edfig1d_conversion.csv: summary of the transition scenarios (the 9 combinations of four classes x trend: achieved+Up -> sustained, achieved+Down -> risk of loss, not achieved+Up -> moving towards achievement, not achieved+Down -> widening gap, etc.)"))

# =====================================================================
# Supplementary Fig. 2 - geographic and development-status decomposition of the achievement gap in the projection period
#   the country-level convention is the same as A.5 (reusing pc above); the grouping convention is the same as Sup Fig. 1 (5 continent groups /
#   6 status groups DD/DL/LDC/LLDC/SIDS/OECD, a country can belong to several groups)
#   (a,b) pooled four classes (A.5) within the groups for the targets/goals; (c) 27 layers x group %Unfavourable;
#   (d) per-layer KW (achievement order Unfav=0/Mixed=1/Fav=2/Full=3)
# =====================================================================
message("Supplementary Fig. 2 ...")
de6 <- new_folder("sup_fig2")
if (!exists("pc") || !exists("STATUS_LV") || !exists("safe_kw_e1"))
  stop("Sup Fig. 2 requires the pc of the Extended Data Fig. 1 section and the STATUS_LV/safe_kw_e1 of the Sup Fig. 1 section - please run the whole script in order")
geo5 <- cmeta %>% select(country, continent, development, sids, ldc, llc, oecd)
pc5 <- pc %>%
  left_join(geo5, by = "country") %>%
  mutate(ach_ord = case_when(
    class == "Unfavourable"    ~ 0,
    class == "Mixed"           ~ 1,
    class == "Favourable"      ~ 2,
    class == "Full achievement" ~ 3))
cont5 <- pc5 %>% filter(!is.na(continent)) %>%
  transmute(layer, layer_id, country, class, ach_ord,
            grouping = "continent", grp = continent)
stat5 <- pc5 %>%
  transmute(layer, layer_id, country, class, ach_ord,
            DD = development == "developed", DL = development == "developing",
            LDC = ldc %in% TRUE, LLDC = llc %in% TRUE,
            SIDS = sids %in% TRUE, OECD = oecd %in% TRUE) %>%
  pivot_longer(all_of(STATUS_LV), names_to = "grp", values_to = "member") %>%
  filter(member %in% TRUE) %>%
  transmute(layer, layer_id, country, class, ach_ord,
            grouping = "status", grp)
p5g <- bind_rows(cont5, stat5)

# (a/b) Pooled within-group shares for the targets/goals (unit: country-layer pair; n_countries = number of countries in the group)
e6ab <- p5g %>%
  mutate(scope = ifelse(layer == "target", "target", "goal")) %>%
  count(scope, grouping, grp, class, name = "n_pairs") %>%
  complete(nesting(scope, grouping, grp),
           class = c("Full achievement", "Favourable", "Mixed", "Unfavourable"),
           fill = list(n_pairs = 0)) %>%
  group_by(scope, grouping, grp) %>%
  mutate(n_tot = sum(n_pairs)) %>%
  ungroup() %>%
  pivot_wider(names_from = class, values_from = n_pairs) %>%
  left_join(p5g %>%
              mutate(scope = ifelse(layer == "target", "target", "goal")) %>%
              distinct(scope, grouping, grp, country) %>%
              count(scope, grouping, grp, name = "n_countries"),
            by = c("scope", "grouping", "grp")) %>%
  mutate(share_full = `Full achievement` / n_tot, share_fav = Favourable / n_tot,
         share_mixed = Mixed / n_tot, share_unfav = Unfavourable / n_tot) %>%
  select(scope, grouping, grp, n_countries, n_pairs = n_tot,
         `Full achievement`, Favourable, Mixed, Unfavourable,
         share_full, share_fav, share_mixed, share_unfav)
save_csv(e6ab, de6, "supfig2_achv_by_group.csv")

# (c) Status counts and the Unfavourable share per layer x group (the n_countries inside a cell = number of countries in that layer and group)
e6c <- p5g %>%
  count(layer, layer_id, grouping, grp, class, name = "n_pairs") %>%
  complete(nesting(layer, layer_id, grouping, grp),
           class = c("Full achievement", "Favourable", "Mixed", "Unfavourable"),
           fill = list(n_pairs = 0)) %>%
  group_by(layer, layer_id, grouping, grp) %>%
  mutate(n_countries = sum(n_pairs)) %>%
  ungroup() %>%
  pivot_wider(names_from = class, values_from = n_pairs) %>%
  mutate(share_unfav = ifelse(n_countries > 0, Unfavourable / n_countries, NA_real_)) %>%
  select(layer, layer_id, grouping, grp, n_countries,
         `Full achievement`, Favourable, Mixed, Unfavourable, share_unfav)
save_csv(e6c, de6, "supfig2_unac_by_layer_group.csv")

# (d) Per-layer KW: achievement order 0/1/2/3 ~ group (safe_kw_e1 follows the convention of Sup Fig. 1; the status groups overlap -> approximation)
e6d <- p5g %>%
  group_by(layer, layer_id, grouping) %>%
  group_modify(~ {
    kw <- safe_kw_e1(.x$ach_ord, .x$grp)
    tibble(kw_H = kw$H, kw_p = kw$p,
           n = n_distinct(.x$country), n_groups = n_distinct(.x$grp))
  }) %>% ungroup()
save_csv(e6d, de6, "supfig2_kw_achv.csv")
save_readme(de6, c(
  "# Supplementary Fig. 2 figure source data (projected achievement gaps by geography and development status)",
  "",
  "- Definition: the country-level achievement status follows A.5 / Extended Data Fig. 1 (reusing pc; for T6 a zero requirement with non-negative growth is recorded as Ac; for T17 countries with a zero baseline are excluded according to Formula 12);",
  "  the grouping follows Sup Fig. 1 (5 continent groups; 6 status groups DD/DL/LDC/LLDC/SIDS/OECD, a country can belong to several groups -> the status-group tests are an approximation)",
  "- supfig2_achv_by_group.csv: Ac/Unc/Unac counts and shares for the pooled targets/goals x group (unit: country-layer pair; n_countries = number of countries in the group)",
  "- supfig2_unac_by_layer_group.csv: status counts and the Unac share for the 27 layers x groups (the n_countries inside a cell = number of countries in that layer and group)",
  "- supfig2_kw_achv.csv: per-layer Kruskal-Wallis (achievement order Unac=0/Unc=1/Ac=2) grouped by continent/status, kw_H/kw_p/n/n_groups"))

# =====================================================================
# Sup Figs. 5-11 - one panel per uncertainty analysis
# =====================================================================
message("Supplementary Figs. 5-11 ...")
ds1 <- new_folder("sup_fig5")
save_csv(a1, ds1, "supfig5_indicator_subset.csv")
save_readme(ds1, "# Sup Fig. 5: distribution of the target/goal values combined over the k-subsets (one composite value score_c per combination, k = number of indicators in the subset; for k = n there is a single combination shown as a single horizontal line, with the medians connected by a line, Algorithm 1)")

ds2 <- new_folder("sup_fig6")
save_csv(a2,   ds2, "supfig6_gap_length.csv")
save_csv(a2pi, ds2, "supfig6_per_indicator.csv")
save_csv(a2s,  ds2, "supfig6_sensitivity.csv")
save_readme(ds2, c(
  "# Sup Fig. 6: bootstrap imputation uncertainty (Algorithm 2)",
  "",
  "- supfig6_gap_length.csv: relative CI width (width/IQR), stability and entropy grouped by gap length g (baseline Linear/LOO, B=1000, alpha=0.05)",
  "- supfig6_per_indicator.csv: relative and absolute CI width, stability, calibration coverage and entropy per indicator",
  "- supfig6_sensitivity.csv: relative width, stability and coverage for the sensitivity analyses (B / alpha / imputer / rolling-window CV)"))

ds3 <- new_folder("sup_fig7")
save_csv(a3o, ds3, "supfig7_window_length.csv")
save_csv(a3c, ds3, "supfig7_sigma_curves.csv")
save_csv(a3x, ds3, "supfig7_cross_model.csv")
save_readme(ds3, c(
  "# Sup Fig. 7: time-series length uncertainty (Algorithm 3, block-bootstrap + ETS/ARIMA/RW)",
  "",
  "- supfig7_window_length.csv: optimal window L_opt per indicator (Kneedle knee; for flat or non-convex curves the fallback L_der is used, flagged by knee_found), the derivative-based L_der, the series length n_years, sigma(L_opt), the L_opt per model and model_spread",
  "- supfig7_sigma_curves.csv: per-indicator sigma_i(L) curve (pooled + ETS/ARIMA/RW; ci95_halfwidth = 1.96*sigma)",
  "- supfig7_cross_model.csv: cross-model consistency table (L_opt_kneedle / L_opt_der / sigma(L_opt) for indicator x model)"))

ds4 <- new_folder("sup_fig9")
save_csv(a5 %>% arrange(desc(Sx_mean_abs)), ds4, "supfig9_perturbation_sensitivity.csv")
save_readme(ds4, "# Sup Fig. 9: sensitivity to a +/-10% perturbation of a single indicator, |Sx| (Algorithm 5)")

message("All figure source data have been written to ", FIG_DATA)
