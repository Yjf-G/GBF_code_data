# =====================================================================
# Sup Fig. 3 - Global trends of MIS and EIS, 2011-2024: annual means, year-on-year changes and Wilcoxon significance
#   Layout: 2 rows x 2 columns
#     (a) global annual means of MIS (target level / goal level)  (b) global annual means of EIS
#     (c) bar chart of the year-on-year changes of MIS + Wilcoxon stars  (d) year-on-year changes of EIS + stars
#   Data source: data/intermediate/mis_eis_annual.csv
#   Wilcoxon: paired within country, two-sided, normal approximation with continuity correction (identical to SI B.2.7 / Supplementary Table 6)
#   Source Data: data/fig_data/sup_fig3/supfig3_global_paths.csv,
#                supfig3_yoy_changes.csv, supfig3_yoy_tests.csv
# =====================================================================
source("fig_code/fig_theme.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

me <- read_csv("data/intermediate/mis_eis_annual.csv", show_col_types = FALSE)
MET <- c("MIS_target", "MIS_goal", "EIS_target", "EIS_goal")
LV_LAB <- c(target = "Target level", goal = "Goal level")

# ---------- Global annual means ----------
glob <- me %>%
  pivot_longer(all_of(MET), names_to = "series", values_to = "value") %>%
  group_by(series, year) %>%
  summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(metric = ifelse(grepl("^MIS", series), "MIS", "EIS"),
         level = ifelse(grepl("target$", series), "target", "goal"))

# ---------- Year-on-year changes ----------
yoy <- glob %>%
  arrange(series, year) %>%
  group_by(series) %>%
  mutate(delta = mean_val - lag(mean_val)) %>%
  filter(!is.na(delta)) %>%
  ungroup() %>%
  mutate(pair = factor(paste0(year - 1, "–", year),
                       levels = paste0(2011:2023, "–", 2012:2024)))

# ---------- Year-on-year paired Wilcoxon (paired within country, identical convention to Table 6) ----------
wtest <- function(s, y0) {
  w <- me %>%
    filter(year %in% c(y0, y0 + 1)) %>%
    select(country, year, value = all_of(s)) %>%
    pivot_wider(names_from = year, values_from = value) %>%
    drop_na()
  d <- w[[as.character(y0 + 1)]] - w[[as.character(y0)]]
  d <- d[abs(d) > 1e-12]
  wt <- wilcox.test(d, exact = FALSE, correct = TRUE)
  tibble(series = s, y0 = y0, n = length(d),
         mean_delta = mean(d), W = unname(wt$statistic), p = wt$p.value)
}
tests <- bind_rows(lapply(MET, function(s)
  bind_rows(lapply(2011:2023, function(y) wtest(s, y)))))
tests <- tests %>%
  mutate(sig = case_when(p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns"))

yoy <- yoy %>% left_join(tests %>% transmute(series, year = y0 + 1, sig),
                         by = c("series", "year"))

# ---------- Source Data ----------
d9 <- "data/fig_data/sup_fig3"
if (!dir.exists(d9)) dir.create(d9, recursive = TRUE)
write_csv(glob, file.path(d9, "supfig3_global_paths.csv"))
write_csv(yoy %>% select(series, pair, delta), file.path(d9, "supfig3_yoy_changes.csv"))
write_csv(tests, file.path(d9, "supfig3_yoy_tests.csv"))
writeLines(c(
  "# Sup Fig. 3 figure source data (global trends of MIS and EIS, 2011-2024)",
  "# supfig3_global_paths.csv : global annual means (equal-weight average over countries)",
  "# supfig3_yoy_changes.csv  : year-on-year changes of the global means",
  "# supfig3_yoy_tests.csv    : year-on-year paired Wilcoxon (= the values of Supplementary Table 6)"),
  file.path(d9, "README.txt"))

# ---------- Plotting ----------
MK <- c(target = PAL8[1], goal = PAL8[2])

p_path <- function(metric, tag) {
  d <- glob %>% filter(.data$metric == !!metric) %>%
    mutate(level = factor(level, levels = c("target", "goal")))
  rng <- range(d$mean_val)
  pad <- diff(rng) * 0.12
  ggplot(d, aes(year, mean_val, colour = level, shape = level)) +
    geom_line(linewidth = 0.4) +
    geom_point(size = 1.1) +
    scale_colour_manual(values = MK, labels = LV_LAB, name = NULL) +
    scale_shape_manual(values = c(target = 16, goal = 17), labels = LV_LAB, name = NULL) +
    scale_x_continuous(breaks = seq(2011, 2024, by = 3),
                       expand = expansion(mult = c(0.03, 0.03))) +
    scale_y_continuous(limits = c(rng[1] - pad, rng[2] + pad),
                       breaks = scales::pretty_breaks(4)) +
    labs(tag = tag, x = "Year", y = paste("Global mean", metric)) +
    theme_nature() +
    theme(axis.text = element_text(size = 5),
          axis.title = element_text(size = 6),
          plot.tag = element_text(size = 8, face = "bold"),
          legend.position = "bottom",
          legend.background = element_blank(),
          legend.key.size = unit(2.6, "mm"),
          legend.text = element_text(size = 4.5),
          plot.margin = margin(2, 2, 1, 1))
}

p_bars <- function(metric, tag) {
  d <- yoy %>% filter(.data$metric == !!metric) %>%
    mutate(level = factor(level, levels = c("target", "goal")))
  rng <- range(d$delta)
  pad <- diff(rng) * 0.22
  ggplot(d, aes(pair, delta, fill = level)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.68) +
    geom_hline(yintercept = 0, linewidth = 0.25) +
    geom_text(aes(label = sig, group = level,
                  vjust = ifelse(delta >= 0, -0.45, 1.45)),
              position = position_dodge(width = 0.72), size = 1.7) +
    scale_fill_manual(values = MK, labels = LV_LAB, name = NULL) +
    scale_y_continuous(limits = c(rng[1] - pad, rng[2] + pad),
                       breaks = scales::pretty_breaks(4)) +
    labs(tag = tag, x = "Year pair", y = paste("Change in global mean", metric)) +
    theme_nature() +
    theme(axis.text = element_text(size = 4.2),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(size = 6),
          plot.tag = element_text(size = 8, face = "bold"),
          legend.position = "none",
          plot.margin = margin(2, 2, 1, 1))
}

fig <- wrap_plots(list(p_path("MIS", "a"), p_path("EIS", "b"),
                       p_bars("MIS", "c"), p_bars("EIS", "d")),
                  ncol = 2, heights = c(1, 1.15)) +
  plot_annotation(theme = theme(plot.margin = margin(1, 2, 1, 1)))
save_nature(fig, "supfig3_global_trends", 183, 118)
