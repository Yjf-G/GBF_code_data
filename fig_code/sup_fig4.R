# =====================================================================
# Sup Fig. 4 - dot-and-line plot of group-mean trajectories: MIS/EIS 2011-2024 (regions / development-status groups)
#   Layout: 4 rows x 2 columns (row order as in Fig 3: MIS target / MIS goal / EIS target / EIS goal)
#   Columns: the six regions (same abbreviations as Fig 3) | the six development-status groups
#   Points = annual group means (equal-weight average over the countries in the group); data source data/intermediate/mis_eis_annual.csv
#   Source Data: data/fig_data/sup_fig4/supfig4_group_trajectories.csv
# =====================================================================
source("fig_code/fig_theme.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

me <- read_csv("data/intermediate/mis_eis_annual.csv", show_col_types = FALSE)
cm <- read_csv("data/intermediate/country_meta.csv", show_col_types = FALSE)

REG_LV <- c("Africa", "Asia", "Europe", "Oceania", "LAC", "North A.")
GRP_LV <- c("DD", "DL", "LDC", "LLDC", "OECD", "SIDS")
MET_LV <- c("MIS_target", "MIS_goal", "EIS_target", "EIS_goal")

mem <- cm %>%
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

val_long <- me %>%
  pivot_longer(all_of(MET_LV), names_to = "metric", values_to = "value") %>%
  left_join(mem, by = "country")

reg_means <- val_long %>%
  filter(!is.na(region6)) %>%
  group_by(metric, region = region6, year) %>%
  summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(group_type = "region")

grp_means <- val_long %>%
  pivot_longer(starts_with("g_"), names_to = "gname", values_to = "in_grp") %>%
  filter(in_grp %in% TRUE) %>%
  mutate(group = sub("^g_", "", gname)) %>%
  group_by(metric, group, year) %>%
  summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(group_type = "group")

traj <- bind_rows(
  reg_means %>% rename(group = region) %>% select(metric, group_type, group, year, mean_val),
  grp_means %>% select(metric, group_type, group, year, mean_val))

# ---------- Source Data ----------
d5 <- "data/fig_data/sup_fig4"
if (!dir.exists(d5)) dir.create(d5, recursive = TRUE)
write_csv(traj, file.path(d5, "supfig4_group_trajectories.csv"))
writeLines(c(
  "# Sup Fig. 4 figure source data (group-mean trajectories 2011-2024)",
  "",
  "- supfig4_group_trajectories.csv: annual group means (equal-weight average over the countries in the group)",
  "  metric: MIS_target / MIS_goal / EIS_target / EIS_goal; group_type: region (6 regions) / group (6 development-status groups)",
  "  group: Africa/Asia/Europe/Oceania/LAC/North A. and DD/DL/LDC/LLDC/OECD/SIDS (same as Fig 3)",
  "  mean_val = annual group mean; same convention as the ks column of Fig 3 (country-level MIS/EIS, 130 countries passing the coverage rule)",
  "",
  "Related: Supplementary Table 7 is the goal-level (MIS/EIS = MIS_goal/EIS_goal) grouped version of this table (it also carries KS p and the monotonic share)."),
  file.path(d5, "README.md"))

# ---------- Panels ----------
ylabs <- c(MIS_target = "MIS (target level)", MIS_goal = "MIS (goal level)",
           EIS_target = "EIS (target level)", EIS_goal = "EIS (goal level)")
yr_rng <- traj %>% group_by(metric) %>%
  summarise(lo = min(mean_val), hi = max(mean_val), .groups = "drop") %>%
  mutate(pad = (hi - lo) * 0.08)

mk_panel <- function(metric, gtype, tag, title = NULL) {
  levels <- if (gtype == "region") REG_LV else GRP_LV
  rng <- yr_rng %>% filter(metric == !!metric)
  d <- traj %>% filter(metric == !!metric, group_type == !!gtype) %>%
    mutate(group = factor(group, levels = levels))
  ggplot(d, aes(year, mean_val, colour = group)) +
    geom_line(linewidth = 0.35) +
    geom_point(size = 0.6) +
    scale_colour_manual(values = PAL8[seq_along(levels)], breaks = levels,
                        drop = FALSE,
                        name = if (gtype == "region") "Region" else "Development group") +
    scale_x_continuous(breaks = seq(2011, 2024, by = 3),
                       labels = as.character(seq(2011, 2024, by = 3)),
                       expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = c(rng$lo - rng$pad, rng$hi + rng$pad),
                       breaks = scales::pretty_breaks(4)) +
    labs(tag = tag, title = title, x = "Year", y = ylabs[[metric]]) +
    theme_nature() +
    theme(axis.text = element_text(size = 5),
          axis.title = element_text(size = 6),
          plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
          plot.tag = element_text(size = 8, face = "bold"),
          legend.key.size = unit(2.8, "mm"),
          legend.text = element_text(size = 4.5),
          legend.title = element_text(size = 5),
          plot.margin = margin(2, 2, 1, 1))
}

rows <- lapply(seq_along(MET_LV), function(i) {
  m <- MET_LV[i]
  tags <- c(letters[(i - 1) * 2 + 1], letters[(i - 1) * 2 + 2])
  ttl <- if (i == 1) c("Geographic regions", "Development-status groups") else c(NULL, NULL)
  (mk_panel(m, "region", tags[1], ttl[1]) | mk_panel(m, "group", tags[2], ttl[2])) +
    plot_layout(widths = c(1, 1))
})

fig <- wrap_plots(rows, ncol = 1, heights = c(1, 1, 1, 1)) +
  plot_layout(guides = "collect") +
  plot_annotation(theme = theme(plot.margin = margin(1, 2, 1, 1)))
save_nature(fig, "supfig4_group_trajectories", 183, 175)
