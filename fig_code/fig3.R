# =====================================================================
# Fig 3 — Spatiotemporal dynamics of MIS and EIS, 2011–2024
#   Design reference: figures/fig_4&5/plot_fig4_fig5_combined_target_first.py
#   Layout: 4 rows x [geographic triangular matrix + ks | development-group triangular matrix + ks | delta map]
#   Row identity (MIS/EIS x target/goal) is given in each small-panel title (no left-hand row labels)
#   Triangular matrix: upper triangle 2011 (blue text), lower triangle 2024 (red text); diagonal corners labelled 2011/2024
#   Values = X - Y (column group minus row group) difference of group means, as stated in the caption
#   Fill = significance (triangle: Dunn-Bonferroni p; ks column: KS p): three red levels + n.s. grey
#   No significance stars; significant colours are clearly separated from grey
#   Square cells (coord_fixed); the identical legend is shown once below each matrix; the map column is widened
#   2026-11-11: countries not included are filled grey (grey80) with a "Not included" legend key at the bottom;
#               matrix numbers 1.3 -> 1.7 mm (reaching the Nature 5 pt lower limit); legend text 4.5 -> 5 pt
# =====================================================================
source("fig_code/fig_theme.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(sf)
  library(rnaturalearth)
})

tri_r <- read_csv("data/fig_data/fig3/fig3_tri_regions.csv", show_col_types = FALSE)
tri_g <- read_csv("data/fig_data/fig3/fig3_tri_groups.csv", show_col_types = FALSE)
ks_df <- read_csv("data/fig_data/fig3/fig3_ks.csv", show_col_types = FALSE)
map_d <- read_csv("data/fig_data/fig3/fig3_map_delta.csv", show_col_types = FALSE)
CAP <- round(as.numeric(quantile(abs(map_d$delta), 0.985, na.rm = TRUE)))

REG_LV <- c("Africa", "Asia", "Europe", "Oceania", "LAC", "North A.")
GRP_LV <- c("DD", "DL", "LDC", "LLDC", "OECD", "SIDS")

# Significance classes (red hues distinguishable from grey; darker = more significant)
SIG_LEVELS <- c("***", "**", "*", "n.s.")
SIG_FILL <- c("***" = "#7F0000", "**" = "#D6604D", "*" = "#F4A582",
              "n.s." = "#E0E0E0")
SIG_LABELS <- c("P < 0.001", "P < 0.01", "P < 0.05", "n.s.")
TXT_2011 <- "#1565C0"   # reference script: 2011 in blue
TXT_2024 <- "#C62828"   # reference script: 2024 in red

level_of <- function(p) {
  factor(ifelse(is.na(p), "n.s.",
                ifelse(p < 0.001, "***",
                       ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "n.s.")))),
         levels = SIG_LEVELS)
}

# ---------- Triangular matrix + ks column ----------
tri_panel <- function(tri, ks, levels, metric, tag, ttl) {
  tri <- tri %>% filter(metric == !!metric) %>%
    mutate(g1 = factor(g1, levels = levels), g2 = factor(g2, levels = levels),
           i1 = as.integer(g1), i2 = as.integer(g2))
  # p_bonf is stored only on the canonical row g1<g2; the 2024 lower triangle takes the p of the mirrored row
  p_canon <- tri %>% filter(i1 < i2) %>%
    transmute(year, j1 = i1, j2 = i2, p_c = p_bonf)
  cells <- bind_rows(
    tri %>% filter(year == 2011, i1 < i2),
    tri %>% filter(year == 2024, i1 > i2)
  ) %>%
    mutate(j1 = pmin(i1, i2), j2 = pmax(i1, i2)) %>%   # diff = X-Y (column minus row), consistent with the figure
    left_join(p_canon, by = c("year", "j1", "j2")) %>%
    mutate(p = ifelse(is.na(p_bonf), p_c, p_bonf),
           lvl = level_of(p),
           label = sprintf("%.1f", diff),
           txt_col = ifelse(as.character(lvl) %in% c("***", "**"), "white",
                            ifelse(year == 2011, TXT_2011, TXT_2024))) %>%
    mutate(g2 = factor(as.character(g2), levels = c(levels, "ks")))
  ks_cells <- ks %>% filter(metric == !!metric) %>%
    mutate(g1 = factor(group, levels = levels),
           g2 = factor("ks", levels = c(levels, "ks")),
           lvl = level_of(ks_p),
           label = ifelse(abs(delta_mean) >= 1,
                          sprintf("%+.1f", delta_mean),
                          sprintf("%+.2f", delta_mean)),
           txt_col = ifelse(as.character(lvl) %in% c("***", "**"), "white", "black"))
  n <- length(levels)
  yr_marks <- tibble(
    g1 = factor(c(levels[1], levels[n]), levels = levels),
    g2 = factor(c(levels[1], levels[n]), levels = c(levels, "ks")),
    label = c("2011", "2024"), col = c(TXT_2011, TXT_2024))
  ggplot() +
    geom_tile(data = cells, aes(g2, g1, fill = lvl), colour = "white",
              linewidth = 0.35) +
    geom_tile(data = ks_cells, aes(g2, g1, fill = lvl), colour = "white",
              linewidth = 0.35) +
    geom_text(data = cells, aes(g2, g1, label = label, colour = txt_col),
              size = 1.77, fontface = "bold") +
    geom_text(data = ks_cells, aes(g2, g1, label = label, colour = txt_col),
              size = 1.77, fontface = "bold") +
    geom_text(data = yr_marks, aes(g2, g1, label = label, colour = col),
              size = 1.77, fontface = "italic") +
    scale_colour_identity(guide = "none") +
    scale_fill_manual(values = SIG_FILL, breaks = SIG_LEVELS,
                      labels = SIG_LABELS, drop = FALSE, guide = "none") +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(limits = rev(levels), drop = FALSE) +
    coord_fixed(expand = FALSE) +
    labs(tag = tag, title = ttl, x = NULL, y = NULL) +
    theme_nature() +
    theme(axis.text.x = element_text(size = 5, face = "bold", angle = 45,
                                     hjust = 1),
          axis.text.y = element_text(size = 5, face = "bold"),
          plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
          plot.tag = element_text(size = 8, face = "bold"),
          panel.grid = element_blank(),
          plot.margin = margin(2, 2, 1, 1))
}

# ---------- Per-country delta map ----------
world <- ne_countries(scale = "small", returnclass = "sf")
map_panel <- function(metric, tag, ttl) {
  w <- world %>% left_join(map_d %>% filter(metric == !!metric),
                           by = c("adm0_a3" = "iso3"))
  ggplot(w) +
    geom_sf(aes(fill = delta), colour = "grey88", linewidth = 0.07) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, na.value = "grey80",
                         limits = c(-CAP, CAP), oob = scales::squish,
                         guide = "none") +
    coord_sf(ylim = c(-58, 84), expand = FALSE) +
    labs(tag = tag, title = ttl, x = NULL, y = NULL) +
    theme_nature() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(), axis.title = element_blank(),
          plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
          plot.tag = element_text(size = 8, face = "bold"),
          plot.margin = margin(2, 1, 1, 1))
}

# ---------- Bottom legend row (exact grid layout; one shared legend band for the three groups) ----------
# The grey "Not included" key sits below the colour bar on the same horizontal line as the tick numbers (2026-11-11 comment)
library(grid)
ramp_cols <- scales::gradient_n_pal(c("#2166AC", "white", "#B2182B"))(seq(0, 1, length.out = 256))
tp <- gpar(fontsize = 5, fontfamily = "Arial")          # legend text 5 pt
tpb <- gpar(fontsize = 5, fontface = "bold", fontfamily = "Arial")
sig_x <- c(26, 41, 55.5, 69)                          # left edge of the significance key squares (mm)
num_x <- 113 + 7 * (0:6)                              # centre of the colour bar tick numbers
num_lab <- sprintf("%.0f", seq(-CAP, CAP, length.out = 7))

leg_grob <- grobTree(
  # Significance (title on the same row as the keys, as in the current version)
  textGrob("Significance", x = unit(15, "mm"), y = unit(6.3, "mm"), gp = tpb),
  rectGrob(x = unit(sig_x, "mm"), y = unit(6.3, "mm"),
           width = unit(3.2, "mm"), height = unit(2.6, "mm"),
           gp = gpar(fill = SIG_FILL, col = NA)),
  textGrob(SIG_LABELS, x = unit(sig_x + 4.4, "mm"), y = unit(6.3, "mm"),
           just = "left", gp = tp),
  # Delta colour bar + title (upper row)
  textGrob(sprintf("Δ 2024–2011 (capped ±%d)", CAP), x = unit(82, "mm"),
           y = unit(6.3, "mm"), just = "left", gp = tpb),
  rasterGrob(as.raster(matrix(ramp_cols, nrow = 1)),
             x = unit(134, "mm"), y = unit(6.3, "mm"),
             width = unit(42, "mm"), height = unit(2.8, "mm")),
  # Tick numbers (lower row)
  textGrob(num_lab, x = unit(num_x, "mm"), y = unit(2.6, "mm"), gp = tp),
  # Not included grey key: below the colour bar, on the same horizontal line as the numbers (right end of the number row)
  rectGrob(x = unit(159.2, "mm"), y = unit(2.6, "mm"),
           width = unit(3.2, "mm"), height = unit(2.6, "mm"),
           gp = gpar(fill = "grey80", col = NA)),
  textGrob("Not included", x = unit(162.2, "mm"), y = unit(2.6, "mm"),
           just = "left", gp = tp)
)

# ---------- Assembly ----------
panel_row <- function(metric, tags, score, level) {
  lv <- tolower(level)
  id <- sprintf("%s %s", score, lv)          # row identity label: MIS target etc.
  p1 <- tri_panel(tri_r, ks_df %>% filter(group_type == "region"),
                  REG_LV, metric, tags[1], sprintf("%s · Geography", id))
  p2 <- tri_panel(tri_g, ks_df %>% filter(group_type == "group"),
                  GRP_LV, metric, tags[2], sprintf("%s · Development", id))
  p3 <- map_panel(metric, tags[3],
                  sprintf("%s change
2011–2024", id))
  (p1 | p2 | p3) +
    plot_layout(widths = c(1, 1, 1.95))
}

rows <- list(
  panel_row("MIS_target", c("a", "b", "c"), "MIS", "Target"),
  panel_row("MIS_goal",   c("d", "e", "f"), "MIS", "Goal"),
  panel_row("EIS_target", c("g", "h", "i"), "EIS", "Target"),
  panel_row("EIS_goal",   c("j", "k", "l"), "EIS", "Goal")
)
leg_grob_shifted <- grobTree(leg_grob,
  vp = viewport(x = unit(-1.94, "mm"), just = "left", clip = "off"))
legrow <- wrap_elements(full = leg_grob_shifted)
fig <- wrap_plots(c(rows, list(legrow)), ncol = 1, heights = c(1, 1, 1, 1, 0.2))
save_nature(fig, "fig3_spatiotemporal_mis_eis", 183, 196)
