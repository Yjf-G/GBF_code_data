# =====================================================================
# Supplementary Fig. 1 — heterogeneity of the baseline summary statistics across continents and development-status groups
#   All cells are equal squares (coord_fixed(ratio=1) plus balanced slot geometry)
#   (a) 23 targets x 8 columns (Continent | Status blocks x 4 parameters):
#       omnibus KW test p (grey = n.s., p<0.05 / p<0.01 / p<0.001)
#   (b) The 4 goals, as in (a)
#   (c) 23 targets x pairwise significance matrix for the 10 continent pairs (4 parameter blocks x 10 pairs;
#       MW p, Bonferroni x 10, omnibus-gated)
#   (d) The 4 goals, as in (c)
#   (e) 23 targets x pairwise significance matrix for the 15 status pairs (4 parameter blocks x 15 pairs; Bonferroni x 15)
#   (f) The 4 goals, as in (e)
# Run from the repository root:  Rscript fig_code/sup_fig1.R
# =====================================================================
source("fig_code/fig_theme.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

e1k  <- read_csv("data/fig_data/sup_fig1/supfig1_kw_tests.csv", show_col_types = FALSE)
e1g  <- read_csv("data/fig_data/sup_fig1/supfig1_pair_gaps.csv", show_col_types = FALSE)
pw_c <- read_csv("data/fig_data/sup_fig1/supfig1_continent_pairwise_tests.csv",
                 show_col_types = FALSE)
pw_s <- read_csv("data/fig_data/sup_fig1/supfig1_status_pairwise_tests.csv",
                 show_col_types = FALSE)

PARAM_LEVELS <- c("mean_v", "sd_v", "slope", "range_v")
P_SHORT <- c(mean_v = "Mean", sd_v = "SD", slope = "Slope", range_v = "Range")
GROUP2 <- c("continent", "status")
GROUP2_LAB <- c(continent = "Continent", status = "Status")

pkat <- function(p) factor(
  ifelse(is.na(p), "n.s.",
  ifelse(p < 0.001, "p < 0.001",
  ifelse(p < 0.01, "p < 0.01",
  ifelse(p < 0.05, "p < 0.05", "n.s.")))),
  levels = c("n.s.", "p < 0.05", "p < 0.01", "p < 0.001"))
P_COLS <- c("n.s." = "#D9D9D9", "p < 0.05" = "#F0B0A8",
            "p < 0.01" = "#D66060", "p < 0.001" = "#B2182B")

# Pair order: 10 continent pairs / 15 status pairs, sorted by descending |dz| gap (auxiliary table supfig1_pair_gaps)
pair_rank <- e1g %>%
  group_by(pair_type, grp1, grp2) %>%
  summarise(m = mean(mean_gap_z), .groups = "drop") %>%
  group_by(pair_type) %>%
  arrange(desc(m), .by_group = TRUE) %>%
  transmute(pair = paste(grp1, grp2, sep = " – "))
PAIR_CONT <- filter(pair_rank, pair_type == "continent")$pair
PAIR_STAT <- filter(pair_rank, pair_type == "status")$pair

# ---- Common: 90 degree x labels and square cells ----
hm_theme <- function(show_y) {
  theme(axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = if (show_y) element_text(size = 5) else element_blank(),
        axis.ticks.y = if (show_y) element_line(linewidth = 0.25) else element_blank(),
        axis.title.y = element_text(size = 6),
        plot.title = element_text(size = 6, hjust = 0),
        legend.text = element_text(size = 5.5),
        legend.title = element_text(size = 6),
        legend.position = "none")
}

# (a/b) Omnibus KW: 2 blocks (Continent | Status) x 4 parameters = 8 columns
X_LEV_O <- unlist(lapply(GROUP2, function(g) paste(g, PARAM_LEVELS, sep = "|")))
hm_omni <- function(df, id_levels, y_title, tag, ttl) {
  n_id <- length(id_levels)
  d <- df %>%
    filter(grouping %in% GROUP2) %>%
    mutate(x = factor(paste(grouping, param, sep = "|"), levels = X_LEV_O),
           id = factor(layer_id, levels = id_levels),
           pk = pkat(kw_p))
  ggplot(d, aes(x, id, fill = pk)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_vline(xintercept = 4.5, colour = "white", linewidth = 0.9) +
    annotate("text", x = c(2.5, 6.5), y = n_id + 0.75,
             label = unname(GROUP2_LAB), size = 1.75, fontface = "bold") +
    scale_fill_manual(values = P_COLS, name = "p-value", drop = FALSE) +
    scale_x_discrete(labels = rep(unname(P_SHORT[PARAM_LEVELS]), 2),
                     expand = expansion(mult = 0, add = c(0.3, 0.3))) +
    scale_y_discrete(expand = expansion(mult = 0, add = c(0.3, 1.5))) +
    coord_fixed(ratio = 1) +
    labs(x = NULL, y = y_title, tag = tag, title = ttl) +
    theme_nature() + hm_theme(show_y = TRUE) +
    theme(plot.margin = margin(1.5, 2, 6.5, 1.5))
}
pa <- hm_omni(filter(e1k, layer == "target"), rev(paste0("T", 1:23)), "Target",
              "a", "Omnibus")
pb <- hm_omni(filter(e1k, layer == "goal"), rev(c("A", "B", "C", "D")), "Goal",
              "b", "Omnibus")

# (c/d/e/f) Pairwise significance matrices: 4 parameter blocks x n pairs
hm_pw <- function(df, pair_order, id_levels, y_title, tag, ttl, show_y = TRUE) {
  n_id <- length(id_levels)
  n_pr <- length(pair_order)
  X_LEV <- unlist(lapply(PARAM_LEVELS,
                         function(p) paste(p, pair_order, sep = "|")))
  d <- df %>%
    mutate(pair = paste(grp1, grp2, sep = " – "),
           x = factor(paste(param, pair, sep = "|"), levels = X_LEV),
           id = factor(layer_id, levels = id_levels),
           pk = pkat(ifelse(!is.na(omnibus_p) & omnibus_p < 0.05, p_bonf, NA_real_)))
  ggplot(d, aes(x, id, fill = pk)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_vline(xintercept = seq(n_pr + 0.5, 3 * n_pr + 0.5, by = n_pr),
               colour = "white", linewidth = 0.9) +
    annotate("text", x = seq(n_pr / 2, 4 * n_pr - n_pr / 2, by = n_pr),
             y = n_id + 0.75, label = unname(P_SHORT[PARAM_LEVELS]),
             size = 1.75, fontface = "bold") +
    scale_fill_manual(values = P_COLS, name = "p-value", drop = FALSE) +
    scale_x_discrete(labels = rep(gsub(" ", "", pair_order), length(PARAM_LEVELS)),
                     expand = expansion(mult = 0, add = c(0.3, 0.3))) +
    scale_y_discrete(expand = expansion(mult = 0, add = c(0.3, 1.5))) +
    coord_fixed(ratio = 1) +
    labs(x = NULL, y = if (show_y) y_title else NULL, tag = tag, title = ttl) +
    theme_nature() + hm_theme(show_y = show_y)
}
pc <- hm_pw(filter(pw_c, layer == "target"), PAIR_CONT, rev(paste0("T", 1:23)),
            "Target", "c", "Continent pairwise", show_y = FALSE)
pd <- hm_pw(filter(pw_c, layer == "goal"), PAIR_CONT, rev(c("A", "B", "C", "D")),
            "Goal", "d", "Continent pairwise", show_y = FALSE)
pe <- hm_pw(filter(pw_s, layer == "target"), PAIR_STAT, rev(paste0("T", 1:23)),
            "Target", "e", "Status pairwise", show_y = TRUE)
pf <- hm_pw(filter(pw_s, layer == "goal"), PAIR_STAT, rev(c("A", "B", "C", "D")),
            "Goal", "f", "Status pairwise", show_y = TRUE)

# ---- Legend: attached below f (full width, ample room); the whole figure uses one shared p colour scale ----
pf <- pf + theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(fill = guide_legend(nrow = 1, title.position = "left"))

# ---- Layout (flat design, controllable proportions; cells are equal squares):
#   rows 1-2: [a|c|placeholder] [b|d|placeholder];  row 3: e full width;  row 4: f full width (legend below)
sp <- ggplot() + theme_void()
fig <- wrap_plots(pa, pb, pc, pd, pe, pf, sp,
                  design = "ACG
BDG
EEE
FFF",
                  widths = c(1, 4.14, 1.07), heights = c(4.2, 1, 4.2, 1))
# 2026-11-14: the explanatory text at the bottom was removed (the information is now in the SI caption); uncomment below to restore it
# fig <- fig + plot_annotation(
#   caption = paste0(
#     "Colour: p-value category of the corresponding test.
# ",
#     "a, b: omnibus Kruskal–Wallis p (across continents or status groups). ",
#     "c–f: pairwise Mann–Whitney U p,
# Bonferroni-corrected within parameter, ",
#     "omnibus-gated (Supplementary Table 3)."),
#   theme = theme(plot.caption = element_text(hjust = 0, size = 5)))
save_nature(fig, "supfig1_geographic_heterogeneity", 183, 270)
