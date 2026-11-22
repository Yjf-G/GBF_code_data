# =====================================================================
# Extended Data Fig. 1 - per-country trend/achievement status distributions for each target/goal and the transitions between them
#   (a) Stacked counts of countries in the four achievement classes (A.5: Full achievement/Favourable/Mixed/Unfavourable)
#       (same colours as Fig 2)
#   (b) Cross-tabulation heatmap of trend (Up/NS/Down) x the four achievement classes (colour depth = within-row share, number = country count)
#   (c) Stacked counts of countries by trend status (Up/NS/Down)
#   2026-11-11: the three-class achievement dimension (Ac/Unc/Unac) is replaced by the four A.5 classes, aligned with Fig 2
#   2026-11-11 v2: the three legend groups (achievement class / Row share colour bar / trend) are collected into a row at the bottom
#     (same convention as Fig 1/4); panel b is widened so that the "Up NS Down" ticks display correctly in every facet
#     (at the old width the 4.5 mm cells were smaller than the 8 mm "Down" label and the labels collided)
# Run from the repository root:  Rscript fig_code/ed_fig1.R
# =====================================================================
source("fig_code/fig_theme.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

e5a <- read_csv("data/fig_data/ed_fig1/edfig1a_achievement_counts.csv", show_col_types = FALSE)
e5b <- read_csv("data/fig_data/ed_fig1/edfig1b_trend_counts.csv", show_col_types = FALSE)
e5c <- read_csv("data/fig_data/ed_fig1/edfig1c_trend_x_achievement.csv", show_col_types = FALSE)

# ---- Four classes (A.5) and colours (identical to Fig 2) ----
cls_lv   <- c("Full achievement", "Favourable", "Mixed", "Unfavourable")
cls_cols <- c(`Full achievement` = "#00669E", `Favourable` = "#56B4E9",
              `Mixed` = "#E69F00", `Unfavourable` = "#C1441E")
trend_cols <- c(Up = "#0072B2", NS = "#999999", Down = "#D55E00")
GLOBAL_Q <- c("T2", "T3", "T7")   # GBF targets with a global scope: country values are a baseline-distance diagnostic, the decision follows the aggregate trajectory
mark_glob <- function(x) ifelse(x %in% GLOBAL_Q, paste0(x, "*"), x)

# Row order: targets T1-T23 on top, goals A-D below
e5c <- e5c %>%
  mutate(lab = factor(layer_id, levels = c(paste0("T", 1:23), "A", "B", "C", "D")),
         trend = factor(trend, levels = c("Up", "NS", "Down")),
         class = factor(class, levels = cls_lv)) %>%
  group_by(lab) %>%
  mutate(share = n_countries / sum(n_countries)) %>%
  ungroup()

# (a) Stacked bars of the four achievement classes
e5a_l <- e5a %>%
  mutate(lab = factor(layer_id, levels = c(paste0("T", 1:23), "A", "B", "C", "D")),
         class = factor(class, levels = cls_lv))
pa <- ggplot(e5a_l, aes(n, lab, fill = class)) +
  geom_col(width = 0.78, colour = "white", linewidth = 0.15) +
  scale_fill_manual(values = cls_cols, name = "Achievement class",
                    labels = cls_lv,
                    guide = guide_legend(order = 1)) +
  scale_x_continuous(expand = expansion(mult = 0.01)) +
  scale_y_discrete(labels = mark_glob) +
  labs(x = "Countries", y = NULL, tag = "a") +
  theme_nature() +
  theme(axis.text.y = element_text(size = 5),
        legend.position = "none")

# (b) Trend x four achievement classes heatmap: class faceted (columns) x trend (within-row columns), 27x3 cells x 4
# Short facet labels (full names in the legend of panel a) to avoid clipping
cls_short <- c(`Full achievement` = "Full", `Favourable` = "Fav.",
               `Mixed` = "Mixed", `Unfavourable` = "Unfav.")
e5c <- e5c %>% mutate(class_lab = factor(cls_short[as.character(class)],
                                         levels = unname(cls_short)))
pb <- ggplot(e5c, aes(trend, lab)) +
  geom_tile(aes(fill = share), colour = "white", linewidth = 0.35) +
  geom_text(aes(label = ifelse(n_countries > 0, n_countries, "")),
            size = 1.7, colour = "grey15") +
  facet_wrap(~class_lab, ncol = 4) +
  scale_fill_gradient(low = "#F2F2F2", high = "#9CC2E5",
                      labels = function(x) paste0(round(100 * x), "%"),
                      name = "Row share",
                      guide = guide_colourbar(order = 2,
                        barwidth = unit(24, "mm"), barheight = unit(2.2, "mm"))) +
  scale_x_discrete(expand = c(0.02, 0.02)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Trend (Up / NS / Down)", y = NULL, tag = "b") +
  theme_nature() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 5),
        strip.text = element_text(size = 6),
        legend.position = "none")

# (c) Stacked bars of trend status
e5b_l <- e5b %>%
  mutate(lab = factor(layer_id, levels = c(paste0("T", 1:23), "A", "B", "C", "D"))) %>%
  pivot_longer(c(Up, NS, Down), names_to = "trend", values_to = "n") %>%
  mutate(trend = factor(trend, levels = c("Up", "NS", "Down")))
pc <- ggplot(e5b_l, aes(n, lab, fill = trend)) +
  geom_col(width = 0.78, colour = "white", linewidth = 0.15) +
  scale_fill_manual(values = trend_cols, name = "Trend",
                    labels = c("Up", "NS", "Down"),
                    guide = guide_legend(order = 3)) +
  scale_x_continuous(expand = expansion(mult = 0.01)) +
  labs(x = "Countries", y = NULL, tag = "c") +
  theme_nature() +
  theme(axis.text.y = element_blank(),
        legend.position = "none")

# Layout: three panels in one row; b is widened (so that the Up/NS/Down ticks of every facet do not collide); the three legend groups are gathered into a row at the bottom
fig <- (pa | pb | pc) +
  plot_layout(widths = c(1.05, 1.5, 0.85), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.key.size = unit(2.6, "mm"))
save_nature(fig, "edfig1_trend_achievement_counts", 183, 175)
