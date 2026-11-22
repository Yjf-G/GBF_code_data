# =====================================================================
# Fig 2 — Trajectory projections fall well short of the GBF ambition
#   (a,b) Stacked bars of the country-level share of the four achievement classes per target/goal
#         (A.5: Full achievement / Favourable / Mixed / Unfavourable); T2/T3/T7 are the GBF
#         global-quantity targets, country values there are distance diagnostics only, and the
#         three x labels carry an * (the global verdict is in Supplementary Table 4).
#         2026-11-20: the top Global-verdict strip is removed (the verdict is folded into the
#         Supplementary Table 4 convention).
#   (c) Achievement rate AR with 50%/95% bootstrap intervals (overall/goals/clusters); the summary
#       rows "All targets/All goals" use a darker, larger point + AR value labels (21%/28%).
#   (d,e) target PGR gap map + country-share histogram (median marked only)
#   (f,g) goal  PGR gap map + country-share histogram (median marked only)
#   2026-11-20 layout: d=target map / e=target bars / f=goal map / g=goal bars;
#     the d/f colour bars share the +/-4 pseudo-log range (=+/-10k, matching the bar x axis),
#     one horizontal colour bar at the bottom; a x labels 5.5pt vertical, b x labels 5.5pt
#     horizontal; all legends collected at the bottom.
# Palette (colour-blind friendly, checked): Full achv #00669E / Fav #56B4E9 /
#   Mixed #E69F00 / Unfav #C1441E (blue-to-orange-red order, echoing the blue(surplus)-
#   red(gap) semantics of the d/f maps); all text black; no formula numbers in labels.
# Layout: abc in one row, d|e in one row, f|g in one row; legends at the bottom
# Run from the repository root:  Rscript fig_code/fig2.R
# =====================================================================
source("fig_code/fig_theme.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(sf); library(rnaturalearth)
})

f2c  <- read_csv("data/fig_data/fig2/fig2c_achievement.csv", show_col_types = FALSE)
f2de <- read_csv("data/fig_data/fig2/fig2de_pgr_by_country.csv", show_col_types = FALSE)
f2f  <- read_csv("data/fig_data/fig2/fig2f_pgr_distribution.csv", show_col_types = FALSE)
e5a4 <- read_csv("data/fig_data/ed_fig1/edfig1a_achievement_counts.csv", show_col_types = FALSE)

# ---- Four classes (A.5) and the colour-blind-safe palette ----
cls_lv   <- c("Full achievement", "Favourable", "Mixed", "Unfavourable")
cls_cols <- c(`Full achievement` = "#00669E", `Favourable` = "#56B4E9",
              `Mixed` = "#E69F00", `Unfavourable` = "#C1441E")
GLOBAL_Q <- c("T2", "T3", "T7")   # GBF global-quantity targets
mark_glob <- function(x) ifelse(x %in% GLOBAL_Q, paste0(x, "*"), x)

# (a,b) stacked bars of the country-level four-class shares (percentages)
mk_achv_stack <- function(df, xlab, lv, ttl, tag, x_angle, show_legend = FALSE) {
  d <- df %>%
    group_by(layer_id) %>%
    mutate(share = n / sum(n)) %>%
    ungroup() %>%
    mutate(layer_id = factor(layer_id, levels = lv),
           class = factor(class, levels = cls_lv))
  ggplot(d, aes(layer_id, share, fill = class)) +
    geom_col(width = 0.78, colour = "white", linewidth = 0.12) +
    scale_fill_manual(values = cls_cols, name = "Achievement class",
                      labels = c("Full", "Fav.", "Mixed", "Unfav.")) +
    scale_x_discrete(labels = mark_glob) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = ttl, x = xlab, y = "Share of countries", tag = tag) +
    theme_nature() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 5.5, angle = x_angle,
                                     hjust = if (x_angle == 90) 1 else 0.5,
                                     vjust = if (x_angle == 90) 0.5 else 1),
          legend.position = if (show_legend) "bottom" else "none")
}

pa <- mk_achv_stack(e5a4 %>% filter(layer == "target"), "Target",
                    paste0("T", 1:23), "Targets (up to 2030)", "a", 90)
pb <- mk_achv_stack(e5a4 %>% filter(layer == "goal"), "Goal",
                    c("A", "B", "C", "D"), "Goals (2050)", "b", 0,
                    show_legend = TRUE)

# (c) achievement rate (readable group labels; thick 50% CI / thin 95% CI;
#     summary rows All targets/All goals with darker, larger points + AR labels)
f2c <- f2c %>%
  mutate(lab = case_when(
    group_type == "overall" & group == "TARGET" ~ "All targets",
    group_type == "overall" & group == "GOAL" ~ "All goals",
    group_type == "goal" ~ paste0("Goal ", group),
    group == "reducing_threats" ~ "Reducing threats",
    group == "meeting_needs" ~ "Meeting people's needs",
    group == "tools_solutions" ~ "Tools and solutions",
    TRUE ~ group),
    lab = factor(lab, levels = rev(c("All targets", "All goals",
                                     "Goal A", "Goal B", "Goal C", "Goal D",
                                     "Reducing threats",
                                     "Meeting people's needs",
                                     "Tools and solutions"))),
    is_summary = as.character(lab) %in% c("All targets", "All goals"))
pc <- ggplot(f2c, aes(AR, lab)) +
  geom_errorbarh(aes(xmin = ci95_lo, xmax = ci95_hi), height = 0,
                 linewidth = 0.35, colour = "grey35") +
  geom_errorbarh(aes(xmin = ci50_lo, xmax = ci50_hi), height = 0,
                 linewidth = 0.7, colour = "#0072B2") +
  geom_point(size = 1.2, colour = "#0072B2") +
  geom_point(data = f2c %>% filter(is_summary), size = 2.2,
             colour = "#003A5D") +
  geom_text(data = f2c %>% filter(is_summary),
            aes(label = paste0(round(100 * AR), "%")),
            hjust = 1.35, size = 2, fontface = "bold", colour = "black") +
  scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     limits = c(0, 0.5), expand = expansion(mult = 0.06)) +
  labs(x = "Achievement rate (AR)", y = NULL,
       subtitle = "50% CI (thick) and 95% CI (thin)",
       tag = "c") +
  theme_nature() +
  theme(axis.text.y = element_text(size = 5.5))

# (d,e,f,g) PGR gap maps + country-share histograms (shared colour bar, median marked only;
# countries not included are grey with a legend key)
world <- ne_countries(scale = "small", returnclass = "sf")
map_df <- f2de %>% filter(!is.na(iso3))
pseudo_log <- function(x) sign(x) * log10(1 + abs(x))
PLIM <- 4   # shared pseudo-log range +/-4 (~+/-10k), matching the bar x axis
map_plot <- function(val_col, ttl, tag, show_legend) {
  w <- world %>%
    left_join(map_df, by = c("adm0_a3" = "iso3"))
  w$val_t <- pseudo_log(w[[val_col]])
  w_na <- w %>% filter(is.na(val_t))
  ggplot(w) +
    geom_sf(aes(fill = val_t), colour = "grey85", linewidth = 0.08) +
    # countries not included: grey fill; the colour mapping only generates the "Not included" key
    geom_sf(data = w_na, aes(colour = "Not included"), fill = "grey80",
            linewidth = 0.08) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         midpoint = 0, na.value = "grey80",
                         limits = c(-PLIM, PLIM), breaks = -PLIM:PLIM,
                         labels = c("-10k", "-1k", "-100", "-10", "0",
                                    "10", "100", "1k", "10k"),
                         name = "PGR (%)",
                         guide = if (show_legend) {
                           guide_colourbar(direction = "horizontal",
                                           barwidth = unit(58, "mm"),
                                           barheight = unit(3, "mm"))
                         } else "none") +
    scale_colour_manual(values = c(`Not included` = "grey85"), name = NULL,
                        guide = if (show_legend) {
                          guide_legend(override.aes =
                            list(fill = "grey80", shape = 22, size = 4))
                        } else "none") +
    coord_sf(ylim = c(-58, 84), expand = FALSE) +
    labs(title = ttl, tag = tag) +
    theme_nature() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank())
}

# country-share histogram: signed pseudo-log equal-width bins + median annotation
dist_plot <- function(layer_key, tag) {
  f <- f2f %>% filter(layer == layer_key)
  md <- f$median_val[1]
  ymax <- max(f$share)
  ggplot(f, aes(xmin = bin_lo, xmax = bin_hi, ymin = 0, ymax = share)) +
    geom_rect(fill = "#9CC2E5", colour = "white", linewidth = 0.15) +
    geom_vline(xintercept = pseudo_log(md), linetype = "dashed",
               linewidth = 0.35, colour = "grey20") +
    annotate("text", x = pseudo_log(md), y = 0.92 * ymax,
             label = paste0("median ", round(md), "%"),
             colour = "black", size = 1.8, hjust = -0.06) +
    scale_x_continuous(breaks = -PLIM:PLIM, limits = c(-PLIM, PLIM),
                       labels = c("-10k", "-1k", "-100", "-10", "0",
                                  "10", "100", "1k", "10k")) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                       breaks = scales::pretty_breaks(n = 3),
                       expand = expansion(mult = c(0, 0.06))) +
    labs(x = "Per-country PGR (%)",
         y = "Share of countries", tag = tag) +
    theme_nature()
}

pd <- (map_plot("pgr_target", "Targets: growth gap to 2030", "d", TRUE) |
         dist_plot("target", "e")) +
  plot_layout(widths = c(1.7, 1))
pe <- (map_plot("pgr_goal", "Goals: growth gap to 2050", "f", FALSE) |
         dist_plot("goal", "g")) +
  plot_layout(widths = c(1.7, 1))

row1 <- (pa | pb | pc) + plot_layout(widths = c(1.5, 0.85, 1.1))
fig <- wrap_plots(row1, pd, pe, ncol = 1,
                  heights = c(1.12, 1.12, 1.12)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal")
save_nature(fig, "fig2_trajectories_achievement", 183, 180)
