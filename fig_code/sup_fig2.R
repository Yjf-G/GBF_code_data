# =====================================================================
# Supplementary Fig. 2 - projected achievement gaps by geography and development status
#   (a, b) Pooled four-class shares (A.5) per group for targets/goals (5 continents + 6 status groups;
#          group order sorted by descending target-level Unfav share, worst class at the top of each block; n = number of countries in the group)
#   (c) Heatmap of %Unfavourable for the 27 layers (T1-T23, A-D; T1 at the top) x 11 groups (each cell labelled with n countries; column order = group order)
#   (d) Per-layer KW significance (achievement order Unfav=0/Mixed=1/Fav=2/Full=3; two columns, Continent | Status)
#   2026-11-11: the three-class achievement dimension is replaced by the four A.5 classes, colours identical to Fig 2 / Extended Data Fig. 1
# Run from the repository root:  Rscript fig_code/sup_fig2.R
# =====================================================================
source("fig_code/fig_theme.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

e6ab <- read_csv("data/fig_data/sup_fig2/supfig2_achv_by_group.csv", show_col_types = FALSE)
e6c  <- read_csv("data/fig_data/sup_fig2/supfig2_unac_by_layer_group.csv", show_col_types = FALSE)
e6d  <- read_csv("data/fig_data/sup_fig2/supfig2_kw_achv.csv", show_col_types = FALSE)

LAB_LV   <- c(paste0("T", 1:23), "A", "B", "C", "D")
cls_lv   <- c("Full achievement", "Favourable", "Mixed", "Unfavourable")
cls_cols <- c(`Full achievement` = "#00669E", `Favourable` = "#56B4E9",
              `Mixed` = "#E69F00", `Unfavourable` = "#C1441E")
GRP_LAB  <- c(continent = "Continent", status = "Status")

# Group order (shared by a-d): the continent block first, then the status block; within a block, descending target-level Unac share
gord <- e6ab %>%
  filter(scope == "target") %>%
  arrange(match(grouping, c("continent", "status")), desc(share_unfav), grp) %>%
  pull(grp)
lab_bar <- e6ab %>% filter(scope == "target") %>%
  { setNames(paste0(.$grp, " (n=", .$n_countries, ")"), .$grp) }

# ---- (a, b) Pooled share stacks ----
mk_bars <- function(sc, xlab, show_y, tag) {
  d <- e6ab %>% filter(scope == sc) %>%
    pivot_longer(c(share_full, share_fav, share_mixed, share_unfav),
                 names_to = "class", values_to = "share") %>%
    mutate(class = factor(
      c(share_full = "Full achievement", share_fav = "Favourable",
        share_mixed = "Mixed", share_unfav = "Unfavourable")[class],
      levels = cls_lv))
  ggplot(d, aes(share, grp, fill = class)) +
    geom_col(width = 0.78, colour = "white", linewidth = 0.15) +
    geom_hline(yintercept = 6.5, colour = "grey75", linewidth = 0.25) +
    scale_y_discrete(limits = rev(gord), labels = lab_bar) +
    scale_fill_manual(values = cls_cols, name = "Achievement class",
                      labels = cls_lv) +
    scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"),
                       expand = expansion(mult = 0.01)) +
    labs(x = xlab, y = NULL, tag = tag) +
    theme_nature() +
    theme(axis.text.y = if (show_y) element_text(size = 5) else element_blank(),
          axis.ticks.y = if (show_y) element_line(linewidth = 0.25) else element_blank(),
          legend.position = if (show_y) "none" else "right")
}
pa <- mk_bars("target", "Share of country–target pairs", TRUE,  "a")
pb <- mk_bars("goal",   "Share of country–goal pairs",   FALSE, "b")

# ---- (c) Per-layer %Unac heatmap ----
pc_hm <- e6c %>%
  mutate(layer_lab = factor(layer_id, levels = rev(LAB_LV))) %>%
  ggplot(aes(grp, layer_lab)) +
  geom_tile(aes(fill = share_unfav), colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(n_countries > 0, n_countries, "")),
            size = 1.5, colour = "grey15") +
  geom_vline(xintercept = 5.5, colour = "grey75", linewidth = 0.25) +
  scale_x_discrete(limits = gord, position = "top") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient(low = "#F7F7F7", high = "#B2182B",
                      labels = function(x) paste0(round(100 * x), "%"),
                      na.value = "white", name = "Unfav. share") +
  labs(x = NULL, y = NULL, tag = "c") +
  theme_nature() +
  theme(axis.text.x = element_text(size = 5, angle = 90, hjust = 0, vjust = 0.5),
        axis.ticks = element_blank(),
        legend.key.width = unit(3, "mm"), legend.key.height = unit(8, "mm"))

# ---- (d) Per-layer KW significance (achievement order 0/1/2/3; same p colour classes as Sup Fig. 1) ----
pkat <- function(p) factor(
  ifelse(is.na(p), "n.s.",
  ifelse(p < 0.001, "p < 0.001",
  ifelse(p < 0.01, "p < 0.01",
  ifelse(p < 0.05, "p < 0.05", "n.s.")))),
  levels = c("n.s.", "p < 0.05", "p < 0.01", "p < 0.001"))
P_COLS <- c("n.s." = "#D9D9D9", "p < 0.05" = "#F0B0A8",
            "p < 0.01" = "#D66060", "p < 0.001" = "#B2182B")

pd <- e6d %>%
  mutate(layer_lab = factor(layer_id, levels = rev(LAB_LV)),
         glab = factor(GRP_LAB[grouping], levels = unname(GRP_LAB))) %>%
  ggplot(aes(glab, layer_lab)) +
  geom_tile(aes(fill = pkat(kw_p)), colour = "white", linewidth = 0.3) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = P_COLS, na.value = "white", guide = "none") +
  labs(x = NULL, y = NULL, tag = "d") +
  theme_nature() +
  theme(axis.text.x = element_text(size = 5, angle = 90, hjust = 0, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

fig <- wrap_plots(
  pa, pb,
  wrap_plots(pc_hm, pd, nrow = 1, widths = c(1, 0.22)),
  nrow = 3, heights = c(1, 1, 2.35))
save_nature(fig, "supfig2_projected_achievement_gaps", 183, 210)
