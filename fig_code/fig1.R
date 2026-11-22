# =====================================================================
# Fig 1 — Baseline across the GBF's goals and targets (2011–2020)
#   (a) Slopes of the 23 targets with 95% CI (point size encodes R2, colour encodes section, filled = significant)
#   (b) Global baseline trend lines for goals A-D (black/grey x solid/dashed, labelled directly in black at the line ends)
#   (c) Slope x R2 quadrant (23 targets only, all point labels black, thin grey leader lines)
# Outliers: x axis clipped to [-1.05, 2.42]; T6 (slope 4.27, CI entirely outside) and T14 (both CI ends outside)
#   are marked with open truncation arrows, and the T6 value 4.27 is printed to the left of the arrow (truncation noted in the caption).
# Legend: the three groups Section / Trend significance / R2(a) are collected together and placed in a row at the bottom.
# Colours: three Okabe-Ito section colours (blue = reducing threats / green = meeting needs / orange = tools and solutions, shared by a and c);
#   the four goal lines use a black/grey x line-type scheme so that goal colours are semantically separated from section colours (final, 2026-11-11).
# Text: all black (Nature: avoid coloured text); labels 5-7 pt, panel letters 8 pt bold.
# Run from the repository root:  Rscript fig_code/fig1.R
# =====================================================================
source("fig_code/fig_theme.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggrepel)
})

f1a <- read_csv("data/fig_data/fig1/fig1a_target_slopes.csv", show_col_types = FALSE)
f1b <- read_csv("data/fig_data/fig1/fig1b_goal_trends.csv", show_col_types = FALSE)
f1c <- read_csv("data/fig_data/fig1/fig1c_quadrant.csv", show_col_types = FALSE)

# ---- Axis truncation parameters (measured on the 2026-11-11 data: the body spans -0.29 to 1.47, only T6/T14 fall outside) ----
CLIP_LO <- -1.0    # truncation bound (where the arrows sit)
CLIP_HI <- 2.25
XLIM    <- c(-1.05, 2.42)

SEC3 <- c("reducing_threats", "meeting_needs", "tools_solutions")
# Both panels share one Section scale definition (identical values/labels/names) so that patchwork can merge the legends
sec_scale <- function(...) {
  scale_colour_manual(values = SEC_COLORS[SEC3], labels = SEC_LABELS[SEC3],
                      name = "Cluster", breaks = SEC3, ...)
}

f1a <- f1a %>%
  mutate(order = as.integer(sub("T", "", layer_id)),
         label = factor(layer_id, levels = rev(paste0("T", 1:23))),
         section = factor(section, levels = SEC3),
         # error-bar ends after truncation; if the whole interval lies outside, set NA and draw no bar
         x_lo = ifelse(ci_lo > CLIP_HI | ci_hi < CLIP_LO, NA_real_,
                       pmax(ci_lo, CLIP_LO)),
         x_hi = ifelse(ci_lo > CLIP_HI | ci_hi < CLIP_LO, NA_real_,
                       pmin(ci_hi, CLIP_HI)),
         trunc_lo = ci_lo < CLIP_LO,
         trunc_hi = ci_hi > CLIP_HI,
         point_out = slope < CLIP_LO | slope > CLIP_HI)

# Truncation-arrow layer: one open arrow per out-of-range end (pointing away from the bounds)
trunc_arrows <- bind_rows(
  f1a %>% filter(trunc_lo) %>% transmute(label, x = CLIP_LO + 0.12, xend = CLIP_LO),
  f1a %>% filter(trunc_hi | point_out) %>%
    transmute(label, x = CLIP_HI - 0.12, xend = CLIP_HI)
) %>% distinct()

# T6 truncation value annotation (data frame plus factor y, so that annotate does not scramble the axis order)
t6_note <- tibble(x = CLIP_HI - 0.20,
                  label = factor("T6", levels = levels(f1a$label)),
                  txt = "4.27")

# (a) slope dot plot
pa <- ggplot(f1a, aes(slope, label, colour = section)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey60") +
  geom_errorbarh(data = f1a %>% filter(!is.na(x_lo)),
                 aes(xmin = x_lo, xmax = x_hi), height = 0, linewidth = 0.3) +
  geom_point(data = f1a %>% filter(!point_out),
             aes(size = r2, shape = significant)) +
  geom_segment(data = trunc_arrows,
               aes(x = x, xend = xend, y = label, yend = label),
               inherit.aes = FALSE, linewidth = 0.3, colour = "grey25",
               arrow = arrow(type = "open", length = unit(1.1, "mm"))) +
  geom_text(data = t6_note, aes(x = x, y = label, label = txt),
            inherit.aes = FALSE, hjust = 1, size = 2, colour = "black") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     name = "Trend significance",
                     labels = c(`TRUE` = "significant", `FALSE` = "not significant"),
                     guide = guide_legend(order = 2,
                                          override.aes = list(size = 2))) +
  scale_size_continuous(range = c(0.7, 4.5), name = quote(R^2 ~ "(a)"),
                        guide = guide_legend(order = 3,
                                             override.aes = list(colour = "black",
                                                                 shape = 16))) +
  sec_scale(guide = guide_legend(order = 1,
                                 override.aes = list(size = 2))) +
  scale_y_discrete(limits = rev(paste0("T", 1:23))) +
  scale_x_continuous(limits = XLIM, breaks = seq(-1, 2, 0.5),
                     expand = expansion(mult = 0.02)) +
  labs(x = "Trend slope", y = "Target", tag = "a") +
  theme_nature() +
  theme(axis.text.y = element_text(size = 5.5))

# (b) goal trend lines - black/grey x solid/dashed, labelled directly in black at the line ends (goal colours separated from section colours)
GOAL_LTY  <- c(A = "solid",  B = "22",     C = "solid",  D = "22")
f1b <- f1b %>% mutate(Goal = factor(layer_id, levels = c("A", "B", "C", "D")))
ends <- f1b %>% group_by(Goal) %>%
  slice_max(year, n = 1) %>%
  mutate(lab = gsub("-", "−", sprintf("%s %+.2f", Goal, slope), fixed = TRUE))
pb <- ggplot(f1b, aes(year, index, colour = Goal, linetype = Goal)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 0.7) +
  geom_text(data = ends, aes(label = lab),
            colour = "black", hjust = -0.65, size = 2, show.legend = FALSE) +
  scale_colour_manual(values = GOAL_COLORS, guide = "none") +
  scale_linetype_manual(values = GOAL_LTY, guide = "none") +
  scale_x_continuous(breaks = seq(2011, 2020, 3), limits = c(2011, 2022.5),
                     expand = c(0.02, 0)) +
  labs(x = "Year", y = "Goal index (0–100)", tag = "b") +
  theme_nature() +
  theme(legend.position = "none")

# (c) quadrant (23 targets only; black labels with thin grey leader lines; T6 truncation arrows follow the panel a convention)
f1c <- f1c %>%
  filter(grepl("^T[0-9]+$", layer_id)) %>%
  mutate(section = factor(section, levels = c("reducing_threats", "meeting_needs",
                                              "tools_solutions")),
         lab = layer_id,
         point_out = slope < CLIP_LO | slope > CLIP_HI)
pc <- ggplot(f1c, aes(slope, r2, colour = section)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.25, colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey70") +
  geom_point(data = f1c %>% filter(!point_out),
             aes(shape = significant), size = 2, alpha = 0.9) +
  geom_segment(data = f1c %>% filter(point_out),
               aes(x = CLIP_HI - 0.12, xend = CLIP_HI, y = r2, yend = r2),
               inherit.aes = FALSE, linewidth = 0.3, colour = "grey25",
               arrow = arrow(type = "open", length = unit(1.1, "mm"))) +
  geom_text(data = f1c %>% filter(point_out),
            aes(y = r2, label = lab), inherit.aes = FALSE,
            x = CLIP_HI - 0.20, hjust = 1, size = 2, colour = "black") +
  geom_text_repel(data = f1c %>% filter(!point_out),
                  aes(label = lab), size = 2, colour = "black",
                  show.legend = FALSE,
                  max.overlaps = Inf, box.padding = 0.5, point.padding = 0.4,
                  max.iter = 3000, force = 3, seed = 42,
                  min.segment.length = 0, segment.colour = "grey70",
                  segment.size = 0.15) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     name = "Trend significance",
                     labels = c(`TRUE` = "significant", `FALSE` = "not significant"),
                     guide = "none") +
  sec_scale(guide = "none") +
  scale_x_continuous(limits = XLIM, breaks = seq(-1, 2, 0.5),
                     expand = expansion(mult = 0.02)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(-0.06, 1.12)) +
  labs(x = "Trend slope", y = quote(R^2), tag = "c") +
  theme_nature()

# Layout: a on the left, b/c on the right; the legend has three groups (Section / Trend significance / R2(a))
# collected into a single row at the bottom, so that panel c is no longer squeezed by the legend.
fig <- (pa | (pb / pc)) +
  plot_layout(widths = c(1, 1.25), guides = "collect") &
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.direction = "horizontal",
        legend.spacing.y = unit(0.4, "mm"))
save_nature(fig, "fig1_baseline_trends", 183, 166)
