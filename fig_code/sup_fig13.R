# =====================================================================
# sup_fig13.R — Supplementary Fig. 13
#   CSC null-model distributions (Null A exchangeability / Null B sign symmetry), observed values as vertical lines + p annotation
#   (a) share_positive (share of countries with a positive contribution, observed 47%)
#   (b) top5_share (share of the positive contributions contributed by the top five, observed 62%)
#   (c) max_share (share of the largest single country, observed 43%)
# Data: run fig_code/csc_null.R first
# Run from the repository root:  Rscript fig_code/sup_fig13.R
# =====================================================================
source("fig_code/fig_theme.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr)
})

d13 <- "data/fig_data/sup_fig13"
drw <- read_csv(file.path(d13, "fig13_null_draws.csv"), show_col_types = FALSE)
smr <- read_csv(file.path(d13, "fig13_null_summary.csv"), show_col_types = FALSE)

MODEL_COLS <- c(A_exchangeability = "#0072B2", B_sign_symmetry = "#D55E00")
MODEL_LABS <- c(A_exchangeability = "Exchangeability",
                B_sign_symmetry = "Sign symmetry",
                C_binomial = "Binomial (n = 130)")

fmt_p <- function(p) ifelse(p < 0.001, "P < 0.001", sprintf("P = %.2f", p))

mk_panel <- function(stat_key, tag, xlim = NULL) {
  d <- drw %>% filter(stat == stat_key) %>%
    mutate(model = factor(model, levels = names(MODEL_COLS)))
  s <- smr %>% filter(stat == stat_key)
  obs_v <- s$observed[1]
  p_lab <- s %>% mutate(lab = sprintf("%s: %s", MODEL_LABS[model], fmt_p(p_value)))
  ann <- paste(p_lab$lab, collapse = "\n")
  p <- ggplot(d, aes(value, fill = model, colour = model)) +
    geom_histogram(bins = 40, linewidth = 0.2,
                   position = "identity", alpha = 0.35) +
    geom_vline(xintercept = obs_v, linewidth = 0.45, colour = "black") +
    annotate("text", x = obs_v, y = Inf, vjust = 1.6, hjust = -0.08,
             size = 2, colour = "black",
             label = sprintf("observed %.0f%%", 100 * obs_v)) +
    annotate("text", x = -Inf, y = Inf, vjust = 1.6, hjust = -0.05,
             size = 1.8, colour = "grey15", label = ann) +
    scale_fill_manual(values = MODEL_COLS, labels = MODEL_LABS, name = "Null model") +
    scale_colour_manual(values = MODEL_COLS, labels = MODEL_LABS, name = "Null model") +
    scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"),
                       limits = xlim) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
    labs(x = NULL, y = "Permutations", tag = tag) +
    theme_nature() +
    theme(legend.position = "none")
  p
}

pa <- mk_panel("share_positive", "a", xlim = c(-0.05, 1.02))
pb <- mk_panel("top5_share", "b", xlim = c(-0.05, 1.02))
pc <- mk_panel("max_share", "c", xlim = c(-0.05, 1.02))

fig <- (pa | pb | pc) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
save_nature(fig, "supfig13_csc_null", 183, 62)
