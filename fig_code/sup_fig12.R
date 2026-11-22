# =====================================================================
# Supplementary Fig. 12 - representative model fits (one example per model form, a-g)
#   Title = indicator number + selected form; subtitle = indicator name
#   Points = observations (global mean series), lines = fits, ribbons = 95% CI,
#   grey dashed line = the 2020 boundary, green dashed line = the target year 2030; a common x axis from 2000 to 2030
#   The legend occupies the 8th cell on its own (to the right of panel g)
# Run from the repository root:  Rscript fig_code/sup_fig12.R
# =====================================================================
source("fig_code/fig_theme.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
})

e4 <- read_csv("data/fig_data/sup_fig12/supfig12_model_fits.csv", show_col_types = FALSE) %>%
  filter(year >= 2000)   # truncated before 2000
FORM_LAB <- c(poly1 = "linear", power_law = "power law", quad = "quadratic",
              asymptotic = "asymptotic", exponential = "exponential",
              segmented = "segmented", DLM = "DLM (high-noise track)")

mk_fit <- function(id, tag, x_axis = FALSE, y_axis = FALSE, legend = FALSE) {
  d <- e4 %>% filter(indicator == id)
  nm <- d$indicator_name[1]
  form <- unname(FORM_LAB[d$model_form[1]])
  ggplot(d, aes(year)) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), fill = "#56B4E9", alpha = 0.25) +
    geom_line(aes(y = fitted, linetype = segment), linewidth = 0.45,
              colour = "#B2182B") +
    geom_point(aes(y = observed), size = 0.5, na.rm = TRUE) +
    geom_vline(xintercept = 2020, linetype = "dashed", linewidth = 0.25,
               colour = "grey55") +
    geom_vline(xintercept = 2030, linetype = "dashed", linewidth = 0.25,
               colour = "#009E73") +
    scale_linetype_manual(values = c(observed = "solid", fitted = "solid",
                                     extrapolation = "22"),
                          breaks = "extrapolation",
                          labels = "Extrapolation", name = NULL) +
    scale_x_continuous(limits = c(2000, 2030), breaks = seq(2000, 2030, 10)) +
    labs(x = if (x_axis) "Year" else NULL,
         y = if (y_axis) "Normalized index (0–100)" else NULL,
         tag = tag,
         title = paste0("Indicator ", id, " — ", form),
         subtitle = nm) +
    theme_nature() +
    theme(axis.text.x = element_text(size = 5),
          axis.ticks.x = element_line(linewidth = 0.25),
          axis.text.y = element_text(size = 5),
          axis.ticks.y = element_line(linewidth = 0.25),
          axis.title = element_text(size = 5.5),
          plot.title = element_text(size = 5.2, hjust = 0),
          plot.subtitle = element_text(size = 3.8, face = "plain",
                                       colour = "grey25", hjust = 0),
          legend.position = if (legend) "bottom" else "none",
          legend.direction = "horizontal",
          legend.text = element_text(size = 5),
          legend.key.width = unit(3.8, "mm"))
}

pa <- mk_fit(11L, "a", y_axis = TRUE)
pb <- mk_fit(13L, "b")
pc <- mk_fit(39L, "c")
pd <- mk_fit(21L, "d")
pe <- mk_fit(49L, "e", x_axis = TRUE, y_axis = TRUE)
pf <- mk_fit(3L,  "f", x_axis = TRUE)
pg <- mk_fit(14L, "g", x_axis = TRUE)

# The legend occupies the 8th cell on its own (to the right of panel g): all graphical elements are drawn manually
pleg <- ggplot() +
  annotate("point", x = 0.16, y = 0.92, size = 0.8) +
  annotate("text", x = 0.28, y = 0.92, hjust = 0, size = 1.9,
           label = "Observations") +
  annotate("segment", x = 0.05, xend = 0.27, y = 0.78, yend = 0.78,
           colour = "#B2182B", linewidth = 0.5) +
  annotate("text", x = 0.28, y = 0.78, hjust = 0, size = 1.9,
           label = "Fitted") +
  annotate("rect", xmin = 0.05, xmax = 0.27, ymin = 0.605, ymax = 0.655,
           fill = "#56B4E9", alpha = 0.35) +
  annotate("text", x = 0.28, y = 0.63, hjust = 0, size = 1.9,
           label = "95% CI") +
  annotate("segment", x = 0.05, xend = 0.27, y = 0.50, yend = 0.50,
           linetype = "22", colour = "#B2182B", linewidth = 0.5) +
  annotate("text", x = 0.28, y = 0.50, hjust = 0, size = 1.9,
           label = "Extrapolation") +
  annotate("segment", x = 0.05, xend = 0.27, y = 0.36, yend = 0.36,
           linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  annotate("text", x = 0.28, y = 0.36, hjust = 0, size = 1.9,
           label = "2020 boundary") +
  annotate("segment", x = 0.05, xend = 0.27, y = 0.22, yend = 0.22,
           linetype = "dashed", colour = "#009E73", linewidth = 0.3) +
  annotate("text", x = 0.28, y = 0.22, hjust = 0, size = 1.9,
           label = "Target year (2030)") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_void()

fig <- wrap_plots(pa, pb, pc, pd, pe, pf, pg, pleg, ncol = 4)
save_nature(fig, "supfig12_model_fits", 183, 100)
