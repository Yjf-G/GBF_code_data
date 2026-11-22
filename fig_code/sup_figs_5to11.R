# =====================================================================
# Supplementary Figs. 5-11 - one panel per uncertainty analysis
#   Sup Fig. 5: distribution of the country-level composite value over k-subsets (Algorithm 1)
#   Sup Fig. 6: bootstrap imputation uncertainty (Algorithm 2)
#   Sup Fig. 7: time-series length uncertainty (Algorithm 3)
#   Sup Fig. 9: single-indicator perturbation sensitivity |Sx| (Algorithm 5)
# =====================================================================
source("fig_code/fig_theme.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

s2 <- read_csv("data/fig_data/sup_fig6/supfig6_gap_length.csv", show_col_types = FALSE)
s2s <- read_csv("data/fig_data/sup_fig6/supfig6_sensitivity.csv", show_col_types = FALSE)
s3 <- read_csv("data/fig_data/sup_fig7/supfig7_window_length.csv", show_col_types = FALSE)
s3c <- read_csv("data/fig_data/sup_fig7/supfig7_sigma_curves.csv", show_col_types = FALSE)
s3x <- read_csv("data/fig_data/sup_fig7/supfig7_cross_model.csv", show_col_types = FALSE)
s4 <- read_csv("data/fig_data/sup_fig9/supfig9_perturbation_sensitivity.csv", show_col_types = FALSE)

# Sup Fig. 5 - targets with >=2 indicators plus the goals (Goal A/B/D; goal C has only 1 indicator and is not shown)
#   one composite value score_c per k-combination (mean across the units); the box plot gives the distribution of the combination values for that k;
#   k = n has only one combination -> a single horizontal line (box width), with a red dashed line connecting the medians across k
s1 <- read_csv("data/fig_data/sup_fig5/supfig5_indicator_subset.csv", show_col_types = FALSE)
s1 <- s1 %>%
  mutate(layer_lab = ifelse(grepl("^T", layer_id), layer_id, paste("Goal", layer_id)),
         layer_lab = factor(layer_lab, levels = c(
           intersect(paste0("T", 1:23), unique(as.character(layer_lab))),
           intersect(paste("Goal", LETTERS[1:4]), unique(as.character(layer_lab)))))) %>%
  group_by(layer_id, k) %>% mutate(n_c = n()) %>% ungroup()
s1m <- s1 %>% group_by(layer_lab, k) %>%
  summarise(med = median(value), .groups = "drop")
p1 <- ggplot(s1, aes(factor(k), value)) +
  geom_boxplot(data = ~subset(.x, n_c >= 2), outlier.size = 0.18,
               linewidth = 0.22, fill = "#9CC2E5") +
  geom_errorbar(data = ~subset(.x, n_c == 1),
                aes(ymin = value, ymax = value), width = 0.75, linewidth = 0.22) +
  geom_line(data = s1m, aes(factor(k), med, group = layer_lab),
            colour = "red", linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~layer_lab, scales = "free_y") +
  labs(x = "Number of indicators combined (k)",
       y = "Composite score (target/goal values)") +
  theme_nature() +
  theme(axis.text.x = element_text(size = 5))
save_nature(p1, "supfig5_indicator_subset", 150, 120)

# Sup Fig. 6 - bootstrap imputation uncertainty (Algorithm 2):
#   (a) relative width of the 95% CI (width/IQR) as a function of the gap length g (filled = mean, hollow = median);
#   (b) sign stability as a function of g - the sign stability of the point estimate decreases for long gaps;
#   (c) sensitivity: relative width under different B / alpha / imputer / rolling-window CV choices
p2a <- ggplot(s2, aes(g, mean_ci_width_rel)) +
  geom_point(size = 1.1, colour = "#0072B2") +
  geom_point(aes(y = median_ci_width_rel), size = 0.9, shape = 21,
             fill = "white", colour = "#0072B2", stroke = 0.3) +
  scale_x_continuous(breaks = function(l) unique(floor(pretty(l)))) +
  labs(x = "Gap length g (years)", y = "95% CI width (relative to IQR)",
       tag = "a") +
  theme_nature()
p2b <- ggplot(s2, aes(g, mean_stability)) +
  geom_point(size = 1.1, colour = "#D55E00") +
  geom_hline(yintercept = 0.95, linetype = "dashed", linewidth = 0.25,
             colour = "grey50") +
  scale_x_continuous(breaks = function(l) unique(floor(pretty(l)))) +
  scale_y_continuous(limits = c(0.4, 1)) +
  labs(x = "Gap length g (years)", y = "Sign stability", tag = "b") +
  theme_nature()
p2c <- s2s %>%
  mutate(f = factor(paste(factor, level),
                    levels = rev(paste(factor, level)))) %>%
  ggplot(aes(f, mean_ci_width_rel)) +
  geom_point(size = 1.1, colour = "#E69F00") +
  coord_flip() +
  labs(x = NULL, y = "95% CI width (relative to IQR)", tag = "c") +
  theme_nature() +
  theme(axis.text.y = element_text(size = 4.5))
p2 <- wrap_plots(p2a, p2b, p2c, ncol = 3, widths = c(1, 1, 1.35))
save_nature(p2, "supfig6_imputation_bootstrap", 183, 62)

# Sup Fig. 7 - time-series length uncertainty (Algorithm 3):
#   (a) distribution of L_opt (Kneedle; histogram + KDE);
#   (b) relative dispersion of sigma_i(L) per indicator (sigma/sigma(L=5), thin grey line) and the median across indicators for each L (thick blue line) - it plateaus for long windows;
#   (c) cross-model consistency: the L_opt of ETS/ARIMA/RW (points, sorted by the pooled L_opt) and the per-indicator spread (vertical lines)
p3a <- ggplot(s3, aes(L_opt)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 2,
                 fill = "#9CC2E5", colour = "white", linewidth = 0.2) +
  geom_density(linewidth = 0.35, colour = "#0072B2") +
  labs(x = expression("Optimal window length L"[opt]*" (years)"),
       y = "Density", tag = "a") +
  theme_nature()
s3rel <- s3c %>%
  group_by(indicator) %>%
  mutate(sigma_rel = sigma / sigma[which.min(L)]) %>%
  ungroup()
s3med <- s3rel %>%
  group_by(L) %>%
  summarise(med = median(sigma_rel, na.rm = TRUE), n = n(), .groups = "drop") %>%
  filter(n >= 10)   # the median line is drawn only up to L values still covered by >=10 indicators, avoiding the single-sample artefact of long series
p3b <- ggplot(s3rel, aes(L, sigma_rel, group = indicator)) +
  geom_line(linewidth = 0.25, colour = "grey65") +
  geom_line(data = s3med, aes(L, med, group = NULL),
            linewidth = 0.6, colour = "#0072B2") +
  labs(x = "Window length L (years)",
       y = expression(sigma(L)*" / "*sigma(5)), tag = "b") +
  theme_nature()
ord <- s3 %>% arrange(L_opt) %>% pull(indicator)
s3span <- s3x %>%
  group_by(indicator) %>%
  summarise(ymin = min(L_opt_kneedle), ymax = max(L_opt_kneedle), .groups = "drop")
p3c <- ggplot(s3x, aes(factor(indicator, levels = ord), L_opt_kneedle)) +
  geom_linerange(data = s3span,
                 aes(factor(indicator, levels = ord), ymin = ymin, ymax = ymax),
                 inherit.aes = FALSE, linewidth = 0.25, colour = "grey55") +
  geom_point(aes(colour = model, shape = model), size = 1.35,
             position = position_dodge(0.7)) +
  scale_shape_manual(values = c(ETS = 16, ARIMA = 17, RW = 15), name = NULL,
                     breaks = c("ETS", "ARIMA", "RW")) +
  scale_colour_manual(values = c(ETS = "#0072B2", ARIMA = "#D55E00", RW = "#009E73"),
                      name = NULL, breaks = c("ETS", "ARIMA", "RW")) +
  labs(x = "Indicator", y = expression("Optimal window length L"[opt]*" (years)"),
       tag = "c") +
  theme_nature() +
  theme(axis.text.x = element_text(size = 3.2, angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "inside", legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1), legend.direction = "vertical")
p3 <- wrap_plots(p3a, p3b, p3c, ncol = 3, widths = c(1, 1, 1.45))
save_nature(p3, "supfig7_window_length", 183, 62)

# Sup Fig. 9
s4 <- s4 %>% group_by(indicator) %>%
  summarise(Sx_mean_abs = mean(Sx_mean_abs), .groups = "drop") %>%
  arrange(desc(Sx_mean_abs)) %>%
  mutate(ind = factor(indicator, levels = unique(indicator)))
p4 <- ggplot(s4, aes(ind, Sx_mean_abs)) +
  geom_col(width = 0.72, fill = "#E69F00") +
  geom_hline(yintercept = 0.2, linetype = "dashed", linewidth = 0.25,
             colour = "#D55E00") +
  labs(x = "Indicator", y = expression("Sensitivity |S"[x]*"|")) +
  theme_nature() +
  theme(axis.text.x = element_text(size = 4.2, angle = 90, hjust = 1, vjust = 0.5))
save_nature(p4, "supfig9_perturbation", 150, 70)

# =====================================================================
# Sup Fig. 8 - country-sample uncertainty (the three phases of Algorithm 4):
#   (a) family of sigma_i(n) curves per indicator + median; vertical line = the actual panel of 130 countries;
#   (b) share of years stable within +/-5% in the full-size bootstrap (the n=30 comparison value is given in the caption);
#   (c) country bootstrap distribution of the net MIS change (goal level)
a6s <- read_csv("data/fig_data/sup_fig8/supfig8_sigma_curves.csv", show_col_types = FALSE)
a6p <- read_csv("data/fig_data/sup_fig8/supfig8_prop_stable.csv", show_col_types = FALSE)
a6h <- read_csv("data/fig_data/sup_fig8/supfig8_net_change_bootstrap.csv", show_col_types = FALSE)
a6hs <- read_csv("data/fig_data/sup_fig8/supfig8_net_change_summary.csv", show_col_types = FALSE)
a6m <- a6s %>% group_by(n) %>% summarise(sigma = median(sigma), .groups = "drop")
p6a <- ggplot(a6s, aes(n, sigma, group = factor(indicator))) +
  geom_line(linewidth = 0.3, colour = "grey55", alpha = 0.7) +
  geom_line(data = a6m, aes(n, sigma, group = 1),
            linewidth = 0.55, colour = "#0072B2") +
  geom_vline(xintercept = 130, linetype = "dashed", linewidth = 0.25,
             colour = "#D55E00") +
  labs(x = "Number of countries sampled", y = "Sampling SD of the indicator mean",
       tag = "a") +
  theme_nature()
p6b <- a6p %>%
  arrange(prop_stable) %>%
  mutate(ind = factor(indicator, levels = indicator)) %>%
  ggplot(aes(ind, prop_stable)) +
  geom_col(width = 0.72, fill = "#E69F00") +
  geom_hline(yintercept = 0.8, linetype = "dashed", linewidth = 0.25,
             colour = "#D55E00") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Indicator", y = "Proportion of stable years", tag = "b") +
  theme_nature() +
  theme(axis.text.x = element_text(size = 3.2, angle = 90, hjust = 1, vjust = 0.5))
p6c <- ggplot(a6h, aes(net_change)) +
  geom_histogram(bins = 30, fill = "#9CC2E5", colour = "white", linewidth = 0.12) +
  geom_vline(xintercept = a6hs$full, linetype = "dashed", linewidth = 0.25,
             colour = "#D55E00") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.4, size = 1.9,
           label = sprintf("%.2f ± %.2f (B = %d)",
                           a6hs$full, a6hs$boot_sd, a6hs$B)) +
  labs(x = "Bootstrap net MIS change",
       y = "Resamples", tag = "c") +
  theme_nature()
p6 <- wrap_plots(p6a, p6b, p6c, ncol = 3, widths = c(1.15, 1.5, 1))
save_nature(p6, "supfig8_country_coverage", 183, 62)

# =====================================================================
# Sup Fig. 10 - weight sensitivity (SI B.2.2; keep_index 130 countries, target/goal two levels)
#   dirichlet_sd is used as the error bar; facet = target level / goal level (= the net change reported in the paper)
s7 <- read_csv("data/fig_data/sup_fig10/supfig10_weight_sensitivity.csv", show_col_types = FALSE)
s7_sd <- s7 %>%
  filter(scheme == "dirichlet_sd") %>%
  transmute(level, sd = net_change)
s7b <- s7 %>%
  filter(scheme != "dirichlet_sd") %>%
  left_join(s7_sd, by = "level") %>%
  mutate(scheme_f = factor(
    recode(scheme,
           "equal" = "Equal",
           "dirichlet_mean" = "Dirichlet (mean ± SD)",
           "PCA_first_component" = "PCA (PC1)",
           "equal_excluding_BCH46" = "Equal excl. BCH 46"),
    levels = c("Equal", "Dirichlet (mean ± SD)", "PCA (PC1)", "Equal excl. BCH 46")),
    level_f = factor(level, levels = c("target", "goal"),
                     labels = c("Target level", "Goal level (headline)")))
p7a <- ggplot(filter(s7b, level == "target"), aes(scheme_f, net_change)) +
  geom_col(width = 0.68, fill = "#9CC2E5") +
  geom_errorbar(data = filter(s7b, level == "target", scheme == "dirichlet_mean"),
                aes(ymin = net_change - sd, ymax = net_change + sd),
                width = 0.16, linewidth = 0.35) +
  labs(x = "Weighting scheme", y = "Target-level MIS change 2011–2024", tag = "a") +
  theme_nature() +
  theme(axis.text.x = element_text(size = 5, angle = 22, hjust = 1))
p7b <- ggplot(filter(s7b, level == "goal"), aes(scheme_f, net_change)) +
  geom_col(width = 0.68, fill = "#9CC2E5") +
  geom_errorbar(data = filter(s7b, level == "goal", scheme == "dirichlet_mean"),
                aes(ymin = net_change - sd, ymax = net_change + sd),
                width = 0.16, linewidth = 0.35) +
  labs(x = "Weighting scheme", y = "Goal-level MIS change 2011–2024 (net change)", tag = "b") +
  theme_nature() +
  theme(axis.text.x = element_text(size = 5, angle = 22, hjust = 1))
p7 <- wrap_plots(p7a, p7b, ncol = 2)
save_nature(p7, "supfig10_weight_sensitivity", 150, 62)

# =====================================================================
# Sup Fig. 11 - endpoint robustness (SI B.2.3; the truncation sensitivity of B.2.4 is computed separately):
#   (a) effect of leaving out one year on the net change, Delta = recomputed - full data (internal years are always 0,
#       only the 2024 endpoint gives -0.12/-0.04); (b) distribution of the absolute extrapolation deviation per country-target
#       (solid line = median 0.65 / dashed line = mean 4.41; the x axis is cut at the 95% quantile)
s8 <- read_csv("data/fig_data/sup_fig11/supfig11_leave_one_year_out.csv", show_col_types = FALSE)
s8e <- read_csv("data/fig_data/sup_fig11/supfig11_extrapolation_test.csv", show_col_types = FALSE)
s8s <- read_csv("data/fig_data/sup_fig11/supfig11_extrapolation_series.csv", show_col_types = FALSE)
s8l <- s8 %>%
  transmute(omitted_year, `MIS` = dMIS, `EIS` = dEIS) %>%
  pivot_longer(-omitted_year, names_to = "index", values_to = "value") %>%
  mutate(index = factor(index, levels = c("MIS", "EIS")))
s8d <- s8l %>% filter(omitted_year == 2024) %>%
  mutate(lab = sprintf("Δ = %.2f", value))
p8a <- ggplot(s8l, aes(omitted_year, value)) +
  geom_hline(yintercept = 0,
             linetype = "dashed", linewidth = 0.25, colour = "grey45") +
  geom_segment(aes(xend = omitted_year, y = 0, yend = value),
               linewidth = 0.5, colour = "#0072B2") +
  geom_point(size = 0.9, colour = "#0072B2") +
  geom_text(data = s8d, aes(label = lab), size = 1.8, vjust = -0.6, hjust = 1.12) +
  facet_wrap(~index, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.30))) +
  labs(x = "Omitted year", y = "Change in net change vs full data (Δ)", tag = "a") +
  theme_nature()
p8b <- ggplot(s8s, aes(abs_deviation)) +
  geom_histogram(bins = 45, fill = "#E69F00", colour = "white", linewidth = 0.12) +
  geom_vline(xintercept = s8e$mean_absolute_deviation,
             linetype = "dashed", linewidth = 0.25, colour = "#D55E00") +
  geom_vline(xintercept = median(s8s$abs_deviation),
             linetype = "solid", linewidth = 0.25, colour = "black") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.4, size = 1.9,
           label = sprintf("median %.2f (solid) · mean %.2f (dashed)\nSD %.2f; n = %s",
                           median(s8s$abs_deviation),
                           s8e$mean_absolute_deviation, s8e$sd,
                           format(s8e$n, big.mark = ","))) +
  scale_x_continuous(limits = c(0, unname(quantile(s8s$abs_deviation, 0.95)))) +
  labs(x = "2024 deviation from 2011–2023\nlinear-trend projection (index points)",
       y = "Country–target series", tag = "b") +
  theme_nature()
p8 <- wrap_plots(p8a, p8b, ncol = 2, widths = c(1.25, 1))
save_nature(p8, "supfig11_endpoint_robustness", 150, 62)

message("Supplementary figures done")
