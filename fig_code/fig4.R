# =====================================================================
# Fig 4 — Trend alignment of 41 indicators globally and country synergy contribution (CSC) during 2011–2024
#   Terminology aligned with the caption and the main text: gross synergy (= SyE = S_t, SI Formula 7;
#   SyD(8)=sgn(S_t), CSC(9)); axis labels use the terms only, without formula numbers (final, 2026-11-20).
#   (a) Annual range of the aligned synergy effect per indicator (dot-and-whisker: min-max bar + mean dot)
#   (b) Percentage of years with positive gross synergy (lollipop: thin stem + endpoint, 50% threshold line)
#   2026-11-11 presentation v2 (axis truncation; row labels carry only the serial number)
#   2026-11-20 presentation v3:
#     the DPSIR-class colouring is dropped -> one common colour #0072B2; row separation uses
#       zebra shading of alternating rows (odd serial numbers on a light-grey band); the row
#       labels remain the serial numbers (full names in Supplementary Table 1)
#     rows are ordered by increasing indicator serial number (top -> bottom)
#     the a axis title gains "Indicator"; the a/b axis labels drop the formulas
#     (a) x axis truncated at +/-3.3, out-of-range whisker ends carry open arrowheads (same
#     convention as Fig 1; 10 of the 41 indicators have 17 whisker ends out of range, the full
#     values are in the source data and are not annotated in the figure);
#     (b) filled = above the 50% line, open = at or below (shape as the second channel)
#   map c: countries not included are grey with a "Not included" legend key (same convention
#     as Figs 2/3)
# Run from the repository root:  Rscript fig_code/fig4.R
# =====================================================================
source("fig_code/fig_theme.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
  library(sf); library(rnaturalearth)
})

COL <- "#0072B2"          # common indicator colour (replacing the DPSIR-class colouring)
BAND <- "grey96"          # zebra banding of alternating rows

f4a <- read_csv("data/fig_data/fig4/fig4a_ga_ranges.csv", show_col_types = FALSE)
f4b <- read_csv("data/fig_data/fig4/fig4b_share_positive_years.csv", show_col_types = FALSE)
f4c <- read_csv("data/fig_data/fig4/fig4c_csc_by_country.csv", show_col_types = FALSE)

# ---- shared row order: increasing indicator serial number (top -> bottom) ----
# y uses numeric positions with a breaks/labels mapping (in ggplot2 4.x geom_tile and a
# discrete y axis are incompatible and would scramble the axis order; a numeric axis
# avoids this entirely; indicator 1 is at the top)
ids_asc <- sort(unique(f4a$indicator))
n_id <- length(ids_asc)
d <- f4a %>% inner_join(f4b, by = "indicator") %>%
  arrange(indicator) %>%
  mutate(ypos = n_id + 1 - match(indicator, ids_asc),
         above50 = share_positive > 0.5)
# banded rows = odd serial numbers (row order ascending -> regular alternation)
band_df <- d %>% filter(indicator %% 2 == 1) %>% distinct(ypos)

# ---- (a) axis truncation (same convention as Fig 1): body within +/-3.3,
#          out-of-range whisker ends marked with open arrowheads ----
CLIP_LO <- -3.3
CLIP_HI <- 3.3
XLIM    <- c(-3.55, 3.62)
da <- d %>%
  mutate(x_lo = pmax(ga_min, CLIP_LO), x_hi = pmin(ga_max, CLIP_HI),
         trunc_lo = ga_min < CLIP_LO, trunc_hi = ga_max > CLIP_HI)
arr_a <- bind_rows(
  da %>% filter(trunc_lo) %>% transmute(ypos, x = CLIP_LO + 0.22, xend = CLIP_LO),
  da %>% filter(trunc_hi) %>% transmute(ypos, x = CLIP_HI - 0.22, xend = CLIP_HI)
)

# identical y expand in both panels keeps the 41 rows aligned
y_exp <- expansion(mult = c(0.012, 0.012))

pa <- ggplot(da, aes(y = ypos)) +
  geom_tile(data = band_df, aes(x = mean(XLIM), y = ypos),
            width = diff(XLIM), height = 1,
            inherit.aes = FALSE, fill = BAND, colour = NA) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey60") +
  geom_segment(aes(x = x_lo, xend = x_hi, yend = ypos),
               colour = COL, linewidth = 0.75, alpha = 0.7, lineend = "butt") +
  geom_segment(data = arr_a,
               aes(x = x, xend = xend, y = ypos, yend = ypos),
               inherit.aes = FALSE, linewidth = 0.3, colour = "grey25",
               arrow = arrow(type = "open", length = unit(1.1, "mm"))) +
  geom_point(aes(x = ga_mean), colour = COL, size = 1.3) +
  scale_x_continuous(limits = XLIM, breaks = seq(-3, 3, 1),
                     expand = expansion(mult = 0.01)) +
  scale_y_continuous(breaks = rev(seq_along(ids_asc)), labels = ids_asc,
                     expand = y_exp) +
  labs(x = "Gross synergy", y = "Indicator", tag = "a") +
  theme_nature() +
  theme(axis.text.y = element_text(size = 5.5),
        legend.position = "none")

# ---------- (b) lollipop + 50% threshold line; shape = above the line or not ----------
pb <- ggplot(d, aes(y = ypos, x = share_positive)) +
  geom_tile(data = band_df, aes(x = 0.51, y = ypos), width = 1.02, height = 1,
            inherit.aes = FALSE, fill = BAND, colour = NA) +
  geom_segment(aes(x = 0, xend = share_positive, yend = ypos),
               colour = COL, linewidth = 0.45, alpha = 0.7) +
  geom_point(aes(shape = above50), colour = COL, fill = "white",
             size = 1.35, stroke = 0.35) +
  geom_vline(xintercept = 0.5, linewidth = 0.35, linetype = "dashed",
             colour = "grey25") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     name = NULL, breaks = c("TRUE", "FALSE"),
                     labels = c(`TRUE` = "Above 50%", `FALSE` = "At or below 50%"),
                     guide = guide_legend(order = 1,
                       override.aes = list(size = 2.2, colour = COL,
                                           fill = "white", stroke = 0.35))) +
  scale_y_continuous(breaks = rev(seq_along(ids_asc)), labels = ids_asc,
                     expand = y_exp) +
  scale_x_continuous(limits = c(0, 1.02),
                     labels = function(x) paste0(round(100 * x), "%")) +
  labs(x = "Percentage of years with positive gross synergy",
       y = NULL, tag = "b") +
  theme_nature() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

# ---------- (c) CSC map (countries not included grey + legend key) ----------
world <- ne_countries(scale = "small", returnclass = "sf")
w <- world %>% left_join(f4c, by = c("adm0_a3" = "iso3"))
w_na <- w %>% filter(is.na(CSC_mean))
pc <- ggplot(w) +
  geom_sf(aes(fill = CSC_mean), colour = "grey85", linewidth = 0.08) +
  geom_sf(data = w_na, aes(colour = "Not included"), fill = "grey80",
          linewidth = 0.08) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, na.value = "grey80",
                       limits = c(-2, 2), oob = scales::squish,
                       name = "Mean CSC\n(capped ±2)") +
  scale_colour_manual(values = c(`Not included` = "grey85"), name = NULL,
                      guide = guide_legend(override.aes =
                        list(fill = "grey80", shape = 22, size = 4))) +
  coord_sf(ylim = c(-58, 84), expand = FALSE) +
  labs(tag = "c", x = NULL, y = NULL) +
  theme_nature() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(), axis.title = element_blank(),
        legend.box = "vertical",
        legend.key.width = unit(3, "mm"), legend.key.height = unit(8, "mm"))

# Layout: top row (a)|(b) share the row order, the shape key is collected at the bottom; (c) full width
pab <- (pa | pb) +
  plot_layout(widths = c(1.02, 0.98), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal")
fig <- wrap_plots(pab, pc, ncol = 1, heights = c(1, 0.82))
save_nature(fig, "fig4_trend_alignment_csc", 183, 160)
