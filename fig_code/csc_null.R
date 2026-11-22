# =====================================================================
# csc_null.R - CSC null model (reviewer comment: are the 61/130 positive contributions merely an artefact of the decomposition?)
#   Null A, exchangeability null: permute the country labels of Delta within each (indicator, year) segment (destroys the persistence of
#           country identity; the marginal distribution within a segment, S_t, g_t and sgn_S are all unchanged)
#   Null B, sign-symmetry null: flip the sign of every Delta independently (Rademacher) (destroys directional synergy -> sgn_S is randomised,
#           |Delta| and S_t are unchanged)
#   Null C, binomial reference: share_positive ~ Binomial(130, 1/2)
#   Check: the unpermuted Delta must reproduce data/intermediate/csc_country.csv country by country
#         (61 positive, Mexico 24.357..., top5 61.6%, max 43.2%)
#   Test statistics: share_positive (share of countries with a positive contribution, 47% in the text), top5_share (62% in the text),
#             max_share (43% in the text); two-sided p for the share, upper-sided p for the concentration measures
#   Written to: data/output/csc_null_draws.csv, csc_null_summary.csv,
#         data/fig_data/sup_fig13/
# Run from the repository root:  Rscript fig_code/csc_null.R   (about 1-3 minutes)
# =====================================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
})
set.seed(20261110)
B_REP <- 2000

# ---------- 0. Data and year-on-year changes (identical formula by formula to gbf_compute_pipeline.R synergy_decompose) ----------
norm <- read_csv("data/intermediate/norm_indicator.csv", show_col_types = FALSE)
cov  <- read_csv("data/intermediate/country_coverage.csv", show_col_types = FALSE)
meta <- read_csv("data/meta.csv", show_col_types = FALSE)
keep <- cov$country[cov$keep_index]

expected_map <- setNames(ifelse(meta$direction == "+", 1, -1), meta$indicator_id)

pairs <- norm %>%
  filter(year %in% 2011:2024, country %in% keep) %>%
  arrange(indicator, country, year) %>%
  group_by(indicator, country) %>%
  mutate(r_next = lead(norm), year_next = lead(year)) %>%
  filter(!is.na(r_next), year_next == year + 1) %>%
  ungroup() %>%
  transmute(indicator, year, country, delta = r_next - norm)

# Segment = (indicator, year); the order of seg_id is fixed
pairs <- pairs %>%
  mutate(seg = paste(indicator, year, sep = "_")) %>%
  arrange(seg, country)
seg_levels <- unique(pairs$seg)
seg_expected <- as.integer(expected_map[sub("_.*$", "", seg_levels)])

# Segment-level quantities (eq. 10/11/12) and the unit base value u = (|Delta| - d_t/C)*100/(S_t*C); CSC = u * sgn_S
seg_stats <- pairs %>%
  group_by(seg) %>%
  summarise(C = n(), d_foreign = sqrt(sum(delta^2)), sum_abs = sum(abs(delta)),
            g_t = mean(delta), .groups = "drop") %>%
  mutate(S_t = (sum_abs - d_foreign) / C,
         sgn_S = ifelse(sign(g_t) == seg_expected[match(seg, seg_levels)], 1, -1),
         w = 100 / (S_t * C),
         ref = d_foreign / C)
stopifnot(setequal(seg_stats$seg, seg_levels))
seg_stats <- seg_stats[match(seg_levels, seg_stats$seg), ]

u_vals <- (abs(pairs$delta) - seg_stats$ref[match(pairs$seg, seg_levels)]) *
  seg_stats$w[match(pairs$seg, seg_levels)]
sgn_cell <- seg_stats$sgn_S[match(pairs$seg, seg_levels)]
cty <- pairs$country
seg_idx <- match(pairs$seg, seg_levels)
n_cty <- n_distinct(cty)

# Country-level CSC mean = mean(u * sgn_S) (segments with S_t==0 give w=Inf/NaN -> excluded, as with na.rm=TRUE)
cty_means <- function(u, sgn_cell) {
  v <- u * sgn_cell
  ok <- is.finite(v)
  m <- tapply(v[ok], cty[ok], mean)
  m
}

# Test statistics
stat_set <- function(m) {
  m <- m[is.finite(m)]
  pos <- sort(m[m > 0], decreasing = TRUE)
  if (length(pos) == 0) return(c(share_positive = 0, top5_share = NA_real_,
                                 max_share = NA_real_))
  c(share_positive = mean(m > 0),
    top5_share = sum(head(pos, 5)) / sum(pos),
    max_share = unname(pos[1]) / sum(pos))
}

# ---------- 1. Check: reproduce csc_country.csv ----------
obs_m <- cty_means(u_vals, sgn_cell)
ref_tab <- read_csv("data/intermediate/csc_country.csv", show_col_types = FALSE)
joined <- data.frame(country = names(obs_m), mine = as.numeric(obs_m),
                     stringsAsFactors = FALSE) %>%
  inner_join(ref_tab %>% select(country, CSC_mean), by = "country")
max_diff <- max(abs(joined$mine - joined$CSC_mean))
n_pos <- sum(joined$mine > 0)
message(sprintf("[check] countries=%d  positive=%d  max|diff| vs csc_country = %.3e",
                nrow(joined), n_pos, max_diff))
stopifnot(max_diff < 1e-6, n_pos == 61)

obs <- stat_set(obs_m)
message(sprintf("[obs] share_positive=%.4f  top5_share=%.4f  max_share=%.4f",
                obs["share_positive"], obs["top5_share"], obs["max_share"]))

# ---------- 2. Null models ----------
cell_by_seg <- split(seq_along(u_vals), seg_idx)
n_seg <- length(seg_levels)

null_run <- function(mode, B) {
  out <- matrix(NA_real_, B, 3,
                dimnames = list(NULL, c("share_positive", "top5_share", "max_share")))
  u0 <- u_vals
  for (b in seq_len(B)) {
    if (mode == "A") {
      # Exchangeability: permute the country assignment of u within the segment
      u <- u0
      for (s in seq_len(n_seg)) {
        ii <- cell_by_seg[[s]]
        u[ii] <- sample(u0[ii], length(ii), replace = FALSE)
      }
      sgn <- sgn_cell
    } else {
      # Sign symmetry: flip the sign independently within each segment (Rademacher)
      u <- u0
      sgn <- sample(c(-1, 1), n_seg, replace = TRUE)[seg_idx]
    }
    out[b, ] <- stat_set(cty_means(u, sgn))
  }
  out
}

message("Null A (exchangeability) ...")
draws_A <- null_run("A", B_REP)
message("Null B (sign symmetry) ...")
draws_B <- null_run("B", B_REP)

# ---------- 3. p values and summary ----------
two_sided_p <- function(null, o) {
  p <- 2 * min(mean(null <= o), mean(null >= o))
  min(1, p)
}
upper_p <- function(null, o) (1 + sum(null >= o)) / (length(null) + 1)

stats_n <- c("share_positive", "top5_share", "max_share")
summarise_null <- function(draws, mode) {
  bind_rows(lapply(stats_n, function(s) {
    null <- draws[, s]; o <- obs[s]
    tibble(stat = s, model = mode, observed = unname(o),
           null_mean = mean(null), null_sd = sd(null),
           q025 = quantile(null, 0.025), q975 = quantile(null, 0.975),
           p_value = if (s == "share_positive") two_sided_p(null, o) else upper_p(null, o),
           p_type = if (s == "share_positive") "two-sided" else "upper-tail",
           B = length(null))
  }))
}
binom_row <- tibble(
  stat = "share_positive", model = "C_binomial",
  observed = unname(obs["share_positive"]),
  null_mean = 0.5, null_sd = sqrt(0.5 * 0.5 / n_cty),
  q025 = qbinom(0.025, n_cty, 0.5) / n_cty, q975 = qbinom(0.975, n_cty, 0.5) / n_cty,
  p_value = binom.test(round(obs["share_positive"] * n_cty), n_cty, 0.5)$p.value,
  p_type = "two-sided", B = NA_integer_)

summary_tab <- bind_rows(summarise_null(draws_A, "A_exchangeability"),
                         summarise_null(draws_B, "B_sign_symmetry"),
                         binom_row)

put_dir <- "data/output"
write_csv(summary_tab, file.path(put_dir, "csc_null_summary.csv"))
draws_long <- bind_rows(
  as_tibble(draws_A) %>% mutate(model = "A_exchangeability", rep = row_number()),
  as_tibble(draws_B) %>% mutate(model = "B_sign_symmetry", rep = row_number())) %>%
  pivot_longer(all_of(stats_n), names_to = "stat", values_to = "value")
write_csv(draws_long, file.path(put_dir, "csc_null_draws.csv"))
message("[out] data/output/csc_null_summary.csv / csc_null_draws.csv")

# ---------- 4. Sup Fig 13 figure source data ----------
d13 <- "data/fig_data/sup_fig13"
if (!dir.exists(d13)) dir.create(d13, recursive = TRUE)
write_csv(draws_long, file.path(d13, "fig13_null_draws.csv"))
write_csv(summary_tab %>% select(stat, model, observed, p_value, p_type),
          file.path(d13, "fig13_null_summary.csv"))
write_lines(c(
  "# Sup Fig 13 figure source data - CSC null model (Null A exchangeability / Null B sign symmetry)",
  "",
  "- fig13_null_draws.csv: share_positive / top5_share / max_share for every permutation",
  "- fig13_null_summary.csv: observed values and p values",
  "",
  "Generated by: Rscript fig_code/csc_null.R; the observed statistics must agree with csc_country.csv (the script contains assertions)"),
  file.path(d13, "README.txt"))

print(summary_tab)
message("CSC NULL DONE")
