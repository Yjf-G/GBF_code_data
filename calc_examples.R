# =====================================================================
# GBF paper - numerical re-computation and cross-check of the Supplementary Information worked examples
# Algorithmic reference: the submitted Supplementary Information (A.1-A.5 methods, B.1-B.2 uncertainty
#   and sensitivity analysis)
# Re-computes every numerical worked example of the Supplementary Information one by one and asserts that
# the numbers match the document (A.2 imputation; A.4.2 Formula 4; A.4.3 Formulas 5-9; A.5 Formulas 10-12;
# A.5 Formulas 10-13); the tables are printed at the same time so that they can be quoted directly.
# Formula numbering of the Supplementary Information: 1 normalized value / 2 goal-target value / 3 MIS /
#   4 EIS / 5 DC / 6 FC / 7 SyE / 8 SyD / 9 CSC / 10 PGRsim / 11 PGRreq / 12 PGRgap / 13 AR.
# The blocks marked "implementation check" are internal consistency checks of the code and are not
#   worked examples of the document.
#
# Run: source("calc_examples.R")
# =====================================================================

suppressPackageStartupMessages(library(dplyr))

ok <- function(label, got, expect, tol = 0.05) {
  pass <- abs(got - expect) < tol
  cat(sprintf("%-58s got=%8.2f expect=%8.2f  %s\n",
              label, got, expect, if (pass) "OK" else "FAIL"))
  if (!pass) stop("value mismatch: ", label)
  invisible(pass)
}

cat("==== A.2 Example (imputation) ====\n")
# Series 2015-2020: (41, NA, NA, 47, 50, NA)
# Missing endpoint 2020 -> forward fill with 50; internal gap 2016-2017 -> linear interpolation between 41 and 47, giving 43.0/45.0
imp <- function(year, value) {           # identical to gbf_compute_pipeline.R::impute_series
  ord <- order(year); y <- year[ord]; v <- value[ord]
  obs <- which(!is.na(v))
  if (obs[1] > 1) v[1:(obs[1] - 1)] <- v[obs[1]]
  if (obs[length(obs)] < length(v)) v[(obs[length(obs)] + 1):length(v)] <- v[obs[length(obs)]]
  nas <- which(is.na(v))
  if (length(nas) > 0) v[nas] <- approx(y[!is.na(v)], v[!is.na(v)], xout = y[nas])$y
  v
}
v <- imp(2015:2020, c(41, NA, NA, 47, 50, NA))
ok("2016 linear interpolation", v[2], 43.0, 1e-9)
ok("2017 linear interpolation", v[3], 45.0, 1e-9)
ok("2020 endpoint forward fill", v[6], 50, 1e-9)
cat("Categorical series (implementation level): 2016-2017 gap -> LOCF carries the 2015 value forward (qualitative rule)\n\n")

cat("==== Formula 1 normalization (implementation check) ====\n")
norm_high <- function(x, worst, best) pmin(100, pmax(0, 100 * (x - worst) / (best - worst)))
norm_low  <- function(x, worst, best) pmin(100, pmax(0, 100 * (worst - x) / (worst - best)))
ok("OHI 64.3 (28.6→100)", norm_high(64.3, 28.6, 100), 50.0)
ok("OHI 25 below worst -> 0", norm_high(25, 28.6, 100), 0)
ok("OHI 110 above best -> 100", norm_high(110, 28.6, 100), 100)
ok("EcoFootprint 15 (18.89→0.61)", norm_low(15, 18.89, 0.61), 21.3)
ok("EcoFootprint 20 → 0", norm_low(20, 18.89, 0.61), 0)
ok("EcoFootprint 0.3 → 100", norm_low(0.3, 18.89, 0.61), 100)
cat("\n")

cat("==== Formula 2/3 aggregation (implementation check) ====\n")
ok("F2 (80+60)/2", mean(c(80, 60)), 70, 1e-9)
ok("F3 MIS (70,30,50,90)/4", mean(c(70, 30, 50, 90)), 60, 1e-9)
cat("Second country has 25% of its indicators missing (>20%) -> excluded from MIS (coverage rule, qualitative)\n\n")

cat("==== A.4.2 Example (Formula 4) ====\n")
eis <- function(r) {
  n <- length(r)
  S <- sum((1 / n) * pi * r^2)
  L <- 2 * (max(r) - min(r)) + (2 * pi / n) * sum(r)
  400 * pi * S / L^2
}
ok("EIS Y=(50,50) = 100.0", eis(c(50, 50)), 100.0)
ok("EIS X=(70,30) = 73.7", eis(c(70, 30)), 73.7)
ok("EIS Z=(30,70) = 73.7 (mirror of X)", eis(c(30, 70)), 73.7)
ok("EIS (50,50,50,50) = 100.0 (implementation check)", eis(c(50, 50, 50, 50)), 100.0)
ok("EIS (80,60,40,20) = 62.8 (implementation check)", eis(c(80, 60, 40, 20)), 62.8)
cat("  Note: X, Z and Y share the same mean (MIS=50); EIS separates evenness - MIS and EIS are complementary by design\n\n")

cat("==== A.5 bootstrap validation (implementation check) ====\n")
cat("2.5%/97.5% quantiles of 10,000 annual rate-of-change estimates = 0.61 / 1.80\n")
cat("The interval does not contain 0 -> the change is significantly positive (qualitative rule; the numbers are the values given in the document)\n\n")

cat("==== A.5 Example (Formulas 10-12) ====\n")
ok("PGR_sim (40→55, 2020–2030)", 100 * (55 - 40) / 40, 37.5, 1e-9)
ok("PGR_req threshold target (Y=100)", 100 * (100 - 40) / 40, 150, 1e-9)
ok("PGR gap = 150 − 37.5", 150 - 37.5, 112.5, 1e-9)
# No-threshold branch: mean historical change +1.0, n=10 -> Y = 40+10 = 50
ok("PGR_req without threshold (Y=50)", 100 * (50 - 40) / 40, 25, 1e-9)
ok("PGR gap = 25 - 37.5 = -12.5 (surplus)", 25 - 37.5, -12.5, 1e-9)
cat("\n")

cat("==== A.5 Example (Formulas 10-13), Target 3 with three countries ====\n")
cat("Country | v2020 | Y2030 | p    | PGRsim | PGRreq | PGR gap | I\n")
dat <- tibble::tibble(
  country = c("A", "B", "C"),
  v2020 = c(50, 40, 30),
  y2030 = c(110, 85, 60),
  p = c(0.02, 0.03, 0.07)
)
dat <- dat %>%
  dplyr::mutate(
    pgr_sim = 100 * (y2030 - v2020) / v2020,
    pgr_req = 100 * (100 - v2020) / v2020,
    pgr_gap = pgr_req - pgr_sim,
    significant = p < 0.05,
    I = as.integer(pgr_sim >= pgr_req & significant)
  )
for (i in 1:3) {
  cat(sprintf("%-7s | %5.0f | %5.0f | %.2f | %6.1f | %6.1f | %7.1f | %d\n",
              dat$country[i], dat$v2020[i], dat$y2030[i], dat$p[i],
              dat$pgr_sim[i], dat$pgr_req[i], dat$pgr_gap[i], dat$I[i]))
}
ok("A: PGR_sim = 120%", dat$pgr_sim[1], 120, 1e-9)
ok("A: PGR_req = 100%", dat$pgr_req[1], 100, 1e-9)
ok("A: PGR gap = -20 (surplus)", dat$pgr_gap[1], -20, 1e-9)
ok("B: PGR_sim = 112.5%", dat$pgr_sim[2], 112.5, 1e-9)
ok("B: PGR_req = 150%", dat$pgr_req[2], 150, 1e-9)
ok("B: PGR gap = 37.5 (shortfall)", dat$pgr_gap[2], 37.5, 1e-9)
ok("C: PGR_sim = 100%", dat$pgr_sim[3], 100, 1e-9)
ok("C: PGR_req ≈ 233.3%", dat$pgr_req[3], 233.3, 0.05)
ok("C: PGR gap = 133.3 (shortfall)", dat$pgr_gap[3], 133.3, 0.05)
stopifnot(identical(dat$I, c(1L, 0L, 0L)))
ok("AR(Target 3) = (1+0+0)/3", mean(dat$I) * 100, 33.3, 0.1)
cat("  I decision: A=1 (sim>=req and significant); B=0 (growth too slow); C=0 (p=0.07, not significant)\n\n")

cat("==== A.4.3 Example (Formulas 5-9) ====\n")
# SI Formulas 5-9: 5 d=|Dr| (DC); 6 d_t=sqrt(Sum Dr^2) (FC); 7 S_t=(1/C)(Sum d - d_t) (SyE);
# 8 sgn(S_t)=+1 if sgn(g_t) matches the expected direction (SyD); 9 CSC=(d - d_t/C)/(S_t*sgnS*C)*100%
syn_case <- function(d, expected_dir = 1) {
  DC <- abs(d)
  FC <- sqrt(sum(d^2))
  C <- length(d)
  S_t <- (sum(DC) - FC) / C
  g <- mean(d)
  sgnS <- ifelse(sign(g) == expected_dir, 1, -1)
  CSC <- (DC - FC / C) / (S_t * sgnS * C) * 100
  list(DC = DC, FC = FC, S_t = S_t, g = g, sgnS = sgnS, CSC = CSC)
}
p1 <- syn_case(c(10, -5, 10))
cat("Period 1: A 40->50, B 50->45, C 60->70; signed changes (+10, -5, +10)\n")
cat(sprintf("  Formula 5  DC = (%.0f, %.0f, %.0f)\n", p1$DC[1], p1$DC[2], p1$DC[3]))
ok("  Formula 6 d_t = √(10²+5²+10²) = 15", p1$FC, 15, 1e-9)
ok("  Formula 7 S_t = (25−15)/3 = 3.33", p1$S_t, 3.33, 0.01)
ok("  Formula 8 g_t = +5", p1$g, 5, 1e-9)
stopifnot(p1$sgnS == 1)
cat("  Formula 8 sgn(S_t) = +1 (consistent with the expected direction)\n")
ok("  Formula 9 CSC_A = (10−5)/(3.33×3) = 50%", p1$CSC[1], 50, 1e-9)
ok("  Formula 9 CSC_B = 0% (exactly at the per-capita reference d_t/C=5)", p1$CSC[2], 0, 1e-9)
ok("  Formula 9 CSC_C = 50%", p1$CSC[3], 50, 1e-9)
ok("  Sum of CSC = +100%", sum(p1$CSC), 100, 1e-6)
p2 <- syn_case(c(-10, -5, -10))
cat("Period 2: all countries decline by the same amount; signed changes (-10, -5, -10)\n")
ok("  d_t identical to Period 1 = 15", p2$FC, 15, 1e-9)
ok("  S_t also = 3.33", p2$S_t, 3.33, 0.01)
ok("  g_t = −25/3 ≈ −8.33", p2$g, -8.33, 0.01)
stopifnot(p2$sgnS == -1)
cat("  Formula 8 sgn(S_t) = -1 (opposite to the expected direction, deterioration is recorded as -1)\n")
ok("  CSC_A = −50%", p2$CSC[1], -50, 1e-9)
ok("  CSC_B = 0%", p2$CSC[2], 0, 1e-9)
ok("  CSC_C = −50%", p2$CSC[3], -50, 1e-9)
ok("  Sum of CSC = -100%", sum(p2$CSC), -100, 1e-6)
cat("  Note: a single country contribution is not limited to +/-100%; the sum over countries = +100% or -100% (sign of the gross synergy)\n\n")

cat("==== A.5 four-class rho classification (reference table) ====\n")
rho_class <- function(rho, sig_pos) {
  if (is.na(rho)) return(NA_character_)
  if (rho >= 1 && sig_pos) return("Full achievement")
  if (rho >= 0.5 && sig_pos) return("Favourable")
  if ((sig_pos && rho > 0 && rho < 0.5) || (rho >= 1 && !sig_pos)) return("Mixed")
  if (rho <= 0 || !sig_pos) return("Unfavourable")
  "Mixed"
}
chk <- tibble::tibble(
  rho = c(1.2, 0.7, 0.3, 1.2, 0.3, -0.1),
  sig = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
  expect = c("Full achievement", "Favourable", "Mixed", "Mixed",
             "Unfavourable", "Unfavourable")
)
chk$got <- mapply(rho_class, chk$rho, chk$sig)
print(chk)
stopifnot(all(chk$got == chk$expect))
cat("  rho classification OK\n\n")

cat("All worked example values agree with the Supplementary Information - CALC EXAMPLES PASSED\n")
