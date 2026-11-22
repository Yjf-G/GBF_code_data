# Targeted re-run: Algorithm 3 time-series length uncertainty (computes a3 only, does not re-run the whole pipeline)
# Inputs: chk_06_norm.rds (country-level normalization) + chk_06b_norm_regional.rds (regional indicators 27-29)
# Outputs: data/output/uncertainty__a3_window_length__{curves,optimal,cross_model}.csv
# Run from the repository root:  Rscript scripts/run_a3_bootstrap.R   (about 5-10 minutes)
GBF_NO_AUTORUN <- TRUE
source("gbf_compute_pipeline.R")
source("gbf_compute_remaining.R")
norm_indicator <- readRDS("data/intermediate/chk_06_norm.rds")
norm_regional  <- readRDS("data/intermediate/chk_06b_norm_regional.rds")
t0 <- Sys.time()
r <- algo3_window_length(bind_rows(norm_indicator, norm_regional),
                         B = 1000, seed = 42)
cat(sprintf("[a3 full] %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
write_csv(r$curves,      "data/output/uncertainty__a3_window_length__curves.csv")
write_csv(r$optimal,     "data/output/uncertainty__a3_window_length__optimal.csv")
write_csv(r$cross_model, "data/output/uncertainty__a3_window_length__cross_model.csv")
cat("saved 3 CSVs\n")
o <- r$optimal
cat(sprintf("L_opt (Kneedle, falls back to L_der when no knee is found): mean %.2f sd %.2f median %.1f range %d-%d\n",
            mean(o$L_opt), sd(o$L_opt), median(o$L_opt), min(o$L_opt), max(o$L_opt)))
cat(sprintf("convex knee found: %d/%d; 5-15 y share: %d/%d (%.0f%%)\n",
            sum(o$knee_found), nrow(o),
            sum(o$L_opt >= 5 & o$L_opt <= 15), nrow(o),
            100 * mean(o$L_opt >= 5 & o$L_opt <= 15)))
cat(sprintf("L_der: mean %.2f median %.1f; model spread: mean %.2f, <=1 y in %d/%d\n",
            mean(o$L_der), median(o$L_der), mean(o$model_spread),
            sum(o$model_spread <= 1), nrow(o)))
print(as.data.frame(o), digits = 3)
cat("A3_FULL_DONE\n")
