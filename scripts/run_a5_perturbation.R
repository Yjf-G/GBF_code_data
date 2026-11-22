# Targeted re-run: Algorithm 5 single-indicator perturbation sensitivity (computes a5 only, does not re-run the whole pipeline)
#   Input convention identical to a1/a3: bind_rows(norm_indicator, norm_regional) -- regional indicators 27-29 are merged in
#   One row per indicator (42 target-mapped indicators); no duplicate rows any more
# Usage: Rscript scripts/run_a5_perturbation.R   (about 1-2 minutes)
GBF_NO_AUTORUN <- TRUE
source("gbf_compute_pipeline.R")
source("gbf_compute_remaining.R")
res <- load_res_from_checkpoints()
i2t <- load_i2t()
t0 <- Sys.time()
a5 <- algo5_perturbation(bind_rows(res$norm_indicator, res$norm_regional), i2t)
cat(sprintf("[a5 full] %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
write_csv(a5, "data/output/uncertainty__a5_perturbation.csv")
if (!dir.exists("data/fig_data/sup_fig4")) dir.create("data/fig_data/sup_fig4", recursive = TRUE)
write_csv(a5 %>% arrange(desc(Sx_mean_abs)),
          "data/fig_data/sup_fig4/supfig4_perturbation_sensitivity.csv")
cat("saved 2 CSVs; rows:", nrow(a5), "\n")
print(as.data.frame(a5 %>% arrange(desc(Sx_mean_abs))), digits = 3)
cat("A5_DONE\n")
