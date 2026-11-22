# Code and data for the Kunming–Montreal Global Biodiversity Framework study

This repository contains the complete analysis code and the data needed to reproduce the results of the manuscript on the progress of the Kunming–Montreal Global Biodiversity Framework (GBF): the computation pipeline, the figure and table source data, the raw indicator data and the rendered figures.

**All paths in the scripts are relative to the repository root. Run everything from the repository root.**

---

## Quick start

```r
setwd("<path-to-this-repository>")

# 1. Verify the implementation against the worked examples of the Supplementary Information
source("calc_examples.R")            # prints "CALC EXAMPLES PASSED"

# 2. Run the main pipeline (n_boot = 10000 for the final results)
source("gbf_compute_pipeline.R")
source("gbf_compute_remaining.R")
res  <- run_pipeline()
save_intermediates(res)
rest <- run_remaining(res)
```

```sh
# 3. Rebuild the figure and table source data and render the figures
Rscript fig_code/make_fig_data.R
Rscript fig_code/fig1.R ; Rscript fig_code/fig2.R ; Rscript fig_code/fig3.R ; Rscript fig_code/fig4.R
Rscript fig_code/ed_fig1.R ; Rscript fig_code/sup_fig1.R ; Rscript fig_code/sup_fig2.R
Rscript fig_code/sup_fig3.R ; Rscript fig_code/sup_fig4.R ; Rscript fig_code/sup_figs_5to11.R
Rscript fig_code/sup_fig12.R
Rscript fig_code/csc_null.R ; Rscript fig_code/sup_fig13.R
Rscript fig_code/endfill_sensitivity.R
Rscript fig_code/inherited_trend_test.R ; Rscript fig_code/lag_vs_ineffectiveness.R
Rscript table_code/make_tables.R
```

The exact run order, the dependencies and the conventions of every input file are described in **[RUN_GUIDE.md](RUN_GUIDE.md)**.

---

## What is in this repository

| Path | Content |
|---|---|
| `gbf_compute_pipeline.R` | The main computation pipeline: the normalization (Formula 1), the target/goal values (Formula 2), MIS (Formula 3), EIS (Formula 4), the synergy family and the country contributions (Formulas 5–9), the PGR family (Formulas 10–12), the achievement rate (Formula 13), the baseline parameters (Supplementary Table 3), the four-class classification (A.5) and the projection with the bootstrap (A.5) |
| `gbf_compute_remaining.R` | The derived analyses behind Supplementary Tables 3–10 (internal products `table_S1`–`table_S6`), the multiple-comparison corrections (the note of Supplementary Table 3), the section/goal summaries (Fig. 2c), the trajectory analyses (B.2.7), the achievement-rate intervals (A.5), the five uncertainty algorithms (B.1/B.2.1), the weight sensitivity (B.2.2), the endpoint robustness (B.2.3/B.2.4) and the figure source data of Figs 1/3/4 |
| `calc_examples.R` | Asserts every numerical worked example of the Supplementary Information |
| `scripts/` | Targeted re-runs of the long-running modules (Algorithms 2, 3 and 5, and the Algorithm 4 / weight-sensitivity / endpoint-robustness refresh) from the checkpoints, without repeating the whole pipeline |
| `fig_code/` | Builds the figure source data and renders the main, Extended Data and Supplementary figures (see [fig_code/README.md](fig_code/README.md)) |
| `table_code/` | Regenerates the source data of the Supplementary Tables (see [table_code/README.md](table_code/README.md)) |
| `data/` | The raw indicator data, the metadata, the intermediate products and the source data of every figure and table (see [data/README.md](data/README.md)) |

---

## The method in one paragraph

Fifty-three candidate indicators (44 used) are harmonized per country and year, imputed inside the 2011–2024 window with the anchors taken from the complete history, and normalized to 0–100 against the best/worst bounds. The indicators are aggregated into the 23 targets and the 4 goals of the GBF (equal-weight means), and combined into the Mean Index Score (MIS, the arithmetic mean) and the Evenness Index Score (EIS, the area of the radar plot relative to the circle with the same perimeter). The trajectories are projected to the target years (2030 for the indicators and targets, 2050 for the goals) with a multi-model average (low-noise series) or a dynamic linear model (high-noise series, selected by the signal-to-noise ratio). The synergy effects and the country contributions follow Formulas 5–9; the projected growth is compared with the growth required by the GBF thresholds (or by the historical trend when no threshold exists) in Formulas 10–12; and the resulting achievement rate and the classification into four trajectory classes follow Formula 13 and section A.5 of the Supplementary Information.

---

## Requirements

- R (developed with R 4.5.3)
- Required packages: `tidyverse` (`dplyr`, `readr`, `tidyr`, `purrr`, `stringr`, `tibble`), `readxl`, `ggplot2`, `scales`, `patchwork`, `sf`, `rnaturalearth`, `ggrepel`
- Optional packages: `dlm` (the DLM projection of the high-noise series; the code falls back to an equivalent built-in Kalman filter), `forecast` (ETS / auto-ARIMA in Algorithm 3; falls back to a linear model), `ragg` (the TIFF export of the figures)

---

## Reproducibility notes

- The main pipeline writes a checkpoint (`data/intermediate/chk_*.rds`) after every stage, so an interrupted run can be resumed; the checkpoints carry a version number and are invalidated automatically when the computation logic changes.
- The final results require the full 10,000-fold bootstrap; `run_pipeline(n_boot = 500)` is only meant for a trial run.
- The intermediate and output products are committed, so the figures and the tables can be rebuilt without repeating the bootstrap. Delete `data/intermediate/` and `data/output/` to recompute everything from the raw data.
- The figures are written to `data/figures/` as PDF (vector), TIFF and PNG (300 dpi, when the `ragg` package is available).

## Citation

Please cite the archived code and data:

> Yi, J., Gregory, R. D., Cao, M., Cao, Y., Dubois, G., Gonzalez, A., Lin, H., Liu, W., Delli, G., Wu, J., Zhang, R., Yu, D. & Xu, H. Code and data for: Urgent and balanced action needed to achieve global biodiversity framework. Zenodo https://doi.org/10.5281/zenodo.24430436 (2026).

The bibliographic details of the associated manuscript (authors, title, journal, DOI) are to be added here at the time of publication.

## Data sources

The raw indicator data in `data/origin_data/` are harmonized extracts of the official sources of the Global Biodiversity Framework indicators and of the related SDG series. The country classification is the UN M49 country classification (`data/geoclass.xlsx`). Please refer to the Supplementary Information of the manuscript for the full list of sources and for the procedures applied to each indicator.
