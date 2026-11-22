# GBF paper — data layout and run guide

**Algorithm reference:** the Supplementary Information of the submission (sections A.1–A.5 methods, B.1–B.2 uncertainty and sensitivity analysis; formula numbering 1–13 as given in A.3–A.5)
**Companion R code:** `gbf_compute_pipeline.R` (main computation) + `gbf_compute_remaining.R` (derived analyses) + `calc_examples.R` (verification of the worked examples)

All paths below are **relative to the repository root**. Run every script from the repository root, e.g.

```r
setwd("<path-to-this-repository>")
source("gbf_compute_pipeline.R")
source("gbf_compute_remaining.R")
```

---

## 1. How the code maps onto the manuscript

| Code module | What it implements |
|---|---|
| `gbf_compute_pipeline.R` | Formula 1 normalization; Formula 2 target/goal values; Formula 3 MIS; Formula 4 EIS; Formulas 5–9 DC/FC/SyE/SyD/CSC; Formulas 10–12 PGR_sim / PGR_req / PGR gap; Formula 13 AR (I = PGR_sim ≥ PGR_req and a rate of change significantly ≠ 0, bootstrap p < 0.05); the four baseline parameters of Supplementary Table 3; the four-class classification of A.5; the projection of A.5 (DLM/PMA) + 10,000 bootstrap resamples |
| `gbf_compute_remaining.R` | Supplementary Tables 3–10 (internal products `table_S1`–`table_S6`); the Bonferroni + BH-FDR corrections of the Supplementary Table 3 note; the section/goal summaries of Fig. 2c; the trajectory analyses of B.2.7 (monotonicity / acceleration / endpoint bias / COVID / regional paths); the achievement-rate intervals of A.5; Algorithms 1–5 of B.1/B.2.1; the weight sensitivity of B.2.2; the endpoint robustness of B.2.3/B.2.4; the figure source data of Figs 1/3/4 |
| `calc_examples.R` | Recomputes and asserts every numerical worked example of the Supplementary Information (A.2 imputation 43/45/50; A.4.2 EIS 100.0/73.7/73.7; A.4.3 the synergy S_t = 3.33 and the CSC ±100%; A.5 the full PGR set and AR = 33.3% for the three countries of Target 3), plus implementation checks of Formulas 1–3 |
| `fig_code/` | Builds the figure source data (`make_fig_data.R`) and renders every main, Extended Data and Supplementary figure |
| `table_code/make_tables.R` | Regenerates the source data of the Supplementary Tables from the pipeline outputs |
| `scripts/` | Targeted re-runs of the long-running modules (Algorithms 2, 3, 5 and the Algorithm 4 / weight-sensitivity / endpoint-robustness refresh) without re-running the whole pipeline |

---

## 2. Directory layout

```
repository root/
├── README.md                       # this guide's short companion: overview and quick start
├── RUN_GUIDE.md                    # this document
├── gbf_compute_pipeline.R          # main pipeline
├── gbf_compute_remaining.R         # derived analyses
├── calc_examples.R                 # worked-example verification (runs on its own)
├── scripts/                        # targeted re-runs of the long modules
│   ├── run_a2_bootstrap.R
│   ├── run_a3_bootstrap.R
│   ├── run_a5_perturbation.R
│   └── run_c3c5_refresh.R
├── fig_code/                       # figure source data + rendering
├── table_code/                     # table source data
└── data/
    ├── meta.csv                    # 53 candidate indicators
    ├── indicator_to_target.csv     # indicator -> target/goal mapping
    ├── target_section.csv          # target <-> section <-> goal
    ├── indicator_dpsir.csv         # indicator typing (Process/Response/Pressure/State)
    ├── geoclass.xlsx               # country classification by region and group
    ├── origin_data/                # the 53 raw indicator CSVs
    ├── intermediate/               # written by save_intermediates()
    ├── output/                     # written by run_remaining()
    ├── fig_data/                   # figure source data, one folder per figure
    ├── tables/                     # source data of the Supplementary Tables
    └── figures/                    # the rendered figures (PDF, TIFF and PNG)
```

---

## 3. Input files

### 3.1 `data/geoclass.xlsx` (country classification, 288 rows × 19 columns)

The UN M49 country classification table and the **single authoritative source** of the country → region/group information. `build_country_meta()` reads it automatically (requires the R package `readxl`); if the file is missing it falls back to aggregating the geographic columns carried by the indicator files in `origin_data/`.

| Column | Use |
|---|---|
| `Geographical region` | Country/territory name (aligned with `GeoAreaName` in the indicator data) |
| `Continent` / `Subcontinent` / `Region` | Continent / subregion / region (Fig 3, Supplementary Table 7, KW test grouping) |
| `country groups1` / `country groups2` | Country groups (OECD and others) |
| `UN M.49 CODE` / `ISO CODES_2` / `ISO CODES_3` | Codes |
| `European Union` / `Eurozone` | EU membership flag (value `EU`) |
| `United Nations` | Value `nation` = sovereign state (193 countries) |
| `Developed and developing regions-developing/-developed` | Development status (`dl` / `dd`) |
| `Least developed countries` | LDC flag (`ldc`, 49 countries) |
| `Landlocked developing countries` | LLDC flag (`llc`, 31 countries) |
| `Small island developing States` | SIDS flag (`sid`, 51 countries) |
| `Transition countries_*` | Transition economies |

`build_country_meta()` outputs the columns `country, iso3, continent, subcontinent, region, group1, group2, is_nation, development, sids, ldc, llc, oecd, eu`.

**Country universe:** when counting countries, only the 193 states with `United Nations == "nation"` are counted. Non-country rows such as continents/regions/`Global average` enter no country-level statistic (coverage denominators, the C of Formulas 13/7, `n_countries`, CSC rankings, etc.).

**Order of processing:** all 44 indicators with `used=Y` are **imputed first** (filling only the missing values of the 2011–2024 window), and the filtering is applied **after** the imputation — (1) the country panel is filtered to the 193-country universe; (2) the **regional indicators 27–29** (RLI food/medicinal species, RLI traded species, RLI impact of invasive alien species, whose raw data are continent/region series) are routed out and loaded **only** into the global/regional analyses (the regional rows of Supplementary Table 10, `regional_indicator_summary.csv`); they take no part in any country-level analysis.

**Grouping convention** (`group_membership()`, used for Supplementary Table 7 / the regional paths of B.2.7 / Fig 3): continent, subcontinent, region, development (developed/developing), SIDS, LDC, LLDC, OECD, EU. The China row (`China` / CHN) must exist — the China row of the CSC table (Supplementary Table 8 source data) depends on it.

### 3.2 `data/meta.csv` (53 rows)

`indicator_id, indicator_name, used, period, countries, worst, best, direction`
- The 44 indicators with `used = Y` enter the main analysis; the 9 candidates with `used = N` are excluded (the reasons are given in Supplementary Table 1 of the Supplementary Information).
- `worst`/`best` are the 0-point and 100-point normalization bounds (annotated forms such as `0.6*` or `100 +` are allowed and cleaned by `parse_bound()`); when they are missing the A.3 decision tree provides the fallback (lower bound = the 2.5 quantile, upper bound = the mean of the top five countries).
- `direction`: `+` higher is better / `−` lower is better.

### 3.3 `data/indicator_to_target.csv`

`indicator_id, indicator_code, target, direction, section, goal`
- The `target` column mixes two kinds of values: `A`–`D` are goal mappings, `1`–`23` are target mappings (the code converts them to `T1`–`T23`).
- An indicator can be mapped to several targets/goals (many-to-many).

### 3.4 `data/target_section.csv`

`id, target, section, goal, theme` — the target ↔ section (reducing_threats / meeting_needs / tools_solutions) ↔ goal ↔ theme mapping.

### 3.5 `data/indicator_dpsir.csv`

`indicator_id, indicator_name, dpsir_class, dpsir_class_alt` — the indicator typing (Process / Response / Pressure / State) with an alternative coding. The typing defines the action-type (Process/Response) versus state-type (Pressure/State) diagnostics of SI B.2.4–B.2.6.

### 3.6 `data/origin_data/indicator*.csv` (53 files)

**All historical years are kept** — the imputation only fills the missing values inside the 2011–2024 window, but the projection fit needs the complete history to identify the pattern. Only the 44 indicators with `used=Y` in `meta.csv` enter the analysis.

Both layouts are accepted (`read_indicator_raw()` recognizes them automatically):
- **wide**: the first column is `GeoAreaName`, followed by four-digit year columns (e.g. `2000`…`2024`)
- **long**: it has `year` and `value` columns

The geographic group columns (`Continent` etc.) are optional — the country grouping follows `geoclass.xlsx`. The regional indicators (27–29, whose `GeoAreaName` is a continent/region name) enter as they are and are not subject to the 193-country universe.

---

## 4. Output files

### 4.1 `data/intermediate/` (`save_intermediates()`)

| File | Content |
|---|---|
| `country_meta.csv` | Country grouping (derived from geoclass.xlsx) |
| `norm_indicator.csv` | Normalized values at the indicator level (Formula 1, country panel) |
| `imputed_all.csv` / `imputed_regional.csv` / `norm_regional.csv` | The full imputation record and the imputed/normalized results for the regional indicators 27–29 |
| `imputation_detail.csv` | ⭐ Cell-by-cell imputation detail: `status` (observed/imputed/missing_no_anchor/discarded_outside_window), `method` (linear_interpolation/forward_fill/backward_fill), the anchor years and values, `value_raw`/`value_imputed`/`value_final`, `grid_added` (flag for the rows created for a missing year) |
| `norm_layer.csv` | Target/goal values (Formula 2, country level, historical segment) |
| `goal_target_ts_country.csv` / `goal_target_ts_region.csv` / `goal_target_ts_global.csv` | ⭐ Complete target/goal **time series** (country / region / global): from the existing data (observed + imputed) to the end of the projection (T1–T23 to 2030, Goals A–D to 2050). The columns include `value_type` (historical/projected), the year-by-year 95% CI of the projected segment (`ci_lo`/`ci_hi`), `model_type`/`model_form`, `n_ind`/`n_countries` |
| `mis_eis_annual.csv` | Annual MIS/EIS series 2011–2024 (Formulas 3/4, at both the goal and the target resolution, country level) |
| `mis_eis_annual_region.csv` / `mis_eis_annual_global.csv` | The annual MIS/EIS series (regional/global level, the equal-weight mean over the countries of the group, both resolutions) |
| `mis_eis_ts_country.csv` / `mis_eis_ts_region.csv` / `mis_eis_ts_global.csv` | ⭐ Complete MIS/EIS **time series** (country / region / global): the 2011–2024 history + the projected segment up to the target year (MIS/EIS_target to 2030, MIS/EIS_goal to 2050, computed from the projected target/goal values with Formulas 3/4); `value_type` marks historical/projected |
| `projections.csv` | Projection to the target year + CI + model form + the bootstrapped annual rate of change |
| `pgr.csv` | PGR_sim / PGR_req / PGR gap (Formulas 10–12) |
| `ar.csv` / `ar_summary.csv` / `ar_ci.csv` | I, AR and the 50%/95% intervals (Formula 13 / A.5) |
| `trajectory_class.csv` | The ρ classification + the confidence flag (A.5) |
| `indicator_projection.csv` | Projection and AR at the indicator level |
| `csc_country.csv` | The country synergy contribution CSC (Formula 9, mean across indicators and years) |
| `baseline_target_goal.csv` / `baseline_per_country.csv` | The four baseline parameters (Supplementary Table 3) |
| `country_coverage.csv` | Indicator coverage per country (more than 20% missing excludes it from MIS/EIS/rankings) |

### 4.2 `data/output/` (`run_remaining()`)

| File | Filled into |
|---|---|
| `imputation_detail.csv` | Cell-by-cell imputation detail (identical to the file of the same name in `intermediate`, kept in `output` as a backup) |
| `table_S6_imputation.csv` | Supplementary Table 9 (per-indicator imputation record, aggregated from the imputation detail: observed cells / imputed cells (linear/forward/backward) / cells that cannot be imputed / rows created for missing years / fills computed outside the window that are not retained) |
| `goal_target_ts_country.csv` / `goal_target_ts_region.csv` / `goal_target_ts_global.csv` | Complete target/goal time series (country / region / global, existing data → end of the projection 2030/2050) |
| `mis_eis_annual.csv` / `mis_eis_annual_region.csv` / `mis_eis_annual_global.csv` | Annual MIS/EIS series (country / region / global, 2011–2024, goal and target resolutions) |
| `mis_eis_ts_country.csv` / `mis_eis_ts_region.csv` / `mis_eis_ts_global.csv` | Complete MIS/EIS time series (country / region / global, history + projected segment up to 2030/2050) |
| `table_S1_projection_models.csv` | Supplementary Table 10 (series type/CV/SNR / the model with the highest weight / AICc / the mean number of models / top weight / the breakpoint year; including the rows of the regional series 27–29) |
| `regional_indicator_summary.csv` | The global/regional annual summary of the indicators 27–29 (used only by the global/regional analyses) |
| `table_S2_baseline_regression.csv` / `_kw.csv` / `_posthoc.csv` | Supplementary Table 3 (slope/R²/p + KW + MWU-Bonferroni + BH-FDR) |
| `table_S3_trajectories.csv` | Supplementary Table 4 (mean(SD) of PGR / ρ / bootstrap p / class / flag) |
| `table_S4_mis_eis_groups.csv` | Supplementary Table 7 (annual MIS/EIS by group + KS + share of monotonicity) |
| `table_S5_csc.csv` | Supplementary Table 8 (country rows incl. China) |
| `table_S6_imputation.csv` | Supplementary Table 9 |
| `b6_section_summary.csv` / `b6_goal_summary.csv` / `b6_kw_section.csv` | section/goal summaries (Fig. 2c) |
| `b7_monotonicity_*.csv` / `b7_acceleration*.csv` / `b7_endpoint_vs_trajectory.csv` / `b7_covid_perturbation.csv` / `b7_regional_paths.csv` | B.2.7 trajectory diagnostics |
| `c1_achievement_ci__*.csv` | The achievement-rate intervals of A.5 (the CI of the abstract/main text) |
| `c5_extrapolation_test.csv` / `c5_leave_one_year_out.csv` / `c5_loo_max.csv` | B.2.3 endpoint robustness |
| `fig1_slopes.csv` / `fig1c_quadrant.csv` / `fig3_*.csv` / `fig4_synergy_binomial.csv` | The figure source data of Figs 1/3/4 |
| `uncertainty__a1..a5__*.csv` | Algorithms 1–5 of B.1/B.2.1 |
| `uncertainty__c3_weight_sensitivity.csv` | B.2.2 (Dirichlet/PCA/equal weights/dropping BCH-46) |
| `csc_null_draws.csv` / `csc_null_summary.csv` | B.2.5 CSC null model (feeds Supplementary Fig. 13) |
| `endfill_*.csv` | B.2.4 end-of-series truncation sensitivity (feeds the new columns of Supplementary Table 9) |
| `inherited_trend_*.csv` | B.2.6 continuation of the pre-GBF trend test |
| `lag_test_*.csv` | the lag-versus-ineffectiveness diagnostics of the Discussion |

> **Sensitivity to adding the dropped candidates back is not implemented:** the 9 candidates with `used=N` were all dropped because they did not satisfy the A.1 selection criteria (too few years, limited country coverage, etc.), so there is no robustness question about the composite results after adding them back.

---

## 5. How to run

### 5.1 Everything in one go (recommended)

```r
setwd("<path-to-this-repository>")
source("gbf_compute_pipeline.R")
source("gbf_compute_remaining.R")

res  <- run_pipeline()          # the main computation (default n_boot = 10000)
save_intermediates(res)         # writes the full export (CSV + RDS) to data/intermediate/

rest <- run_remaining(res)      # the derived analyses (each block writes to data/output/)
```

**Incremental writing and resuming**:
- `run_pipeline()` writes `data/intermediate/chk_*.rds` (01_meta … 13_ind_projection, 13 checkpoints) as soon as a stage is finished. An error in the middle does not lose the completed stages.
- The checkpoints carry a version number (`chk_00_version.rds`): when the computation logic changes, the old checkpoints are **invalidated automatically** and recomputed, so stale intermediate results can never be loaded.
- `resume = TRUE` is the default: on a rerun the existing checkpoints are loaded and the corresponding stages are skipped. To force the recomputation of a stage, delete the corresponding `chk_*.rds` (or the whole directory); alternatively call `run_pipeline(resume = FALSE)`.
- `run_remaining()` writes `data/output/*.csv` as soon as each table/figure dataset is computed instead of writing everything at the end. If an error occurs, the tables that were already computed are on disk.

> **Runtime warning:** the full `run_pipeline()` includes a 10,000-fold bootstrap for every country and every layer and can take tens of minutes. For a quick trial run use `run_pipeline(n_boot = 500)` (the final results must use 10,000).
> The **dlm** package is installed automatically at run time (if it is missing, `install.packages("dlm")` is attempted; if the installation fails, the code falls back to the built-in equivalent DLM). When the **forecast** package is present, Algorithm 3 uses ETS, otherwise it falls back to a linear model.

### 5.2 Verifying the worked examples (no data needed, seconds)

```r
source("calc_examples.R")       # should print "CALC EXAMPLES PASSED"
```

### 5.3 Internal smoke test (no data needed)

```r
source("gbf_compute_pipeline.R")
source("gbf_compute_remaining.R")
smoke_test()                    # should print "SMOKE TEST PASSED"
```

### 5.4 Targeted re-runs of the long modules

These scripts recompute a single long-running module from the checkpoints already on disk, without repeating the whole pipeline:

```sh
Rscript scripts/run_a2_bootstrap.R     # Algorithm 2: bootstrap imputation uncertainty (about 1-3 minutes)
Rscript scripts/run_a3_bootstrap.R     # Algorithm 3: time-series length uncertainty (about 5-10 minutes)
Rscript scripts/run_a5_perturbation.R  # Algorithm 5: single-indicator perturbation (about 1-2 minutes)
Rscript scripts/run_c3c5_refresh.R     # Algorithm 4 + B.2.2 + B.2.3 (a4 about 2-5 minutes)
```

### 5.5 Figures and tables

```sh
Rscript fig_code/make_fig_data.R       # the figure source data -> data/fig_data/
Rscript fig_code/fig1.R                # ... fig4.R, ed_fig1.R, sup_fig1.R ... sup_fig4.R,
                                       #    sup_figs_5to11.R, sup_fig12.R
Rscript fig_code/csc_null.R            # the CSC null model (about 1-3 minutes; contains an assertion that reproduces csc_country)
Rscript fig_code/sup_fig13.R           # Supplementary Fig. 13 (depends on the products of csc_null.R)
Rscript fig_code/endfill_sensitivity.R # the endpoint-imputation truncation sensitivity (contains an assertion that reproduces Table 6;
                                       #   the new columns of Table 9 depend on this script, so it has to run before make_tables.R)
Rscript fig_code/inherited_trend_test.R # the inherited-trend continuation test (the empirical support of the counterfactual projection framing)
Rscript fig_code/lag_vs_ineffectiveness.R # lag vs ineffective action (the legacy stock / the age of the gains / the within-target signature)
Rscript table_code/make_tables.R        # the source data of the Supplementary Tables -> data/tables/
```

The rendered figures are written to `data/figures/` as PDF, TIFF and PNG (TIFF/PNG when the `ragg` package is available).

---

## 6. Quick reference of the methods and their code symbols

| Item of the Supplementary Information | Where in the code | Key point |
|---|---|---|
| A.2 Imputation | `expand_impute_grid()` + `impute_series()` / `imputation_detail()` | **The rows for the missing years are created first** (the wide-table indicators 18, 39, 26, 30 etc. have missing year columns) → **only the missing values inside the 2011–2024 window are filled** (the values outside the window stay NA); forward/backward fill at the endpoints; linear interpolation inside; **the anchors are all the existing observations (including the history outside the window)**, not only the points inside the period to be imputed; categorical LOCF; the cell-by-cell detail is written to `imputation_detail.csv`; a country is kept when its coverage is ≥ 80% |
| Unit conventions | `UNIT_SCALE` (`load_all_raw()`) | Indicator 18 (Extent of natural ecosystems) was in proportions (0–1) while the meta bounds are percentages (13.57025/99.35725) → **on 2026-10-06 the raw data of `indicator18_ExtNatEco.csv` were corrected directly by a factor of 100**; the `UNIT_SCALE` mechanism is kept in case a unit mismatch has to be registered again |
| Formula 1 | `build_normalized()` | `100×(x−worst)/(best−worst)`, truncated to 0/100; for indicators where lower is better the sign is flipped naturally by the definition of the bounds |
| Formulas 3/4 | `compute_mis()` / `eis_from_values()` | MIS is the arithmetic mean year by year; EIS = 400πS/L² (the area of the radar plot relative to the area of the circle with the same perimeter × 100) |
| A.5 Projection | `fit_pma()` / `project_entity()` | **The fit uses the complete historical series**; **the noise criterion is the SNR** (process variance / observation variance, MLE of a random walk plus noise model, Tittensor et al.): SNR < 1 → **DLM** (local linear trend state space, ref 10; when the `dlm` package is missing the built-in equivalent Kalman filter is used and it is no longer silently skipped); SNR ≥ 1 → **PMA, a genuine multi-model averaging** (Tittensor et al. 2014, Science): a family of 18 parametric models + quad + segmented, multi-start nls, **AICc weights, a weighted average prediction over the confidence set with ΔAICc < 2**, the unconditional variance Σ w_i(σ_i² + (pred_i − pred_avg)²) (nls uses the delta method, which includes the parameter uncertainty); the target years are 2030 / 2050 |
| A.5 Bootstrap | `bootstrap_rate_of_change()` | The model path after 2020 is resampled 10,000 times; the 2.5%/97.5% quantiles form the 95% CI; a change is significant if the interval does not contain 0 |
| Formula 11 | `required_target_value()` | For the threshold targets (T2/T3/T6/T7/T16/T18/T19, see `THRESHOLD_TARGETS`) Y = 100; otherwise Y = v_2020 + |m|·Δt with Y ≥ v_2020 |
| Formula 13 | `compute_ar()` | I = 1{PGR_sim ≥ PGR_req and the bootstrap p < 0.05 (the rate of change is significantly ≠ 0)}; AR is averaged over the layer |
| A.5 | `classify_rho()` | ρ = PGR_sim/PGR_req; Full (ρ ≥ 1) / Favourable (0.5 ≤ ρ < 1) / Mixed / Unfavourable; a ±0.1 perturbation of the thresholds gives the confidence flag; PGR_req = 0 with non-negative growth is recorded as achieved (the T6 case) |
| Formulas 5–9 | `synergy_decompose()` | d = |Δr|; d_t = √(Σ Δr²); S_t = (1/C)(Σ d − d_t); the sign of S_t is compared with the expected direction; CSC = (d − d_t/C)/(S_t·sgnS·C)×100%, and the sum over countries is ±100%; **the country contributions CSC follow the same convention as MIS/EIS/rankings: countries with more than 20% of the indicators missing are removed (A.2)** |
| Baseline parameters | `baseline_params()` | The four parameters over 2011–2020 (mean/SD/slope/range; Supplementary Table 3); KW + MWU-Bonferroni (the post-hoc is only run when the omnibus p < 0.05) |
| Multiple comparisons | `bh_fdr()` | Bonferroni within each family + a global BH-FDR (the note of Supplementary Table 3) |

---

## 7. Dependencies

| Package | Required | Use |
|---|---|---|
| `tidyverse` (dplyr, readr, tidyr, purrr, stringr, tibble) | ✅ | all data processing |
| `readxl` | ✅ | reads `geoclass.xlsx` (the code falls back when it is missing) |
| `ggplot2`, `scales` | ✅ | the figures |
| `patchwork` | ✅ | the multi-panel figures |
| `sf`, `rnaturalearth` | ✅ | the maps |
| `ggrepel` | ✅ | the label placement in the figures |
| `dlm` | optional | the DLM projection of the high-noise series (falls back to the built-in equivalent Kalman filter) |
| `forecast` | optional | ETS / auto-ARIMA in Algorithm 3 (falls back to a linear model) |
| `ragg` | optional | the TIFF export of the figures |

The code was developed with R 4.5.3.

---

## 8. FAQ

**Q1. `find_indicator_file()` reports that it cannot find an indicator file.**
Answer: the files under `origin_data/` must be named `indicator<id>_*.csv` or `indicator<id>.csv` (e.g. `indicator1_RLI.csv`).

**Q2. The continent column of an indicator file disagrees with `geoclass.xlsx`.**
Answer: `geoclass.xlsx` is authoritative (it is the single source of the country → region/group information). The geographic columns inside the indicator files are only a fallback.

**Q3. The China row of the CSC table is empty.**
Answer: check that the row `Geographical region = China` (CHN) exists in `geoclass.xlsx` and that the indicator data contain the country name `China`.

**Q4. The bootstrap is too slow.**
Answer: use `n_boot = 500` for a trial run; the final results must use 10,000 (the Supplementary Information A.5 requires it). `run_remaining(n_boot = 2000)` in `gbf_compute_remaining.R` controls the number of resamples of the derived analyses.

**Q5. How can I check that the code agrees with the formulas of the Supplementary Information?**
Answer: `calc_examples.R` asserts every numerical worked example, and the header comment of each script maps the implemented quantities onto the numbered formulas and the sections of the Supplementary Information (A.2–A.5, B.1–B.2).
