# fig_code — figure source data and rendering

`make_fig_data.R` prepares the source data for every figure, and the individual scripts render the figures. The products are written to `data/fig_data/` (the source data, one folder per figure, each with a README describing its CSVs) and to `data/figures/` (the rendered figures as PDF, TIFF and 300 dpi PNG; the TIFF and PNG files are written when the `ragg` package is available).

All paths are relative to the repository root; run every script from there.

## Figures and scripts

| Figure | Script | Content |
|---|---|---|
| Fig. 1 | `fig1.R` | Baseline across the GBF's goals and targets (2011–2020): (a) the slopes of the 23 targets with 95% CI, (b) the global baseline trend lines of goals A–D, (c) the slope × R² quadrant |
| Fig. 2 | `fig2.R` | Trajectory projections fall well short of the GBF ambition: (a,b) the country-level four-class shares per target/goal, (c) the achievement rate with 50%/95% intervals, (d–g) the PGR gap maps and country-share distributions at the target (d,e) and goal (f,g) levels |
| Fig. 3 | `fig3.R` | Spatiotemporal dynamics of MIS and EIS, 2011–2024: the triangular matrices of the pairwise differences, the KS tests and the per-country Δ map |
| Fig. 4 | `fig4.R` | Trend alignment of 41 indicators and country synergy contribution (CSC), 2011–2024: (a) the annual range of the aligned synergy effect per indicator, (b) the share of years with a positive synergy, (c) the map of the country contributions |
| Extended Data Fig. 1 | `ed_fig1.R` | Country counts by trend and achievement class for each target/goal and the transitions between them |
| Supplementary Fig. 1 | `sup_fig1.R` | Heterogeneity of the baseline summary statistics across continents and development-status groups: the omnibus tests and the pairwise significance matrices |
| Supplementary Fig. 2 | `sup_fig2.R` | Projected achievement gaps by geography and development status |
| Supplementary Fig. 3 | `sup_fig3.R` | Global trends of MIS and EIS, 2011–2024: annual means, year-on-year changes and Wilcoxon significance |
| Supplementary Fig. 4 | `sup_fig4.R` | Group mean trajectories of MIS and EIS, 2011–2024 |
| Supplementary Fig. 5 | `sup_figs_5to11.R` | The uncertainty of the indicator subset (Algorithm 1) |
| Supplementary Fig. 6 | `sup_figs_5to11.R` | The uncertainty of the imputation (Algorithm 2) |
| Supplementary Fig. 7 | `sup_figs_5to11.R` | The uncertainty of the time-series length (Algorithm 3) |
| Supplementary Fig. 8 | `sup_figs_5to11.R` | The uncertainty of the country sample (Algorithm 4) |
| Supplementary Fig. 9 | `sup_figs_5to11.R` | The perturbation of a single indicator (Algorithm 5) |
| Supplementary Fig. 10 | `sup_figs_5to11.R` | The weight sensitivity (SI B.2.2) |
| Supplementary Fig. 11 | `sup_figs_5to11.R` | The endpoint robustness (SI B.2.3) |
| Supplementary Fig. 12 | `sup_fig12.R` | The representative model fits |
| Supplementary Fig. 13 | `csc_null.R` + `sup_fig13.R` | The CSC null-model distributions (exchangeability / sign symmetry) |

`fig_theme.R` holds the shared ggplot theme and the `save_nature()` helper that writes the PDF/TIFF files into `data/figures/`.

## Run order

```sh
Rscript fig_code/make_fig_data.R        # the figure source data -> data/fig_data/
Rscript fig_code/fig1.R                 # ... fig4.R, ed_fig1.R, sup_fig1.R ... sup_fig4.R,
                                        #    sup_figs_5to11.R, sup_fig12.R
Rscript fig_code/csc_null.R             # the CSC null model (about 1-3 minutes; asserts that it reproduces csc_country)
Rscript fig_code/sup_fig13.R            # Supplementary Fig. 13 (depends on the products of csc_null.R)
Rscript fig_code/endfill_sensitivity.R  # the endpoint-imputation truncation sensitivity (asserts that it reproduces Table 6;
                                        #   the new columns of Table 9 depend on this script, so it must run before table_code/make_tables.R)
Rscript fig_code/inherited_trend_test.R # the inherited-trend continuation test (the empirical support of the counterfactual projection framing)
Rscript fig_code/lag_vs_ineffectiveness.R # lag vs ineffective action (the legacy stock / the age of the gains / the within-target signature)
```

## Notes

- The folders and file names follow the final figure numbering and the wording of the final figure titles.
- The pipeline output **`uncertainty__c3_weight_sensitivity`** keeps its internal product name (its figure-source-data copy is `supfig10_weight_sensitivity.csv`).
- The indicator typing in `data/indicator_dpsir.csv` (Process / Response / Pressure / State) separates the action-type from the state-type diagnostics of SI B.2.4–B.2.6; `lag_vs_ineffectiveness.R` recomputes its within-target slopes from this typing. (Panels a and b of Fig. 4 were originally coloured by this typing; the final figure uses one common colour with zebra row bands.)
- The old Extended Data Fig. 2 (the China case study) and the DPSIR typology figure were removed after the manual review; no script in this folder generates them any more.
