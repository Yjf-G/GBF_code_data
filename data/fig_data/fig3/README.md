# Fig 3 figure source data (Spatiotemporal dynamics of MIS and EIS, 2011-2024)

- fig3_tri_regions.csv: triangular matrix of the pairwise differences between regions (6 regions: Africa/Asia/Europe/Oceania/LAC/North A.);
  diff = column-group mean - row-group mean (X-Y, identical to the values in the figure; 2011 in blue in the upper triangle, 2024 in red in the lower triangle); p_bonf = Dunn post-hoc test (Bonferroni, coloured by significance in the figure)
- fig3_tri_groups.csv: triangular matrix of the same format for the development-status groups (DD/DL/LDC/LLDC/OECD/SIDS)
- fig3_ks.csv: KS test per group (2011 vs 2024 country distributions) and delta_mean = 2024 mean - 2011 mean (the "ks" column; the fill colour of the figure gives the significance of the KS p)
- fig3_map_delta.csv: change of each country in 2024 relative to 2011 (blue->red = negative->positive)
- fig3_paired_year_on_year.csv: year-on-year paired tests (paired Wilcoxon for 4 metrics x 13 year pairs: metric/year pair/mean_diff/p_paired)

Definition: country-level MIS/EIS (both the goal and the target resolution); group means = equal-weight mean over the countries in the group.
