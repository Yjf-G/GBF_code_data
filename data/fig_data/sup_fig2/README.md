# Supplementary Fig. 2 figure source data (projected achievement gaps by geography and development status)

- Definition: the country-level achievement status follows A.5 / Extended Data Fig. 1 (reusing pc; for T6 a zero requirement with non-negative growth is recorded as Ac; for T17 countries with a zero baseline are excluded according to Formula 12);
  the grouping follows Sup Fig. 1 (5 continent groups; 6 status groups DD/DL/LDC/LLDC/SIDS/OECD, a country can belong to several groups -> the status-group tests are an approximation)
- supfig2_achv_by_group.csv: Ac/Unc/Unac counts and shares for the pooled targets/goals x group (unit: country-layer pair; n_countries = number of countries in the group)
- supfig2_unac_by_layer_group.csv: status counts and the Unac share for the 27 layers x groups (the n_countries inside a cell = number of countries in that layer and group)
- supfig2_kw_achv.csv: per-layer Kruskal-Wallis (achievement order Unac=0/Unc=1/Ac=2) grouped by continent/status, kw_H/kw_p/n/n_groups
