# Sup Fig. 7: time-series length uncertainty (Algorithm 3, block-bootstrap + ETS/ARIMA/RW)

- supfig7_window_length.csv: optimal window L_opt per indicator (Kneedle knee; for flat or non-convex curves the fallback L_der is used, flagged by knee_found), the derivative-based L_der, the series length n_years, sigma(L_opt), the L_opt per model and model_spread
- supfig7_sigma_curves.csv: per-indicator sigma_i(L) curve (pooled + ETS/ARIMA/RW; ci95_halfwidth = 1.96*sigma)
- supfig7_cross_model.csv: cross-model consistency table (L_opt_kneedle / L_opt_der / sigma(L_opt) for indicator x model)
