# data — inputs, intermediate results and source data

Everything the R pipeline reads and writes lives here. All the paths used by the scripts are relative to the repository root.

## Input files

| File | Content |
|---|---|
| `meta.csv` | The metadata of the 53 candidate indicators: `indicator_id, indicator_name, used, period, countries, worst, best, direction`. The 44 indicators with `used = Y` enter the main analysis; `worst`/`best` are the 0-point and 100-point normalization bounds; `direction` is `+` (higher is better) or `−` (lower is better) |
| `indicator_to_target.csv` | The indicator → target/goal mapping: `indicator_id, indicator_code, target, direction, section, goal`. The `target` column mixes the goal mappings `A`–`D` with the target mappings `1`–`23`; an indicator may be mapped to several targets |
| `target_section.csv` | The target ↔ section (`reducing_threats` / `meeting_needs` / `tools_solutions`) ↔ goal ↔ theme mapping for the 4 goals and the 23 targets |
| `indicator_dpsir.csv` | The indicator typing of each indicator (`Process` / `Response` / `Pressure` / `State`, with an alternative `dpsir_class_alt` coding); it defines the action-type/state-type diagnostics of SI B.2.4–B.2.6 |
| `geoclass.xlsx` | The UN M49 country classification (288 rows × 19 columns). It is the authoritative source of the country → continent/subcontinent/region/group assignments, of the 193-country universe (`United Nations == "nation"`), of the development status and of the SIDS/LDC/LLDC/OECD/EU flags |
| `origin_data/` | The 53 raw indicator files, one CSV per indicator (`indicator<id>_<code>.csv`). They are the harmonized extracts of the official sources (GBF, SDG, FAOSTAT, …) |

### Conventions of the raw indicator files

- Both layouts are accepted: **wide** (the first column is `GeoAreaName`, followed by four-digit year columns) and **long** (with `year` and `value` columns).
- All historical years are kept: the imputation only fills the missing values inside the 2011–2024 window, but the projection fit needs the complete history to identify the pattern.
- The geographic group columns inside the indicator files are optional; the grouping follows `geoclass.xlsx`.
- The regional indicators 27–29 (RLI food/medicinal species, RLI traded species, RLI impact of invasive alien species) have continent/region series rather than country series. They enter the global/regional analyses only and are not subject to the 193-country universe.

## Products of the pipeline

| Folder | Content |
|---|---|
| `intermediate/` | The intermediate products of `run_pipeline()` (75 files): the country metadata, the normalized values at the indicator level, the full imputation record and the cell-by-cell imputation detail, the target/goal values (Formula 2), the complete target/goal and MIS/EIS time series at the country, regional and global level, the baseline parameters (Supplementary Table 3), the projections with their 95% CI, the PGR family (Formulas 10–12), the achievement rate (Formula 13), the trajectory classification (A.5), the country synergy contributions (Formula 9), the indicator-level projections and the country coverage. It also holds the `chk_*.rds` checkpoints, which allow the pipeline to be resumed without repeating the completed stages |
| `output/` | The products of `run_remaining()` (84 files): the source data of Supplementary Tables 3–10 (internal products `table_S1`–`table_S6`), the Bonferroni/FDR tables, the section/goal summaries (Fig. 2c), the trajectory analyses (B.2.7), the achievement-rate intervals (A.5), the extrapolation and leave-one-year-out outputs (B.2.3), the five uncertainty algorithms (B.1/B.2.1), the weight sensitivity (B.2.2), the CSC null model, the endpoint-imputation sensitivity (B.2.4), the inherited-trend and lag diagnostics (B.2.6), and the figure source data of Figs 1/3/4 |
| `fig_data/` | The figure source data, one folder per figure, each with a `README.md` describing the columns of its CSVs: `fig1`–`fig4`, `ed_fig1`, `sup_fig1`–`sup_fig13` |
| `tables/` | The source data of the Supplementary Tables 3–11, one CSV per table |
| `figures/` | The rendered figures (PDF for the vector version; TIFF and PNG at 300 dpi, written by `fig_code/fig_theme.R` when the `ragg` package is available) |

The `intermediate/` and `output/` folders are committed so that the tables and the figures can be regenerated without repeating the 10,000-fold bootstrap of the main pipeline; delete them if you prefer to recompute everything from the raw data.
