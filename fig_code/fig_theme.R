# =====================================================================
# fig_theme.R — Shared styling and output utilities for Nature submission figures
#   Specification highlights (Nature artwork guide):
#   - Width: two-column 183 mm / single-column 83–89 mm; Arial font 5–7 pt
#   - Output: vector PDF + 300 dpi TIFF + 300 dpi PNG; white background; thin axis lines; no heavy boxes or grids
#   - Colours: colour-blind friendly (Okabe–Ito); GBF section colours match the figure legend
#     Blue = reducing threats, green = meeting people's needs, orange = tools and solutions
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
})

SEC_COLORS <- c(reducing_threats = "#0072B2", meeting_needs = "#009E73",
                tools_solutions = "#E69F00", Goals = "#4D4D4D")
SEC_LABELS <- c(reducing_threats = "Reducing threats",
                meeting_needs = "Meeting people's needs",
                tools_solutions = "Tools and solutions", Goals = "Goals")
GOAL_COLORS <- c(A = "#0072B2", B = "#D55E00", C = "#009E73", D = "#CC79A7")
PAL8 <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9",
          "#D55E00", "#999999", "#000000")

theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(family = "Arial", size = base_size, colour = "black"),
      axis.title = element_text(size = base_size, colour = "black"),
      axis.text = element_text(size = base_size - 1, colour = "black"),
      axis.line = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "black"),
      legend.title = element_text(size = base_size - 0.5),
      legend.text = element_text(size = base_size - 1),
      legend.key.size = unit(3.2, "mm"),
      legend.background = element_blank(),
      legend.margin = margin(0, 0, 0, 0.5),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.title = element_text(size = base_size, face = "bold"),
      plot.tag = element_text(size = base_size + 1, face = "bold"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.background = element_blank(),
      panel.background = element_blank(),
      plot.margin = margin(1.5, 2, 1.5, 1.5)
    )
}

# Save PDF (vector) + TIFF (300 dpi) + PNG (300 dpi)
save_nature <- function(p, name, width_mm, height_mm, dir = "data/figures") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  ggsave(file.path(dir, paste0(name, ".pdf")), p,
         width = w, height = h, units = "in", device = cairo_pdf)
  ragg::agg_tiff(file.path(dir, paste0(name, ".tiff")),
                 width = w, height = h, units = "in", res = 300,
                 background = "white")
  print(p)
  dev.off()
  ragg::agg_png(file.path(dir, paste0(name, ".png")),
                width = w, height = h, units = "in", res = 300,
                background = "white")
  print(p)
  dev.off()
  message("[figure] ", name, "  ", width_mm, "x", height_mm, " mm -> ", dir)
}

# Common small utilities
scale_sec_fill <- function(...) {
  scale_fill_manual(values = SEC_COLORS, labels = SEC_LABELS, name = "Section", ...)
}
scale_sec_color <- function(...) {
  scale_colour_manual(values = SEC_COLORS, labels = SEC_LABELS, name = "Section", ...)
}
