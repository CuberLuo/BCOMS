setwd("plot")
library(ggplot2)
library(dplyr)
library(patchwork)

df <- data.frame(
  Scenario = factor(rep(c(4, 5, 8, 9), each = 14)),
  Method = factor(
    rep(c(
      "samplesize30-high", "samplesize30-low",
      "samplesize36-high", "samplesize36-low",
      "samplesize42-high", "samplesize42-low",
      "samplesize48-high", "samplesize48-low",
      "samplesize54-high", "samplesize54-low",
      "samplesize60-high", "samplesize60-low",
      "samplesize66-high", "samplesize66-low"
    ), 4),
    levels = c(
      "samplesize30-high", "samplesize30-low",
      "samplesize36-high", "samplesize36-low",
      "samplesize42-high", "samplesize42-low",
      "samplesize48-high", "samplesize48-low",
      "samplesize54-high", "samplesize54-low",
      "samplesize60-high", "samplesize60-low",
      "samplesize66-high", "samplesize66-low"
    )
  ),
  Sel = c(
    73.6, 71, 77.4, 74.9, 78.5, 76.4, 82.3, 78.3, 77.9, 78.8, 81.4, 81.1, 82.2, 79,
    75.4, 77.6, 82.5, 77.7, 84, 78.9, 81.7, 81.9, 85.3, 82.1, 86.2, 81.3, 85.8, 81.9,
    84.4, 58.7, 82.8, 57.5, 83.8, 57.6, 82.5, 56.9, 85, 58.3, 83.6, 60.2, 84.4, 55.6,
    62.8, 61.5, 70.1, 66.4, 68.5, 69.9, 72.5, 69, 71.2, 70.9, 73.7, 70.2, 76.3, 71.6
  ),
  Pat = c(
    59.89, 47.18, 60.26, 49.75, 62.69, 51.08, 62.16, 50.72, 64.26, 51.32, 64.66, 53.73, 65.2, 53.79,
    52.82, 40.99, 53.46, 40.69, 55.8, 42.46, 56.04, 43.17, 57.94, 44.19, 60.2, 45.12, 58.72, 45,
    56.1, 42.74, 57.44, 42.67, 57.86, 42.55, 58.85, 43.99, 59.52, 44.6, 59.15, 45.09, 59.96, 44.7,
    38.39, 29.5, 40.39, 31.12, 41.39, 32.33, 43.26, 33.48, 44.26, 34.41, 46.14, 35.02, 47.05, 35.48
  )
) %>%
  mutate(
    Group = ifelse(grepl("high$", Method), "High", "Low"),
    SampleSize = factor(gsub("-.*", "", Method),
                        levels = c(
                          "samplesize30", "samplesize36", "samplesize42",
                          "samplesize48", "samplesize54", "samplesize60", "samplesize66"
                        ),
                        labels = c("30", "36", "42", "48", "54", "60", "66")
    ),
    Scenario = factor(paste0("Scenario", Scenario),
                      levels = c("Scenario4", "Scenario5", "Scenario8", "Scenario9")
    )
  )

line_cols <- c(
  "Scenario4" = "#E63946",
  "Scenario5" = "#457B9D",
  "Scenario8" = "#F1A208",
  "Scenario9" = "#2A9D8F"
)

point_shapes <- c(
  "Scenario4" = 16,   # ● 实心圆
  "Scenario5" = 17,   # ▲ 实心三角
  "Scenario8" = 15,   # ■ 实心方块
  "Scenario9" = 23    # ◆ 实心菱形（fill填色）
)

plot_fun <- function(dat, yvar, ylab, ylim_vals = NULL) {
  p <- ggplot(dat, aes(
    x = SampleSize,
    y = .data[[yvar]],
    group = Scenario,
    color = Scenario,
    shape = Scenario,
    fill  = Scenario
  )) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 5, stroke = 0.8) +
    scale_color_manual(values = line_cols) +
    scale_fill_manual(values = line_cols) +
    scale_shape_manual(values = point_shapes) +
    guides(
      fill  = "none",
      color = guide_legend(override.aes = list(
        size  = 4,
        shape = c(16, 17, 15, 23),
        fill  = c("#E63946", "#457B9D", "#F1A208", "#2A9D8F")
      )),
      shape = "none"   # 只保留color图例，shape图例合并进去
    ) +
    labs(
      x = "Sample Size", y = ylab,
      color = "Scenario", shape = "Scenario"
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position    = "bottom",
      legend.key.width   = unit(1.5, "cm"),
      legend.key.size    = unit(0.6, "cm"),
      legend.text        = element_text(size = 11),
      axis.title         = element_text(size = 12),
      axis.text          = element_text(size = 11),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.4)
    )
  if (!is.null(ylim_vals)) {
    p <- p + coord_cartesian(ylim = ylim_vals)
  }
  p
}

p_sel_high <- plot_fun(filter(df, Group == "High"), "Sel", "(BCOMS-H) Selection of OBDC (%)", c(50, 100))
p_sel_low  <- plot_fun(filter(df, Group == "Low"),  "Sel", "(BCOMS-L) Selection of OBDC (%)", c(50, 100))
p_pat_high <- plot_fun(filter(df, Group == "High"), "Pat", "(BCOMS-H) Patients allocated to OBDC (%)", c(30, 100))
p_pat_low  <- plot_fun(filter(df, Group == "Low"),  "Pat", "(BCOMS-L) Patients allocated to OBDC (%)", c(30, 100))

p_all <- (p_sel_high | p_sel_low) / (p_pat_high | p_pat_low) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

pdf("output/sensitivity_samplesize.pdf", width = 12, height = 9)
print(p_all)
dev.off()