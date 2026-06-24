setwd("plot")
library(ggplot2)
library(dplyr)
library(patchwork)

# ===================== 数据：4 个 Scenario + 4 种方法 =====================
df <- data.frame(
  Scenario = factor(rep(c(4,5,8,9), each = 4)), # 4 个场景，每个 4 种方法
  Method = factor(
    rep(c("Uniform(TITE-BCOMS-L)", "Exponential(TITE-BCOMS-L)", "Uniform(TITE-BCOMS-H)", "Exponential(TITE-BCOMS-H)"), 4),
    levels = c("Uniform(TITE-BCOMS-L)", "Exponential(TITE-BCOMS-L)", "Uniform(TITE-BCOMS-H)", "Exponential(TITE-BCOMS-H)")
  ),
  # 4 个 Scenario × 4 个方法 的 Sel 数值
  Sel = c(
    80.9, 79.1, 81.0, 81.4, # Scenario 4
    73.3, 73.4, 82.2, 82.6, # Scenario 5
    49.6, 53.0, 74.0, 73.2, # Scenario 8
    65.7, 67.4, 70.2, 73.3  # Scenario 9
  ),
  # 4 个 Scenario × 4 个方法 的 Pat 数值
  Pat = c(
    53.7, 52.0, 62.83, 61.93, # Scenario 4
    38.08, 37.49, 47.79, 48.27, # Scenario 5
    39.96, 40.65, 51.28, 52.12, # Scenario 8
    28.66, 28.76, 35.64, 36.50  # Scenario 9
  )
)

# ===================== 4 种方法配色（专业蓝黄配色，区分度高） =====================
colors <- c(
  "Uniform(TITE-BCOMS-L)"  = "#3d6b8a",
  "Exponential(TITE-BCOMS-L)"  = "#a8c4d4",
  "Uniform(TITE-BCOMS-H)" = "#9a6b1f",
  "Exponential(TITE-BCOMS-H)" = "#e8c97a"
)

# ===================== 绘图函数 =====================
make_plot <- function(data, y_var, y_label) {
  ggplot(data, aes(x = Scenario, y = .data[[y_var]], fill = Method)) +
    geom_col(position = position_dodge(width = 0.9), color = "black", na.rm = TRUE) +
    scale_fill_manual(values = colors) +
    labs(x = "Scenario", y = y_label) +
    theme_classic() +
    theme(
      legend.position = "top",
      legend.key.width = unit(1.5, "cm"),
      plot.title = element_text(hjust = 0)
    )+
    # 强制图例2行，按顺序填充：前3个第一行，后3个第二行
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
}

# ===================== 出图 =====================
p1 <- make_plot(df, "Sel", "Selection of OBDC (%)")
p2 <- make_plot(df, "Pat", "Patients allocated to OBDC (%)")

p_final <- p1 / p2
pdf("output/sensitivity_arrivaldist.pdf", width = 12, height = 9)
print(p_final)
dev.off()
