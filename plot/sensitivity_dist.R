setwd("plot")
library(ggplot2)
library(dplyr)
library(patchwork)

# ===================== 数据：4 个 Scenario + 6 种方法 =====================
df <- data.frame(
    Scenario = factor(rep(c(4,5,8,9), each = 6)), # 4 个场景，每个 6 种方法
    Method = factor(
        rep(c(
            "Uniform(TITE-BCOMS-L)",
            "Weibull(TITE-BCOMS-L)",
            "Log-Logistic(TITE-BCOMS-L)",
            "Uniform(TITE-BCOMS-H)",
            "Weibull(TITE-BCOMS-H)",
            "Log-Logistic(TITE-BCOMS-H)"
        ), 4),
        levels = c(
            "Uniform(TITE-BCOMS-L)",
            "Weibull(TITE-BCOMS-L)",
            "Log-Logistic(TITE-BCOMS-L)",
            "Uniform(TITE-BCOMS-H)",
            "Weibull(TITE-BCOMS-H)",
            "Log-Logistic(TITE-BCOMS-H)"
        )
    ),
    # 4 个 Scenario × 6 个方法 的 Sel 数值
    Sel = c(
        79.1, 80.9, 78.1, 81.3, 81, 81.3, # Scenario 4
        71.9, 73.3, 75.3, 81.5, 82.2, 81, # Scenario 5
        48.5, 49.6, 52.7, 72.6, 74, 72.5, # Scenario 8
        64.8, 65.7, 65.1, 74.3, 70.2, 72.2 # Scenario 9
    ),
    # 4 个 Scenario × 6 个方法 的 Pat 数值
    Pat = c(
        52.24, 53.7, 52.22, 61.64, 62.83, 62.52, # Scenario 4
        37.37, 38.08, 37.75, 47.14, 47.79, 47.23, # Scenario 5
        40.03, 39.96, 40.78, 51.02, 51.28, 50.92, # Scenario 8
        28.43, 28.66, 28.32, 36.01, 35.64, 37.08 # Scenario 9
    )
)
# ===================== 6 种方法配色（专业蓝黄配色，区分度高） =====================
colors <- c(
    "Uniform(TITE-BCOMS-L)"  = "#3d6b8a",
    "Weibull(TITE-BCOMS-L)"  = "#74a0c0",
    "Log-Logistic(TITE-BCOMS-L)"  = "#a8c4d4",
    "Uniform(TITE-BCOMS-H)" = "#9a6b1f",
    "Weibull(TITE-BCOMS-H)" = "#c49b4d",
    "Log-Logistic(TITE-BCOMS-H)" = "#e8c97a"
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
            legend.key.width = unit(1.2, "cm"), # 加宽图例
            plot.title = element_text(hjust = 0)
        )+
        # 强制图例2行，按顺序填充：前3个第一行，后3个第二行
        guides(fill = guide_legend(nrow = 2, byrow = TRUE))
}

# ===================== 出图 =====================
p1 <- make_plot(df, "Sel", "Selection of OBDC (%)")
p2 <- make_plot(df, "Pat", "Patients allocated to OBDC (%)")

p_final <- p1 / p2
pdf("output/sensitivity_dist.pdf", width = 12, height = 9)
print(p_final)
dev.off()
