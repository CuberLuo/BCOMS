setwd("plot")
library(ggplot2)
library(dplyr)
library(patchwork)

# ===================== 数据：4 个 Scenario + 6 种方法 =====================
df <- data.frame(
    Scenario = factor(rep(c(4,5,8,9), each = 6)), # 4 个场景 × 6 种方法
    Method = factor(
        rep(c(
            "Preset Prior(BCOMS-L)",
            "Strong Prior(BCOMS-L)",
            "Weak Prior(BCOMS-L)",
            "Preset Prior(BCOMS-H)",
            "Strong Prior(BCOMS-H)",
            "Weak Prior(BCOMS-H)"
        ), 4),
        levels = c(
            "Preset Prior(BCOMS-L)",
            "Strong Prior(BCOMS-L)",
            "Weak Prior(BCOMS-L)",
            "Preset Prior(BCOMS-H)",
            "Strong Prior(BCOMS-H)",
            "Weak Prior(BCOMS-H)"
        )
    ),
    # 4×6=24 个 Sel 数值
    Sel = c(
        77.1, 76.0, 76.3, 76.6, 78.8, 77.2, # Scenario 4
        77.8, 75.9, 74.9, 78.6, 83.9, 83.3, # Scenario 5
        56.6, 55.2, 54.8, 84.5, 83.2, 85.0, # Scenario 8
        69.2, 69.8, 68.8, 71.7, 71.2, 71.4  # Scenario 9
    ),
    # 4×6=24 个 Pat 数值
    Pat = c(
        49.17, 49.25, 49.16, 60.84, 60.34, 60.41, # Scenario 4
        40.75, 39.56, 39.51, 52.68, 53.93, 53.53, # Scenario 5
        42.98, 42.69, 42.55, 57.89, 57.34, 57.43, # Scenario 8
        31.54, 31.17, 31.04, 41.30, 40.43, 41.02  # Scenario 9
    )
)

# ===================== 6 种方法配色（专业蓝黄配色，区分度高） =====================
colors <- c(
    "Preset Prior(BCOMS-L)"    = "#3d6b8a",
    "Strong Prior(BCOMS-L)"      = "#74a0c0",
    "Weak Prior(BCOMS-L)"   = "#a8c4d4",
    "Preset Prior(BCOMS-H)"   = "#9a6b1f",
    "Strong Prior(BCOMS-H)"     = "#c49b4d",
    "Weak Prior(BCOMS-H)"  = "#e8c97a"
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
            legend.key.width = unit(1.3, "cm"),
            plot.title = element_text(hjust = 0)
        )+
        # 强制图例2行，按顺序填充：前3个第一行，后3个第二行
        guides(fill = guide_legend(nrow = 2, byrow = TRUE))
}

# ===================== 出图 =====================
p1 <- make_plot(df, "Sel", "Selection of OBDC (%)")
p2 <- make_plot(df, "Pat", "Patients allocated to OBDC (%)")

p_final <- p1 / p2
pdf("output/sensitivity_prior.pdf", width = 12, height = 9)
print(p_final)
dev.off()
