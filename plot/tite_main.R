setwd("plot")
library(ggplot2)
library(dplyr)
library(patchwork)

df <- data.frame(
  Scenario = factor(rep(1:12, each = 4)),
  Method = factor(rep(c("BCOMS-L", "TITE-BCOMS-L", "BCOMS-H", "TITE-BCOMS-H"), 12),
    levels = c("BCOMS-L", "TITE-BCOMS-L", "BCOMS-H", "TITE-BCOMS-H")
  ),
  # Sel = select OBDC correctly (%)
  Sel = c(
    86.4, 87.7, 64.1, 82.3,   # Scenario 1  
    84.8, 81.9, 90.6, 87,     # Scenario 2  
    78.7, 76.2, 71.7, 73.9,   # Scenario 3  
    77.1, 80.9, 76.6, 81,     # Scenario 4  
    77.8, 73.3, 78.6, 82.2,   # Scenario 5  
    74.9, 72.2, 77.1, 74.4,   # Scenario 6  
    72, 66.3, 87.5, 83.8,     # Scenario 7  
    56.6, 49.6, 84.5, 74,     # Scenario 8  
    69.2, 65.7, 71.7, 70.2,   # Scenario 9  
    70.1, 73.1, 55.3, 66,     # Scenario 10 
    NA, NA, NA, NA,           # Scenario 11
    NA, NA, NA, NA            # Scenario 12
  ),
  # Pat = average percent of patients treated at OBDC correctly (%)
  Pat = c(
    56, 59.44, 48.58, 56.29,   # Scenario 1
    62.65, 60.4, 80.38, 72.76, # Scenario 2
    69.31, 67.51, 65.83, 65.49,# Scenario 3
    49.17, 53.7, 60.84, 62.83, # Scenario 4
    40.75, 38.08, 52.68, 47.79,# Scenario 5
    41.06, 39.9, 49.94, 43.75, # Scenario 6
    41.38, 40.75, 51.04, 49.37,# Scenario 7
    42.98, 39.96, 57.89, 51.28,# Scenario 8
    31.54, 28.66, 41.3, 35.64, # Scenario 9
    33.46, 33.89, 30, 31.77,   # Scenario 10
    NA, NA, NA, NA,            # Scenario 11
    NA, NA, NA, NA             # Scenario 12
  ),
  # Duration = average trial duration
  Duration = c(
    16.5, 14.4, 13.6, 12,     # Scenario 1  
    16.3, 13.8, 13.1, 10.8,   # Scenario 2  
    16.1, 13.6, 13.3, 11.4,   # Scenario 3  
    16.4, 14.2, 13.4, 11.6,   # Scenario 4  
    16.5, 14.3, 13.4, 11.3,   # Scenario 5  
    16.2, 13.5, 13.2, 11.3,   # Scenario 6  
    23.7, 20.8, 20.1, 17.7,   # Scenario 7  
    19.7, 16.8, 22.2, 19,     # Scenario 8  
    22.2, 19, 18.8, 16.4,     # Scenario 9  
    20.4, 17.4, 19.2, 17.2,   # Scenario 10 
    7.3,  7.5,  12.4, 10.6,   # Scenario 11 
    23.7, 20.2, 20.2, 17.4    # Scenario 12 
  )
)

colors <- c(
  "BCOMS-L"       = "#a8c4d4",
  "TITE-BCOMS-L"  = "#3d6b8a",
  "BCOMS-H"      = "#e8c97a",
  "TITE-BCOMS-H" = "#9a6b1f"
)

make_plot <- function(data, y_var, y_label) {
  ggplot(data, aes(x = Scenario, y = .data[[y_var]], fill = Method)) +
    geom_col(position = position_dodge(width = 0.9), color = "black", na.rm = TRUE) +
    scale_fill_manual(values = colors) +
    labs(x = "Scenario", y = y_label) +
    theme_classic() +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0)
    )
}

p1 <- make_plot(df, "Sel", "Selection of OBDC (%)")
p2 <- make_plot(df, "Pat", "Patients allocated to OBDC (%)")
p3 <- make_plot(df, "Duration", "Trial duration (months)")

pdf("output/tite_main.pdf", width = 12, height = 9)
print(p1 / p2 / p3)
dev.off()
