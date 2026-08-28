setwd("D:/Secondary Metabolites")

library(ggplot2)
library(ggbreak)
library(tidyr)

# 1. Read raw data
my_data <- read.csv("Individual_phenolics.csv", check.names = FALSE)
colnames(my_data)[1] <- "Compounds"

# 2. Reshape to long format 
long_data <- pivot_longer(
  my_data,
  cols = -Compounds,
  names_to = "Sample",
  values_to = "Value"
)

# 3. Assign Group based on first letter of Sample name
long_data$Group <- ifelse(substr(long_data$Sample, 1, 1) == "C", "Control", "Treatment")
long_data$Group <- factor(long_data$Group, levels = c("Control", "Treatment"))

# 4. Calculate Mean and SE per Compound x Group
mean_data <- aggregate(Value ~ Compounds + Group, data = long_data, FUN = mean)
sd_data   <- aggregate(Value ~ Compounds + Group, data = long_data, FUN = sd)
n_data    <- aggregate(Value ~ Compounds + Group, data = long_data, FUN = length)

summary_data <- mean_data
names(summary_data)[3] <- "Mean"
summary_data$SD <- sd_data$Value
summary_data$N  <- n_data$Value
summary_data$SE <- summary_data$SD / sqrt(summary_data$N)

# 5. Plot 
ggplot(summary_data, aes(x = Compounds, y = Mean, fill = Group)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.7,
           color = "black") +
  geom_errorbar(
    aes(ymin = Mean - SE, ymax = Mean + SE),
    position = position_dodge(width = 0.7),
    width = 0.25,
    linewidth = 0.4
  ) +
  
  scale_y_break(c(300, 800, 2000, 55000), scales = 0.25, space = 0.03, expand = FALSE) +
  scale_y_continuous(
    limits = c(0, 120000),
    breaks = c(0, 50, 100, 150, 200, 250, 1000, 1500, 60000, 90000, 120000),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  scale_fill_manual(
    values = c(
      "Control" = "#00B050",
      "Treatment" = "#002060"
    ),
    labels = c(
      "Control" = "Control",
      "Treatment" = "0.1 mA Electrical Treatment"
    )
  ) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black", hjust = 0.5),
    axis.title.y = element_text(color = "black", hjust = 0.5),
    
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(color = "black", size = 9),
    legend.key.size = unit(0.4, "cm"),      
    legend.spacing.x = unit(0.4, "cm"),
    
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right = element_blank(),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  
  labs(
    x = "Compounds",
    y = "Phenolic content (ng/g)"
  )

# 6. Save figure 
ggsave("Phenolics_BrokenAxis.png", width = 8, height = 6, dpi = 600)
ggsave("Phenolics_BrokenAxis.pdf", width = 8, height = 6)
