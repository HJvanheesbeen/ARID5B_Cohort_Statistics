# ============================================================
#  Total Social Time (Summed Across Experiments)
#  Comparison Between WT and Arid5bQ522* Mice
#  Two-tailed t-test with outlier detection
# ============================================================


library(dplyr)
library(ggplot2)

# --- 1. Load data ---
data_file <- "**Path**/TotalSocialTime.csv"
data_raw <- read.csv(data_file, stringsAsFactors = FALSE)

# --- 2. Rename for clarity ---
data <- data_raw %>%
  rename(
    Genotype = 1,
    Exp1 = 2,
    Familiar_Ex2 = 3,
    Novel_Ex2 = 4
  ) %>%
  mutate(Genotype = factor(Genotype, levels = c("WT", "Q522*")))

# --- 3. Compute total social time (sum across all experiments) ---
data <- data %>%
  mutate(Total_social_time = Exp1 + Familiar_Ex2 + Novel_Ex2)

# --- 4. Identify and remove outliers (1.5x IQR rule per genotype) ---
detect_outliers <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower <- q1 - 1.5 * iqr
  upper <- q3 + 1.5 * iqr
  (x < lower) | (x > upper)
}

data <- data %>%
  group_by(Genotype) %>%
  mutate(is_outlier = detect_outliers(Total_social_time)) %>%
  ungroup()

outliers <- data %>% filter(is_outlier)
data_clean <- data %>% filter(!is_outlier)

# --- 5. Perform two-tailed independent t-test ---
ttest_result <- t.test(Total_social_time ~ Genotype, data = data_clean, var.equal = FALSE)

# --- 6. Extract results ---
df <- ttest_result$parameter
t_value <- ttest_result$statistic
p_value <- ttest_result$p.value
means <- data_clean %>%
  group_by(Genotype) %>%
  summarise(
    mean = mean(Total_social_time),
    sd = sd(Total_social_time),
    n = n()
  )

# --- 7. Save text summary---
output_text <- sprintf(
  "Two-sample (two-tailed) t-test comparing total social interaction time between genotypes (WT vs Arid5bQ522*):\n
WT (mean ± SD, n) = %.2f ± %.2f (n = %d)
Q522* (mean ± SD, n) = %.2f ± %.2f (n = %d)
t(%.1f) = %.3f, p = %.4f

Outliers removed: %d
Outlier values (if any): %s
",
  means$mean[1], means$sd[1], means$n[1],
  means$mean[2], means$sd[2], means$n[2],
  df, t_value, p_value,
  nrow(outliers),
  ifelse(nrow(outliers) > 0,
         paste(round(outliers$Total_social_time, 2), collapse = ", "),
         "None")
)

output_file <- "**Path**/TotalSocialTime_ttest_summary_outliers.txt"
writeLines(output_text, output_file)

cat("✅ Statistical summary saved to:\n", output_file, "\n")

# --- 8. Plot and save SVG ---
plot_file <- "**Path**/TotalSocialTime_boxplot_outliers_removed.svg"

p <- ggplot(data_clean, aes(x = Genotype, y = Total_social_time, fill = Genotype)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.8, size = 2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Total social interaction time (summed across experiments)",
    y = "Total social time (s)",
    x = ""
  ) +
  scale_fill_manual(values = c("WT" = "#F6BE00", "Q522*" = "#E97451")) +
  coord_cartesian(ylim = c(0, 400)) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

ggsave(plot_file, plot = p, width = 5, height = 5, device = "svg")

cat("✅ SVG plot saved to:\n", plot_file, "\n")
