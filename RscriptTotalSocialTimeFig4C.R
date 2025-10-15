# ============================================================
#  Total Social Time (Summed Across Experiments)
#  Comparison Between WT and Arid5bQ522* Mice
#  Two-tailed t-test with Grubbs outlier detection
# ============================================================

library(dplyr)
library(ggplot2)
library(outliers)  # Make sure this is installed

# --- 1. Load data ---
data_file <- "***Path***/NumericalData/TotalSocialTime.csv"
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

# --- 4. Identify and remove outliers per genotype using Grubbs' test ---
remove_grubbs_outliers <- function(x, alpha = 0.05) {
  x_clean <- x
  outlier_vals <- c()
  
  repeat {
    if(length(x_clean) < 3) break  # Grubbs requires >=3 values
    grubbs_test <- try(grubbs.test(x_clean, two.sided = TRUE), silent = TRUE)
    if(inherits(grubbs_test, "try-error")) break  # no more outliers
    p_val <- grubbs_test$p.value
    if(p_val < alpha) {
      # extract outlier value from test result
      outlier_val <- as.numeric(strsplit(grubbs_test$alternative, " ")[[1]][3])
      x_clean <- x_clean[x_clean != outlier_val]
      outlier_vals <- c(outlier_vals, outlier_val)
    } else {
      break
    }
  }
  
  list(clean = x_clean, outliers = outlier_vals)
}

# --- 5. Apply Grubbs per genotype ---
cleaned_data <- data.frame()
outlier_data <- data.frame()

for(gt in levels(data$Genotype)) {
  subset_gt <- data$Total_social_time[data$Genotype == gt]
  res <- remove_grubbs_outliers(subset_gt)
  
  # Append cleaned data
  cleaned_data <- rbind(cleaned_data,
                        data.frame(Genotype = gt, Total_social_time = res$clean))
  
  # Append outliers
  if(length(res$outliers) > 0) {
    outlier_data <- rbind(outlier_data,
                          data.frame(Genotype = gt, Total_social_time = res$outliers))
  }
  
  # Print summary for this genotype
  cat(sprintf("Genotype '%s': %d outliers removed: %s\n",
              gt,
              length(res$outliers),
              if(length(res$outliers) > 0) paste(round(res$outliers, 2), collapse = ", ") else "None"))
}


# --- 6. Perform two-tailed independent t-test on cleaned data ---
ttest_result <- t.test(Total_social_time ~ Genotype, data = cleaned_data, var.equal = FALSE)

# --- 7. Extract results ---
df <- ttest_result$parameter
t_value <- ttest_result$statistic
p_value <- ttest_result$p.value

means <- cleaned_data %>%
  group_by(Genotype) %>%
  summarise(
    mean = mean(Total_social_time),
    sd = sd(Total_social_time),
    n = n()
  )

# --- 8. Save text summary ---
output_file <- dirname(data_file)  # Save in same folder as input
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
  nrow(outlier_data),
  ifelse(nrow(outlier_data) > 0,
         paste(round(outlier_data$Total_social_time, 2), collapse = ", "),
         "None")
)

summary_file <- file.path(output_file, "TotalSocialTime_ttest_summary_grubbs.txt")
writeLines(output_text, summary_file)
cat("Statistical summary saved to:\n", summary_file, "\n")

# --- 9. Plot and save SVG ---
plot_file <- file.path(output_file, "TotalSocialTime_boxplot_grubbs_removed.svg")

p <- ggplot(cleaned_data, aes(x = Genotype, y = Total_social_time, fill = Genotype)) +
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
cat("SVG plot saved to:\n", plot_file, "\n")
