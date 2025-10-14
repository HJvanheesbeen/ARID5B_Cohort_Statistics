# ============================================================
#  WT vs HET ARID5B mouse weights across development
#  Handles unequal n, independent and repeated measures
#  Includes raw values for validation
# ============================================================

library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(tidyr)

# --- 1. Load data ---
data_file <- "**Path**/NumericalDataWeight.csv"
data_raw <- read.csv(data_file, stringsAsFactors = FALSE, fill = TRUE)

# Clean column names
colnames(data_raw) <- gsub("\\.+", "_", colnames(data_raw))
colnames(data_raw) <- trimws(colnames(data_raw))

# Convert numeric columns
numeric_cols <- setdiff(colnames(data_raw), "Genotype")
data <- data_raw %>%
  mutate(across(all_of(numeric_cols), ~ as.numeric(ifelse(. == "", NA, .))))

# --- 2. Independent ages (P1, P6, P60) ---
independent_ages <- c("P1","P6","P60")
independent_results <- lapply(independent_ages, function(age) {
  df_age <- data %>% select(Genotype, all_of(age)) %>% filter(!is.na(.data[[age]]))
  
  x <- df_age %>% filter(Genotype=="WT") %>% pull(all_of(age))
  y <- df_age %>% filter(Genotype=="HET") %>% pull(all_of(age))
  
  if(length(x) > 0 & length(y) > 0) {
    ttest <- t.test(y, x, var.equal = TRUE, alternative = "less")
    data.frame(
      Age = age,
      Mean_WT = mean(x),
      Mean_HET = mean(y),
      Diff = mean(y) - mean(x),
      t_stat = ttest$statistic,
      df = ttest$parameter,
      p_value = ttest$p.value,
      Values_WT = paste(round(x,2), collapse=", "),
      Values_HET = paste(round(y,2), collapse=", "),
      Test = "Independent t-test"
    )
  }
}) %>% bind_rows()

# --- 3. Repeated ages (P8, P21) ---
df_repeated <- data %>%
  select(Genotype, P8, P21) %>%
  mutate(mouse_id = row_number()) %>%
  pivot_longer(cols=c(P8,P21), names_to="Age", values_to="Weight") %>%
  filter(!is.na(Weight)) %>%
  mutate(Age = factor(Age, levels=c("P8","P21")),
         Genotype = factor(Genotype))

m <- lmer(Weight ~ Genotype*Age + (1|mouse_id), data=df_repeated)
anova_res <- anova(m)

# emmeans for per-age comparisons
emm <- emmeans(m, ~ Genotype | Age)
contrast_res <- contrast(emm, method="pairwise", adjust="mvt")
contrast_df <- as.data.frame(summary(contrast_res))
contrast_df$p_one_tailed <- ifelse(contrast_df$estimate < 0, contrast_df$p.value/2, 1)

# Build summary table for repeated ages
emm_df <- as.data.frame(emm)

repeated_summary <- contrast_df %>%
  transmute(
    Age = Age,
    Mean_WT = emm_df %>% filter(Genotype=="WT") %>% pull(emmean),
    Mean_HET = emm_df %>% filter(Genotype=="HET") %>% pull(emmean),
    Diff = estimate,
    t_stat = t.ratio,
    df = df,
    p_value = p_one_tailed,
    # Output raw values for validation
    Values_WT = sapply(levels(df_repeated$Age), function(a) paste(round(df_repeated$Weight[df_repeated$Age==a & df_repeated$Genotype=="WT"],2), collapse=", "))[Age],
    Values_HET = sapply(levels(df_repeated$Age), function(a) paste(round(df_repeated$Weight[df_repeated$Age==a & df_repeated$Genotype=="HET"],2), collapse=", "))[Age],
    Test = "Repeated mixed model"
  )

# --- 4. Combine all results ---
full_summary <- bind_rows(independent_results, repeated_summary)

# --- 5. Save final table ---
output_file <- "**Path**/Figure_4c_full_summary_with_values.txt"
write.table(full_summary, file=output_file, sep="\t", row.names=FALSE, quote=FALSE)

cat("Combined per-age summary table with raw values ready at:\n", output_file, "\n")
