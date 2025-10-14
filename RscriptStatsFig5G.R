# ============================================================
    ###Experimental design###
#  To asses increases in ARID5B expression area per cell of a set of designed and identified ARID5B variants == automated assessment of nuclear localization defect ==
#  Data represent average expression areas from ImageJ-quantified ROI boundaries (ARID5B expression per cell);
#  Each dot corresponds to the mean from non overlapping single cells on coverslips obtained from 3 different experimental plates.
# ============================================================
   ###Statistical approach###
#  Processes quantitative imaging data to compute n, mean, and SEM, and performs one-way ANOVA
#  (Variants that impair nuclear localization are expected to increase expression area per cell),
#  Along with pairwise t-tests vs control (ARID5B Isoform I) and Bonferroni correction for multiple comparisons.
# ============================================================


# Load necessary library
library(dplyr)
library(reshape2)

# -----------------------------
# 1. Load raw data
# -----------------------------
data_file <- "**Path**/Figure_5g_rawnumericaldata.txt"
data_raw <- read.table(data_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)

# Convert all columns to numeric (non-numeric become NA)
data <- data_raw %>% mutate(across(everything(), ~as.numeric(.)))

# -----------------------------
# 2. Summary statistics table (n, mean, SEM)
# -----------------------------
summary_table <- data.frame(
  Group = colnames(data),
  n = sapply(data, function(col) sum(!is.na(col))),
  mean = sapply(data, function(col) mean(col, na.rm=TRUE)),
  SEM  = sapply(data, function(col) sd(col, na.rm=TRUE)/sqrt(sum(!is.na(col))))
)

# Save summary table
write.table(summary_table,
            file = "**Path**/Figure_5g_summary_table.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# -----------------------------
# 3. Pairwise t-tests vs control (one-tailed, greater)
# -----------------------------
control_name <- colnames(data)[1]
other_groups <- setdiff(colnames(data), control_name)

pairwise_results <- lapply(other_groups, function(g) {
  x <- data[[control_name]][!is.na(data[[control_name]])]
  y <- data[[g]][!is.na(data[[g]])]
  
  ttest <- t.test(y, x, var.equal = TRUE, alternative = "greater")
  
  data.frame(
    Comparison = paste(control_name, "-", g),
    diff = mean(y) - mean(x),
    t_stat = ttest$statistic,
    df = ttest$parameter,
    p.value = ttest$p.value,
    CI_lower = ttest$conf.int[1],
    CI_upper = ttest$conf.int[2]
  )
})

pairwise_results_df <- do.call(rbind, pairwise_results)

# Bonferroni correction
pairwise_results_df$p.adj <- p.adjust(pairwise_results_df$p.value, method = "bonferroni")

# Significance stars
pairwise_results_df$Significance <- cut(
  pairwise_results_df$p.adj,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "ns")
)

# Save t-test table
write.table(pairwise_results_df,
            file = "**Path**/Figure_5g_ttest_vs_control.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# -----------------------------
# 4. ANOVA (one-way)
# -----------------------------
anova_long <- reshape2::melt(data, variable.name="Group", value.name="Value")
anova_result <- aov(Value ~ Group, data = anova_long)
anova_summary <- summary(anova_result)[[1]]

F_value <- anova_summary["Group","F value"]
df1 <- anova_summary["Group","Df"]
df2 <- anova_summary["Residuals","Df"]
p_value <- anova_summary["Group","Pr(>F)"]

anova_text <- paste0(
  "A one-way ANOVA was performed across the groups. ",
  "This revealed a significant effect of group (F(", df1, ", ", df2, ") = ", 
  round(F_value,2), ", P = ", signif(p_value,3), ", one-tailed, alpha = 0.05)."
)

# Save ANOVA summary as a text file
writeLines(anova_text,
           con = "**Path**/Figure_5g_anova_summary.txt")

cat("Tables and ANOVA summary saved successfully.\n")
