# Alpha diversity (Shannon index)
#
# Computes Shannon diversity from a phyloseq object and visualises it
# by treatment arm across time points or at a selected time point.
#
# Inputs:
# - physeq: phyloseq object with sample_data
# - group_var: grouping variable (e.g. intervention ARM)
# - time_var: time variable
#
# Output:
# - ggplot object

library(phyloseq)
library(ggplot2)

plot_alpha_shannon <- function(physeq,
                               group_var = "ARM",
                               time_var  = "Time",
                               time_level = NULL,     # e.g. "Week1"; if NULL -> all time points
                               facet = TRUE,
                               palette = NULL,
                               base_size = 12) {
  stopifnot(inherits(physeq, "phyloseq"))

  # 1) Alpha diversity
  alpha_df <- estimate_richness(physeq, measures = "Shannon")

  # 2) Metadata
  meta_df <- as(sample_data(physeq), "data.frame")

  # 3) Add grouping variables
  if (!group_var %in% colnames(meta_df)) stop("group_var not found in sample_data(): ", group_var)
  if (!time_var  %in% colnames(meta_df)) stop("time_var not found in sample_data(): ", time_var)

  alpha_df[[group_var]] <- meta_df[[group_var]]
  alpha_df[[time_var]]  <- meta_df[[time_var]]

  # 4) Optional: select one time point
  if (!is.null(time_level)) {
    alpha_df <- alpha_df[alpha_df[[time_var]] == time_level, , drop = FALSE]
    if (nrow(alpha_df) == 0) stop("No samples found for ", time_var, " == ", time_level)
    facet <- FALSE
  }

  # 5) Plot
  p <- ggplot(alpha_df, aes_string(x = group_var, y = "Shannon", fill = group_var)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.65) +
    geom_jitter(width = 0.15, size = 1.3, alpha = 0.6) +
    theme_classic(base_size = base_size) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(size = 10),
      axis.text.x = element_text(angle = 30, hjust = 1)
    ) +
    labs(
      title = if (is.null(time_level)) "Alpha diversity across time points"
              else paste0("Alpha diversity at ", time_var, " = ", time_level),
      x = NULL,
      y = "Shannon index"
    )

  if (facet) {
    p <- p + facet_wrap(stats::as.formula(paste("~", time_var)), nrow = 1, scales = "free_y")
  }

  if (!is.null(palette)) {
    p <- p + scale_fill_manual(values = palette)
  }

  p
}

# Example: your original “all time points faceted”
arm_cols <- c(
  "Prolact+6 H^2MF arm" = "#0072B2",
  "Bovine protein arm" = "#D55E00"
)

p_all <- plot_alpha_shannon(physeq, group_var = "ARM", time_var = "Time",
                            time_level = NULL, facet = TRUE, palette = arm_cols)
p_all

# Example: single time point (replaces per-time scripts)
p_day1_to_7 <- plot_alpha_shannon(physeq, group_var = "ARM", time_var = "Time",
                            time_level = "day_1_to_7", palette = arm_cols)
p_day1_to_7
