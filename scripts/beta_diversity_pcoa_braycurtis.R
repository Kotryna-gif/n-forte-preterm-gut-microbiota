# Beta diversity (Bray-Curtis) + PCoA
#
#
# Computes Bray-Curtis distances from a phyloseq object (relative abundance)
# and produces a PCoA plot coloured by a chosen variable, optionally faceted
# by time point.
#
# Inputs:
# - physeq: phyloseq object with sample_data
# - color_var: sample_data column mapped to colour (e.g. "ARM")
# - time_var: sample_data column used for faceting (e.g. "Time"); set NULL to disable
# - time_level: optionally select a single time point (e.g. "Week1"); NULL keeps all
# - palette: optional named vector of colours for discrete variables
#
# Output:
# - list with $plot (ggplot), $ordination, $distance

library(phyloseq)
library(ggplot2)

plot_pcoa_bray <- function(physeq,
                           color_var = "ARM",
                           time_var  = "Time",
                           time_level = NULL,
                           facet = TRUE,
                           transform = c("relative_abundance", "none"),
                           palette = NULL,
                           seed = 42,
                           point_size = 1.8,
                           point_alpha = 0.85,
                           base_size = 12) {
  stopifnot(inherits(physeq, "phyloseq"))
  transform <- match.arg(transform)

  meta_df <- as(sample_data(physeq), "data.frame")

  if (!color_var %in% colnames(meta_df)) stop("color_var not found in sample_data(): ", color_var)
  if (!is.null(time_var) && !time_var %in% colnames(meta_df)) stop("time_var not found in sample_data(): ", time_var)

  # Optional: subset to one time point
  if (!is.null(time_level)) {
    if (is.null(time_var)) stop("time_level provided but time_var is NULL.")
    keep <- meta_df[[time_var]] == time_level
    if (!any(keep, na.rm = TRUE)) stop("No samples found for ", time_var, " == ", time_level)
    physeq <- prune_samples(keep, physeq)
    facet <- FALSE
  }

  # Transform counts if requested
  ps_use <- physeq
  if (transform == "relative_abundance") {
    ps_use <- transform_sample_counts(ps_use, function(x) x / sum(x))
  }

  set.seed(seed)

  # Bray-Curtis distance + PCoA
  bray <- phyloseq::distance(ps_use, method = "bray")
  ord  <- ordinate(ps_use, method = "PCoA", distance = bray)

  # PCoA plot
  p <- plot_ordination(ps_use, ord, color = color_var) +
    geom_point(size = point_size, alpha = point_alpha) +
    coord_equal() +
    theme_classic(base_size = base_size) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8),
      strip.background = element_blank(),
      strip.text = element_text(size = 10)
    ) +
    labs(
      title = "PCoA (Bray-Curtis)",
      color = color_var
    )

  if (!is.null(time_var) && facet) {
    p <- p + facet_wrap(stats::as.formula(paste("~", time_var)), nrow = 1)
  }

  if (!is.null(palette)) {
    p <- p + scale_colour_manual(values = palette)
  }

  list(plot = p, ordination = ord, distance = bray)
}

# Example usage (intervention arm, faceted by time)
arm_cols <- c(
  "Prolact+6 H^2MF arm" = "#0072B2",
  "Bovine protein arm" = "#D55E00"
)

# res <- plot_pcoa_bray(physeq, color_var = "ARM", time_var = "Time", palette = arm_cols)
# res$plot
