# Longitudinal Bray-Curtis distance to each infant's baseline
# N-FORTE; related to Figure S1B
# Input: a phyloseq object with id, ARM, and Time in sample_data.
# Output: data, fitted model, ANOVA, BH-adjusted contrasts, ggplot.
# No input data are modified and no files are written by the script.

library(phyloseq)
library(ggplot2)
library(lmerTest)
library(emmeans)

distance_to_baseline_lmm <- function(
    physeq,
    id_var = "id", arm_var = "ARM", time_var = "Time",
    time_levels = c("Day 1-7", "Day 8-14", "Day 15-21",
                    "Day 22-28", "GW 33+0-33+6")) {

  stopifnot(inherits(physeq, "phyloseq"))
  meta <- as(sample_data(physeq), "data.frame")
  vars <- c(id_var, arm_var, time_var)
  if (!all(vars %in% names(meta))) {
    stop("Missing metadata variables: ", paste(setdiff(vars, names(meta)), collapse = ", "))
  }

  keep <- complete.cases(meta[, vars, drop = FALSE]) &
    as.character(meta[[time_var]]) %in% time_levels &
    sample_sums(physeq)[rownames(meta)] > 0
  ps <- prune_samples(rownames(meta)[keep], physeq)
  if (nsamples(ps) == 0) stop("No eligible samples.")

  ps <- transform_sample_counts(ps, function(x) x / sum(x))
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu <- t(otu)
  meta <- as(sample_data(ps), "data.frame")[rownames(otu), , drop = FALSE]

  dat <- data.frame(
    Sample = rownames(otu),
    ID = factor(as.character(meta[[id_var]])),
    ARM = factor(as.character(meta[[arm_var]])),
    Time = factor(as.character(meta[[time_var]]),
                  levels = time_levels, ordered = TRUE)
  )

  # Select the first Day 1-7 sample for each infant, as in the source analysis.
  baseline_rows <- which(dat$Time == time_levels[1])
  baseline_idx <- baseline_rows[
    match(as.character(dat$ID), as.character(dat$ID[baseline_rows]))
  ]
  keep <- !is.na(baseline_idx)
  if (!any(keep)) stop("No infants have a Day 1-7 baseline sample.")

  x <- otu[keep, , drop = FALSE]
  x0 <- otu[baseline_idx[keep], , drop = FALSE]
  dat <- dat[keep, , drop = FALSE]

  # Pairwise Bray-Curtis distance, mathematically equivalent to vegdist.
  dat$bray_to_baseline <- rowSums(abs(x - x0)) / rowSums(x + x0)

  model <- lmerTest::lmer(bray_to_baseline ~ Time * ARM + (1 | ID), data = dat)
  tests <- anova(model)

  emm <- emmeans::emmeans(model, ~ ARM | Time)
  pairs <- emmeans::contrast(emm, method = "pairwise")
  # Adjust across the five time-window contrasts, not within each window.
  contrasts <- as.data.frame(summary(pairs, by = NULL, adjust = "BH"))

  means <- aggregate(
    bray_to_baseline ~ Time + ARM, data = dat,
    FUN = function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x)))
  )
  means <- data.frame(
    Time = means$Time, ARM = means$ARM,
    mean = means$bray_to_baseline[, "mean"],
    se = means$bray_to_baseline[, "se"]
  )

  p <- ggplot(dat, aes(Time, bray_to_baseline, group = ID, colour = ARM)) +
    geom_line(alpha = 0.18, linewidth = 0.4, show.legend = FALSE) +
    geom_point(size = 1.2, alpha = 0.5, show.legend = FALSE) +
    geom_line(
      data = means,
      aes(x = Time, y = mean, group = ARM, colour = ARM),
      linewidth = 1.1, inherit.aes = FALSE
    ) +
    geom_errorbar(
      data = means,
      aes(x = Time, ymin = mean - se, ymax = mean + se,
          group = ARM, colour = ARM),
      width = 0.12, inherit.aes = FALSE
    ) +
    labs(x = NULL, y = "Bray-Curtis distance to baseline", colour = NULL) +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 30, hjust = 1))

  list(data = dat, model = model, anova = tests,
       contrasts = contrasts, plot = p)
}

# Example usage:
# result <- distance_to_baseline_lmm(physeq)
# result$anova
# result$contrasts
# result$plot
