# PERMANOVA (adonis2) on Bray-Curtis distances (by = "margin")
#
# - Prunes to complete cases for specified variables
# - Coerces factor/numeric variables safely
# - Computes Bray-Curtis on relative abundance
# - Runs adonis2 with partial R2 (by = "margin")
# - Provides a compact barplot of partial R2 with permutation p-values


# ---- Example usage (adjust variable names to your sample_data) ----
# fac_vars <- c("Run_name", "Site_name", "Delivery_mode", "Time", "ARM")
# num_vars <- c("Gestational_age_wks")


library(phyloseq)
library(vegan)
library(ggplot2)
library(scales)

make_rel_bray <- function(ps) {
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
  phyloseq::distance(ps_rel, method = "bray")
}

coerce_meta_types <- function(meta, fac_vars = character(), num_vars = character()) {
  meta <- data.frame(meta, check.names = FALSE)
  if (length(fac_vars)) meta[fac_vars] <- lapply(meta[fac_vars], \(x) factor(as.character(x)))
  if (length(num_vars)) meta[num_vars] <- lapply(meta[num_vars], \(x) suppressWarnings(as.numeric(as.character(x))))
  meta
}

run_permanova_margin <- function(physeq,
                                 fac_vars,
                                 num_vars = character(),
                                 permutations = 999,
                                 seed = 42) {
  stopifnot(inherits(physeq, "phyloseq"))

  vars <- c(fac_vars, num_vars)

  meta <- coerce_meta_types(sample_data(physeq), fac_vars = fac_vars, num_vars = num_vars)

  keep <- complete.cases(meta[, vars, drop = FALSE])
  if (!any(keep)) stop("No complete cases after filtering on: ", paste(vars, collapse = ", "))

  physeq_cc <- prune_samples(rownames(meta)[keep], physeq)
  meta_cc   <- coerce_meta_types(sample_data(physeq_cc), fac_vars = fac_vars, num_vars = num_vars)

  set.seed(seed)
  bray <- make_rel_bray(physeq_cc)

  # Align metadata row order to distance labels
  lab <- attr(bray, "Labels")
  meta_cc <- meta_cc[lab, , drop = FALSE]

  # Build formula: bray ~ var1 + var2 + ...
  form <- stats::as.formula(paste("bray ~", paste(vars, collapse = " + ")))

  permanova <- adonis2(
    form,
    data = meta_cc,
    permutations = permutations,
    by = "margin"
  )

  list(result = permanova, meta = meta_cc, distance = bray, formula = form)
}

permanova_to_table <- function(permanova) {
  tab <- as.data.frame(permanova)
  tab$term <- rownames(tab)
  tab <- subset(tab, !(term %in% c("Residual", "Total")))

  # Robustly pick R2 column
  r2_col <- if ("R2" %in% names(tab)) "R2" else names(tab)[grep("^R2", names(tab))[1]]
  names(tab)[names(tab) == r2_col] <- "R2"

  tab
}

plot_permanova_r2 <- function(permanova,
                              highlight_term = NULL,
                              title = "PERMANOVA effect sizes",
                              subtitle = "Bars = partial R2 (by = 'margin'); labels = permutation p-values") {
  tab <- permanova_to_table(permanova)

  tab$plab <- paste0("p=", signif(tab[["Pr(>F)"]], 2))
  tab$term <- factor(tab$term, levels = tab$term[order(tab$R2)])

  if (!is.null(highlight_term)) {
    tab$highlight <- ifelse(tab$term == highlight_term, "highlight", "other")
  } else {
    tab$highlight <- "other"
  }

  ggplot(tab, aes(x = term, y = R2)) +
    geom_col(aes(fill = highlight), width = 0.75, colour = "black", linewidth = 0.2) +
    geom_text(aes(label = plab), hjust = -0.05, size = 3) +
    coord_flip(ylim = c(0, max(tab$R2, na.rm = TRUE) * 1.15)) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_fill_manual(values = c("highlight" = "grey70", "other" = "white")) +
    labs(x = NULL, y = "Partial R2", title = title, subtitle = subtitle) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none")
}


# out <- run_permanova_margin(physeq, fac_vars = fac_vars, num_vars = num_vars, permutations = 999, seed = 42)
# print(out$result)
# p <- plot_permanova_r2(out$result, highlight_term = "ARM", title = NULL, subtitle = NULL)
# p

