# N-FORTE: Mother's own milk vs donor milk genus associations (Figure 4B)
# Extracted from the final pseudocount=1 CLR analysis in N-FORTE R Markdown.
# Input: analysis-ready phyloseq object 'physeq' containing BM_group,
#        with labels 'Mother' (>90% MOM) and 'Donor' (<10% MOM).
# Intermediate milk-source samples must already have BM_group missing/other.
# No beta-diversity analysis, PERMANOVA, or automatic file exports.

library(phyloseq)
library(dplyr)
library(tidyr)
library(pheatmap)

bm_ps_in <- physeq
bm_time_var <- "Time"
bm_group_var <- "BM_group"
bm_batch_var <- "Run_name"
bm_top_n <- 20
bm_pseudo <- 1
bm_time_levels <- c("Day 1-7", "Day 8-14", "Day 15-21", "Day 22-28")

bm_to_chr <- function(x) {
  if (is.list(x)) {
    vapply(x, function(z) if (length(z) == 0) NA_character_ else as.character(z[[1]]),
           FUN.VALUE = character(1))
  } else {
    as.character(x)
  }
}

bm_get_meta <- function(ps) {
  m <- as.data.frame(sample_data(ps), stringsAsFactors = FALSE)
  m[[bm_time_var]] <- bm_to_chr(m[[bm_time_var]])
  m[[bm_group_var]] <- bm_to_chr(m[[bm_group_var]])
  m[[bm_batch_var]] <- bm_to_chr(m[[bm_batch_var]])
  m
}

bm_otu_mat_samples_by_taxa <- function(ps) {
  X <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) X <- t(X)
  X
}

bm_clr <- function(X, pseudo = 1) {
  logX <- log(X + pseudo)
  logX - rowMeans(logX)
}

bm_fit_lm_one_time <- function(ps_time) {
  meta <- bm_get_meta(ps_time)
  keep <- !is.na(meta[[bm_group_var]]) &
    meta[[bm_group_var]] %in% c("Mother", "Donor") &
    !is.na(meta[[bm_batch_var]])
  ps_time <- prune_samples(rownames(meta)[keep], ps_time)
  meta <- bm_get_meta(ps_time)
  if (nrow(meta) < 3 || length(unique(meta[[bm_group_var]])) < 2) return(NULL)

  meta[[bm_group_var]] <- factor(meta[[bm_group_var]], levels = c("Donor", "Mother"))
  meta[[bm_batch_var]] <- factor(meta[[bm_batch_var]])

  X <- bm_otu_mat_samples_by_taxa(ps_time)
  X <- X[rowSums(X) > 0, , drop = FALSE]
  meta <- meta[rownames(X), , drop = FALSE]
  X_clr <- bm_clr(X, pseudo = bm_pseudo)

  res <- lapply(colnames(X_clr), function(genus) {
    df <- data.frame(
      y = X_clr[, genus],
      BM_group = meta[[bm_group_var]],
      Run_name = meta[[bm_batch_var]]
    )
    fit <- lm(y ~ BM_group + Run_name, data = df)
    cf <- summary(fit)$coefficients
    term <- "BM_groupMother"  # Donor is the reference group
    data.frame(
      Genus = genus,
      beta = if (term %in% rownames(cf)) unname(cf[term, "Estimate"]) else NA_real_,
      p = if (term %in% rownames(cf)) unname(cf[term, "Pr(>|t|)"]) else NA_real_,
      n = nrow(df),
      n_mother = sum(df$BM_group == "Mother"),
      n_donor = sum(df$BM_group == "Donor"),
      stringsAsFactors = FALSE
    )
  })
  res <- bind_rows(res)
  res$padj <- p.adjust(res$p, method = "BH") # within time window
  res
}

# Aggregate once and retain unclassified genera, as in the original script.
bm_ps_genus <- tax_glom(bm_ps_in, taxrank = "Genus", NArm = FALSE)
tax_df <- as.data.frame(tax_table(bm_ps_genus), stringsAsFactors = FALSE)
genus_label <- tax_df$Genus
genus_label[is.na(genus_label) | genus_label == ""] <- "Genus_unclassified"
genus_label <- make.unique(genus_label)
tax_map <- data.frame(
  TaxonID = taxa_names(bm_ps_genus),
  GenusLabel = genus_label,
  stringsAsFactors = FALSE
)

bm_meta0 <- bm_get_meta(bm_ps_genus)
bm_lm_res <- lapply(bm_time_levels, function(tt) {
  ps_tt <- prune_samples(rownames(bm_meta0)[bm_meta0[[bm_time_var]] == tt], bm_ps_genus)
  if (nsamples(ps_tt) == 0) return(NULL)
  ps_tt <- prune_taxa(taxa_sums(ps_tt) > 0, ps_tt)
  out <- bm_fit_lm_one_time(ps_tt)
  if (is.null(out)) return(NULL)
  out$Time <- tt
  out
})

bm_lm_res_df <- bind_rows(bm_lm_res) %>%
  relocate(Time, Genus, beta, p, padj, n, n_mother, n_donor) %>%
  left_join(tax_map, by = c("Genus" = "TaxonID")) %>%
  mutate(Genus = ifelse(is.na(GenusLabel), Genus, GenusLabel)) %>%
  select(-GenusLabel)

# Original top-genus ranking: lowest q across windows; tie-break max |beta|.
bm_top_genera <- bm_lm_res_df %>%
  group_by(Genus) %>%
  summarise(
    min_padj = suppressWarnings(min(padj, na.rm = TRUE)),
    max_abs_beta = suppressWarnings(max(abs(beta), na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(min_padj, desc(max_abs_beta)) %>%
  slice_head(n = bm_top_n) %>%
  pull(Genus)

bm_lm_res_df <- bm_lm_res_df %>%
  mutate(Time = factor(Time, levels = bm_time_levels)) %>%
  arrange(Time)

bm_beta_mat <- bm_lm_res_df %>%
  filter(Genus %in% bm_top_genera) %>%
  select(Time, Genus, beta) %>%
  pivot_wider(names_from = Time, values_from = beta) %>%
  tibble::column_to_rownames("Genus") %>%
  as.matrix()
bm_beta_mat <- bm_beta_mat[, intersect(bm_time_levels, colnames(bm_beta_mat)), drop = FALSE]

bm_padj_mat <- bm_lm_res_df %>%
  filter(Genus %in% bm_top_genera) %>%
  select(Time, Genus, padj) %>%
  pivot_wider(names_from = Time, values_from = padj) %>%
  tibble::column_to_rownames("Genus") %>%
  as.matrix()
bm_padj_mat <- bm_padj_mat[rownames(bm_beta_mat), colnames(bm_beta_mat), drop = FALSE]

# Plot is returned to the active graphics device; nothing is saved.
bm_sig_mat <- ifelse(!is.na(bm_padj_mat) & bm_padj_mat < 0.05, "*", "")
bm_col_fun <- colorRampPalette(c("#D7191C", "white", "#2C7BB6"))
bm_heatmap <- pheatmap(
  bm_beta_mat,
  color = bm_col_fun(101),
  breaks = seq(-max(abs(bm_beta_mat), na.rm = TRUE),
                max(abs(bm_beta_mat), na.rm = TRUE), length.out = 102),
  cluster_rows = TRUE, cluster_cols = FALSE,
  display_numbers = bm_sig_mat,
  number_color = "black", fontsize_number = 12,
  silent = TRUE
)

# Inspect bm_lm_res_df (all model results), bm_beta_mat, bm_padj_mat.
# To display: grid::grid.newpage(); grid::grid.draw(bm_heatmap$gtable)
