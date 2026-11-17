# N-FORTE: HMO–genus associations (Figure 3B; corresponding tables)
# Extracted from the HMO heatmap analysis in N-FORTE R Markdown.
# Input: analysis-ready phyloseq object named physeq.
# No PERMANOVA, output files, or automatic exports.

library(phyloseq)
library(dplyr)
library(tidyr)
library(purrr)
library(compositions)
library(broom)
library(ggplot2)

TIME_LEVEL <- "Day 8-14"
BATCH_VAR <- "Run_name"

# Final analysis list: excludes the duplicate variable lnfi_img_kg_day.
HMO_VARS <- c(
  "x3_f_lmg_kg_day", "lnd_img_kg_day", "lnfii_img_kg_day",
  "ldf_tmg_kg_day", "x2_f_lmg_kg_day", "lnf_img_kg_day",
  "l_nn_tmg_kg_day", "ln_tmg_kg_day", "ls_tcmg_kg_day",
  "x6_s_lmg_kg_day", "x3_s_lmg_kg_day", "ls_tamg_kg_day",
  "ls_tbmg_kg_day", "dsln_tmg_kg_day"
)

hmo_pretty <- c(
  x3_f_lmg_kg_day = "3-FL", x2_f_lmg_kg_day = "2-FL",
  x3_s_lmg_kg_day = "3-SL", x6_s_lmg_kg_day = "6-SL",
  lnf_img_kg_day = "LNFP I", lnfii_img_kg_day = "LNFP III",
  lnd_img_kg_day = "LNDFH I", ln_tmg_kg_day = "LNT",
  l_nn_tmg_kg_day = "LNnT", ldf_tmg_kg_day = "LDFT",
  ls_tcmg_kg_day = "LSTc", ls_tamg_kg_day = "LSTa",
  ls_tbmg_kg_day = "LSTb", dsln_tmg_kg_day = "DSLNT"
)

TOP_N <- 25
SIG_Q <- 0.10

# 1. Subset to day 8–14, aggregate counts at genus rank.
ps_tt <- subset_samples(physeq, Time == TIME_LEVEL)
ps_tt <- prune_samples(sample_sums(ps_tt) > 0, ps_tt)
ps_genus <- tax_glom(ps_tt, taxrank = "Genus")

# Match the original taxonomy fallbacks for the heatmap.
tt <- as.data.frame(tax_table(ps_genus), stringsAsFactors = FALSE)
tax_label <- function(x) {
  genus <- x[["Genus"]]
  family <- x[["Family"]]
  order <- x[["Order"]]
  genus <- ifelse(is.na(genus) | genus == "" | genus == "NA", NA, genus)
  family <- ifelse(is.na(family) | family == "" | family == "NA", NA, family)
  order <- ifelse(is.na(order) | order == "" | order == "NA", NA, order)
  dplyr::case_when(
    !is.na(genus) & !is.na(family) ~ paste0(genus, " (", family, ")"),
    !is.na(genus) ~ genus,
    !is.na(family) ~ paste0("Unclassified (", family, ")"),
    !is.na(order) ~ paste0("Unclassified (", order, ")"),
    TRUE ~ "Unclassified"
  )
}
tax_map <- data.frame(
  taxa_id = rownames(tt),
  taxa_label = make.unique(apply(tt, 1, tax_label)),
  stringsAsFactors = FALSE
)

# 2. CLR on genus counts, pseudocount = 1.
otu <- as(otu_table(ps_genus), "matrix")
if (taxa_are_rows(ps_genus)) otu <- t(otu)
otu_clr <- as.data.frame(compositions::clr(otu + 1))
otu_clr$SampleID <- rownames(otu_clr)

meta <- data.frame(sample_data(ps_genus), check.names = FALSE)
meta$SampleID <- rownames(meta)
stopifnot(all(c(HMO_VARS, BATCH_VAR) %in% names(meta)))
meta_sub <- meta %>%
  select(SampleID, all_of(HMO_VARS), all_of(BATCH_VAR)) %>%
  mutate(
    across(all_of(HMO_VARS), ~ suppressWarnings(as.numeric(as.character(.)))),
    batch = factor(.data[[BATCH_VAR]])
  ) %>%
  select(-all_of(BATCH_VAR))

dat <- left_join(meta_sub, otu_clr, by = "SampleID")
genus_cols <- setdiff(names(otu_clr), "SampleID")

# 3. Original genus-by-HMO linear models, adjusted for sequencing run.
fit_one <- function(g, h, dat) {
  df <- data.frame(batch = dat$batch, hmo = dat[[h]], clr = dat[[g]])
  df <- df[complete.cases(df), , drop = FALSE]
  if (nrow(df) < 20) return(NULL)
  if (is.na(sd(df$hmo)) || sd(df$hmo) == 0) return(NULL)
  fit <- lm(clr ~ hmo + batch, data = df)
  broom::tidy(fit) %>%
    filter(term == "hmo") %>%
    transmute(genus = g, hmo = h, estimate, p.value)
}

results <- expand_grid(genus = genus_cols, hmo = HMO_VARS) %>%
  mutate(res = map2(genus, hmo, ~ fit_one(.x, .y, dat))) %>%
  select(-genus, -hmo) %>%
  unnest(res) %>%
  left_join(tax_map, by = c("genus" = "taxa_id")) %>%
  mutate(genus_label = coalesce(taxa_label, genus))

# Original heatmap: BH within each HMO across genera.
results2 <- results %>%
  group_by(hmo) %>%
  mutate(q.value = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    q.matrix = p.adjust(p.value, method = "BH"), # also report across all tests
    stars = case_when(
      q.value < 0.001 ~ "***",
      q.value < 0.01 ~ "**",
      q.value < 0.05 ~ "*",
      q.value < SIG_Q ~ "·",
      TRUE ~ ""
    )
  )

# 4. Same top-25 selection and heatmap structure as original R Markdown.
top_taxa <- results2 %>%
  group_by(genus_label) %>%
  summarise(max_abs_beta = max(abs(estimate), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_abs_beta)) %>%
  slice_head(n = TOP_N) %>%
  pull(genus_label)

plot_df <- results2 %>%
  filter(genus_label %in% top_taxa) %>%
  mutate(
    hmo = factor(hmo, levels = names(hmo_pretty)),
    hmo_pretty = factor(hmo_pretty[as.character(hmo)], levels = unname(hmo_pretty))
  )
row_order <- plot_df %>%
  group_by(genus_label) %>%
  summarise(max_abs = max(abs(estimate), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_abs)) %>%
  pull(genus_label)
plot_df$genus_label <- factor(plot_df$genus_label, levels = rev(row_order))

p_hmo <- ggplot(plot_df, aes(x = hmo_pretty, y = genus_label, fill = estimate)) +
  geom_tile() +
  geom_text(aes(label = stars), size = 5, na.rm = TRUE) +
  scale_fill_gradient2(
    low = "#5E3C99", mid = "white", high = "#2C7FB8",
    midpoint = 0, name = "β", limits = c(-0.015, 0.015),
    breaks = c(-0.015, 0, 0.015), oob = scales::squish
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 10)
  ) +
  labs(x = NULL, y = NULL)

# Inspect: results2 (model coefficients and q values), p_hmo (heatmap).
# LSTa outlier sensitivity is documented separately in the original R Markdown.
