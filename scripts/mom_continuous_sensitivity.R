# N-FORTE: continuous mother's own milk proportion — genus sensitivity analysis
# Extracted from the original R Markdown sensitivity section (Table S7).
# NOTE: The original model used log10(relative abundance + 1e-5), NOT CLR.
# Input: analysis-ready phyloseq object 'physeq'. No automatic exports.

library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)

PS_OBJ <- physeq
TIME_LEVELS <- c("Day 1-7", "Day 8-14", "Day 15-21", "Day 22-28")
SUBJECT_VAR <- "id"
MOM_VAR <- "Proportion_of_MM_of_total_BM"
PSEUDOCOUNT <- 1e-5
COVARS <- character(0)  # unadjusted original model; optional additional variables

# Aggregate to genus, transform to relative abundance, select top 30 by mean RA.
physeq_genus <- tax_glom(PS_OBJ, taxrank = "Genus", NArm = FALSE)
physeq_genus_ra <- transform_sample_counts(physeq_genus, function(x) x / sum(x))

mom_genus_all <- psmelt(physeq_genus_ra) %>%
  mutate(
    Genus = as.character(Genus),
    Time = as.character(Time),
    MOM_percent = suppressWarnings(as.numeric(.data[[MOM_VAR]])),
    MOM_10pct = MOM_percent / 10,
    RA_percent = Abundance * 100,
    log10_RA = log10(Abundance + PSEUDOCOUNT)
  ) %>%
  filter(
    Time %in% TIME_LEVELS,
    !is.na(Genus), Genus != "",
    !is.na(MOM_percent), !is.na(.data[[SUBJECT_VAR]])
  )

top30_table <- mom_genus_all %>%
  group_by(Genus) %>%
  summarise(
    mean_RA = mean(Abundance, na.rm = TRUE),
    prevalence = mean(Abundance > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_RA)) %>%
  slice_head(n = 30)
top30_genera <- top30_table$Genus

mom_genus_df <- mom_genus_all %>%
  filter(Genus %in% top30_genera) %>%
  mutate(
    Time = factor(Time, levels = TIME_LEVELS),
    Genus = factor(Genus, levels = top30_genera)
  ) %>%
  droplevels()

# Same per-genus/per-window linear model and effect scale as original.
fit_mom_slope <- function(data, covars = character(0)) {
  required_columns <- c("log10_RA", "MOM_10pct", SUBJECT_VAR, covars)
  data_model <- data %>%
    filter(complete.cases(across(all_of(required_columns))))
  n_samples <- nrow(data_model)
  n_infants <- n_distinct(data_model[[SUBJECT_VAR]])
  n_mom_values <- n_distinct(data_model$MOM_10pct)
  prevalence <- mean(data_model$Abundance > 0, na.rm = TRUE)

  empty_result <- function() tibble(
    n_samples = n_samples, n_infants = n_infants,
    n_MOM_values = n_mom_values, prevalence = prevalence,
    beta_per_10pct = NA_real_, SE = NA_real_,
    lower_CI = NA_real_, upper_CI = NA_real_,
    t_value = NA_real_, p_nominal = NA_real_,
    fold_change_per_10pct = NA_real_
  )
  if (n_samples < 10 || n_mom_values < 2) return(empty_result())

  model <- lm(
    reformulate(c("MOM_10pct", covars), response = "log10_RA"),
    data = data_model
  )
  cf <- summary(model)$coefficients
  if (!"MOM_10pct" %in% rownames(cf)) return(empty_result())
  ci <- confint(model, parm = "MOM_10pct", level = 0.95)
  beta <- cf["MOM_10pct", "Estimate"]

  tibble(
    n_samples = n_samples, n_infants = n_infants,
    n_MOM_values = n_mom_values, prevalence = prevalence,
    beta_per_10pct = beta,
    SE = cf["MOM_10pct", "Std. Error"],
    lower_CI = ci[1], upper_CI = ci[2],
    t_value = cf["MOM_10pct", "t value"],
    p_nominal = cf["MOM_10pct", "Pr(>|t|)"],
    fold_change_per_10pct = 10^beta
  )
}

mom_results <- mom_genus_df %>%
  group_by(Time, Genus) %>%
  group_modify(~ fit_mom_slope(.x, covars = COVARS)) %>%
  ungroup() %>%
  group_by(Time) %>%
  mutate(q_BH_time = p.adjust(p_nominal, method = "BH")) %>%
  ungroup() %>%
  mutate(
    q_BH_all = p.adjust(p_nominal, method = "BH"),
    significant_nominal = p_nominal < 0.05,
    significant_BH_time = q_BH_time < 0.05,
    significant_BH_all = q_BH_all < 0.05
  ) %>%
  arrange(Time, q_BH_time)

# Inspect: top30_table and mom_results. No file is saved.
