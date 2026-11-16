# ============================================================
# BETA DIVERSITY LONGITUDINAL: DISTANCE TO OWN BASELINE
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(vegan)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

# -----------------------------
# USER SETTINGS
# -----------------------------
PS_OBJ   <- physeq
TIME_VAR <- "Time"
ARM_VAR  <- "ARM"
ID_VAR   <- "id"
BASELINE_LEVEL <- "Day 1-7"

time_levels <- c("Day 1-7","Day 8-14","Day 15-21","Day 22-28","GW 33+0-33+6")

arm_cols <- c(
  "HMF" = "#0072B2",
  "BMF" = "#D55E00"
)

# -----------------------------
# 1) Prepare phyloseq object
# -----------------------------
sample_data(PS_OBJ)[[TIME_VAR]] <- factor(
  as.character(sample_data(PS_OBJ)[[TIME_VAR]]),
  levels = time_levels,
  ordered = TRUE
)

ps0 <- prune_samples(
  !is.na(sample_data(PS_OBJ)[[TIME_VAR]]) &
    !is.na(sample_data(PS_OBJ)[[ARM_VAR]]) &
    !is.na(sample_data(PS_OBJ)[[ID_VAR]]),
  PS_OBJ
)

ps0 <- prune_samples(sample_sums(ps0) > 0, ps0)
ps_rel <- transform_sample_counts(ps0, function(x) x / sum(x))

meta_df <- as(sample_data(ps_rel), "data.frame") %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(
    Time = factor(as.character(.data[[TIME_VAR]]), levels = time_levels, ordered = TRUE),
    ARM  = factor(as.character(.data[[ARM_VAR]])),
    ID   = factor(as.character(.data[[ID_VAR]]))
  )

# -----------------------------
# 2) Build sample-by-taxa matrix
# -----------------------------
otu_mat <- as(otu_table(ps_rel), "matrix")

if (taxa_are_rows(ps_rel)) {
  otu_mat <- t(otu_mat)
}

otu_df <- as.data.frame(otu_mat) %>%
  tibble::rownames_to_column("Sample")

dat <- meta_df %>%
  left_join(otu_df, by = "Sample")

taxa_cols <- setdiff(colnames(dat), c("Sample", colnames(meta_df)))

# -----------------------------
# 3) Keep only IDs with a baseline sample
# -----------------------------
baseline_samples <- dat %>%
  filter(Time == BASELINE_LEVEL) %>%
  distinct(ID)

dat2 <- dat %>%
  semi_join(baseline_samples, by = "ID")

# -----------------------------
# 4) Distance to own baseline
# -----------------------------
get_dist_to_baseline <- function(df_one_id, taxa_cols, baseline_level) {
  baseline_row <- df_one_id %>% filter(Time == baseline_level)

  # If multiple baseline samples exist, keep first
  baseline_row <- baseline_row[1, , drop = FALSE]

  baseline_vec <- as.numeric(baseline_row[, taxa_cols])

  df_one_id %>%
    rowwise() %>%
    mutate(
      bray_to_baseline = vegan::vegdist(
        rbind(baseline_vec, as.numeric(c_across(all_of(taxa_cols)))),
        method = "bray"
      )[1]
    ) %>%
    ungroup()
}

beta_long <- dat2 %>%
  group_by(ID) %>%
  group_modify(~ get_dist_to_baseline(.x, taxa_cols, BASELINE_LEVEL)) %>%
  ungroup() %>%
  select(Sample, ID, ARM, Time, bray_to_baseline)

# -----------------------------
# 5) Summary for mean ± SE
# -----------------------------
beta_sum <- beta_long %>%
  group_by(Time, ARM) %>%
  summarise(
    n = n(),
    mean = mean(bray_to_baseline, na.rm = TRUE),
    sd   = sd(bray_to_baseline, na.rm = TRUE),
    se   = sd / sqrt(n),
    .groups = "drop"
  )

# -----------------------------
# 6) Mixed model
# -----------------------------
fit_beta <- lmer(bray_to_baseline ~ Time * ARM + (1 | ID), data = beta_long)

cat("\n================ BETA LMM ================\n")
print(anova(fit_beta))

# Optional post hoc: ARM within each time point
beta_emm <- emmeans(fit_beta, ~ ARM | Time)
beta_pairs <- contrast(beta_emm, method = "pairwise", adjust = "BH") %>%
  as.data.frame()

cat("\n========== BETA: ARM WITHIN TIME ==========\n")
print(beta_pairs)

# Extract global p-values for annotation
beta_anova <- as.data.frame(anova(fit_beta))
p_time_beta <- beta_anova["Time", "Pr(>F)"]
p_arm_beta  <- beta_anova["ARM", "Pr(>F)"]
p_int_beta  <- beta_anova["Time:ARM", "Pr(>F)"]

beta_label <- paste0(
  "LMM: Time p=", signif(p_time_beta, 3),
  " | ARM p=", signif(p_arm_beta, 3),
  " | Time×ARM p=", signif(p_int_beta, 3)
)

# -----------------------------
# 7) Plot
# -----------------------------
p_beta_longitudinal <-
  ggplot(beta_long, aes(x = Time, y = bray_to_baseline, group = ID, colour = ARM)) +
  geom_line(alpha = 0.18, linewidth = 0.4, show.legend = FALSE) +
  geom_point(
    shape = 21, fill = "white", color = "grey30",
    stroke = 0.35, size = 1.5, alpha = 1,
    position = position_jitter(width = 0.05, height = 0)
  ) +
  geom_line(
    data = beta_sum,
    aes(x = Time, y = mean, group = ARM, colour = ARM),
    linewidth = 1.1,
    show.legend = TRUE
  ) +
  geom_point(
    data = beta_sum,
    aes(x = Time, y = mean, fill = ARM),
    shape = 21, size = 2.7, color = "black", stroke = 0.35,
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    data = beta_sum,
    aes(x = Time, ymin = mean - se, ymax = mean + se, colour = ARM),
    width = 0.12, linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  scale_colour_manual(values = arm_cols, name = NULL) +
  scale_fill_manual(values = arm_cols, guide = "none") +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.subtitle = element_text(size = 10)
  ) +
  labs(
    x = NULL,
    y = "Bray–Curtis distance to baseline",
    subtitle = beta_label
  )

p_beta_longitudinal
