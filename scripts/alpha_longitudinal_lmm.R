# ============================================================
# ALPHA DIVERSITY LONGITUDINAL
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
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

time_levels <- c("Day 1-7","Day 8-14","Day 15-21","Day 22-28","GW 33+0-33+6")

arm_cols <- c(
  "HMF" = "#0072B2",
  "BMF" = "#D55E00"
)

# -----------------------------
# 1) Alpha diversity + metadata
# -----------------------------
alpha_df <- estimate_richness(PS_OBJ, measures = "Shannon")
meta_df  <- as(sample_data(PS_OBJ), "data.frame")

alpha_df <- alpha_df %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(
    !!TIME_VAR := meta_df[[TIME_VAR]][match(Sample, rownames(meta_df))],
    !!ARM_VAR  := meta_df[[ARM_VAR]][match(Sample, rownames(meta_df))],
    !!ID_VAR   := meta_df[[ID_VAR]][match(Sample, rownames(meta_df))]
  ) %>%
  filter(
    !is.na(.data[[TIME_VAR]]),
    !is.na(.data[[ARM_VAR]]),
    !is.na(.data[[ID_VAR]]),
    !is.na(Shannon)
  ) %>%
  mutate(
    Time = factor(as.character(.data[[TIME_VAR]]), levels = time_levels, ordered = TRUE),
    ARM  = factor(as.character(.data[[ARM_VAR]])),
    ID   = factor(as.character(.data[[ID_VAR]]))
  )

# -----------------------------
# 2) Summary for mean ± SE
# -----------------------------
alpha_sum <- alpha_df %>%
  group_by(Time, ARM) %>%
  summarise(
    n = n(),
    mean = mean(Shannon, na.rm = TRUE),
    sd   = sd(Shannon, na.rm = TRUE),
    se   = sd / sqrt(n),
    .groups = "drop"
  )

# -----------------------------
# 3) Mixed model
# -----------------------------
fit_alpha <- lmer(Shannon ~ Time * ARM + (1 | ID), data = alpha_df)

cat("\n================ ALPHA LMM ================\n")
print(anova(fit_alpha))

# Optional post hoc: ARM within each time point
alpha_emm <- emmeans(fit_alpha, ~ ARM | Time)
alpha_pairs <- contrast(alpha_emm, method = "pairwise", adjust = "BH") %>%
  as.data.frame()

cat("\n========== ALPHA: ARM WITHIN TIME ==========\n")
print(alpha_pairs)

# Extract global p-values for annotation
alpha_anova <- as.data.frame(anova(fit_alpha))
p_time_alpha <- alpha_anova["Time", "Pr(>F)"]
p_arm_alpha  <- alpha_anova["ARM", "Pr(>F)"]
p_int_alpha  <- alpha_anova["Time:ARM", "Pr(>F)"]

alpha_label <- paste0(
  "LMM: Time p=", signif(p_time_alpha, 3),
  " | ARM p=", signif(p_arm_alpha, 3),
  " | Time×ARM p=", signif(p_int_alpha, 3)
)

# -----------------------------
# 4) Plot
# -----------------------------
p_alpha_longitudinal <-
  ggplot(alpha_df, aes(x = Time, y = Shannon, group = ID, colour = ARM)) +
  geom_line(alpha = 0.18, linewidth = 0.4, show.legend = FALSE) +
  geom_point(
    shape = 21, fill = "white", color = "grey30",
    stroke = 0.35, size = 1.5, alpha = 1,
    position = position_jitter(width = 0.05, height = 0)
  ) +
  geom_line(
    data = alpha_sum,
    aes(x = Time, y = mean, group = ARM, colour = ARM),
    linewidth = 1.1,
    show.legend = TRUE
  ) +
  geom_point(
    data = alpha_sum,
    aes(x = Time, y = mean, fill = ARM),
    shape = 21, size = 2.7, color = "black", stroke = 0.35,
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    data = alpha_sum,
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
    y = "Shannon index",
    subtitle = alpha_label
  )

p_alpha_longitudinal
