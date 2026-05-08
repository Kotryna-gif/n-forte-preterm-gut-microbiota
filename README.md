# n-forte-preterm-gut-microbiota
Analysis code for longitudinal gut microbiota profiling in extremely preterm infants (N-forte trial)


This repository contains analysis code associated with a multicentre randomized controlled trial
investigating longitudinal gut microbiota development in extremely preterm infants (N-forte trial).

Infants were randomized to receive either a bovine-protein based or a human milk–derived (Prolacta) fortifier
and followed longitudinally during the neonatal period. Gut microbiota samples were collected at
five weekly time points.

Raw sequencing data are deposited in NCBI SRA BioProject:PRJNA1415696.

## Repository contents
This repository provides scripts and documentation used for microbiota analyses reported in the
associated manuscript. No raw sequencing data or identifiable clinical data are included.

Structure:
- scripts/   : analysis scripts (R)
- figures/   : code used to generate manuscript figures
- metadata/  : variable descriptions and templates

## How to use the scripts (high level)
The scripts in `scripts/` are parameterised: the same analysis code is reused across different
explanatory variables (e.g. intervention group, mother's own milk intake, macronutrient exposures)
by changing the input arguments (e.g. `group_var` / `color_var`).


## Alpha diversity
- `alpha_diversity_shannon.R`
  - Main function: `plot_alpha_shannon(physeq, group_var, time_var, time_level = NULL, ...)`

## Beta diversity (PCoA)
- `beta_diversity_pcoa_braycurtis.R`
  - Main function: `plot_pcoa_bray(physeq, color_var, time_var, time_level = NULL, ...)`

## PERMANOVA
- `permanova_bray_margin.R`
  - Main function: `run_permanova_margin(physeq, fac_vars, num_vars, ...)`
  - Plot helper: `plot_permanova_r2(permanova_result, highlight_term = NULL, ...)`

## Notes
This repository is intended to support transparency and reproducibility of the published analyses.
