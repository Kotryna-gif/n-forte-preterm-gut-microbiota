# n-forte-preterm-gut-microbiota
Analysis code for longitudinal gut microbiota profiling in extremely preterm infants (N-Forte trial)

# N-Forte preterm gut microbiota

This repository contains analysis code associated with a multicentre randomized controlled trial
investigating longitudinal gut microbiota development in extremely preterm infants (N-Forte trial).

Infants were randomized to receive either a bovine-protein based or a human milk–derived (Prolacta) fortifier
and followed longitudinally during the neonatal period. Gut microbiota samples were collected at
five weekly time points.

Raw sequencing data are deposited in NCBI SRA (accession pending).

## Repository contents
This repository provides scripts and documentation used for microbiota analyses reported in the
associated manuscript. No raw sequencing data or identifiable clinical data are included.

Planned structure:
- scripts/   : analysis scripts (R)
- figures/   : code used to generate manuscript figures
- metadata/  : variable descriptions and templates

## How to use the scripts (high level)
The scripts in `scripts/` are parameterised: the same analysis code is reused across different
explanatory variables (e.g. intervention arm, mother's own milk intake, macronutrient exposures)
by changing the input arguments (e.g. `group_var` / `color_var`).

## Notes
This repository is intended to support transparency and reproducibility of the published analyses.
