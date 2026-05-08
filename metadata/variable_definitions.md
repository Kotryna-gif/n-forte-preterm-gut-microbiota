# Variable definitions

This document describes the key sample-level variables referenced in the analysis
scripts. No sample-level metadata or values are provided here.

## Core variables

- **Intervention group**  
  Intervention group in the N-forte randomized controlled trial.  
  Infants were randomized to receive either a bovine-based fortifier or a
  human milk–derived fortifier (Prolacta).

- **Time**  
  Postnatal sampling time point.  
  Samples were collected at five weekly time points during the neonatal period:
  - day 1–7  
  - day 8–14  
  - day 15–21  
  - day 22–28  
  - gw 33+0–33+6


- **Delivery_mode**  
  Mode of delivery at birth.  
  Categorical variable (e.g. Vaginal delivery or Caesarean section).

- **Gestational_age_wks**  
  Gestational age at birth, expressed in completed weeks.

## Technical / study design variables

- **Run_name**  
  Sequencing run identifier used to account for technical batch effects.

- **Site_name**  
  Study site / neonatal intensive care unit where the infant was enrolled.

## Notes

- All variables are stored in the `sample_data` component of the phyloseq object.
- Full sample metadata are deposited via the corresponding SRA/BioSample records.
- This repository does not contain individual-level clinical or identifying data.
