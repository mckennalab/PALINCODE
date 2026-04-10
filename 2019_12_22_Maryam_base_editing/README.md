# 2019_12_22_Maryam_base_editing

Detailed two-target study (Pal1 and Pal3) across multiple base editor
variants (ABE, CBE) with and without UMIs. This is the first directory
that introduces the **directional LEFT / RIGHT / BOTH / NONE** outcome
classification that the rest of the project uses.

## Inputs

- `data/all_samples.txt` — aggregated per-position base counts for 12
  samples covering both palindromic targets and four conditions
  (e.g. `MF_ABPal1`, `MF_BEPal1`, `MF_Pal1_Ct`, …). Generated upstream
  from the alignment pipeline.

## Notebooks and scripts

- `notebooks/target_1_3_round2_analysis.ipynb` — main analysis. Loads the
  per-sample stats, computes 5×20 (`ACGTN` × position) pileup matrices,
  identifies A→G (ABE) or C→T (CBE) conversions inside windowed regions,
  and assigns each read to a LEFT / RIGHT / BOTH / NONE class.
- `notebooks/create_per_position_plots.R` — generates the faceted bar
  plots of per-base editing rates, split by sample and palindrome.

## Outputs

In `data/`:

- `palindrome_table.txt` — summary of editing directionality, indexed by
  window size, position, and outcome class.

In `plots/`: 14+ PNGs covering per-base editing rates (Pal1, Pal3, and
combined), control-subtracted editing, and the LEFT / RIGHT / BOTH
distributions.

## Notes

This is where the LRBN ("LEFT / RIGHT / BOTH / NEITHER") nomenclature
that runs through the rest of the repo originates.
