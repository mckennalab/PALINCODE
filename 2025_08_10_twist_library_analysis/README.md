# 2025_08_10_twist_library_analysis

Comparative editing-rate analyses across the Twist v1 library, broken down
by **guide length (17 / 18 / 20 bp)**, **tracrRNA scaffold (old vs. new)**,
and **base editor variant (ABEmax / 8e / control)**. This directory
produced the first drafts of the per-target layout and balance plots that
appear in the manuscript figures.

## Inputs

- Processed data from `aidan/Data/all_length{17,18,20}.csv` (upstream of
  this directory).
- `18_feature_selection_2025_08_11.tsv` — feature-selection output used
  by the highlight plots.

## Outputs (PDFs)

Per-length target layouts:

- `target_17_layout.pdf`
- `target_18_layout.pdf`
- `target_20_layout.pdf`

Per-scaffold / per-editor 18 bp comparisons:

- `target_18_layout_8e_new_18_.pdf`, `target_18_layout_8e_old_18_.pdf`
- `target_18_layout_abemax_new_18_.pdf`,
  `target_18_layout_abemax_old_18_.pdf`
- `target_18_layout_control_new_18_.pdf`,
  `target_18_layout_control_old_18_.pdf`

Library balance and summary plots:

- `balance_old_twist.pdf`, `balance_plot_libraries.pdf`
- `target_tbl_1_018.pdf`, `target_tbl_1_018_editing_per_target.pdf`
- `highlight_editing_rates.pdf`

## Scripts

- `figure_3/editing_rates.R` — position-difference mapping and
  target-specific editing comparisons used by the figure-3 panels.
- `figure_3/sup_figure_3_feature_plot.R` — feature-importance plotting
  used in the supplementary figure.
