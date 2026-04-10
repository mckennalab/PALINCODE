# 2025_08_11_twist2

Per-target / per-position base editing rate processing for the **three
PAM-manipulated Twist pools** (PAM_TWIST_1, PAM_TWIST_2, PAM_TWIST_3).
This is the conversion + feature-importance step that feeds the modeling
work in `../2025_11_30_twist_analysis/`.

## Inputs

Large stats files from the three pools (17–31 MB each):

- `PAM_TWIST_1_017.txt`, `PAM_TWIST_1_018.txt` — pool 1, 17 bp and 18 bp.
- `PAM_TWIST_2_018.txt`, `PAM_TWIST_2_N18.txt` — pool 2, 18 bp and N18
  variant.
- `PAM_TWIST_3_018.txt` — pool 3, 18 bp.

## Notebook

- `PAM manipulated twist pool 3-18-21.ipynb` — design rationale and
  oligo construction logic for the PAM-manipulated pools.

## Conversion script

- `convert_stats_to_editing_table.py` — same converter as in
  `../2023_02_05_stats_to_base_editing/`, vendored here so the
  directory is self-contained. Handles palindrome-specific
  post-processing (LEFT / RIGHT / BOTH / NEITHER classification) and
  outputs position-by-base editing matrices.

## Outputs

- `total_editing_sites.txt` (~9.9 MB) — aggregated per-position editing
  table across all targets and pools. This is the file consumed
  downstream by `../2025_11_30_twist_analysis/`.
- `feature_importance_2025_08_12.pdf`,
  `feature_importance_2025_08_12.png` — XGBoost feature-importance
  ranking from an early ML pass on the editing-rate prediction task.
- `Screenshot 2025-08-11 at 9.29.42 PM.png` — exploratory snapshot.
