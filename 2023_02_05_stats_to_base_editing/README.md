# 2023_02_05_stats_to_base_editing

**Pipeline directory** that converts raw aligner stats files into per-target
base editing event tables. This is the "stats → tidy tables" hub that
downstream analyses (`../2023_07_24_base_editing_heatmap_counts/`,
`../2025_08_*`, `../2025_11_30_twist_analysis/`) consume.

## What it does

For each stats file:

1. Filter to PASS records.
2. Extract the edited positions per read.
3. Categorize each read into one of LEFT / RIGHT / BOTH / NEITHER based on
   which side(s) of the palindrome were edited.
4. Aggregate per-target counts and write a tidy event table.

## Key scripts and notebooks

- `convert_stats_to_editing_table.py` — Python implementation. Use this if
  you don't have a Scala toolchain.
- `convert_full_editing_output_to_per_target_file.scala` — equivalent
  Scala implementation that aggregates UMI-corrected events across
  positions and LRBN outcomes.
- `2023_02_06_convert_to_event_table.ipynb` — full notebook walkthrough of
  the pipeline (stats load → filter → classify → write table).
- `target_1_3_round2_analysis.ipynb` — single-target analysis example
  using the output tables.
- `2023_07_19_quartet.R` and `2023_02_10_test.R` — statistical /
  simulation work used to validate the editing-rate quartet plots.

## Inputs

- `data/2021_04_13_MF_Lib_12/NpT7_g17_N8.stats.gz` — example input stats
  file (full data sets are not stored in this repo).
- `test_rates.txt.gz`, `test.output.gz` — small test inputs used during
  development.

## Outputs

- TSV / CSV event tables containing `(sequence, proportion,
  editing_sites, left_right_both)` columns per target.
- Quartet PDFs validating the statistical properties of the rate
  estimator (`example_quartet.png.pdf`,
  `quartet_plot_balance_*_editing_rate_*_integrations_*_replicates_a.pdf`,
  `simulation_by_rate_parameters*.pdf`).

## Reproducing

Run the Python script directly on a `.stats(.gz)` file:

```
python convert_stats_to_editing_table.py <input.stats[.gz]> <output.tsv>
```

The notebook walks through the same logic step-by-step if you want to
adapt the categorization rules.
