# figures

Final paper figures and the scripts that generate them. Each figure has
its own subdirectory containing R scripts, input data, and PDF / PNG /
Affinity-Designer outputs.

## Layout

```
figures/
├── figure1/        Early base editing efficiency by guide length
├── figure2/        Non-palindromic vs. palindromic target comparison (T7, RNF2)
├── figure3/        Editing outcome distributions across editor variants
├── figure4/        Clone-depth analysis from 293T lineage trees
├── figure5/        Single-cell editing outcomes and clonal trees in A375
├── sup_figure1/    Supplementary figure 1
├── sup_figure2/    Supplementary figure 2
├── sup_figure3/    Supplementary figure 3
├── sup_figure4/    Supplementary figure 4
└── sup_tables/     Supplementary tables
```

## Per-figure summary

| Figure | What it shows | Main inputs |
|---|---|---|
| **figure1** | Base editing efficiency vs. guide length (MF_Lib_06, MF_Lib_07, ABEmax) — bar plots and entropy analyses. | gzipped early editing tables. |
| **figure2** | Editing of non-palindromic (T7, RNF2) vs. palindromic (PalT7, PalRNF2) targets — heatmaps and point plots. | editing summary tables. |
| **figure3** | Per-target editing rates across three base editors (ABEmax, 8e, control). | `../2025_11_30_twist_analysis/Data/full_twist2_data.csv` |
| **figure4** | Distribution of cell clone depths from 293T lineage trees. | MIX phylogenetic `.contree` file. |
| **figure5** | Single-cell editing outcomes and per-clone phylogenies in A375 cells. | scRNA-seq + lineage data from `../2025_11_07_clone_analysis/` |

## Top-level scripts

- `2025_10_15_293T_clone_9_clone_analysis.R` — 293T clone-9 analysis
  driving the figure-4 panels.

## Reproducing

Each figure subdirectory's R scripts can be run on their own (set the
working directory to the subdir first). The most cross-cutting input is
[`../2025_11_30_twist_analysis/Data/full_twist2_data.csv`](../2025_11_30_twist_analysis/Data/),
which is referenced by several figure-generation scripts.
