# 2025_11_07_clone_analysis

Single-cell **clonal analysis** linking lineage recordings to scRNA-seq
observations. The goal is to connect cell-level transcriptomic identity
(scVI / leiden clustering) to per-cell editing genotypes at two
palindromic targets (PalT7 and PalRNF2 variants).

## Inputs

- `2025_11_07_obs.tsv` — single-cell observation matrix: cell barcodes,
  scVI clustering, leiden clusters, QC metrics (`n_genes`, coverage,
  mitochondrial %).
- `2025_11_07_lineage_recordings.tsv` — per-read lineage data: cell
  barcode, target coordinates, alignment metrics, and the LEFT / RIGHT /
  BOTH / NEITHER outcome at each of two palindromic targets.

## Notebooks

- `Untitled.ipynb` — placeholder notebook (currently empty). Analysis is
  in active development.

## Notes

This directory provides the **input tables** used by the single-cell
panels in `../figures/figure5/`. If you are looking for the actual figure
generation, see the figure5 R scripts there.
