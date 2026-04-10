# palindrome_analysis

**Original 2020 proof-of-concept analysis directory.** Establishes the
core LEFT / RIGHT / BOTH / NEITHER classification scheme on the MF_8_T18
sample and produces the first interactive D3 viewers for browsing
per-position editing patterns. Treat this as the historical foundation
for the rest of the repository — newer analyses live in the dated
`2023_*` and `2025_*` directories.

## Core notebooks

- `base_editing_statistics.ipynb` — main analysis. Walks the
  `MF_8_T18.stats` file, drops insertions, classifies each read as
  LEFT / RIGHT / BOTH / NEITHER based on A→G (forward) and T→C (reverse)
  conversions, and aggregates per-position base counts.
- `2020_11_25_base_editing.ipynb` — earlier per-sample analysis.
- `2020_12_17_analysis/` — extended per-base / per-guide analyses
  (Dec 2020 – Jan 2021).

## Reference and per-sample inputs

- `palt7.fa` — palindromic T7 reference, with companion files
  `palt7.fa.cutSites`, `palt7.fa.primers`, `palt7.fa.sites`.
- `MF_8_T18.stats` — primary aligner stats input.
- `MF_8_T18.perBase`, `MF_8_T18.topReadCounts`, `MF_8_T18.topReadEvents`,
  `MF_8_T18.topReadEventsNew` — per-base and per-read summaries derived
  from the stats file.
- `2020_11_25_data/` — 10 merged stats files from 4 samples (AB8 and ct
  guides, old/new constructs, rounds 17–20).
- `2020_11_25_editing_stats.txt` and `_old.txt` — aggregated editing
  statistics tables.
- `mouse_fish_human_fly_filtered_Hsu90_9-13-20.txt` — cross-species
  filtered guide table used as a comparison set.

## Library design helpers

- `array_order/` — guide design and array-ordering notebooks for the
  pre-Twist arrayed library.
- `old_scaffold_18/` — references for the older 18 bp scaffold.

## R analysis

- `2020_11_29_twist_editing_analysis.R` — first R-side visualization of
  the editing outcome distributions.

## Interactive viewers (HTML)

- `heatmap.html` — D3 v3 heatmap driven by `bases.txt` and
  `left_right.txt`.
- `read_editing_mutlihistogram.html` — multi-histogram visualization of
  per-position editing.
- `JS_files.js` — supporting JavaScript.

Open either HTML file in a browser. If your browser blocks local file
loads, serve the directory with `python -m http.server` and visit
`http://localhost:8000/`.

## Newer follow-ups

The interactive viewer concept is reworked (with newer data) in
[`../2023_07_24_base_editing_heatmap_counts/`](../2023_07_24_base_editing_heatmap_counts/).
The LRBN classification logic is reimplemented in
[`../2023_02_05_stats_to_base_editing/`](../2023_02_05_stats_to_base_editing/).
