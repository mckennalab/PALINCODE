# 2020_02_2016

Larger pooled palindrome library screen comparing wild-type (control), ABE,
and CBE base editor outcomes on the full 314-target palindrome library
(rather than the two-target deep-dive in
`../2019_12_22_Maryam_base_editing/`).

## Inputs

- `data/stats.tar.gz` — compressed bundle of three stats files from the
  custom aligner:
  - `pigpallib_S7_R1.stats` — WT / control
  - `ABpigpallib_S8_R1.stats` — ABE
  - `BEpigpallib_S9_R1.stats` — CBE
- The palindrome reference FASTA from
  [`../2019_11_27_Maryam_palendromes/`](../2019_11_27_Maryam_palendromes/)
  (314 known palindromic targets).

## Notebook

- `notebooks/palindrome_library_analysis_v2.ipynb` — loads the GESTALT
  stats files, walks each merged read, and tallies editing outcomes
  (`NONE` / `LEFT` / `RIGHT` / `BOTH` / `BAD`) using window-based
  conversion rules:
  - **ABE:** A→G at protospacer positions 2–8 and 12–18
  - **CBE:** C→T at the same positions
  Outputs contingency tables and per-guide tallies.

## Outputs

- `data/tallied_outcomes.tsv` — per-guide summary with an `is_known` flag
  and counts in each outcome class (`none`, `left`, `right`, `both`,
  `bad`).
- `plots/` — 4 PNGs showing library occurrence and base editing
  distributions per editor type.
- `slides/` — slide deck assets used in lab presentations.

## Headline numbers (from notebook)

- ~5.4k WT UMIs, ~7.4k ABE, ~6.2k CBE.
- ~30–40% of reads from known guides match the expected on-target.
