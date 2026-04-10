# 2021_08_22_twist_v2 — Twist library v2 design and references

Refined library design and reference files for the **second Twist library
(v2)**. The main goals of v2 were:

- Compare **old vs. new tracrRNA scaffolds** head-to-head.
- Test **17 bp and 18 bp guide variants** in a controlled background.
- Provide reference FASTAs and YAML configs for the alignment pipeline.
- Include MARC1 references for downstream single-cell barcoding work.

The downstream sequencing analysis of this library lives in
`../2025_08_10_twist_library_analysis/`,
`../2025_08_11_twist2/`, and
`../2025_11_30_twist_analysis/`.

## Reference FASTAs

- `17guide1.fasta` — 17 bp guide reference (~115 KB).
- `18guide1.fa`, `18guide1.fa.fai` — 18 bp guide reference + index.
- `all_MARC1_references.fa` — MARC1 reference set used for single-cell
  barcoding experiments (~17 KB).
- `all_MARC1_references_upper.fa`, `*_clipped.fa` — variant
  case-normalized / clipped versions.

## Pipeline configs (YAML)

The `.yml` files describe how the alignment pipeline should treat each
guide variant. The relevant ones are:

- `17guide.yml`
- `18guide1.yml`, `18guide2.yml`, `18guide3.yml`
- `marc1_example_biopsy_1.yml` — example MARC1 biopsy config.

## Test data and visualization

- `17guide_test.bam` (~5.7 MB) and `18guide_test.bam` (~9.9 MB) — small
  aligned read sets used to sanity-check the references.
- `visualize.py` — CIGAR / mismatch parser for inspecting BAM alignments.
- `2021_08_03_editing_rates.R` — early editing-rate calculation by
  sample / outcome / editing type.
- `Untitled.ipynb` — exploratory notebook.

## Other contents

- `MF/` — additional MF-related references.
- `src/` — source / helper scripts.
- `images/` — supporting images for design notes.
- `twist_v1_oligos/` — copies of the v1 oligos for cross-comparison.
