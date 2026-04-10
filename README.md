# palindrome_base_editing

Analyses, library designs, and figures for a study of CRISPR base editing on
**palindromic target sequences**. The repository spans the full project
lifecycle, from initial pilot experiments (2019) through synthesized Twist
oligo libraries, machine-learning models of edit outcomes, and the final
publication figures.

> **Status:** research code accompanying an in-preparation manuscript
> (McKenna lab). Citation will be added once the preprint is posted.

## What this project is about

Palindromic protospacers can in principle be cut/edited from either strand,
and we observed that base editing of these sites produces a characteristic
mixture of **LEFT / RIGHT / BOTH / NEITHER** outcomes depending on the guide,
the editor variant, and the local sequence context. The work in this
repository asks:

1. How do **target length and base position** within the protospacer change
   editing rates and the directionality of editing?
2. How does **flanking sequence context** modulate the same outcomes?
3. Can we **predict per-target editing efficiency and direction** from
   sequence features alone?
4. How do these patterns play out in single-cell **lineage recording**
   contexts (clonal experiments)?

The project uses several large pooled oligo libraries synthesized by Twist
Bioscience, paired-end deep sequencing of the edited libraries, and a custom
alignment / event-classification pipeline that produces per-target "stats"
files which are the input to most of the analyses here.

## Repository layout

The repo is organized as a flat collection of dated experiment directories.
Each directory is largely **self-contained** — it has its own data,
notebooks, and outputs, and depends on other directories only through shared
reference FASTAs and intermediate stats files.

| Directory | What's in it |
|---|---|
| [`palindrome_analysis/`](palindrome_analysis/) | Original 2020 proof-of-concept analysis: extracts LEFT/RIGHT/BOTH/NEITHER outcomes from MF_8_T18 stats and produces the first interactive D3 heatmap viewers. Foundational reference for the rest of the repo. |
| [`2019_11_27_Maryam_palendromes/`](2019_11_27_Maryam_palendromes/) | First-pass processing of the 314-target synthetic palindrome library; produces per-position base counts. |
| [`2019_12_12_control_data/`](2019_12_12_control_data/) | Control vs. edited comparison with UMI-aware deduplication. |
| [`2019_12_22_Maryam_base_editing/`](2019_12_22_Maryam_base_editing/) | Detailed two-target (Pal1, Pal3) study across ABE/CBE editor variants with directional LEFT/RIGHT/BOTH/NONE calls. |
| [`2020_02_2016/`](2020_02_2016/) | Larger pooled palindrome library screen across WT / ABE / CBE samples. |
| [`2021_03_19_twist_library_generation/`](2021_03_19_twist_library_generation/) | **Design** of the Twist v1 oligo library; produces ordered FASTAs. |
| [`2021_08_22_twist_v2/`](2021_08_22_twist_v2/) | **Design + reference setup** for Twist v2 (refined tracrRNA scaffolds, 17/18bp guides, MARC1 references). |
| [`2022_07_18_base_editing_t7_rnf2/`](2022_07_18_base_editing_t7_rnf2/) | Per-position editing heatmaps for the T7 and RNF2 palindromic targets across four pig-UMI samples. |
| [`2023_02_05_stats_to_base_editing/`](2023_02_05_stats_to_base_editing/) | **Pipeline:** converts aligner stats files into per-target editing event tables (Python + Scala + R). The "from stats to tables" hub. |
| [`2023_07_24_base_editing_heatmap_counts/`](2023_07_24_base_editing_heatmap_counts/) | Browser-based interactive D3.js heatmap viewer (`index.html`) for exploring editing patterns by sample and position. |
| [`2025_08_10_twist_library_analysis/`](2025_08_10_twist_library_analysis/) | Comparative editing-rate analyses across 17/18/20bp guides and tracrRNA variants; first-pass paper figure drafts. |
| [`2025_08_11_twist2/`](2025_08_11_twist2/) | Stats → editing table processing for the three PAM-manipulated Twist pools; generates feature-importance plots used in modeling. |
| [`2025_11_07_clone_analysis/`](2025_11_07_clone_analysis/) | Single-cell clone analysis: lineage recordings + scVI-clustered observations linking cell identities to editing genotypes. |
| [`2025_11_30_twist_analysis/`](2025_11_30_twist_analysis/) | **Most recent and most thoroughly documented analysis dir.** Reproducibility + ML modeling (XGBoost / SMOTE) of edit outcomes from sequence features. See its own README. |
| [`figures/`](figures/) | Final paper figures: `figure1`–`figure5` and `sup_figure1`–`sup_figure4`, each with R scripts, input data, and PDF/PNG/AI outputs. |
| [`data_deposits/`](data_deposits/) | GEO submission staging area (sequencing-run metadata template). |

Each subdirectory contains its own `README.md` describing inputs, the
notebooks/scripts that drive the analysis, and what outputs to expect.

## Conceptual data flow

```
   pooled palindrome library  ──►  paired-end sequencing
                                          │
                                          ▼
                            custom aligner ──► .stats files
                                          │
                                          ▼
   2023_02_05_stats_to_base_editing/  (stats → per-target event tables)
                                          │
                  ┌───────────────────────┼───────────────────────┐
                  ▼                       ▼                       ▼
       per-position editing      LEFT/RIGHT/BOTH/NEITHER     interactive viewer
       heatmaps                  outcome tables              (2023_07_24_…)
                  │                       │
                  └───────────┬───────────┘
                              ▼
                  2025_11_30_twist_analysis/
                  (correlations, sensitivity, ML modeling)
                              │
                              ▼
                          figures/
```

## Reproducing the analyses

The notebooks and scripts in this repo were developed across several years
and several environments; there is no single pinned environment file. The
guidance below is what we recommend if you want to rerun things end-to-end.

### Software you will need

- **Python 3.9+** with: `numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`,
  `scikit-learn`, `xgboost`, `imbalanced-learn`, `logomaker`, `seqfold`,
  `jupyter`. (See `2025_11_30_twist_analysis/README.md` for the most
  complete list — that directory is the most recent and most reproducible.)
- **R 4.1+** with `tidyverse`, `ggplot2`, `cowplot`, and `ape` (used by the
  figure scripts and the editing-rate plots).
- **Scala / sbt** — only if you want to rerun
  `2023_02_05_stats_to_base_editing/convert_full_editing_output_to_per_target_file.scala`.
  The Python equivalent (`convert_stats_to_editing_table.py`) is sufficient
  for most use cases.
- A web browser if you want to use the interactive D3 viewers in
  `2023_07_24_base_editing_heatmap_counts/` and `palindrome_analysis/`.

### Where the raw data lives

Large raw sequencing files are **not stored in this repo**. Most analysis
directories contain processed `.stats`, `.txt`, `.csv`, or `.fasta.gz` files
that are sufficient to reproduce the figures from. Raw FASTQ data will be
deposited in GEO; the submission template is in
[`data_deposits/`](data_deposits/).

If a notebook references a file you don't see locally, check:

1. The same filename inside `2025_11_30_twist_analysis/Data/` (the
   consolidated dataset for the most recent analyses), or
2. `2023_07_24_base_editing_heatmap_counts/full_data.csv` (the consolidated
   dataset used by earlier figures).

### Suggested order to read things in

If you are new to the project and want to understand how the analyses build
on each other, we suggest:

1. Skim [`palindrome_analysis/`](palindrome_analysis/) to see the original
   LEFT/RIGHT/BOTH/NEITHER classification on a single sample.
2. Read [`2023_02_05_stats_to_base_editing/`](2023_02_05_stats_to_base_editing/)
   to understand how raw aligner stats become tidy event tables.
3. Open the interactive viewer in
   [`2023_07_24_base_editing_heatmap_counts/`](2023_07_24_base_editing_heatmap_counts/)
   (`index.html`) to develop intuition for the data.
4. Read the README in
   [`2025_11_30_twist_analysis/`](2025_11_30_twist_analysis/) for the
   modeling and statistical analyses that go into the manuscript.
5. Browse [`figures/`](figures/) for the final publication figures and
   their generating R scripts.

## Key terminology

| Term | Meaning |
|---|---|
| **LRBN** | LEFT / RIGHT / BOTH / NEITHER — categorization of which side(s) of a palindromic protospacer were base-edited in a given read. |
| **NN** | "Non-NEITHER" — analyses restricted to reads with at least one edit. |
| **FE** | Folding energy (kcal/mol) for the predicted RNA secondary structure of the guide. |
| **Twist v1 / v2** | Two generations of pooled oligo libraries synthesized by Twist Bioscience to test sequence-context effects. |
| **MARC1** | Reference set used for single-cell barcoding experiments. |
| **stats file** | Per-read output from the custom aligner: alignment quality, UMI, target match, and per-position events. |

## License

This project is released under the [MIT License](LICENSE).

## Citation

A citation will be added here once the manuscript / preprint is posted.

## Contact

Aaron McKenna — [@aaronmck](https://github.com/aaronmck). For analysis
questions, please open an issue on this repository.
