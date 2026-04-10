# 2019_12_12_control_data

Control vs. edited comparison of the palindrome library, with **UMI-aware
deduplication** added on top of the pileup pipeline from
`2019_11_27_Maryam_palendromes/`. The goal is to quantify the background
editing rate in unedited control libraries and contrast it with the edited
samples.

## Inputs

- `data/merged_reads.fastq.aligned.gz` — control / baseline alignment
  (~16 MB).
- `data/merged_umi_controls.fastq.aligned.gz` — UMI-deduplicated control
  sample (~24 MB).
- `data/merged_umi_edited.fastq.aligned.gz` — UMI-deduplicated edited
  sample (~19 MB).
- The palindrome reference FASTA from
  [`../2019_11_27_Maryam_palendromes/`](../2019_11_27_Maryam_palendromes/).

## Notebooks and scripts

- `notebooks/palindrome_processing.ipynb` — same per-position pileup as in
  `2019_11_27`, but with stricter indel filtering (drops reads with more
  than 30 dashes).
- `notebooks/UMI_palindrome_processing.ipynb` — UMI-aware variant that
  collapses duplicate reads sharing a UMI before tabulating editing.
- `notebooks/plot_editing.R`, `notebooks/umi_processing_rnd2.R` — R
  visualization scripts.

## Outputs

In `data/`:

- `palindromes_base_counts_control.txt`,
  `palindromes_base_counts_edited.txt` — per-position editing rates for
  the control and edited libraries.
- `palindromes_base_counts_control_UMI.txt`,
  `palindromes_base_counts_edited_UMI.txt` — UMI-deduplicated versions of
  the above.
- `palindrome_sequence_controls.txt`, `palindrome_sequence_cases.txt` —
  target-matching statistics (fraction of reads hitting the expected vs.
  unexpected target).

In `plots/`: PNG plots of editing proportions and per-position editing
rates split by library and UMI status.
