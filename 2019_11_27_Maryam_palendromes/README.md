# 2019_11_27_Maryam_palendromes

First-pass processing of the synthetic palindromic target library. This is
the **foundational reference dataset** that downstream directories build on
— the FASTA here defines the canonical set of palindromic targets used
throughout the project.

## Inputs

- `synth_constructs_2018_08_28.fasta` — synthetic palindrome reference
  library (~314 unique palindromic targets, ~420 KB).
- `aligned_sequences.fasta.gz` — paired guide–target reads aligned to the
  reference (~50 MB compressed; ~405k aligned read pairs covering ~310 of
  the 314 palindromes).

## Notebook

- `palindrome_processing.ipynb` — loads the palindrome FASTA, parses the
  aligned reads, extracts each (guide, target) pair, performs a base
  editing pileup, and writes per-position editing rates.

## Outputs

- `palindromes_base_counts.txt` — tabular per-target output with the guide
  sequence, the per-target read count, and the editing frequency at each
  of the 20 protospacer positions (`base1`–`base20`).

## Reproducing

Open the notebook and run all cells. No external data dependencies beyond
the two files in this directory.
