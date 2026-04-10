# 2022_07_18_base_editing_t7_rnf2

Per-position base editing analysis at the **T7** and **RNF2** loci across
four pig-UMI samples. Outputs heatmaps showing editing frequency at every
protospacer position for each sample.

## Inputs

Four gzipped FASTAs of aligned paired sequences (each read carries a 26 bp
target for both T7 and RNF2):

- `p7_PigUMI_8_CTCTTCTGCT_sequences.fasta.gz`
- `p7_PigUMI_ct_GGCGCCTTAA_sequences.fasta.gz`
- `p7_PigUMI_pRDA_CGAGGTATAA_sequences.fasta.gz`
- `p7_PigUMI_pTpR8_GCTGCTGGTA_sequences.fasta.gz`

## Notebook

- `2022_07_18_process_guides.ipynb` — parses each FASTA, extracts the T7
  and RNF2 target windows, computes per-position base counts, normalizes
  by read depth, and renders PNG heatmaps for each sample × locus
  combination.

## Outputs

PNG heatmaps written next to the notebook, one per (sample, locus)
combination, showing per-position base editing rates at T7 and RNF2.
