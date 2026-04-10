# Twist Analysis - Palindrome Base Editing Study

This directory contains Jupyter notebooks analyzing CRISPR base editing outcomes from Twist library experiments. The analysis focuses on palindromic target sequences and their editing patterns.

## Notebooks

### 1. Highest_EditRate_Targets.ipynb
**Purpose**: Identify and characterize high-efficiency guide RNAs through sequence analysis and free energy correlation studies.

**Key Analyses**:
- Nucleotide composition analysis across all 208 unique targets (20bp sequences)
- Logo plot generation for sequence motif visualization
- Top 20 targets extraction for different edit categories:
  - Any edits
  - LEFT, RIGHT, BOTH, NEITHER edits
  - NN (non-NEITHER) variants
- Folding energy (FE) correlation analysis:
  - Scatter plots: edit rates vs target free energy
  - Statistical testing: Wilcoxon rank-sum tests comparing top 20 vs bottom targets
  - Bar plots: FE distribution comparisons
- Coverage sensitivity analysis across experiments 1-3

**Main Datasets**:
- `./Data/match_only_twist2_data.csv` - Complete Twist2 data
- `./Data/LOGO_match_only_twist2_data.csv` - Aggregated experiments 1-3 for any_edits
- `./Data/LRBN_match_only_twist2_data.csv` - Left/Right/Both/Neither categorized data

**Key Outputs**:
- Logo plots for sequence motifs
- Correlation plots between free energy and edit rates
- Statistical comparisons (p-values from Wilcoxon tests)
- Ruleset data exports for different edit categories

**Key Findings**:
- Free energy correlation with edit rates varies by edit type
- No significant FE difference for any_edits (p=0.198)
- Marginal significance for RIGHT edits (p=0.063)

---

### 2. Twist2_LRBN_Corr.ipynb
**Purpose**: Examine correlation of edit outcomes (Left, Right, Both, Neither) across replicate experiments.

**Key Analyses**:
- Cross-experiment correlation analysis for LRBN edit patterns
- Replicate reproducibility assessment
- Data quality filtering (removing NaN columns)
- Both inclusive and exclusive analyses of NEITHER edits
- High-correlation target identification

**Main Datasets**:
- `../Data/LRBN_match_only_expt0_data.csv`
- `../Data/LRBN_match_only_twist2_expt1_data.csv`
- `../Data/LRBN_match_only_twist2_expt2_data.csv`
- `../Data/LRBN_match_only_twist2_expt3_data.csv`

**Analysis Types**:
- NN (non-NEITHER) fraction correlations
- LEFT, RIGHT, BOTH edit correlations with/without NEITHER
- Index consistency checks across experiments

---

### 3. Twist2_Sensitivity.ipynb
**Purpose**: Assess robustness of edit rate correlations across different coverage thresholds.

**Key Analyses**:
- Coverage threshold sensitivity (0-39) for correlation stability
- Any_edits correlation across experiments at varying coverage
- LRBN edit correlations:
  - LEFT, RIGHT, BOTH, NEITHER
  - NN variants (excluding NEITHER)
- Line plots showing correlation trends vs coverage threshold

**Main Datasets**:
- `../Data/Sensitivity/twist2_expt123_anyedit_corrs.csv`
- `../Data/Sensitivity/twist2_expt123_LRBN_LEFT_corrs.csv`
- `../Data/Sensitivity/twist2_expt123_LRBN_RIGHT_corrs.csv`
- `../Data/Sensitivity/twist2_expt123_LRBN_BOTH_corrs.csv`
- `../Data/Sensitivity/twist2_expt123_LRBN_NEITHER_corrs.csv`
- NN variants for LEFT, RIGHT, BOTH

**Key Outputs**:
- Coverage sensitivity plots for each edit type
- Correlation stability assessment across thresholds
- Experiment-to-experiment correlation comparisons

**Key Insights**:
- Determines optimal coverage thresholds for reliable correlation measurements
- Identifies which edit types show stable correlations vs coverage-dependent patterns

---

### 4. MaryamSMOTEClean.ipynb
**Purpose**: Machine learning prediction of CRISPR edit outcomes using SMOTE (Synthetic Minority Oversampling Technique) to handle class imbalance.

**Key Analyses**:
- **XGBoost Regression**: Predicting fractional editing rates (`frac_NN_edits`)
  - Hyperparameter grid search (learning_rate, n_estimators, max_depth)
  - Performance: correlation and RMSE on test sets
- **SMOTE Classification**: Multi-class edit outcome prediction (LEFT, RIGHT, BOTH)
  - SMOTENC for mixed continuous/categorical features
  - Handles imbalanced classes in edit outcomes
- **Binary Classification**: High vs low efficiency guides
  - 90th and 95th percentile thresholds
  - Feature importance analysis (XGBoost and Logistic Regression)
- Cross-length analysis: 17, 18, and 20 nucleotide guide RNAs

**Main Datasets**:
- `./Data/full_target_table_expts123_dummies_len20.csv`
- `./Data/full_target_table_expts123_dummies_len17.csv`
- `./Data/full_target_table_expts123_dummies_len18.csv`

**Features**:
- **Sequence features**: One-hot encoded positional nucleotides (Pos0_A, Pos1_T, etc.)
- **Compositional features**: GC%, A%, folding energy (FE)
- **Target variables**:
  - `frac_NN_edits` - regression target
  - `NN_most` - multi-class edit type
  - `NN_frac_over90/95` - binary high-efficiency classification

**Machine Learning Pipeline**:
1. Train-test split (80/20 for regression, 75/25 for classification)
2. SMOTENC oversampling on training data only
3. Model training (XGBoost or ElasticNet Logistic Regression)
4. Evaluation with precision, recall, accuracy, F1-score
5. Feature importance visualization

**Key Outputs**:
- Regression performance plots (true vs predicted)
- Confusion matrices for multi-class and binary classification
- Feature importance rankings (XGBoost gain scores)
- Coefficient visualizations for logistic regression (positive/negative)
- Export: `*_feature_selection_2025_08_11.tsv` files

**Key Findings**:
- XGBoost achieves strong correlation for edit rate prediction
- SMOTE successfully balances minority edit outcome classes
- Feature importance reveals critical sequence positions and motifs
- Both tree-based (XGBoost) and linear (LogReg) approaches provide complementary insights

**Model Comparison**:
- **XGBoost**: Better performance, non-linear relationships, feature interactions
- **Logistic Regression**: Interpretable coefficients, linear decision boundaries

---

## Data Directory Structure

```
./Data/
├── match_only_twist2_data.csv              # Main Twist2 dataset
├── LOGO_match_only_twist2_data.csv         # Aggregated any_edits data
├── LRBN_match_only_twist2_data.csv         # LRBN categorized data
├── ruleset/                                 # Nucleotide fraction data for rulesets
│   ├── any_edit_bar_melt.csv
│   ├── LEFT_edit_bar_melt.csv
│   ├── RIGHT_edit_bar_melt.csv
│   ├── BOTH_edit_bar_melt.csv
│   └── NEITHER_edit_bar_melt.csv
├── Sensitivity/                             # Coverage sensitivity analysis
│   ├── twist2_expt123_anyedit_corrs.csv
│   ├── twist2_expt123_LRBN_LEFT_corrs.csv
│   ├── twist2_expt123_LRBN_RIGHT_corrs.csv
│   ├── twist2_expt123_LRBN_BOTH_corrs.csv
│   ├── twist2_expt123_LRBN_NEITHER_corrs.csv
│   └── *_NN_corrs.csv variants
└── full_target_table_expts123_dummies_len*.csv  # ML feature matrices
```

## Figures Directory Structure

```
./Figures/Twist2/
├── LogoPlots/                              # Sequence logo visualizations
├── FE/                                     # Folding energy analysis
│   ├── FE_anyedit_corr.png
│   ├── FE_LEFTedit_corr.png
│   ├── FE_RIGHTedit_corr.png
│   ├── FE_BOTHedit_corr.png
│   └── NN_* variants
└── Correlations/                           # Coverage sensitivity plots
    ├── any_edit_coverage_sensitivity.png
    ├── LEFT_coverage_sensitivity.png
    ├── RIGHT_coverage_sensitivity.png
    ├── BOTH_coverage_sensitivity.png
    ├── NEITHER_coverage_sensitivity.png
    └── *_NN_* variants

./revisit_figs/SMOTE_diff_lengths/
├── outcomes_w_neither.png                  # Edit outcome distributions
├── outcomes_no_neither.png
├── SMOTE_most_editpos_NN_CM_len*.png      # Multi-class confusion matrices
├── SMOTE_xgb_top_90_NN_FI_len*.png        # XGBoost feature importance
└── SMOTE_*_top_90_perc_guides_NN_CM_len*.png  # Binary classification results
```

## Key Terminology

- **LRBN**: Left, Right, Both, Neither - categorization of edit positions in palindromes
- **NN**: Non-NEITHER - analyses excluding targets with no edits
- **FE**: Folding Energy (kcal/mol) - predicted RNA secondary structure stability
- **SMOTE/SMOTENC**: Synthetic minority oversampling for imbalanced classification
- **Target sequence**: 20bp palindromic guide RNA sequence
- **Coverage threshold**: Minimum read depth for including a target in analysis

## Dependencies

- pandas, numpy - data manipulation
- matplotlib, seaborn - visualization
- logomaker - sequence logo generation
- seqfold - RNA folding energy calculation
- scikit-learn - machine learning utilities
- xgboost - gradient boosting models
- imbalanced-learn (imblearn) - SMOTE implementation
- scipy.stats - statistical testing

## Analysis Workflow

1. **Quality Control**: Sensitivity analysis to determine optimal coverage thresholds
2. **Correlation Analysis**: Assess replicate reproducibility (Twist2_LRBN_Corr)
3. **Feature Analysis**: Identify sequence determinants (Highest_EditRate_Targets)
4. **Predictive Modeling**: Train ML models for outcome prediction (MaryamSMOTEClean)

## Contact

For questions about these analyses, please refer to the original research publication or contact the study authors.
