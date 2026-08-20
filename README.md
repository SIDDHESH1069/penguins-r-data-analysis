# Penguins R Data Analysis — Week 3

Statistical analysis and predictive modeling on the **palmerpenguins** dataset in R.

## Contents
- `week3_model.R` — full reproducible analysis script
- `data/penguins.csv` — raw dataset
- `outputs/week3_console_output.txt` — complete console log of the run
- `docs/Week3_Statistical_Analysis_and_Predictive_Modeling_R.docx` — written report
- `figures/` — 13 exported PNG figures

## Analysis
1. **Exploratory statistics** — Shapiro-Wilk normality, Bartlett variance homogeneity, Welch t-tests / Wilcoxon, two-way ANOVA with interaction, chi-square independence, Pearson/Spearman correlations.
2. **Modeling** — four binary logistic regression models predicting penguin sex, compared with stratified 10-fold cross-validation and stepwise AIC selection.
3. **Diagnostics** — ROC curves, confusion matrix, threshold sweep, calibration, binned residuals, Cook's distance / leverage, linear model diagnostics.

## Reproduce
```bash
Rscript week3_model.R
```
Requires R (>= 4.3) with `ggplot2` and `MASS`.
