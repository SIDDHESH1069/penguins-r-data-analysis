# Palmer Penguins — Data Analysis & Visualization in R

Internship project (Weeks 1–2) exploring the **Palmer Penguins** dataset (344 penguins, 3 species, 3 islands in the Palmer Archipelago, Antarctica) using R.

## What this project is about

The goal is to turn a raw, slightly messy ecological dataset into a clear story: what distinguishes Adelie, Chinstrap and Gentoo penguins, how their body measurements relate to each other, and how those patterns differ by island, sex and year.

## What I am doing

**Week 1 — Data Cleaning & Preliminary Analysis**
- Loaded the raw dataset and inspected structure, types and ranges
- Handled missing values and invalid/ambiguous records (e.g. unknown sex)
- Converted variables to correct types, checked outliers
- Produced descriptive statistics and group summaries

**Week 2 — Data Visualization & Insight Communication**
- Designed 10 visualizations with `ggplot2`, `lattice` and base R:
  1. Bar chart — species count per island
  2. Histogram — flipper length distribution by species
  3. Density curves — body mass by species
  4. Boxplot + jitter — body mass by species and sex
  5. Scatter + linear fits — flipper length vs body mass
  6. Scatter — bill length vs bill depth (Simpson's paradox)
  7. Line chart — mean body mass by year (2007–2009)
  8. Heatmap — correlation matrix of body measurements
  9. Scatterplot matrix (lattice SPLOM)
  10. Mosaic plot — island × species composition
- Backed the visuals with statistics: one-way ANOVA on body mass by species and a linear regression of body mass on flipper length
- Wrote an interpretation of every figure for a non-technical reader

## Key findings

- Gentoo penguins are clearly the heaviest and longest-flippered; Adelie and Chinstrap overlap in mass but differ in bill shape
- Flipper length and body mass are strongly positively correlated (r ≈ 0.87)
- Bill length vs bill depth shows a negative trend overall but a positive trend within each species — a textbook Simpson's paradox
- Islands are not species-neutral: Gentoo are found only on Biscoe, Chinstrap only on Dream, Adelie on all three
- Males are consistently heavier than females within every species

## Files

| File | Description |
|---|---|
| `week2_viz.R` | Self-contained R script that downloads the data and produces all 10 figures plus the statistics |
| `penguins.csv` | The dataset used |
| `Week2_Data_Visualization_and_Insight_Communication_R.docx` | Full written report with narrative, code, figures and interpretation |

## How to reproduce

```r
install.packages(c("ggplot2", "lattice"))
source("week2_viz.R")
```

Figures are written to a `figs/` folder.

## Tools

R 4.5 · ggplot2 · lattice · base R graphics

## Data source

Horst A.M., Hill A.P., Gorman K.B. (2020). *palmerpenguins: Palmer Archipelago (Antarctica) penguin data.*
