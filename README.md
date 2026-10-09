# BIOL 520/820 – Assignment 1: Data Visualization

Group report and R code for Assignment 1 (*Statistical Methods in Life Sciences*, Nazarbayev University, Fall 2026).
The assignment asks for four essential graphical methods of the discipline, two biomedical datasets from Kaggle,
at least two plots per dataset with a justification of the chosen graphical method, and a discussion of the findings.

## Datasets

| | Dataset | Kaggle page | Size |
|---|---|---|---|
| 1 | **Breast cancer gene expression – CuMiDa (GSE45827)** | [brunogrisci/breast-cancer-gene-expression-cumida](https://www.kaggle.com/datasets/brunogrisci/breast-cancer-gene-expression-cumida) | 151 samples × 54,675 probes, 6 classes |
| 2 | **Heart failure clinical records** | [andrewmvd/heart-failure-clinical-data](https://www.kaggle.com/datasets/andrewmvd/heart-failure-clinical-data) | 299 patients × 13 variables |

## Figures

| File (in `figures/pdf` and `figures/png`) | Graphical method | Dataset |
|---|---|---|
| `fig1_1_pca` | Scatter plot (PCA score plot) + scree bars | 1 |
| `fig1_2_heatmap` | Clustered heatmap with annotation tracks | 1 |
| `fig1_3_volcano` | Scatter plot (volcano plot, basal-like vs luminal A) | 1 |
| `fig1_4_marker_violins` | Violin + box plots of marker genes by subtype | 1 |
| `fig2_1_clinical_violins` | Violin + box plots of clinical variables by outcome | 2 |
| `fig2_2_scatter_ef_creat` | Scatter plot with marginal densities | 2 |
| `fig2_3_km_ejection` | Kaplan–Meier survival curves with risk table | 2 |
| `fig2_4_corr_heatmap` | Correlation heatmap | 2 |

Every figure is saved both as a vector PDF (used in the report) and a 300-dpi PNG (for slides).

## Repository layout

```
scripts/00_download_data.sh      download both datasets (+ GPL570 probe annotation)
scripts/01_cumida_breast_plots.R dataset 1: PCA, heatmap, volcano, marker violins
scripts/02_heart_failure_plots.R dataset 2: violins, scatter, Kaplan–Meier, correlation heatmap
scripts/03_render_report.R       knit report/report.Rmd -> report/Group_X.pdf
scripts/utils/theme.R            shared ggplot2 theme, colour palette, save helpers
scripts/run_all.sh               run everything in order
figures/pdf, figures/png         the plots
results/                         tables behind the plots (class counts, PCA variance, DE results, tests, ...)
report/report.Rmd                R Markdown source of the report; report/Group_X.pdf is the submission
data/raw/                        downloaded data (git-ignored; see data/README.md)
```

## Reproduce

```bash
Rscript scripts/install_packages.R   # once: CRAN packages (and tinytex::install_tinytex() for PDF output)
bash scripts/run_all.sh              # download data, make figures, knit report
```

Developed with R 4.4.2; packages: ggplot2, dplyr, tidyr, data.table, ggrepel, patchwork, pheatmap,
matrixStats, survival, rmarkdown/knitr.

## Submission notes

* `report/Group_X.pdf` is the report to upload to Moodle: rename it to your group number and fill in the
  member names in the YAML header of `report/report.Rmd` (`author:`) before the final knit.
* Dataset 2 (heart failure clinical records) was chosen to complement the gene-expression dataset with a
  low-dimensional clinical dataset; swap it for another Kaggle biomedical dataset if the group prefers.
* The PNG versions of the figures in `figures/png/` are intended for the 20-minute presentation slides.
