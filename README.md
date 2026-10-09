# BIOL 520/820 – Assignment 1: Data Visualization

Group report and R code for Assignment 1 of *Statistical Methods in Life Sciences* (Nazarbayev University, Fall 2026).
The assignment asks for four essential graphical methods, two biomedical datasets from Kaggle, at least two plots per
dataset with a justification, and a discussion of the findings. The graphs use only what the course has covered so far:
descriptive statistics, principles of graph construction, ggplot2, checking the normality assumption, and
correlation / simple linear regression.

## Datasets

| | Dataset | Kaggle page | Size |
|---|---|---|---|
| 1 | Breast cancer gene expression – CuMiDa (GSE45827) | [brunogrisci/breast-cancer-gene-expression-cumida](https://www.kaggle.com/datasets/brunogrisci/breast-cancer-gene-expression-cumida) | 151 samples × 54,675 probes, 6 classes |
| 2 | Indian Liver Patient Records (ILPD) | "Indian Liver Patient Records" (uploaded by Md. Faysal Mahmud); same data as [uciml/indian-liver-patient-records](https://www.kaggle.com/datasets/uciml/indian-liver-patient-records) | 583 patients × 11 columns |

## The four graphical methods and the eight figures

| Method | Dataset 1 | Dataset 2 |
|---|---|---|
| Histogram + density | `fig1_histogram_esr1` | `fig5_histogram_bilirubin` |
| Box plot | `fig2_boxplot_subtype` | `fig6_boxplot_liver_tests` |
| Scatter plot + regression line | `fig3_scatter_esr1_gata3` | `fig7_scatter_protein_albumin` |
| Normal Q–Q plot | `fig4_qq_esr1` | `fig8_qq_alt` |

Each figure is saved as a vector PDF in `figures/pdf/` and a 300-dpi PNG in `figures/png/` (for the slides).

## Repository layout

```
scripts/00_download_data.sh   download both datasets (Kaggle CLI if configured, otherwise public mirrors)
scripts/01_cumida_breast.R    dataset 1: histogram, box plots, scatter + regression, Q-Q plots
scripts/02_liver_patients.R   dataset 2: histogram, box plots, scatter + regression, Q-Q plots
scripts/03_render_report.R    knit report/report.Rmd -> report/Group_X.docx and report/Group_X.pdf
scripts/utils/theme.R         colours, ggplot2 theme, save_fig()
scripts/run_all.sh            run everything in order
figures/pdf, figures/png      the figures
results/                      descriptive tables, test output (Shapiro-Wilk, correlation, regression)
report/report.Rmd             report source; Group_X.docx (editable) and Group_X.pdf (submission)
data/raw/                     downloaded data (git-ignored; see data/README.md)
```

## Reproduce

```bash
Rscript scripts/install_packages.R   # once
bash scripts/run_all.sh
```

## Submission notes

* Fill in the group number and member names in the YAML header of `report/report.Rmd`, re-run
  `scripts/03_render_report.R`, and rename `Group_X.pdf` to your group number before uploading to Moodle.
* `Group_X.docx` is the same report in Word format for editing.
