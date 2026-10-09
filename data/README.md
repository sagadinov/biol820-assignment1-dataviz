# Data

Run `bash scripts/00_download_data.sh` from the project root to populate `data/raw/`.
Raw files are not committed (the gene-expression matrix is ~140 MB).

| File | Source | Size |
|---|---|---|
| `Breast_GSE45827.csv` | Kaggle: [brunogrisci/breast-cancer-gene-expression-cumida](https://www.kaggle.com/datasets/brunogrisci/breast-cancer-gene-expression-cumida) (CuMiDa, GEO GSE45827); mirror: <https://sbcb.inf.ufrgs.br/cumida> | 151 samples x 54,677 columns (`samples`, `type`, 54,675 probes) |
| `heart_failure_clinical_records_dataset.csv` | Kaggle: [andrewmvd/heart-failure-clinical-data](https://www.kaggle.com/datasets/andrewmvd/heart-failure-clinical-data); mirror: UCI ML Repository (id 519) | 299 patients x 13 variables |
| `GPL570.annot.gz` | NCBI GEO platform annotation (Affymetrix HG-U133 Plus 2.0) used to map probe IDs to gene symbols | 54,675 probes |

The script uses the Kaggle CLI when it is installed and configured (`~/.kaggle/kaggle.json`);
otherwise it downloads the identical files from the public mirrors listed above.
