# Data

Run `bash scripts/00_download_data.sh` from the project root to populate `data/raw/`.
Raw files are not committed (the gene-expression matrix is ~140 MB).

| File | Source | Size |
|---|---|---|
| `Breast_GSE45827.csv` | Kaggle: [brunogrisci/breast-cancer-gene-expression-cumida](https://www.kaggle.com/datasets/brunogrisci/breast-cancer-gene-expression-cumida) (CuMiDa, GEO GSE45827); mirror <https://sbcb.inf.ufrgs.br/cumida> | 151 rows × 54,677 columns (`samples`, `type`, 54,675 probes) |
| `indian_liver_patient.csv` | Kaggle: "Indian Liver Patient Records" (Md. Faysal Mahmud), identical to [uciml/indian-liver-patient-records](https://www.kaggle.com/datasets/uciml/indian-liver-patient-records); mirror: UCI ML Repository id 225 (header line added by the script) | 583 rows × 11 columns |

The script uses the Kaggle CLI when it is configured (`~/.kaggle/kaggle.json`); otherwise it downloads the same
files from the public mirrors.
