#!/usr/bin/env bash
# Download the two Kaggle datasets used in BIOL 520/820 Assignment 1.
#
# Dataset 1: Breast cancer gene expression - CuMiDa (GSE45827)
#   Kaggle : https://www.kaggle.com/datasets/brunogrisci/breast-cancer-gene-expression-cumida
#   Mirror : https://sbcb.inf.ufrgs.br/cumida  (original CuMiDa repository, identical CSV)
# Dataset 2: Indian Liver Patient Records (ILPD)
#   Kaggle : "Indian Liver Patient Records" (uploaded by Md. Faysal Mahmud); same data as
#            https://www.kaggle.com/datasets/uciml/indian-liver-patient-records
#   Mirror : https://archive.ics.uci.edu/dataset/225  (UCI, same 583 rows; header added below)
#
# If the Kaggle CLI is configured (~/.kaggle/kaggle.json) the Kaggle copies are used,
# otherwise the public mirrors.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW="$ROOT/data/raw"; mkdir -p "$RAW"; cd "$RAW"
have_kaggle=0
if command -v kaggle >/dev/null 2>&1 && [ -f "$HOME/.kaggle/kaggle.json" ]; then have_kaggle=1; fi

# ---- Dataset 1: CuMiDa Breast_GSE45827.csv (~140 MB) -----------------------------------
if [ ! -f Breast_GSE45827.csv ]; then
  if [ $have_kaggle -eq 1 ]; then
    kaggle datasets download -d brunogrisci/breast-cancer-gene-expression-cumida --unzip -p .
  else
    curl -L --retry 3 -o Breast_GSE45827.csv \
      "https://sbcb.inf.ufrgs.br/data/cumida/Genes/Breast/GSE45827/Breast_GSE45827.csv"
  fi
fi

# ---- Dataset 2: indian_liver_patient.csv (583 rows x 11 columns) -----------------------
if [ ! -f indian_liver_patient.csv ]; then
  if [ $have_kaggle -eq 1 ]; then
    kaggle datasets download -d uciml/indian-liver-patient-records --unzip -p .
  else
    curl -L --retry 3 -o ilpd.zip \
      "https://archive.ics.uci.edu/static/public/225/ilpd+indian+liver+patient+dataset.zip"
    unzip -o -q ilpd.zip
    # the UCI file has no header line; add the column names used on Kaggle
    { echo "Age,Gender,Total_Bilirubin,Direct_Bilirubin,Alkaline_Phosphotase,Alamine_Aminotransferase,Aspartate_Aminotransferase,Total_Protiens,Albumin,Albumin_and_Globulin_Ratio,Dataset";
      cat "Indian Liver Patient Dataset (ILPD).csv"; } > indian_liver_patient.csv
    rm -f ilpd.zip "Indian Liver Patient Dataset (ILPD).csv"
  fi
fi
echo "Data files in $RAW:"; ls -lh "$RAW"
