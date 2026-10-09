#!/usr/bin/env bash
# Download the two Kaggle datasets used in BIOL 520/820 Assignment 1.
#
# Dataset 1: Breast cancer gene expression (CuMiDa, GSE45827)
#   Kaggle : https://www.kaggle.com/datasets/brunogrisci/breast-cancer-gene-expression-cumida
#   Mirror : https://sbcb.inf.ufrgs.br/cumida  (the original CuMiDa repository; identical CSV)
# Dataset 2: Heart failure clinical records
#   Kaggle : https://www.kaggle.com/datasets/andrewmvd/heart-failure-clinical-data
#   Mirror : https://archive.ics.uci.edu/dataset/519  (UCI; identical CSV)
# Extra  : GPL570 (Affymetrix HG-U133 Plus 2.0) probe annotation from NCBI GEO, used to map
#          probe IDs to gene symbols for labelling.
#
# If the Kaggle CLI is installed and configured (~/.kaggle/kaggle.json) the Kaggle copies are
# downloaded; otherwise the public mirrors of the same files are used.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW="$ROOT/data/raw"
mkdir -p "$RAW"
cd "$RAW"

have_kaggle=0
if command -v kaggle >/dev/null 2>&1 && [ -f "$HOME/.kaggle/kaggle.json" ]; then have_kaggle=1; fi

# ---- Dataset 1: CuMiDa Breast_GSE45827 -------------------------------------------------
if [ ! -f Breast_GSE45827.csv ]; then
  if [ $have_kaggle -eq 1 ]; then
    echo "[kaggle] brunogrisci/breast-cancer-gene-expression-cumida"
    kaggle datasets download -d brunogrisci/breast-cancer-gene-expression-cumida --unzip -p .
  else
    echo "[mirror] CuMiDa Breast_GSE45827.csv (~140 MB)"
    curl -L --retry 3 -o Breast_GSE45827.csv \
      "https://sbcb.inf.ufrgs.br/data/cumida/Genes/Breast/GSE45827/Breast_GSE45827.csv"
  fi
fi

# ---- Dataset 2: Heart failure clinical records ----------------------------------------
if [ ! -f heart_failure_clinical_records_dataset.csv ]; then
  if [ $have_kaggle -eq 1 ]; then
    echo "[kaggle] andrewmvd/heart-failure-clinical-data"
    kaggle datasets download -d andrewmvd/heart-failure-clinical-data --unzip -p .
  else
    echo "[mirror] UCI heart failure clinical records"
    curl -L --retry 3 -o heart_failure.zip \
      "https://archive.ics.uci.edu/static/public/519/heart+failure+clinical+records.zip"
    unzip -o -q heart_failure.zip && rm -f heart_failure.zip
  fi
fi

# ---- GPL570 probe annotation (gene symbols) ---------------------------------------------
if [ ! -f GPL570.annot.gz ]; then
  echo "[GEO] GPL570.annot.gz (~8 MB)"
  curl -L --retry 3 -o GPL570.annot.gz \
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"
fi

echo "Done. Files in $RAW:"
ls -lh "$RAW"
