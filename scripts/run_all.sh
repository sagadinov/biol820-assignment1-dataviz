#!/usr/bin/env bash
# Reproduce everything: download data -> figures -> report (DOCX + PDF).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"
if command -v module >/dev/null 2>&1 && ! command -v Rscript >/dev/null 2>&1; then
  module load R/4.4.2-gfbf-2024a            # NU HPC cluster
fi
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"
bash scripts/00_download_data.sh
Rscript scripts/01_cumida_breast.R
Rscript scripts/02_liver_patients.R
Rscript scripts/03_render_report.R
