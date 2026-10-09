#!/usr/bin/env bash
# Reproduce everything: download data -> figures for dataset 1 -> figures for dataset 2 -> PDF report.
# On the NU HPC cluster R is provided as a module; elsewhere make sure Rscript is on PATH.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
if command -v module >/dev/null 2>&1 && ! command -v Rscript >/dev/null 2>&1; then
  module load R/4.4.2-gfbf-2024a
fi
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"
bash scripts/00_download_data.sh
Rscript scripts/01_cumida_breast_plots.R
Rscript scripts/02_heart_failure_plots.R
Rscript scripts/03_render_report.R
