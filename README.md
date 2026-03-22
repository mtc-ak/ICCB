# ICCB

Microbiome-based colorectal cancer prediction using public metagenomic cohorts.

## Project Summary
This repository contains code and public project materials for a study evaluating gut microbiome-based machine learning models for colorectal cancer prediction. The analysis compares Logistic Regression, Random Forest, and XGBoost using a development cohort and an independent external validation cohort.

## Study Design
- Development cohort: `ZellerG_2014`
- External validation cohort: `YachidaS_2019`
- Outcome: colorectal cancer (`CRC`) versus control
- Preprocessing: rarefaction, prevalence filtering, CLR transformation
- Evaluation: AUC, PR-AUC, recall, F1-score, Brier score, ECE, external validation

## Main Analysis File
- `codex_crc/run_crc_pipeline.R`

## Public Contents
This repository is intended to contain:
- analysis code
- manuscript drafts and documentation
- public figures
- lightweight summary outputs

## Private / Large Data Policy
Large files, local runtime directories, and private resources should be managed separately in:
- `ICCB-data`

Examples include:
- raw or downloaded datasets
- local package libraries
- temporary directories
- large intermediate outputs
- `.rds` artifacts when not needed publicly

## Repositories
- Public homepage: https://mtc-ak.github.io/ICCB
- Public GitHub repository: https://github.com/mtc-ak/ICCB
- Private data repository: https://github.com/mtc-ak/ICCB-data

## Notes
This project is intended for research documentation, reproducible analysis, and public-facing project presentation.
