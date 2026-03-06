# 🧬 Historical Control Data (HCD) Quality Control & Variance Analysis Toolkit

Tools for **quality control**, **process monitoring**, and **regulatory-grade variance analysis** of Historical Control Data (HCD) from transgenic rodent mutation assays (OECD TG 488).

This repository provides two complementary analysis pipelines:

1. **QC Control Charts (I-MR Charts)** — detect out-of-control studies and abnormal animal values
2. **Variance Components Analysis (REML Models)** — quantify between-study and within-study variability

Designed for **genetic toxicology laboratories**, **regulatory submissions**, and **statistical review** workflows.

## Overview

This repository contains two annotated R scripts:

- `hcd_qc_clean.R` — generates Individual–Moving Range (I-MR) control charts for mutation frequency quality control
- `hcd_variance_enhanced.R` — performs REML-based variance components analysis for historical control data

The QC workflow focuses on **study stability and outlier detection**, while the variance workflow focuses on **regulatory-style decomposition of variability** across studies and sampling days.

---

## QC Control Charts (I-MR)

### Purpose

Detect studies or animals that deviate from expected mutation frequency variation using **Individual–Moving Range control charts**, a standard quality control approach for continuous laboratory data. The QC script uses study-mean normalization so that the analysis focuses on **within-study animal-to-animal variation** rather than systematic study-level shifts.

### Key Features

- Individual–Moving Range control charts
- Study-mean normalization
- Moving Range–based sigma estimation using `MR̄ / d2`
- ±2 SD warning limits and ±3 SD action limits
- Automatic flagging of out-of-control observations
- Stability index calculation
- Publication-quality chart output

### Input

Pre-cleaned CSV containing columns such as:

- `Study_No`
- `Sampling day`
- `Animal No`
- `Mut Freq x 10-6`
- `Log10_MF`

The script includes flexible column matching to handle variation in Excel-exported column names.

### Outputs

- QC summary table with mean, moving range, sigma, and control limits
- Flagged dataset with animal-level classifications
- PNG control chart

### Methodology

The QC script:

1. Reads cleaned input data
2. Identifies required study and log-transformed mutation frequency columns
3. Optionally excludes Study 1 from control-limit calculation because of oversized group size
4. Estimates sigma from average moving range
5. Calculates a stability index from direct SD versus MR-based sigma
6. Normalizes animal values to the overall study mean
7. Flags observations as:
   - `In control`
   - `Warning (outside 2σ)`
   - `Action (outside 3σ)`
8. Produces a control chart and summary outputs

This workflow is intended for routine monitoring of historical control mutation frequency performance. 

### Example Command

```bash
Rscript hcd_qc_clean.R input.csv summary.csv flagged_data.csv plot.png 31
