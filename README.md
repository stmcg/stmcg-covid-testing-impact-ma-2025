# Self-Testing for SARS-CoV-2 Using Antigen Detecting Diagnostics

This repository contains the data and code for the manuscript **"Population health and implementation outcomes of self-testing for SARS-CoV-2 using antigen detecting diagnostics: a systematic review and meta-analysis"** by Brümmer et al.

---

# Part 1: Data Processing and Meta-Analyses

The data processing and meta-analyses were performed in R (version 4.4.1) with the following packages:

- **DescTools** (version 0.99.57)
- **dplyr** (version 1.1.4)
- **metafor** (version 4.6-0)
- **openxlsx** (version 4.2.7.1)
- **readxl** (version 1.4.3)

## Files 

- `Data extracted.xlsx` (`data` folder): Extracted data from the primary studies.
- `data-processing.R` (`code` folder): Processes the data. Results are stored in `SRInput1.csv` (`results` folder).
- `data-analysis.R` (`code` folder): Performs meta-analysis for main and subgroup analyses. Results are stored in `meta_analysis_results.xlsx` for the primary analyses and `meta_analysis_results_sensitivity.xlsx` for the sensitivity analyses with inverse-variance weighting (`results` folder).
- `create-forest-plot-data.R` (`code` folder): Formats the output of the meta-analyses for constructing forest plots. Results are stored in `forest_plot_data.xlsx` (`results` folder).

# Part 2: Forest Plots

The forest plots were created in Python (version 3.9.13) with the following libraries:

- **forestplot** (version 0.4.0)
- **pandas** (version 1.5.2)
- **Matplotlib** (version 3.5.1)

## Files

- `forest_plot_data.xlsx` (`results` folder): Input data for forest plots.
- `plots.py` (`code` folder): Generates forest plots. Results are saved in the folder `results/figures`. 