Reproducible material for the article Association of Annual Open Surgical Volume and Short-Term Outcomes After Open Abdominal Aortic Aneurysm Repair: A Report from the Quality Registry of the German Society for Vascular Surgery and Vascular Medicine (DGG)

# Association of Annual Open Surgical Volume and Short-Term Outcomes After Open Abdominal Aortic Aneurysm Repair

**A Report from the Quality Registry of the German Society for Vascular Surgery and Vascular Medicine (DGG)**

---

## Overview

This repository contains the analytic code, supplementary materials, and reproducibility resources accompanying the manuscript:

> **Peters F, et al.** "Association of Annual Open Surgical Volume and Short-Term Outcomes After Open Abdominal Aortic Aneurysm Repair: A Report from the Quality Registry of the German Society for Vascular Surgery and Vascular Medicine (DGG)." *European Journal of Vascular and Endovascular Surgery* (under review / in press).

The study investigates the relationship between annual facility-level open repair volume and short-term outcomes following elective open abdominal aortic aneurysm (AAA) repair, using data from the DIGG registry (Deutsches Institut für Gefäßmedizin und Gefäßchirurgie Quality Registry).

---

## Background

Open repair of abdominal aortic aneurysms (OAR) is a complex, high-stakes procedure. Despite the growing role of endovascular repair (EVAR), open surgery remains the standard of care for a substantial proportion of patients. Evidence from other complex surgical procedures suggests a robust inverse relationship between institutional volume and adverse outcomes — yet this relationship has not been comprehensively characterized for elective OAR in the German healthcare context.

This study addresses this gap using a large, prospective, nationwide registry with granular procedure- and facility-level information.

---

## Data Source

Analyses are based on the **DIGG registry** — the quality registry of the **Deutsche Gesellschaft für Gefäßchirurgie und Gefäßmedizin (DGG)**. This registry captures mandatory quality indicators and clinical variables for vascular procedures performed at participating German hospitals.

> **Data access:** The DIGG registry data used in this study are not publicly available due to data protection regulations. Researchers interested in accessing these data should contact the DGG directly. The code in this repository is provided to ensure full analytical transparency and reproducibility given appropriate data access.

---

## Outcomes

| Outcome | Definition |
|---|---|
| **In-hospital mortality** | Death during index hospital stay |
| **Failure to rescue (FtR)** | Death following a major postoperative complication |
| **Major complications** | Composite of predefined severe adverse events |
| **Prolonged length of stay** | Length of stay exceeding the 75th percentile |

---

## Statistical Methods

The analytic pipeline combines several complementary approaches to ensure robustness:

### Primary Analysis
- **Generalized additive mixed-effects models (GAMM)** with penalised splines for annual open OAR volume and random effects for treating facility, fitted via `mgcv`
- Odds ratios (ORs) and 95% confidence intervals reported at pre-specified volume thresholds

### Doubly Robust Estimation
- Individual-level **propensity scores** estimated via logistic regression for high- vs. low-volume centre assignment
- **Inverse probability weighting (IPW)** combined with outcome regression to yield doubly robust treatment effect estimates
- Facility-level **open procedure share** incorporated as an instrument for sensitivity analysis

### Threshold Detection
- **Segmented regression** with the Davies test to identify statistically supported volume thresholds in continuous volume–outcome relationships

### Robustness Checks
- Multiple model specifications including adjustment for case-mix, emergency medicine capability, and structural facility characteristics
- Sensitivity analyses with alternative volume categorisations and outcome definitions

---

## Repository Structure

```
.
├── README.md
├── LICENSE
│
├── R/
│   ├── 00_packages.R              # Package dependencies
│   ├── 01_data_preparation.R      # Cleaning, variable construction, centre aggregation
│   ├── 02_descriptives.R          # Table 1 (tableone + flextable)
│   ├── 03_gamm_primary.R          # GAMM mixed-effects models (primary analysis)
│   ├── 04_propensity_score.R      # PS estimation and IPW
│   ├── 05_doubly_robust.R         # Doubly robust estimation
│   ├── 06_threshold_detection.R   # Segmented regression and Davies test
│   ├── 07_sensitivity.R           # Robustness and sensitivity checks
│   └── 08_figures.R               # Publication-quality figures (forest plots, density plots, GAM curves)
│
├── output/
│   ├── figures/                   # Publication-ready figure exports (.pdf, .png)
│   └── tables/                    # Formatted table exports
│
└── supplementary/
    ├── supplementary_methods.pdf  # Extended methods description
    └── supplementary_tables.pdf   # Supplementary results tables
```

---

## Software and Reproducibility

All analyses were conducted in **R** (version ≥ 4.3). The following key packages are used:

| Package | Purpose |
|---|---|
| `data.table` | Data manipulation and aggregation |
| `mgcv` | Generalised additive mixed models |
| `MatchIt` / `WeightIt` | Propensity score methods and IPW |
| `segmented` | Segmented regression and Davies test |
| `tableone` | Descriptive statistics (Table 1) |
| `flextable` | Formatted table output |
| `ggplot2` | Data visualisation |
| `patchwork` | Figure composition |

To install all required packages, run:

```r
source("R/00_packages.R")
```

A `renv` lockfile is included to enable exact package version restoration:

```r
renv::restore()
```

---

## Key Findings

- **Failure to rescue** and **in-hospital mortality** both decreased with increasing annual open OAR volume (OR ≈ 0.70 and ≈ 0.73 per volume increment, respectively), robust across all model specifications.
- **Emergency medicine capability** at the treating facility was independently associated with an approximately **50% reduction in failure-to-rescue odds**, underscoring the importance of structural hospital characteristics beyond surgical volume alone.
- Volume effects on major complications and prolonged stay were attenuated after full covariate adjustment.
- Threshold analyses suggested that benefits may plateau beyond a moderate annual case volume, with no single sharp threshold identified.

---

## Citation

If you use the code or methods from this repository, please cite the original publication:

```
Peters F, et al. Association of Annual Open Surgical Volume and Short-Term
Outcomes After Open Abdominal Aortic Aneurysm Repair: A Report from the
Quality Registry of the German Society for Vascular Surgery and Vascular
Medicine (DGG). European Journal of Vascular and Endovascular Surgery.
[Year]; [Volume]([Issue]):[Pages]. DOI: [DOI]
```

---

## Authors and Affiliations

**Corresponding author:**  
PD Dr. med. Frederik Peters  
Hamburgisches Krebsregister / Department of Epidemiology and Health Services Research  
University of Hamburg, Medical Faculty  
Hamburg, Germany

---

## Funding and Conflicts of Interest

Please refer to the published manuscript for full declarations of funding sources and conflicts of interest.

---

## License

This code is released under the [MIT License](LICENSE). Data are not included and remain subject to the data governance framework of the DGG registry.

---

## Contact

For questions regarding the analysis code or methodology, please open a GitHub Issue or contact the corresponding author via the published manuscript.
