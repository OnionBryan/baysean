# Bayesian Moderated Mediation Tutorial

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Companion code repository for:

 Bayesian Moderated Mediation Analysis: A Step-by-Step Guide 

## Overview

This repository contains complete R code for conducting Bayesian moderated mediation analysis (first-stage moderation; Hayes Model 8) using the `brms` package. The tutorial demonstrates:

- Estimation of the Index of Moderated Mediation (IMM = a₃ × b)
- Monte Carlo simulation for power analysis
- Four empirical demonstrations across different research domains
- Direct comparison with frequentist bootstrapping (lavaan)
- Prior sensitivity analysis

## Requirements

- R ≥ 4.1.0
- Stan/RStan (for brms backend)

### R Packages

```r
install.packages(c(
  "brms",
  "tidyverse",
  "posterior",
  "bayesplot",
  "lavaan",
  "knitr",
  "kableExtra",
  "carData",
  "processR",
  "furrr",
  "future"
))
```

### Optional: Jupyter Notebook Support

A Jupyter notebook version of the tutorial is available (`bayesian_moderated_mediation_tutorial.ipynb`). To run it, you need the R kernel for Jupyter:

```r
# In R:
install.packages("IRkernel")
IRkernel::installspec()
```

Then open the notebook with Jupyter:

```bash
jupyter notebook bayesian_moderated_mediation_tutorial.ipynb
```

## Repository Structure

```
baysean_repo/
├── R/
│   ├── 01_data_preparation.R      # Data cleaning and standardization
│   ├── 02_specify_priors.R        # Prior specification examples
│   ├── 03_fit_model.R             # Model fitting with brms
│   ├── 04_diagnostics.R           # Convergence diagnostics
│   ├── 05_compute_indirect.R      # Indirect effects and IMM
│   ├── 06_test_moderation.R       # Simple slopes and J-N analysis
│   ├── 07_reporting.R             # Publication-ready tables
│   └── complete_template.R        # Full analysis template
├── empirical_analysis/
│   ├── garcia_protest.R           # Social psychology example (SIGNIFICANT)
│   ├── slid_labor.R               # Labor economics example (SIGNIFICANT)
│   ├── talor_media.R              # Communication example (NULL)
│   └── facialburns_clinical.R     # Clinical psychology example (NULL)
├── figures/                        # Output figures (generated)
├── bayesian_moderated_mediation_tutorial.ipynb  # Jupyter notebook (optional)
├── simulation_functions.R          # Monte Carlo power simulation
├── prior_sensitivity.R             # Sensitivity analysis code
├── frequentist_comparison.R        # lavaan bootstrap comparison
├── DESCRIPTION                     # R package metadata
├── renv.lock                       # Package version lock file
├── REFERENCES.md                   # Full bibliography with DOI links
├── LICENSE
├── .gitignore
└── README.md
```

## Quick Start

```r
# Load packages
library(brms)
library(tidyverse)
library(posterior)

# Source the complete template
source("R/complete_template.R")

# Or run step-by-step:
source("R/01_data_preparation.R")
source("R/02_specify_priors.R")
source("R/03_fit_model.R")
source("R/04_diagnostics.R")
source("R/05_compute_indirect.R")
```

## Data Sources

| Dataset | Package | Description |
|---------|---------|-------------|
| Garcia | `processR` | Protest and perceived sexism (N=129) |
| SLID | `carData` | Canadian labor survey (N=2000) |
| Tal-Or | `processR` | Media placement effects (N=123) |
| FacialBurns | Simulated | Clinical distress (N=98) |

## Key Findings

- **Power**: N ≥ 250 required for 80% power to detect medium IMM (0.10)
- **Method choice**: Bayesian preferred for N < 300; either method for N > 300
- **Convergence**: Bayesian and frequentist results converge with adequate sample size

## References

See [REFERENCES.md](REFERENCES.md) for the complete bibliography with DOI links for all cited papers, including:

- Hayes (2015) on the Index of Moderated Mediation
- Bürkner (2017, 2018) on the brms package
- Vehtari et al. (2021) on improved R-hat diagnostics
- All empirical example source papers

## Citation

```bibtex
@article{gibson2025bayesian,
  title={Bayesian Moderated Mediation Analysis: A Step-by-Step Guide},
  author={Gibson, Bryan},
  journal={Psychological Methods},
  year={2025},
  publisher={American Psychological Association}
}
```

## Contact

Bryan Gibson
Georgia State University
bgibson@student.gsu.edu
GitHub: [@onionbryan](https://github.com/onionbryan)

## License

MIT License - see [LICENSE](LICENSE) for details.
