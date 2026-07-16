<!-- badges: start -->
[![Development Status](https://img.shields.io/badge/status-active%20development-yellow)](https://github.com/Breeding-Insight/BIGpopA)
[![R-CMD-check](https://github.com/Breeding-Insight/BIGpopA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Breeding-Insight/BIGpopA/actions/workflows/R-CMD-check.yaml)
[![R](https://img.shields.io/badge/R-%3E%3D%204.4-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)
[![GitHub issues](https://img.shields.io/github/issues/Breeding-Insight/BIGpopA)](https://github.com/Breeding-Insight/BIGpopA/issues)
[![GitHub pull requests](https://img.shields.io/github/issues-pr/Breeding-Insight/BIGpopA)](https://github.com/Breeding-Insight/BIGpopA/pulls)
[![GitHub Release](https://img.shields.io/github/v/release/Breeding-Insight/BIGpopA?include_prereleases)](https://github.com/Breeding-Insight/BIGpopA/releases/latest)
<!-- [![codecov](https://app.codecov.io/gh/Breeding-Insight/BIGpopA/graph/badge.svg?token=PJUZMRN1NF)](https://app.codecov.io/gh/Breeding-Insight/BIGpopA) -->
<!-- badges: end -->

<div align="center">
  <img width="250" height="250" alt="BIGpopa_logo" src="https://github.com/user-attachments/assets/3560cb73-6939-4af6-ab41-045b8a1598c3" />
</div>


# Breeding Insight Genomics population Analyses
### Pedigree Validation and Breed/Line Composition Estimation for Diploid and Polyploid Species
BIGpopA is an R package developed by [Breeding Insight](https://breedinginsight.org/) that provides tools for pedigree quality control and genomic breed/line composition estimation in diploid and polyploid breeding populations. It is designed to help researchers and breeders identify pedigree errors, assign parentage from SNP genotype data, and estimate genome-wide breed or line composition.

### Installation
To install the latest CRAN version of BIGpopA:

```R
install.packages("BIGpopA")
library(BIGpopA)
```
To install the development version of BIGpopA, install from GitHub using `remotes`:
```R
install.packages("remotes")
remotes::install_github("Breeding-Insight/BIGpopA", dependencies = TRUE)
library(BIGpopA)
```
##### Note: BIGpopA is currently in development. Please report any bugs or issues on the GitHub Issues page.

### Funding
BIGpopA development is supported by Breeding Insight, a USDA-funded initiative based at the University of Florida - IFAS.

## Citation
If you use BIGpopA in your research, please cite as:

Chinchilla-Vargas, Josue, and Breeding Insight Team. 2026. "BIGpopA: Pedigree Validation and Breed/Line Composition Estimation for Diploid and Polyploid Species." R package version 1.0.6. https://github.com/Breeding-Insight/BIGpopA.
