<!-- badges: start -->
[![Development Status](https://img.shields.io/badge/status-active%20development-yellow)](https://github.com/Breeding-Insight/BIGpopA)
[![R-CMD-check](https://github.com/Breeding-Insight/BIGpopA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Breeding-Insight/BIGpopA/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/BIGpopA)](https://CRAN.R-project.org/package=BIGpopA)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/BIGpopA)](https://cran.r-project.org/package=BIGpopA)
[![CRAN monthly downloads](https://cranlogs.r-pkg.org/badges/BIGpopA)](https://cran.r-project.org/package=BIGpopA)
[![R](https://img.shields.io/badge/R-%3E%3D%204.4-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)
[![GitHub issues](https://img.shields.io/github/issues/Breeding-Insight/BIGpopA)](https://github.com/Breeding-Insight/BIGpopA/issues)
[![GitHub pull requests](https://img.shields.io/github/issues-pr/Breeding-Insight/BIGpopA)](https://github.com/Breeding-Insight/BIGpopA/pulls)
[![GitHub Release](https://img.shields.io/github/v/release/Breeding-Insight/BIGpopA?include_prereleases)](https://github.com/Breeding-Insight/BIGpopA/releases/latest)
<!-- [![codecov](https://app.codecov.io/gh/Breeding-Insight/BIGpopA/graph/badge.svg?token=PJUZMRN1NF)](https://app.codecov.io/gh/Breeding-Insight/BIGpopA) -->
<!-- badges: end -->

<div align="center">
  <img width="250" height="250" alt="BIGpopa_logo" src="man/figures/BIGpopa_logo.png" />
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

### Main functions

| Function | Purpose |
|---|---|
| `check_ped()` | Detect and correct pedigree errors (duplicates, conflicting trios, missing parents, cycles, inconsistent sex roles) |
| `validate_pedigree()` | Validate parent-offspring trios against SNP genotypes using Mendelian error rates; also returns per-marker mismatch counts (`marker_summary`) |
| `find_parentage()` | Assign the most likely parent(s) to progeny from candidate parents; `top_n` reports several candidates, `tied_candidates` flags ambiguity, and `method = "fill_pedigree"` fills missing parents of a pedigree |
| `allele_freq_poly()` | Compute reference population allele frequencies |
| `solve_composition_poly()` | Estimate genome-wide breed/line composition by quadratic programming |
| `assign_composition_line()` | Percentages per line and predicted line, with an assignment threshold (default 50%) |
| `vcf_to_dosage()` | Convert a VCF (`.vcf` / `.vcf.gz`) to allele dosages for any ploidy |
| `ped_to_dosage()` | Convert a PLINK `.ped` (with optional `.map`) to allele dosages |

Pedigree validation and parentage assignment support any ploidy, using a polysomic Mendelian test for even ploidy and a homozygosity-based check for odd ploidy.

### Genotype input formats

`validate_pedigree()`, `find_parentage()`, `allele_freq_poly()` and `solve_composition_poly()` accept genotypes in any of these formats:

| Format | Notes |
|---|---|
| Text file (`.txt`, `.tsv`, `.csv`) | ID column followed by marker columns coded as allele-B dosage (0, 1, ..., ploidy) |
| VCF (`.vcf`, `.vcf.gz`) or `vcfR` object | `GT` calls converted to ALT-allele dosage using the function's `ploidy` |
| PLINK `.ped` (+ optional `.map`) | Diploid only; the `.map` supplies marker names |
| `data.frame` / `data.table` / `matrix` | Already-loaded dosage data |

```R
# Same call, different input formats
find_parentage("genotypes.vcf.gz", "parents.txt", "progeny.txt", ploidy = 4)
find_parentage("genotypes.ped",    "parents.txt", "progeny.txt")

# Breed/line composition from a reference and a validation VCF
freq <- allele_freq_poly("reference.vcf", populations, ploidy = 2)
comp <- solve_composition_poly("validation.vcf", freq, ploidy = 2)
```

### Shiny app
BIGpopA powers the pedigree and composition modules of [Familia](https://github.com/Breeding-Insight/familia), a point-and-click interface for the same analyses.

### Funding
BIGpopA development is supported by Breeding Insight, a USDA-funded initiative based at the University of Florida - IFAS.

## Citation
If you use BIGpopA in your research, please cite as:

Chinchilla-Vargas, Josue, and Breeding Insight Team. 2026. "BIGpopA: Pedigree Validation and Breed/Line Composition Estimation for Diploid and Polyploid Species." R package version 2.1.0. https://github.com/Breeding-Insight/BIGpopA.
