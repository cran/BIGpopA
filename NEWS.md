# BIGpopA 1.0.6
- Added a "BIGpopA Tutorial" vignette demonstrating the pedigree cleaning, validation, parentage assignment, and breed/line composition workflow

# BIGpopA 1.0.5
- Updated package title and example on function descriptions to pass CRAN review

# BIGpopA 1.0.4
- Removed LICENSE file to pass CRAN automated tests

# BIGpopA 1.0.3
- Fixed typo on repo name and patched find_parentage
- Removed /dev from repository
- Added words to WORDLIST


# BIGpopA 1.0.2
- Renamed package to BIGpopA and made repo public


# popR 1.0.2
- Initial release of `BIGpopA` as a standalone package.
- `BIGpopA` contains pedigree validation and breed/line composition functions
  previously found in [BIGr](https://github.com/Breeding-Insight/BIGr),
  where they will no longer be maintained going forward.
- These functions are the backbone of the pedigree and composition modules
  in the [familia](https://github.com/Breeding-Insight/familia) Shiny app.

## Functions included

- `check_ped()` — checks and corrects common pedigree errors (duplicate rows,
  conflicting trios, missing parents, cycles, inconsistent sex roles)
- `find_parentage()` — assigns most likely parent(s) to progeny using
  Mendelian error rates or homozygous mismatch rates
- `validate_pedigree()` — validates parent-offspring trios against SNP
  genotype data and outputs a corrected pedigree
- `allele_freq_poly()` — computes allele frequencies for diploid and
  polyploid reference populations
- `solve_composition_poly()` — estimates genome-wide breed/line composition
  using quadratic programming
