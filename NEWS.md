# BIGpopA 2.1.0
- New `vcf_to_dosage()` converts a VCF (`.vcf` / `.vcf.gz`) or `vcfR` object into allele-B dosages for any ploidy, as an `id` + markers table or an individuals x markers matrix.
- New `ped_to_dosage()` converts a PLINK `.ped` file (diploid) into allele dosages. A `.map` file is optional and only supplies marker names. The counted allele is the alphabetically/numerically last allele at each marker, or can be supplied to code a second file like a first one.
- `find_parentage()`, `validate_pedigree()`, `allele_freq_poly()` and `solve_composition_poly()` now accept genotypes as a TXT/TSV/CSV file, a VCF file or `vcfR` object, a PLINK `.ped` file, a data.frame / data.table, or a matrix. VCF input is converted with `vcf_to_dosage()` using the function's `ploidy`; `.ped` input requires `ploidy = 2`.
- When the reference panel is a `.ped` file, `allele_freq_poly()` stores the counted allele per marker and `solve_composition_poly()` uses it to code a validation `.ped` file consistently.
- `allele_freq_poly()` and `solve_composition_poly()` also accept a data.frame with an `id` / `ID` column; existing matrix and row-named data.frame inputs are unchanged.
- Genotype read errors in `find_parentage()` and `validate_pedigree()` now include the underlying error message.
- Column names `id`, `male_parent`, `female_parent` and `sex` are now matched ignoring case, spaces and dots (e.g. `ID`, `Male_Parent`, `FEMALE PARENT`) in `check_ped()`, `validate_pedigree()`, `find_parentage()` and genotype tables.
- Fixed `solve_composition_poly()` failing when `Y` contains a single animal.
- Removed the `ped`, `groups`, `mia`, `sire` and `dam` arguments from `solve_composition_poly()`. They relied on internal helpers that were not carried over from BIGr and always failed. The function now takes `Y`, `X` and `ploidy`.

# BIGpopA 2.0.0
- `find_parentage()` and `validate_pedigree()` now support any ploidy through a new `ploidy` argument (default 2). Genotypes may be coded as allele-B dosage (0, 1, ..., ploidy).
- Even ploidy uses a generalized polysomic gamete-range Mendelian test; it is exact for autopolyploids and conservative for allopolyploids (correct trios are never wrongly flagged).
- Odd ploidy (e.g. triploid), where balanced gametes are undefined, automatically falls back to a model-free opposite-homozygote exclusion evaluated on homozygous-informative markers only.
- Diploid results (`ploidy = 2`) are unchanged from previous versions.

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
