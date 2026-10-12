# BIGpopA 2.2.0
- `validate_pedigree()` and `find_parentage()`: `founders_file` now also accepts a character vector of ids or a data.frame (first column), in addition to a file path. All inputs can now be supplied from memory.
- `validate_pedigree()`: removed the suggested replacement parents (`best_male_candidate`, `best_male_candidate_error_pct`, `best_female_candidate`, `best_female_candidate_error_pct`) from `full_results` and the per-status tables. Use `find_parentage(method = "fill_pedigree")` on `corrected_pedigree` to fill in missing parents.
- `validate_pedigree()`: new `marker_summary` element in the returned list gives, for every marker, the number of trios tested, the number with a Mendelian mismatch and the mismatch percentage, sorted worst first. It covers all trios with both parents known and genotyped, and helps spot markers (e.g. paralogous or miscalled) that cause many failures.
- `find_parentage()`: new `method = "fill_pedigree"` fills in missing parents of a pedigree. `progeny_file` must have `id`, `male_parent` and `female_parent` columns (e.g. `validate_pedigree()$corrected_pedigree`); parents coded `0`/`NA`/empty are unknown. A known male searches only females, a known female only males, and rows with both unknown get the full pair search (skip founders with the new `founders_file`). Rows with both parents known are skipped. Results gain `search_mode` and `known_parent_error_pct`, and a known parent without genotypes gets `status = "known_parent_no_genotype"`.
- `find_parentage()`: new `top_n` argument reports the best `top_n` candidates per progeny, ranked by Mendelian error (ties broken by markers tested). Rank 1 stays in the base columns and ranks 2 to `top_n` are added as suffix columns (`male_parent_2`, `female_parent_2`, `mendelian_error_pct_2`, `markers_tested_2` for `best_pair`; `best_match_2` and the same error/marker columns for the other methods). The default `top_n = 1` keeps the previous output, including `show_ties`.
- `find_parentage()`: new `tied_candidates` column gives the number of candidates (pairs for `best_pair`) that share each progeny's lowest error. Values above 1 mean the best assignment is ambiguous, including for the single-parent methods, which previously did not report ties.
- `validate_pedigree()`: fixed `corrected_pedigree` for `low_markers_remove_female_parent`, which also set the male parent to `0`. Low-marker trios now remove only the identified parent(s).
- New `assign_composition_line()` turns the output of `solve_composition_poly()` into the composition table used by Familia: percentages per reference population and a `Predicted line`, labelled `Undetermined` when the highest proportion is below `assign_threshold` (default 50). `solve_composition_poly()` is unchanged.
- `validate_pedigree()`: new `marker_trio_table = FALSE` argument. When TRUE the result includes `marker_trio_table`, a marker-by-trio table with one column per trio (`id_male_parent_female_parent`) giving each marker's status: `match`, `mismatch_male`, `mismatch_female`, `mismatch_both`, `mismatch_combination` (each parent compatible alone but not the pair), `untestable` (odd ploidy) or `missing_*` named after the individuals without a call (e.g. `missing_male`, `missing_male_female`, `missing_progeny`).

# BIGpopA 2.1.1
- `validate_pedigree()`: a trio that fails the Mendelian error threshold while neither parent fails the single-parent homozygous check now gets `recommended_correction = "unresolved"` instead of `"keep_both"`. This happens, for example, with a recorded self whose offspring carries alleles the parent lacks: the pedigree is wrong, but the faulty parent cannot be identified. These trios are left unchanged in `corrected_pedigree` for manual review.
- `validate_pedigree()` plot: bars are coloured by trio error first, so no trio above the threshold is shown as Pass. Unresolved trios are shown in grey ("Fail - Parents Unresolved"), and the legend wraps onto two rows so all labels are visible.

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
