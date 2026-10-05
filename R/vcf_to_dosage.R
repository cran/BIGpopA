#' Convert a VCF File to a BIGpopA Genotype Matrix
#'
#' Reads a variant call format (VCF) file, or an existing `vcfR` object, and
#' returns allele dosages in the layout used by the rest of BIGpopA:
#' individuals in rows, markers in columns, coded as the dosage of the
#' alternate (B) allele (0, 1, ..., ploidy) with NA for missing calls.
#'
#' The function is ploidy agnostic. Dosage is obtained by counting the
#' non-reference alleles in each `GT` call, so a diploid `0/1` returns 1 and a
#' tetraploid `0/1/1/1` returns 3, with no assumption about the ploidy level.
#' When `ploidy` is not supplied it is inferred from the number of allele slots
#' in the calls themselves (the most common slot count across the file). Calls
#' whose slot count disagrees with the file ploidy, for instance a diploid call
#' in an otherwise tetraploid file, are set to missing and reported, since the
#' downstream functions require a single ploidy level.
#'
#' Phased (`|`) and unphased (`/`) separators are both accepted. A call is
#' treated as missing when any of its allele slots is `.`.
#'
#' @param vcf Path to a VCF file (`.vcf` or `.vcf.gz`), OR a `vcfR` object.
#' @param ploidy Integer >= 2, or NULL (default) to infer the ploidy from the
#'   `GT` field. When supplied, it is checked against the ploidy observed in the
#'   file and a warning is issued if they disagree.
#' @param format Character. `"data.frame"` (default) returns a data.table with
#'   an `id` column followed by one column per marker, the input expected by
#'   [find_parentage()] and [validate_pedigree()]. `"matrix"` returns a numeric
#'   matrix with individuals in named rows and markers in named columns, the
#'   input expected by [allele_freq_poly()] and [solve_composition_poly()].
#' @param dosage_field Character or NULL (default). Name of a FORMAT field
#'   holding dosages directly (for example `"DS"` or `"UD"`). When supplied,
#'   that field is read instead of counting alleles in `GT`, and `ploidy` must
#'   be given (it cannot be inferred from a dosage field).
#' @param ref_dosage Logical. If TRUE, count reference alleles instead of
#'   alternate alleles (`ploidy` minus the alternate dosage). Default is FALSE,
#'   which matches the "dosage of allele B" coding used throughout BIGpopA.
#' @param biallelic_only Logical. If TRUE (default), markers with more than one
#'   ALT allele are dropped. If FALSE, every non-reference allele is counted
#'   towards the dosage, so a multiallelic marker is collapsed to
#'   reference vs non-reference.
#' @param marker_names Character. `"auto"` (default) uses the VCF `ID` column
#'   when it is complete and unique, and falls back to `CHROM_POS` otherwise.
#'   `"id"` forces the `ID` column and `"position"` forces `CHROM_POS`.
#' @param min_marker_call_rate Numeric between 0 and 1. Markers called in fewer
#'   than this proportion of individuals are dropped. Default is 0 (no
#'   filtering).
#' @param min_sample_call_rate Numeric between 0 and 1. Individuals called at
#'   fewer than this proportion of the retained markers are dropped. Applied
#'   after the marker filter. Default is 0 (no filtering).
#' @param verbose Logical. If TRUE, prints a conversion report. Default is TRUE.
#' @param ... Further arguments passed to [vcfR::read.vcfR()], for example
#'   `nrows` to read only the first records of a large file.
#'
#' @return A data.table with an `id` column followed by marker columns
#'   (`format = "data.frame"`), or a numeric matrix of individuals by markers
#'   (`format = "matrix"`). The attribute `"ploidy"` records the ploidy used.
#'
#' @details
#' `vcfR` is used to read the file and is listed under Suggests, so it must be
#' installed before calling this function
#' (`install.packages("vcfR")`).
#'
#' @examples
#' \donttest{
#' if (requireNamespace("vcfR", quietly = TRUE)) {
#'   # Write a small tetraploid VCF (3 markers x 4 samples) to a temporary file
#'   vcf_path <- tempfile(fileext = ".vcf")
#'   writeLines(c(
#'     "##fileformat=VCFv4.3",
#'     "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
#'     paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
#'           "FORMAT", "S1", "S2", "S3", "S4", sep = "\t"),
#'     paste("chr1", "100", "snp1", "A", "T", ".", ".", ".", "GT",
#'           "0/0/0/0", "0/0/0/1", "0/1/1/1", "1/1/1/1", sep = "\t"),
#'     paste("chr1", "200", "snp2", "C", "G", ".", ".", ".", "GT",
#'           "0/0/1/1", "./././.", "0/0/0/0", "0/1/1/1", sep = "\t"),
#'     paste("chr1", "300", "snp3", "G", "A", ".", ".", ".", "GT",
#'           "1/1/1/1", "1/1/1/1", "0/0/1/1", "0/0/0/1", sep = "\t")
#'   ), vcf_path)
#'
#'   # tetraploid VCF -> id + marker columns for find_parentage()
#'   genotypes <- vcf_to_dosage(vcf_path, ploidy = 4)
#'   print(genotypes)
#'
#'   # same data as an individuals x markers matrix for allele_freq_poly()
#'   geno_mat <- vcf_to_dosage(vcf_path, ploidy = 4, format = "matrix")
#'
#'   unlink(vcf_path)
#' }
#' }
#'
#' @author Josue Chinchilla-Vargas
#'
#' @seealso [find_parentage()], [validate_pedigree()], [allele_freq_poly()]
#'
#' @importFrom data.table as.data.table setattr
#' @export
vcf_to_dosage <- function(vcf,
                          ploidy               = NULL,
                          format               = c("data.frame", "matrix"),
                          dosage_field         = NULL,
                          ref_dosage           = FALSE,
                          biallelic_only       = TRUE,
                          marker_names         = c("auto", "id", "position"),
                          min_marker_call_rate = 0,
                          min_sample_call_rate = 0,
                          verbose              = TRUE,
                          ...) {

  #### Input Validation ####
  format       <- base::match.arg(format)
  marker_names <- base::match.arg(marker_names)

  if (!base::requireNamespace("vcfR", quietly = TRUE))
    stop("Package 'vcfR' is required to read VCF files. Install it with install.packages(\"vcfR\").")
  if (!base::is.null(ploidy)) .check_ploidy(ploidy)
  if (!base::is.null(dosage_field) && base::is.null(ploidy))
    stop("ploidy must be supplied when dosage_field is used; it cannot be inferred from a dosage field.")
  for (arg in c("min_marker_call_rate", "min_sample_call_rate")) {
    val <- base::get(arg)
    if (base::length(val) != 1 || !base::is.numeric(val) || base::is.na(val) ||
        val < 0 || val > 1)
      stop(arg, " must be a single number between 0 and 1.")
  }

  #### Read the VCF ####
  if (base::is.character(vcf) && base::length(vcf) == 1) {
    if (!base::file.exists(vcf))
      stop("VCF file not found: ", vcf)
    vcf_obj <- vcfR::read.vcfR(vcf, verbose = verbose, ...)
  } else if (base::inherits(vcf, "vcfR")) {
    vcf_obj <- vcf
  } else {
    stop("vcf must be a path to a .vcf / .vcf.gz file or a vcfR object.")
  }

  # getFIX() drops to a named vector when the VCF has one variant; keep a matrix
  get_fix <- function(obj) {
    f <- vcfR::getFIX(obj)
    if (!base::is.null(f) && base::is.null(base::dim(f)))
      f <- base::matrix(f, nrow = 1, dimnames = base::list(NULL, base::names(f)))
    f
  }

  fix <- get_fix(vcf_obj)
  if (base::is.null(fix) || base::nrow(fix) == 0)
    stop("The VCF contains no variants.")

  #### Drop multiallelic markers ####
  n_markers_in <- base::nrow(fix)
  if (biallelic_only) {
    keep <- !base::grepl(",", fix[, "ALT"], fixed = TRUE) & !base::is.na(fix[, "ALT"])
    if (base::sum(keep) == 0)
      stop("No biallelic markers remain. Set biallelic_only = FALSE to keep multiallelic markers.")
    if (base::any(!keep)) {
      warning(base::sum(!keep), " multiallelic marker(s) were dropped. ",
              "Set biallelic_only = FALSE to keep them.", call. = FALSE)
      vcf_obj <- vcf_obj[keep, ]
      fix     <- get_fix(vcf_obj)
    }
  }

  #### Marker names ####
  ids       <- fix[, "ID"]
  positions <- base::paste(fix[, "CHROM"], fix[, "POS"], sep = "_")
  id_usable <- !base::any(base::is.na(ids) | ids == "." | ids == "") &&
    !base::any(base::duplicated(ids))

  if (marker_names == "id") {
    if (!id_usable)
      stop("The VCF ID column is missing, empty or duplicated; use marker_names = \"position\" or \"auto\".")
    snp_names <- ids
  } else if (marker_names == "position") {
    snp_names <- positions
  } else {
    snp_names <- if (id_usable) ids else positions
  }
  if (base::any(base::duplicated(snp_names))) {
    warning("Duplicated marker names were made unique with a numeric suffix.", call. = FALSE)
    snp_names <- base::make.unique(snp_names, sep = "_")
  }

  #### Extract genotype calls ####
  if (!base::is.null(dosage_field)) {

    # Dosages are already numeric in the requested FORMAT field
    dose <- vcfR::extract.gt(vcf_obj, element = dosage_field,
                             as.numeric = TRUE, IDtoRowNames = FALSE)
    if (base::all(base::is.na(dose)))
      stop("FORMAT field '", dosage_field, "' is absent or empty in this VCF.")
    obs_max <- base::suppressWarnings(base::max(dose, na.rm = TRUE))
    if (base::is.finite(obs_max) && obs_max > ploidy)
      warning("Dosages above the requested ploidy (", ploidy,
              ") were found in field '", dosage_field, "'.", call. = FALSE)
    n_wrong_ploidy <- 0L

  } else {

    gt <- vcfR::extract.gt(vcf_obj, element = "GT", IDtoRowNames = FALSE)
    if (base::is.null(gt) || base::all(base::is.na(gt)))
      stop("No GT field found in this VCF.")

    # Genotype strings repeat heavily, so the parsing is done once per unique
    #   call and mapped back to the full matrix. This keeps large files cheap.
    calls <- base::unique(base::as.vector(gt))
    idx   <- base::match(gt, calls)
    n_obs <- base::tabulate(idx, nbins = base::length(calls))

    alleles <- base::strsplit(calls, "[/|]")                   # phased or unphased
    n_slots <- base::lengths(alleles)
    is_miss <- base::is.na(calls) | n_slots == 0 |
      base::vapply(alleles,
                   function(a) base::any(base::is.na(a) | a == "."),
                   base::logical(1))

    # dosage = allele slots that are not the reference allele (0)
    dose_u <- base::rep(NA_real_, base::length(calls))
    dose_u[!is_miss] <- base::vapply(alleles[!is_miss],
                                     function(a) base::sum(a != "0"),
                                     base::numeric(1))

    #### Ploidy ####
    obs_slots <- n_slots[!is_miss]
    obs_count <- n_obs[!is_miss]
    if (base::length(obs_slots) == 0)
      stop("Every genotype call is missing; nothing to convert.")
    slot_totals <- base::tapply(obs_count, obs_slots, base::sum)
    obs_ploidy  <- base::as.integer(base::names(slot_totals)[base::which.max(slot_totals)])

    ploidy_mismatch <- FALSE
    if (base::is.null(ploidy)) {
      ploidy <- obs_ploidy
      if (verbose)
        base::cat("Inferred ploidy from the GT field:", ploidy, "\n")
    } else {
      ploidy_mismatch <- obs_ploidy != ploidy
    }

    off_ploidy     <- !is_miss & n_slots != ploidy
    n_wrong_ploidy <- base::sum(n_obs[off_ploidy])

    # One warning per situation: a contradicting ploidy, or a few stray calls
    if (ploidy_mismatch) {
      warning("Requested ploidy (", ploidy, ") does not match the ploidy most ",
              "common in the VCF (", obs_ploidy, "). ", n_wrong_ploidy,
              " call(s) that do not match the requested ploidy were set to missing.",
              call. = FALSE)
    } else if (n_wrong_ploidy > 0) {
      warning(n_wrong_ploidy, " call(s) with a ploidy other than ", ploidy,
              " were set to missing.", call. = FALSE)
    }

    dose_u[is_miss | off_ploidy] <- NA_real_
    dose <- base::matrix(dose_u[idx],
                         nrow = base::nrow(gt), ncol = base::ncol(gt),
                         dimnames = base::list(NULL, base::colnames(gt)))
    base::storage.mode(dose) <- "integer"   # dosages are whole numbers
  }

  #### Reference instead of alternate dosage ####
  ploidy <- base::as.integer(ploidy)
  if (ref_dosage) dose <- ploidy - dose

  #### Individuals in rows, markers in columns ####
  base::rownames(dose) <- snp_names
  geno <- base::t(dose)
  if (base::is.null(base::rownames(geno)) || base::any(base::rownames(geno) == ""))
    stop("Sample names could not be read from the VCF.")

  #### Call-rate filtering ####
  n_samples_in <- base::nrow(geno)
  n_snps_in    <- base::ncol(geno)

  if (min_marker_call_rate > 0) {
    marker_cr <- base::colMeans(!base::is.na(geno))
    geno      <- geno[, marker_cr >= min_marker_call_rate, drop = FALSE]
    if (base::ncol(geno) == 0)
      stop("No markers passed min_marker_call_rate = ", min_marker_call_rate, ".")
  }
  if (min_sample_call_rate > 0) {
    sample_cr <- base::rowMeans(!base::is.na(geno))
    geno      <- geno[sample_cr >= min_sample_call_rate, , drop = FALSE]
    if (base::nrow(geno) == 0)
      stop("No individuals passed min_sample_call_rate = ", min_sample_call_rate, ".")
  }

  #### Verbose report ####
  if (verbose) {
    base::cat("\n=== VCF Conversion Report ===\n")
    base::cat("Ploidy:                ", ploidy, "\n")
    base::cat("Markers read:          ", n_markers_in, "\n")
    base::cat("Markers returned:      ", base::ncol(geno), "\n")
    base::cat("Individuals returned:  ", base::nrow(geno), " of ", n_samples_in, "\n", sep = "")
    base::cat(base::sprintf("Missing calls:          %.2f%%\n",
                            100 * base::mean(base::is.na(geno))))
    if (n_wrong_ploidy > 0)
      base::cat("Calls set to missing for wrong ploidy:", n_wrong_ploidy, "\n")
    base::cat("Coding:                 dosage of the ",
              if (ref_dosage) "reference" else "alternate", " allele (0..", ploidy, ")\n",
              sep = "")
  }

  #### Return in the requested layout ####
  if (format == "matrix") {
    base::attr(geno, "ploidy") <- ploidy
    return(geno)
  }

  out <- data.table::as.data.table(geno, keep.rownames = "id")
  data.table::setattr(out, "ploidy", ploidy)   # by reference, keeps the data.table valid
  return(out)
}
