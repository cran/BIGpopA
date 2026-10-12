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
  format       <- match.arg(format)
  marker_names <- match.arg(marker_names)

  if (!requireNamespace("vcfR", quietly = TRUE))
    stop("Package 'vcfR' is required to read VCF files. Install it with install.packages(\"vcfR\").")
  if (!is.null(ploidy)) .check_ploidy(ploidy)
  if (!is.null(dosage_field) && is.null(ploidy))
    stop("ploidy must be supplied when dosage_field is used; it cannot be inferred from a dosage field.")
  for (arg in c("min_marker_call_rate", "min_sample_call_rate")) {
    val <- get(arg)
    if (length(val) != 1 || !is.numeric(val) || is.na(val) ||
        val < 0 || val > 1)
      stop(arg, " must be a single number between 0 and 1.")
  }

  #### Read the VCF ####
  if (is.character(vcf) && length(vcf) == 1) {
    if (!file.exists(vcf))
      stop("VCF file not found: ", vcf)
    vcf_obj <- vcfR::read.vcfR(vcf, verbose = verbose, ...)
  } else if (inherits(vcf, "vcfR")) {
    vcf_obj <- vcf
  } else {
    stop("vcf must be a path to a .vcf / .vcf.gz file or a vcfR object.")
  }

  # getFIX() drops to a named vector when the VCF has one variant; keep a matrix
  get_fix <- function(obj) {
    f <- vcfR::getFIX(obj)
    if (!is.null(f) && is.null(dim(f)))
      f <- matrix(f, nrow = 1, dimnames = list(NULL, names(f)))
    f
  }

  fix <- get_fix(vcf_obj)
  if (is.null(fix) || nrow(fix) == 0)
    stop("The VCF contains no variants.")

  #### Drop multiallelic markers ####
  n_markers_in <- nrow(fix)
  if (biallelic_only) {
    keep <- !grepl(",", fix[, "ALT"], fixed = TRUE) & !is.na(fix[, "ALT"])
    if (sum(keep) == 0)
      stop("No biallelic markers remain. Set biallelic_only = FALSE to keep multiallelic markers.")
    if (any(!keep)) {
      warning(sum(!keep), " multiallelic marker(s) were dropped. ",
              "Set biallelic_only = FALSE to keep them.", call. = FALSE)
      vcf_obj <- vcf_obj[keep, ]
      fix     <- get_fix(vcf_obj)
    }
  }

  #### Marker names ####
  ids       <- fix[, "ID"]
  positions <- paste(fix[, "CHROM"], fix[, "POS"], sep = "_")
  id_usable <- !any(is.na(ids) | ids == "." | ids == "") &&
    !any(duplicated(ids))

  if (marker_names == "id") {
    if (!id_usable)
      stop("The VCF ID column is missing, empty or duplicated; use marker_names = \"position\" or \"auto\".")
    snp_names <- ids
  } else if (marker_names == "position") {
    snp_names <- positions
  } else {
    snp_names <- if (id_usable) ids else positions
  }
  if (any(duplicated(snp_names))) {
    warning("Duplicated marker names were made unique with a numeric suffix.", call. = FALSE)
    snp_names <- make.unique(snp_names, sep = "_")
  }

  #### Extract genotype calls ####
  if (!is.null(dosage_field)) {

    # Dosages are already numeric in the requested FORMAT field
    dose <- vcfR::extract.gt(vcf_obj, element = dosage_field,
                             as.numeric = TRUE, IDtoRowNames = FALSE)
    if (all(is.na(dose)))
      stop("FORMAT field '", dosage_field, "' is absent or empty in this VCF.")
    obs_max <- suppressWarnings(max(dose, na.rm = TRUE))
    if (is.finite(obs_max) && obs_max > ploidy)
      warning("Dosages above the requested ploidy (", ploidy,
              ") were found in field '", dosage_field, "'.", call. = FALSE)
    n_wrong_ploidy <- 0L

  } else {

    gt <- vcfR::extract.gt(vcf_obj, element = "GT", IDtoRowNames = FALSE)
    if (is.null(gt) || all(is.na(gt)))
      stop("No GT field found in this VCF.")

    # Genotype strings repeat heavily, so the parsing is done once per unique
    #   call and mapped back to the full matrix. This keeps large files cheap.
    calls <- unique(as.vector(gt))
    idx   <- match(gt, calls)
    n_obs <- tabulate(idx, nbins = length(calls))

    alleles <- strsplit(calls, "[/|]")                   # phased or unphased
    n_slots <- lengths(alleles)
    is_miss <- is.na(calls) | n_slots == 0 |
      vapply(alleles,
                   function(a) any(is.na(a) | a == "."),
                   logical(1))

    # dosage = allele slots that are not the reference allele (0)
    dose_u <- rep(NA_real_, length(calls))
    dose_u[!is_miss] <- vapply(alleles[!is_miss],
                                     function(a) sum(a != "0"),
                                     numeric(1))

    #### Ploidy ####
    obs_slots <- n_slots[!is_miss]
    obs_count <- n_obs[!is_miss]
    if (length(obs_slots) == 0)
      stop("Every genotype call is missing; nothing to convert.")
    slot_totals <- tapply(obs_count, obs_slots, sum)
    obs_ploidy  <- as.integer(names(slot_totals)[which.max(slot_totals)])

    ploidy_mismatch <- FALSE
    if (is.null(ploidy)) {
      ploidy <- obs_ploidy
      if (verbose)
        cat(sprintf("Inferred ploidy from the GT field: %s\n", ploidy))
    } else {
      ploidy_mismatch <- obs_ploidy != ploidy
    }

    off_ploidy     <- !is_miss & n_slots != ploidy
    n_wrong_ploidy <- sum(n_obs[off_ploidy])

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
    dose <- matrix(dose_u[idx],
                         nrow = nrow(gt), ncol = ncol(gt),
                         dimnames = list(NULL, colnames(gt)))
    storage.mode(dose) <- "integer"   # dosages are whole numbers
  }

  #### Reference instead of alternate dosage ####
  ploidy <- as.integer(ploidy)
  if (ref_dosage) dose <- ploidy - dose

  #### Individuals in rows, markers in columns ####
  rownames(dose) <- snp_names
  geno <- t(dose)
  if (is.null(rownames(geno)) || any(rownames(geno) == ""))
    stop("Sample names could not be read from the VCF.")

  #### Call-rate filtering ####
  n_samples_in <- nrow(geno)
  n_snps_in    <- ncol(geno)

  if (min_marker_call_rate > 0) {
    marker_cr <- colMeans(!is.na(geno))
    geno      <- geno[, marker_cr >= min_marker_call_rate, drop = FALSE]
    if (ncol(geno) == 0)
      stop("No markers passed min_marker_call_rate = ", min_marker_call_rate, ".")
  }
  if (min_sample_call_rate > 0) {
    sample_cr <- rowMeans(!is.na(geno))
    geno      <- geno[sample_cr >= min_sample_call_rate, , drop = FALSE]
    if (nrow(geno) == 0)
      stop("No individuals passed min_sample_call_rate = ", min_sample_call_rate, ".")
  }

  #### Verbose report ####
  if (verbose) {
    cat("\n=== VCF Conversion Report ===\n")
    cat(sprintf("Ploidy:                 %s\n", ploidy))
    cat(sprintf("Markers read:           %s\n", n_markers_in))
    cat(sprintf("Markers returned:       %d\n", ncol(geno)))
    cat(sprintf("Individuals returned:   %d of %s\n", nrow(geno), n_samples_in))
    cat(sprintf("Missing calls:          %.2f%%\n", 100 * mean(is.na(geno))))
    if (n_wrong_ploidy > 0)
      cat(sprintf("Calls set to missing for wrong ploidy: %s\n", n_wrong_ploidy))
    cat(sprintf("Coding:                 dosage of the %s allele (0..%s)\n",
                if (ref_dosage) "reference" else "alternate", ploidy))
  }

  #### Return in the requested layout ####
  if (format == "matrix") {
    attr(geno, "ploidy") <- ploidy
    return(geno)
  }

  out <- data.table::as.data.table(geno, keep.rownames = "id")
  data.table::setattr(out, "ploidy", ploidy)   # by reference, keeps the data.table valid
  return(out)
}
