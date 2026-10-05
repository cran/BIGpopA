#' Convert a PLINK .ped File to a BIGpopA Genotype Matrix
#'
#' Reads a PLINK 1 text pedigree file (`.ped`) and returns allele dosages in
#' the layout used by the rest of BIGpopA: individuals in rows, markers in
#' columns, coded as the number of copies of the counted allele (0, 1, 2)
#' with NA for missing calls. PLINK `.ped` files are diploid only.
#'
#' The `.ped` format stores allele letters rather than a reference/alternate
#' pair, so the counted allele (allele B) is chosen per marker. By default it
#' is the alphabetically / numerically last allele observed at the marker, so
#' `1/2` coded files count allele `2` and an `A/G` marker counts `G`. To code a
#' second file the same way as a first one (for example validation genotypes
#' against a reference panel), pass the first file's `counted_allele`
#' attribute through the `counted_allele` argument. Otherwise a marker that is
#' monomorphic in one file can be coded in the opposite direction.
#'
#' Only the individual ID (column 2) is used from the six leading pedigree
#' columns. A `.map` file is optional and only supplies marker names; it is
#' not used for allele coding. Without it, markers are named `SNP1`, `SNP2`,
#' ... in file order.
#'
#' @param ped Path to a PLINK `.ped` file (whitespace separated, two allele
#'   columns per marker).
#' @param map `"auto"` (default) uses a `.map` file with the same name next to
#'   the `.ped` when one exists. A path uses that `.map` file; NULL ignores any
#'   `.map` and names markers by position.
#' @param format Character. `"data.frame"` (default) returns a data.table with
#'   an `id` column followed by one column per marker, the input expected by
#'   [find_parentage()] and [validate_pedigree()]. `"matrix"` returns a numeric
#'   matrix with individuals in named rows and markers in named columns, the
#'   input expected by [allele_freq_poly()] and [solve_composition_poly()].
#' @param counted_allele Named character vector, or NULL (default). Allele to
#'   count at each marker, named by marker. Markers not listed use the default
#'   rule. Typically the `counted_allele` attribute of a previous result.
#' @param missing_codes Character vector of allele codes treated as missing.
#'   Default is `c("0", "N", "-", ".")`. A call is missing when either allele
#'   is missing.
#' @param verbose Logical. If TRUE, prints a conversion report. Default is TRUE.
#'
#' @return A data.table with an `id` column followed by marker columns
#'   (`format = "data.frame"`), or an integer matrix of individuals by markers
#'   (`format = "matrix"`). The attribute `"ploidy"` is 2 and the attribute
#'   `"counted_allele"` holds the allele counted at each returned marker.
#'   Markers with more than two alleles are dropped with a warning.
#'
#' @examples
#' ped_path <- tempfile(fileext = ".ped")
#' writeLines(c(
#'   "F1 P1 0 0 1 -9 A A G G A G",
#'   "F1 P2 0 0 2 -9 A G G G G G",
#'   "F1 O1 P1 P2 1 -9 A G G G A G"
#' ), ped_path)
#'
#' genotypes <- ped_to_dosage(ped_path, map = NULL, verbose = FALSE)
#' print(genotypes)
#' attr(genotypes, "counted_allele")
#'
#' unlink(ped_path)
#'
#' @author Josue Chinchilla-Vargas
#'
#' @seealso [vcf_to_dosage()], [find_parentage()], [allele_freq_poly()]
#'
#' @importFrom data.table fread as.data.table setattr
#' @export
ped_to_dosage <- function(ped,
                          map            = "auto",
                          format         = c("data.frame", "matrix"),
                          counted_allele = NULL,
                          missing_codes  = c("0", "N", "-", "."),
                          verbose        = TRUE) {

  #### Input Validation ####
  format <- base::match.arg(format)
  if (!base::is.character(ped) || base::length(ped) != 1)
    stop("ped must be a path to a PLINK .ped file.")
  if (!base::file.exists(ped))
    stop("PED file not found: ", ped)
  if (!base::is.null(counted_allele) &&
      (!base::is.character(counted_allele) || base::is.null(base::names(counted_allele))))
    stop("counted_allele must be a named character vector (names = marker names).")

  #### Read the .ped ####
  raw <- data.table::fread(ped, header = FALSE, colClasses = "character",
                           data.table = FALSE)
  n_allele_cols <- base::ncol(raw) - 6L
  if (n_allele_cols < 2 || n_allele_cols %% 2 != 0)
    stop("The .ped file must have 6 pedigree columns followed by two allele ",
         "columns per marker (compound genotypes are not supported).")
  n_markers <- n_allele_cols %/% 2L

  ids <- raw[[2]]
  dup <- base::unique(ids[base::duplicated(ids)])
  if (base::length(dup) > 0)
    stop("Duplicated individual IDs (column 2) in the .ped file: ",
         base::paste(dup, collapse = ", "))

  alleles <- base::as.matrix(raw[, -(1:6), drop = FALSE])
  alleles[alleles %in% missing_codes] <- NA
  a1   <- alleles[, base::seq(1, n_allele_cols, by = 2), drop = FALSE]
  a2   <- alleles[, base::seq(2, n_allele_cols, by = 2), drop = FALSE]
  miss <- base::is.na(a1) | base::is.na(a2)
  a1[miss] <- NA
  a2[miss] <- NA

  #### Marker names ####
  map_path <- if (base::identical(map, "auto")) {
    cand <- base::sub("\\.ped$", ".map", ped, ignore.case = TRUE)
    if (cand != ped && base::file.exists(cand)) cand else NULL
  } else {
    map
  }

  if (!base::is.null(map_path)) {
    if (!base::file.exists(map_path))
      stop("MAP file not found: ", map_path)
    map_df <- data.table::fread(map_path, header = FALSE, colClasses = "character",
                                data.table = FALSE)
    if (base::nrow(map_df) != n_markers)
      stop("The .map file lists ", base::nrow(map_df), " markers but the .ped file has ",
           n_markers, ".")
    snp_names <- map_df[[2]]
    no_name   <- base::is.na(snp_names) | snp_names %in% c("", ".")
    snp_names[no_name] <- base::paste(map_df[[1]][no_name],
                                      map_df[[base::ncol(map_df)]][no_name], sep = "_")
  } else {
    snp_names <- base::paste0("SNP", base::seq_len(n_markers))
  }
  if (base::any(base::duplicated(snp_names))) {
    warning("Duplicated marker names were made unique with a numeric suffix.", call. = FALSE)
    snp_names <- base::make.unique(snp_names, sep = "_")
  }

  #### Counted allele per marker ####
  both      <- base::rbind(a1, a2)
  n_alleles <- base::apply(both, 2, function(v) base::length(base::unique(v[!base::is.na(v)])))
  counted   <- base::apply(both, 2, function(v) {
    v <- v[!base::is.na(v)]
    if (base::length(v) == 0) NA_character_ else base::max(v)
  })
  base::names(counted) <- snp_names

  if (!base::is.null(counted_allele)) {
    use <- snp_names[snp_names %in% base::names(counted_allele)]
    counted[use] <- counted_allele[use]
  }

  #### Drop multiallelic markers ####
  keep <- n_alleles <= 2
  if (base::any(!keep))
    warning(base::sum(!keep), " marker(s) with more than two alleles were dropped.",
            call. = FALSE)
  if (base::sum(keep) == 0)
    stop("No biallelic markers remain in the .ped file.")

  #### Dosage = copies of the counted allele ####
  cnt  <- base::matrix(counted, nrow = base::nrow(a1), ncol = n_markers, byrow = TRUE)
  geno <- (a1 == cnt) + (a2 == cnt)
  base::storage.mode(geno) <- "integer"
  base::dimnames(geno) <- base::list(ids, snp_names)
  geno    <- geno[, keep, drop = FALSE]
  counted <- counted[keep]

  #### Verbose report ####
  if (verbose) {
    base::cat("\n=== PED Conversion Report ===\n")
    base::cat("Ploidy:                 2 (PLINK .ped)\n")
    base::cat("Marker names from:     ", if (base::is.null(map_path)) "position (no .map)" else map_path, "\n")
    base::cat("Markers read:          ", n_markers, "\n")
    base::cat("Markers returned:      ", base::ncol(geno), "\n")
    base::cat("Individuals returned:  ", base::nrow(geno), "\n")
    base::cat(base::sprintf("Missing calls:          %.2f%%\n",
                            100 * base::mean(base::is.na(geno))))
  }

  #### Return in the requested layout ####
  if (format == "matrix") {
    base::attr(geno, "ploidy")         <- 2L
    base::attr(geno, "counted_allele") <- counted
    return(geno)
  }

  out <- data.table::as.data.table(geno, keep.rownames = "id")
  data.table::setattr(out, "ploidy", 2L)
  data.table::setattr(out, "counted_allele", counted)
  return(out)
}
