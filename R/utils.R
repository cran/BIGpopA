utils::globalVariables(c(
  # data.table internals
  ":=", ".SD",
  
  # find_parentage.R
  "id", "sex", "male_parent", "female_parent",
  "mendelian_error_pct", "plot_status", "status",
  
  # validate_pedigree.R
  "trio_mendelian_error_pct", "recommended_correction"
))

#### Ploidy-general Mendelian consistency helpers ####

#' Gamete allele-dosage bounds
#'
#' Lower/upper bound on the number of B alleles a balanced gamete can carry
#' for a parent of dosage \code{g} at an even \code{ploidy}. Assumes polysomic
#' (autopolyploid) inheritance: random chromosome segregation, no double
#' reduction. For allopolyploids (disomic inheritance) these bounds are a
#' conservative superset, so correct trios are never wrongly flagged, but some
#' true errors may go undetected.
#'
#' @param g numeric vector or matrix of parental dosages (0..ploidy).
#' @param ploidy even integer ploidy level.
#' @return Object matching \code{g} holding the gamete dosage bound.
#' @noRd
.gamete_lo <- function(g, ploidy) base::pmax(0L, g - ploidy / 2L)
.gamete_hi <- function(g, ploidy) base::pmin(ploidy / 2L, g)

#' Flag impossible offspring dosages
#'
#' Returns TRUE where an offspring dosage cannot arise from the two parental
#' dosages under polysomic inheritance. Operates elementwise on vectors or
#' matrices and reduces exactly to the diploid 0/1/2 rules when ploidy = 2.
#'
#' @param male,female,offspring numeric vectors/matrices of dosages, aligned.
#' @param ploidy even integer ploidy level.
#' @return Logical object matching the inputs (TRUE = Mendelian error).
#' @noRd
mendelian_error <- function(male, female, offspring, ploidy) {
  lo <- .gamete_lo(male, ploidy) + .gamete_lo(female, ploidy)
  hi <- .gamete_hi(male, ploidy) + .gamete_hi(female, ploidy)
  (offspring < lo) | (offspring > hi)
}

#' Per-marker Mendelian mismatch indicator
#'
#' Dispatches on ploidy parity. Even ploidy uses the polysomic gamete-range
#' test (\code{mendelian_error}), which draws on all co-genotyped markers and
#' both parents jointly. Odd ploidy (e.g. triploid), where balanced gametes are
#' undefined, falls back to a model-free opposite-homozygote exclusion: a marker
#' is a mismatch only when the offspring is homozygous and one parent is the
#' opposite homozygote. Reduces exactly to the even-ploidy test when ploidy is
#' even.
#'
#' @param male,female,offspring dosage vectors/matrices (0..ploidy), aligned.
#' @param ploidy integer ploidy level.
#' @return Logical object matching the inputs (TRUE = mismatch).
#' @noRd
.mend_mismatch <- function(male, female, offspring, ploidy) {
  if (ploidy %% 2 == 0)
    return(mendelian_error(male, female, offspring, ploidy))
  o_hom <- offspring == 0 | offspring == ploidy
  ((male   == 0 | male   == ploidy) & male   != offspring & o_hom) |
    ((female == 0 | female == ploidy) & female != offspring & o_hom)
}

#' Per-marker testability indicator
#'
#' Companion to \code{.mend_mismatch} giving the markers that can return a
#' verdict. Even ploidy counts every co-genotyped marker; odd ploidy counts only
#' markers where the offspring is homozygous and at least one parent is
#' homozygous (the homozygosity-informative set).
#'
#' @param male,female,offspring dosage vectors/matrices (0..ploidy), aligned.
#' @param ploidy integer ploidy level.
#' @return Logical object matching the inputs (TRUE = testable).
#' @noRd
.mend_testable <- function(male, female, offspring, ploidy) {
  if (ploidy %% 2 == 0)
    return(!base::is.na(male) & !base::is.na(female) & !base::is.na(offspring))
  o_hom <- !base::is.na(offspring) & (offspring == 0 | offspring == ploidy)
  m_hom <- !base::is.na(male)   & (male   == 0 | male   == ploidy)
  f_hom <- !base::is.na(female) & (female == 0 | female == ploidy)
  o_hom & (m_hom | f_hom)
}

#' Validate a ploidy argument
#'
#' Stops if \code{ploidy} is not an integer >= 2.
#'
#' @param ploidy value supplied by the user.
#' @return Invisibly TRUE if valid; otherwise an error is thrown.
#' @noRd
.check_ploidy <- function(ploidy) {
  if (base::length(ploidy) != 1 || !base::is.numeric(ploidy) ||
      base::is.na(ploidy) || ploidy < 2 || ploidy != base::round(ploidy))
    base::stop("ploidy must be an integer >= 2.")
  base::invisible(TRUE)
}

#### Column name helpers ####

#' Normalize a column name for matching
#'
#' Lower-cases and trims a name and turns spaces, dots and dashes into
#' underscores, so "ID", "Male Parent" and "female.parent" match "id",
#' "male_parent" and "female_parent".
#'
#' @param x character vector of column names.
#' @return Normalized character vector.
#' @noRd
.normalize_name <- function(x) {
  base::gsub("[ .-]+", "_", base::tolower(base::trimws(x)))
}

#' Rename columns to their standard names, ignoring case
#'
#' For each standard name in \code{cols} that is not already present, the first
#' column whose normalized name matches it is renamed. Other columns (e.g.
#' marker names) are left untouched. data.tables are copied first so the
#' caller's object is never modified by reference.
#'
#' @param x data.frame or data.table.
#' @param cols character vector of standard column names.
#' @return \code{x} with matching columns renamed.
#' @noRd
.standardize_names <- function(x, cols) {
  nm   <- base::names(x)
  norm <- .normalize_name(nm)
  old  <- new <- base::character(0)
  for (col in cols) {
    if (col %in% nm) next
    hit <- base::which(norm == col & !(nm %in% old))
    if (base::length(hit) >= 1) {
      old <- c(old, nm[hit[1]])
      new <- c(new, col)
    }
  }
  if (base::length(old) == 0) return(x)
  if (data.table::is.data.table(x)) {
    x <- data.table::copy(x)
    data.table::setnames(x, old, new)
  } else {
    base::names(x)[base::match(old, base::names(x))] <- new
  }
  x
}

#### Genotype input helpers ####

#' Detect the source format of a genotype input
#'
#' Classifies a genotype input so every exported function can accept the same
#' set of formats. New file formats are added here and in \code{.read_genotypes}.
#'
#' @param x genotype input: file path, vcfR object, data.frame / data.table,
#'   or matrix.
#' @return One of "vcf", "plink", "text", "data.frame", or "matrix".
#' @noRd
.genotype_source <- function(x) {
  if (base::inherits(x, "vcfR"))  return("vcf")
  if (base::is.matrix(x))         return("matrix")
  if (base::is.data.frame(x))     return("data.frame")
  if (base::is.character(x) && base::length(x) == 1) {
    if (base::grepl("\\.vcf(\\.gz)?$", x, ignore.case = TRUE)) return("vcf")
    if (base::grepl("\\.ped$",         x, ignore.case = TRUE)) return("plink")
    return("text")
  }
  base::stop("Genotypes must be a file path (TXT/TSV/CSV, VCF or PLINK .ped), ",
             "a vcfR object, a data.frame / data.table, or a matrix.")
}

#' Read genotypes from any supported input
#'
#' Single entry point used by find_parentage(), validate_pedigree(),
#' allele_freq_poly() and solve_composition_poly(). Text files are read with
#' data.table::fread(); VCF files (.vcf / .vcf.gz) and vcfR objects are
#' converted with vcf_to_dosage(); PLINK .ped files with ped_to_dosage();
#' in-memory objects are passed through.
#'
#' @param x genotype input (see \code{.genotype_source}).
#' @param ploidy integer ploidy passed to vcf_to_dosage(); NULL infers it.
#'   Must be 2 (or NULL) for PLINK .ped files.
#' @param format "data.frame" returns a data.table with an id column followed
#'   by marker columns. "matrix" returns individuals in named rows and markers
#'   in columns; a data.frame without an id / ID column is returned unchanged
#'   (row names are assumed to hold the IDs).
#' @param verbose logical, passed to vcf_to_dosage() / ped_to_dosage().
#' @param counted_allele named character vector passed to ped_to_dosage() so a
#'   .ped file is coded like a previous one; ignored for other formats.
#' @return Genotypes in the requested layout.
#' @noRd
.read_genotypes <- function(x, ploidy = NULL,
                            format         = c("data.frame", "matrix"),
                            verbose        = FALSE,
                            counted_allele = NULL) {
  format <- base::match.arg(format)
  src    <- .genotype_source(x)

  # VCF: convert straight to the requested layout
  if (src == "vcf")
    return(vcf_to_dosage(x, ploidy = ploidy, format = format, verbose = verbose))

  # PLINK .ped: diploid only; a .map next to the file supplies marker names
  if (src == "plink") {
    if (!base::is.null(ploidy) && ploidy != 2)
      base::stop("PLINK .ped files are diploid; use ploidy = 2.")
    return(ped_to_dosage(x, format = format, counted_allele = counted_allele,
                         verbose = verbose))
  }

  # Text file: read to a table, then fall through to the in-memory handling
  from_text <- src == "text"
  if (from_text) {
    if (!base::file.exists(x))
      base::stop("Genotype file not found: ", x)
    x   <- data.table::fread(x, sep = "auto", data.table = FALSE)
    src <- "data.frame"
  }

  # ID column: "id" in any case (id, ID, Id, ...); exact "id" preferred
  id_col <- if (src == "data.frame") {
    nm  <- base::names(x)
    hit <- if ("id" %in% nm) "id" else nm[.normalize_name(nm) == "id"]
    if (base::length(hit) >= 1) hit[1] else NA_character_
  } else {
    NA_character_
  }
  # Text files have no row names, so the first column holds the IDs
  if (from_text && format == "matrix" && base::is.na(id_col))
    id_col <- base::names(x)[1]

  if (format == "data.frame") {
    if (src == "matrix")
      return(data.table::as.data.table(x, keep.rownames = "id"))
    out <- data.table::as.data.table(x)
    if (!base::is.na(id_col) && id_col != "id")
      data.table::setnames(out, id_col, "id")
    return(out)
  }

  # format == "matrix": individuals in named rows, markers in columns
  if (src == "matrix" || base::is.na(id_col)) return(x)
  x   <- base::as.data.frame(x)   # data.table subsetting differs from data.frame
  ids <- base::as.character(x[[id_col]])
  x   <- base::as.matrix(x[, base::names(x) != id_col, drop = FALSE])
  base::rownames(x) <- ids
  x
}

#'
#' Performs whole genome breed composition prediction.
#'
#' @param Y numeric vector of genotypes (with names as SNPs) from a single animal.
#'   coded as dosage of allele B \code{{0, 1, 2, ..., ploidy}}
#' @param X numeric matrix of allele frequencies from reference animals
#' @param p numeric indicating number of breeds represented in X
#' @param names character names of breeds
#' @return data.frame of breed composition estimates
#' @import quadprog
#' @importFrom stats cor
#' @references Funkhouser SA, Bates RO, Ernst CW, Newcom D, Steibel JP. Estimation of genome-wide and locus-specific
#' breed composition in pigs. Transl Anim Sci. 2017 Feb 1;1(1):36-44.
#'
#' @noRd
QPsolve <- function(Y, X) {
  
  # Remove NAs from Y and remove corresponding
  #   SNPs from X. Ensure Y is numeric
  Ymod <- Y[!is.na(Y)]
  Xmod <- X[names(Ymod), ]
  
  # Determine properties from X matrix - the number of parameters (breeds) p
  #   and the names of those parameters.
  p <- ncol(X)
  names <- colnames(X)
  
  # perfom steps needed to solve OLS by framing
  # as a QP problem
  # Rinv - matrix should be of dimensions px(p+1) where p is the number of variables in X
  Rinv <- solve(chol(t(Xmod) %*% Xmod))
  
  # C - the first column is a sum restriction (all equal to 1) and the rest of the columns an identity matrix
  C <- cbind(rep(1, p), diag(p))
  
  # b2 - This should be a vector of length p+1 the first element is the value of the sum (1)
  #   the other elements are the restriction of individual coefficients (>)
  #   so a value 0 produces positive coefficients
  b2 <- c(1, rep(0, p))
  
  # dd - this should be a matrix NOT a vector
  dd <- (t(Ymod) %*% Xmod)
  
  qp <- solve.QP(Dmat = Rinv, factorized = TRUE, dvec = dd, Amat = C, bvec = b2, meq = 1)
  beta <- qp$solution
  rr <- cor(Ymod, Xmod %*% beta)^2
  result <- c(beta, rr)
  names(result) <- c(names, "R2")
  return(result)
}

