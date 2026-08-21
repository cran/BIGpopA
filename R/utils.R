utils::globalVariables(c(
  # data.table internals
  ":=", ".SD",
  
  # find_parentage.R
  "id", "sex", "male_parent", "female_parent",
  "mendelian_error_pct", "plot_status", "status",
  
  # validate_pedigree.R
  "trio_mendelian_error_pct", "recommended_correction",
  
  # breedtools internal helpers
  "QPseparate", "QPsolve_par"
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

