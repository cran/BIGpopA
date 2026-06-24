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

