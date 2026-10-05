#' Compute Genome-Wide Breed Composition
#'
#' Computes genome-wide breed/ancestry composition using quadratic programming
#' on a batch of animals.
#'
#' @param Y Genotypes (columns) from all animals (rows) in the population,
#'   coded as dosage of allele B (0, 1, 2, ..., ploidy), as any of: a numeric
#'   matrix or data.frame with named rows; a data.frame with an `id` / `ID`
#'   column followed by SNP columns; a path to a TSV/CSV/TXT file whose first
#'   column holds the IDs; a path to a VCF file (`.vcf` or `.vcf.gz`) or a
#'   `vcfR` object, converted with [vcf_to_dosage()] using `ploidy`; or a path
#'   to a PLINK `.ped` file (diploid only), converted with [ped_to_dosage()].
#'   A `.ped` file is coded with the `counted_allele` attribute of `X` when
#'   present (see [allele_freq_poly()]), so reference and validation `.ped`
#'   files count the same allele at each marker.
#' @param X numeric matrix of allele frequencies (rows) from each reference
#'   panel (columns). Frequencies are relative to allele B.
#' @param ploidy integer. The ploidy level of the species (e.g., 2 for diploid,
#'   3 for triploid).
#'
#' @return A matrix with one row per animal, one column per reference
#'   population (estimated proportions, summing to 1), and an `R2` column.
#'
#' @references Funkhouser SA, Bates RO, Ernst CW, Newcom D, Steibel JP.
#'   Estimation of genome-wide and locus-specific breed composition in pigs.
#'   Transl Anim Sci. 2017 Feb 1;1(1):36-44.
#'
#' @importFrom quadprog solve.QP
#'
#' @examples
#' allele_freqs_matrix <- matrix(
#'   c(0.625, 0.500,
#'     0.500, 0.500,
#'     0.500, 0.500,
#'     0.750, 0.500,
#'     0.625, 0.625),
#'   nrow = 5, ncol = 2, byrow = TRUE,
#'   dimnames = list(paste0("SNP", 1:5), c("VarA", "VarB"))
#' )
#'
#' val_geno_matrix <- matrix(
#'   c(2, 1, 2, 3, 4,
#'     3, 4, 2, 3, 0),
#'   nrow = 2, ncol = 5, byrow = TRUE,
#'   dimnames = list(paste0("Test", 1:2), paste0("SNP", 1:5))
#' )
#'
#' composition <- solve_composition_poly(Y = val_geno_matrix,
#'                                       X = allele_freqs_matrix,
#'                                       ploidy = 4)
#' print(composition)
#'
#' @export
solve_composition_poly <- function(Y,
                                   X,
                                   ploidy = 2) {

  # Accept text files, VCFs, PLINK .ped, and in-memory tables (animals x SNPs).
  #   A .ped file reuses the reference panel's allele coding when available.
  Y <- .read_genotypes(Y, ploidy = ploidy, format = "matrix",
                       counted_allele = attr(X, "counted_allele"))

  # Functions require Y to be animals x SNPs. Transpose
  Y <- t(Y)
  
  # SNPs in Y should only be the ones present in X
  Y <- Y[rownames(Y) %in% rownames(X), , drop = FALSE]   # keep a 1-animal Y as a matrix
  
  # Adjust dosage based on ploidy (default is 2)
  Y <- Y / ploidy

  # Solve the composition of each animal (column of Y)
  results <- t(apply(Y, 2, QPsolve, X))
  return (results)
}