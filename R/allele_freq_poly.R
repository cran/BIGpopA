#' Compute Allele Frequencies for Populations
#'
#' Computes allele frequencies for specified populations given SNP array data.
#'
#' @param geno Genotypes coded as the dosage of allele B (0, 1, 2, ..., ploidy)
#'   as any of: a matrix or data.frame with individuals in rows (named) and
#'   SNPs in columns (named); a data.frame with an `id` / `ID` column followed
#'   by SNP columns; a path to a TSV/CSV/TXT file whose first column holds the
#'   IDs; a path to a VCF file (`.vcf` or `.vcf.gz`) or a `vcfR` object,
#'   converted with [vcf_to_dosage()] using `ploidy`; or a path to a PLINK
#'   `.ped` file (diploid only), converted with [ped_to_dosage()].
#'
#' @details When `geno` is a `.ped` file, the returned matrix carries a
#'   `counted_allele` attribute. [solve_composition_poly()] uses it to code a
#'   validation `.ped` file with the same counted allele at each marker.
#' @param populations list of named populations. Each population has a vector
#'   of IDs that belong to the population. Allele frequencies will be derived
#'   from all animals in each population.
#' @param ploidy integer indicating the ploidy level (default is 2 for diploid).
#'
#' @return A matrix of allele frequencies with SNPs in rows and populations in
#'   columns.
#'
#' @references Funkhouser SA, Bates RO, Ernst CW, Newcom D, Steibel JP.
#'   Estimation of genome-wide and locus-specific breed composition in pigs.
#'   Transl Anim Sci. 2017 Feb 1;1(1):36-44.
#'
#' @examples
#' geno_matrix <- matrix(
#'   c(4, 1, 4, 0,
#'     2, 2, 1, 3,
#'     0, 4, 0, 4,
#'     3, 3, 2, 2,
#'     1, 4, 2, 3),
#'   nrow = 4, ncol = 5, byrow = FALSE,
#'   dimnames = list(paste0("Ind", 1:4), paste0("S", 1:5))
#' )
#'
#' pop_list <- list(
#'   PopA = c("Ind1", "Ind2"),
#'   PopB = c("Ind3", "Ind4")
#' )
#'
#' allele_freqs <- allele_freq_poly(geno = geno_matrix,
#'                                  populations = pop_list,
#'                                  ploidy = 4)
#' print(allele_freqs)
#'
#' @export
allele_freq_poly <- function(geno, populations, ploidy = 2) {

  # Accept text files, VCFs, PLINK .ped, and in-memory tables (individuals x SNPs)
  geno    <- .read_genotypes(geno, ploidy = ploidy, format = "matrix")
  counted <- base::attr(geno, "counted_allele")   # set only for .ped input

  # Initialize returned df
  X <- matrix(NA, nrow = ncol(geno), ncol = length(populations))
  
  # Subset geno into different populations
  for (i in 1:length(populations)) {
    
    # Get name of ith item in the list (population name)
    pop_name <- names(populations[i])
    
    # Subset geno to only include genotypes of IDs in pop
    pop_geno <- geno[rownames(geno) %in% populations[[i]], ]
    
    # Calculate allele frequencies
    al_freq <- colMeans(pop_geno, na.rm = TRUE) / ploidy
    
    # Add to X
    X[, i] <- al_freq
  }
  
  # Label X with populations and SNPs
  colnames(X) <- names(populations)
  rownames(X) <- colnames(geno)

  # Keep the .ped allele coding so solve_composition_poly() codes validation
  #   .ped files the same way
  if (!is.null(counted)) attr(X, "counted_allele") <- counted

  return(X)
}