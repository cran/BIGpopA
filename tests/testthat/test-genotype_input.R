#### Helpers: write a dosage matrix as a VCF or a text file ####

# Dosage matrix (individuals x markers) -> VCF with GT calls of the given ploidy
write_dosage_vcf <- function(geno, ploidy) {
  path <- base::tempfile(fileext = ".vcf")
  to_gt <- function(d) {
    if (base::is.na(d)) return(base::paste(base::rep(".", ploidy), collapse = "/"))
    base::paste(c(base::rep("0", ploidy - d), base::rep("1", d)), collapse = "/")
  }
  head <- c(
    "##fileformat=VCFv4.3",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    base::paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                  "INFO", "FORMAT", base::rownames(geno)), collapse = "\t")
  )
  body <- base::vapply(base::seq_len(base::ncol(geno)), function(j) {
    calls <- base::vapply(geno[, j], to_gt, character(1))
    base::paste(c("chr1", 1000 + j, base::colnames(geno)[j], "A", "T",
                  ".", ".", ".", "GT", calls), collapse = "\t")
  }, character(1))
  base::writeLines(c(head, body), path)
  path
}

# Dosage matrix -> tab-separated text file with an ID column first
write_dosage_txt <- function(geno, id_name = "id") {
  path <- base::tempfile(fileext = ".txt")
  df   <- base::data.frame(base::rownames(geno), geno, check.names = FALSE)
  base::names(df)[1] <- id_name
  utils::write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE)
  path
}

# Diploid dosage matrix -> PLINK .ped (dosage = copies of `alt`), optional .map
write_dosage_ped <- function(geno, ref = "A", alt = "G", with_map = FALSE) {
  path <- base::tempfile(fileext = ".ped")
  to_alleles <- function(d) {
    if (base::is.na(d)) return(c("0", "0"))
    c(base::rep(ref, 2 - d), base::rep(alt, d))
  }
  lines <- base::vapply(base::seq_len(base::nrow(geno)), function(i) {
    calls <- base::unlist(base::lapply(geno[i, ], to_alleles))
    base::paste(c("FAM", base::rownames(geno)[i], "0", "0", "0", "-9", calls),
                collapse = " ")
  }, character(1))
  base::writeLines(lines, path)
  if (with_map) {
    map <- base::paste("1", base::colnames(geno), "0",
                       1000 + base::seq_len(base::ncol(geno)), sep = "\t")
    base::writeLines(map, base::sub("\\.ped$", ".map", path))
  }
  path
}

#### Fixtures ####

# Diploid trio data: Off1 and Off2 are consistent with P1 x P2
base::set.seed(42)
n_snp <- 30
p1 <- base::sample(0:2, n_snp, replace = TRUE)
p2 <- base::sample(0:2, n_snp, replace = TRUE)
gamete <- function(g) base::ifelse(g == 1, base::sample(0:1, length(g), replace = TRUE), g / 2)
off1 <- gamete(p1) + gamete(p2)
off2 <- gamete(p1) + gamete(p2)
p3   <- base::sample(0:2, n_snp, replace = TRUE)

dip_geno <- base::rbind(P1 = p1, P2 = p2, P3 = p3, Off1 = off1, Off2 = off2)
base::colnames(dip_geno) <- base::paste0("SNP", base::seq_len(n_snp))
base::storage.mode(dip_geno) <- "integer"

dip_df <- data.table::as.data.table(dip_geno, keep.rownames = "id")

parents <- base::data.frame(id = c("P1", "P2", "P3"), sex = c("M", "F", "F"))
progeny <- base::data.frame(id = c("Off1", "Off2"))
ped     <- base::data.frame(id = c("Off1", "Off2"),
                            male_parent = c("P1", "P1"),
                            female_parent = c("P2", "P3"))

# Tetraploid reference / validation data for the breedtools functions
tet_ref <- base::rbind(
  base::matrix(base::sample(2:4, 4 * n_snp, replace = TRUE), nrow = 4),
  base::matrix(base::sample(0:2, 4 * n_snp, replace = TRUE), nrow = 4)
)
base::dimnames(tet_ref) <- base::list(base::paste0("Ref", 1:8),
                                      base::paste0("SNP", base::seq_len(n_snp)))
base::storage.mode(tet_ref) <- "integer"
pops <- base::list(PopA = base::paste0("Ref", 1:4), PopB = base::paste0("Ref", 5:8))

tet_val <- tet_ref[c(1, 6), , drop = FALSE]
base::rownames(tet_val) <- c("Val1", "Val2")

#### .genotype_source ####

test_that(".genotype_source classifies every supported input", {
  expect_identical(BIGpopA:::.genotype_source("a.vcf"),    "vcf")
  expect_identical(BIGpopA:::.genotype_source("a.VCF.gz"), "vcf")
  expect_identical(BIGpopA:::.genotype_source("a.txt"),    "text")
  expect_identical(BIGpopA:::.genotype_source(dip_df),     "data.frame")
  expect_identical(BIGpopA:::.genotype_source(dip_geno),   "matrix")
  expect_error(BIGpopA:::.genotype_source(42), "Genotypes must be")
})

#### find_parentage / validate_pedigree ####

test_that("find_parentage gives identical results from VCF, text and data.frame", {
  skip_if_not_installed("vcfR")
  run <- function(g) find_parentage(g, parents, progeny, verbose = FALSE,
                                    plot_results = FALSE)$full_results

  from_df  <- run(dip_df)
  from_vcf <- run(write_dosage_vcf(dip_geno, ploidy = 2))
  from_txt <- run(write_dosage_txt(dip_geno))

  expect_equal(from_vcf, from_df, ignore_attr = TRUE)
  expect_equal(from_txt, from_df, ignore_attr = TRUE)
})

test_that("validate_pedigree gives identical results from VCF and data.frame", {
  skip_if_not_installed("vcfR")
  run <- function(g) validate_pedigree(ped, g, verbose = FALSE,
                                       plot_results = FALSE)$full_results

  expect_equal(run(write_dosage_vcf(dip_geno, ploidy = 2)), run(dip_df),
               ignore_attr = TRUE)
})

test_that("an 'ID' column in a data.frame is accepted as the id column", {
  df_upper <- data.table::copy(dip_df)
  data.table::setnames(df_upper, "id", "ID")
  res <- find_parentage(df_upper, parents, progeny, verbose = FALSE,
                        plot_results = FALSE)$full_results
  expect_equal(base::nrow(res), 2)
})

test_that("a VCF read error is reported with its cause", {
  expect_error(
    find_parentage("missing_file.vcf", parents, progeny, verbose = FALSE,
                   plot_results = FALSE),
    "Error reading input files"
  )
})

#### allele_freq_poly / solve_composition_poly ####

test_that("allele_freq_poly matches across matrix, VCF, text and ID data.frame", {
  skip_if_not_installed("vcfR")
  ref_mat <- allele_freq_poly(tet_ref, pops, ploidy = 4)

  ref_vcf <- allele_freq_poly(write_dosage_vcf(tet_ref, ploidy = 4), pops, ploidy = 4)
  ref_txt <- allele_freq_poly(write_dosage_txt(tet_ref, id_name = "Sample"), pops, ploidy = 4)
  ref_df  <- allele_freq_poly(
    base::data.frame(ID = base::rownames(tet_ref), tet_ref, check.names = FALSE),
    pops, ploidy = 4
  )

  expect_equal(ref_vcf, ref_mat)
  expect_equal(ref_txt, ref_mat)
  expect_equal(ref_df,  ref_mat)
})

test_that("solve_composition_poly matches between matrix and VCF input", {
  skip_if_not_installed("vcfR")
  freq     <- allele_freq_poly(tet_ref, pops, ploidy = 4)
  from_mat <- solve_composition_poly(tet_val, freq, ploidy = 4)
  from_vcf <- solve_composition_poly(write_dosage_vcf(tet_val, ploidy = 4),
                                     freq, ploidy = 4)
  expect_equal(from_vcf, from_mat)
})

#### PLINK .ped ####

test_that("ped_to_dosage counts the last allele and treats 0 as missing", {
  path <- base::tempfile(fileext = ".ped")
  base::writeLines(c(
    "F P1 0 0 1 -9 1 1 2 2 1 2",
    "F P2 0 0 2 -9 1 2 0 0 2 2"
  ), path)
  geno <- ped_to_dosage(path, map = NULL, verbose = FALSE)

  expect_identical(geno$id, c("P1", "P2"))
  expect_identical(base::names(geno)[-1], c("SNP1", "SNP2", "SNP3"))
  expect_equal(geno$SNP1, c(0L, 1L))
  expect_equal(geno$SNP2, c(2L, NA))
  expect_equal(geno$SNP3, c(1L, 2L))
  expect_equal(base::unname(base::attr(geno, "counted_allele")), c("2", "2", "2"))
})

test_that("a sibling .map supplies marker names; map = NULL ignores it", {
  path <- write_dosage_ped(dip_geno, with_map = TRUE)
  expect_identical(base::names(ped_to_dosage(path, verbose = FALSE))[-1],
                   base::colnames(dip_geno))
  expect_identical(base::names(ped_to_dosage(path, map = NULL, verbose = FALSE))[2],
                   "SNP1")
})

test_that("counted_allele forces the coding of a monomorphic marker", {
  path <- base::tempfile(fileext = ".ped")
  base::writeLines("F V1 0 0 0 -9 A A G G", path)
  default <- ped_to_dosage(path, map = NULL, verbose = FALSE)
  forced  <- ped_to_dosage(path, map = NULL, verbose = FALSE,
                           counted_allele = c(SNP1 = "G", SNP2 = "G"))
  expect_equal(default$SNP1, 2L)   # only A observed -> A counted
  expect_equal(forced$SNP1,  0L)   # reference counts G -> 0 copies
  expect_equal(forced$SNP2,  2L)
})

test_that("markers with more than two alleles are dropped with a warning", {
  path <- base::tempfile(fileext = ".ped")
  base::writeLines(c("F P1 0 0 0 -9 A C A G", "F P2 0 0 0 -9 T T A A"), path)
  expect_warning(geno <- ped_to_dosage(path, map = NULL, verbose = FALSE),
                 "more than two alleles")
  expect_identical(base::names(geno)[-1], "SNP2")
})

test_that("find_parentage gives identical results from .ped and data.frame", {
  run <- function(g) find_parentage(g, parents, progeny, verbose = FALSE,
                                    plot_results = FALSE)$full_results
  expect_equal(run(write_dosage_ped(dip_geno, with_map = TRUE)), run(dip_df),
               ignore_attr = TRUE)
})

test_that(".ped input requires ploidy = 2", {
  expect_error(
    find_parentage(write_dosage_ped(dip_geno), parents, progeny, ploidy = 4,
                   verbose = FALSE, plot_results = FALSE),
    "diploid"
  )
})

test_that("a single-sample validation .ped is coded like the reference .ped", {
  # Diploid reference with two populations; validation = one homozygous-rich sample
  dip_ref <- base::rbind(
    base::matrix(base::sample(1:2, 4 * n_snp, replace = TRUE), nrow = 4),
    base::matrix(base::sample(0:1, 4 * n_snp, replace = TRUE), nrow = 4)
  )
  base::dimnames(dip_ref) <- base::dimnames(tet_ref)
  # A PopB sample (dosages 0/1): its 0-dosage markers show only allele A, so
  #   without the reference coding they would be flipped
  dip_val <- dip_ref[6, , drop = FALSE]
  base::rownames(dip_val) <- "Val1"
  val_ped <- write_dosage_ped(dip_val, with_map = TRUE)

  freq_mat <- allele_freq_poly(dip_ref, pops)
  from_mat <- solve_composition_poly(dip_val, freq_mat)

  freq_ped <- allele_freq_poly(write_dosage_ped(dip_ref, with_map = TRUE), pops)
  from_ped <- solve_composition_poly(val_ped, freq_ped)

  expect_false(base::is.null(base::attr(freq_ped, "counted_allele")))
  expect_equal(from_ped, from_mat, tolerance = 1e-8, ignore_attr = TRUE)

  # Without the reference coding the single sample is coded differently
  uncoded <- ped_to_dosage(val_ped, format = "matrix", verbose = FALSE)
  expect_false(base::isTRUE(base::all.equal(uncoded, dip_val, check.attributes = FALSE)))
})

test_that("row-named data.frames are still accepted unchanged", {
  ref_df <- base::as.data.frame(tet_ref)
  expect_equal(allele_freq_poly(ref_df, pops, ploidy = 4),
               allele_freq_poly(tet_ref, pops, ploidy = 4))
})
