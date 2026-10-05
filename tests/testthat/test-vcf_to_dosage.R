#### Helper: write a small VCF to a temporary file ####

write_test_vcf <- function(gt_rows, samples, ids = NULL, alts = NULL,
                           format = "GT:AD:DP") {
  path <- base::tempfile(fileext = ".vcf")
  head <- c(
    "##fileformat=VCFv4.3",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage\">",
    base::paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                  "INFO", "FORMAT", samples), collapse = "\t")
  )
  body <- base::vapply(base::seq_along(gt_rows), function(i) {
    calls <- base::vapply(gt_rows[[i]], function(g) {
      if (format == "GT:DS") {
        ds <- if (base::grepl(".", g, fixed = TRUE)) "." else
          base::as.character(base::sum(base::strsplit(g, "[/|]")[[1]] != "0"))
        base::paste(g, ds, sep = ":")
      } else {
        base::paste0(g, ":10,10:20")
      }
    }, character(1))
    base::paste(c("chr1", 1000 + i * 10,
                  if (base::is.null(ids)) "." else ids[i],
                  "A", if (base::is.null(alts)) "T" else alts[i],
                  ".", ".", "NS=4", format, calls), collapse = "\t")
  }, character(1))
  base::writeLines(c(head, body), path)
  path
}

samples <- c("S1", "S2", "S3", "S4")

tetraploid <- list(
  c("0/0/0/0", "0/0/0/1", "0/1/1/1", "1/1/1/1"),
  c("0/0/1/1", "./././.", "0/0/0/0", "0/1/1/1"),
  c("1/1/1/1", "1/1/1/1", "0/0/1/1", "0/0/0/1")
)

#### Tests ####

test_that("a tetraploid VCF converts to id + marker dosage columns", {
  skip_if_not_installed("vcfR")
  path <- write_test_vcf(tetraploid, samples)
  geno <- vcf_to_dosage(path, verbose = FALSE)

  expect_true(data.table::is.data.table(geno))
  expect_identical(base::names(geno)[1], "id")
  expect_identical(geno$id, samples)
  expect_identical(base::dim(geno), c(4L, 4L))
  expect_identical(base::names(geno)[-1],
                   c("chr1_1010", "chr1_1020", "chr1_1030"))
  expect_equal(base::as.numeric(base::unlist(geno[1, -1])), c(0, 2, 4))
  expect_equal(base::as.numeric(geno[[2]]), c(0, 1, 3, 4))
  expect_true(base::is.na(geno[["chr1_1020"]][2]))     # ./././. is missing
})

test_that("ploidy is inferred from the GT field for any ploidy level", {
  skip_if_not_installed("vcfR")
  levels_by_ploidy <- list(
    "2" = list(c("0/0", "0/1", "1/1", "./.")),
    "3" = list(c("0/0/0", "0/0/1", "0/1/1", "1/1/1")),
    "4" = list(c("0/0/0/0", "0/0/0/1", "0/0/1/1", "1/1/1/1")),
    "6" = list(c("0/0/0/0/0/0", "0/0/0/1/1/1", "1/1/1/1/1/1", "0/1/1/1/1/1"))
  )
  expected <- list("2" = c(0, 1, 2, NA), "3" = c(0, 1, 2, 3),
                   "4" = c(0, 1, 2, 4),  "6" = c(0, 3, 6, 5))

  for (p in base::names(levels_by_ploidy)) {
    path <- write_test_vcf(levels_by_ploidy[[p]], samples)
    geno <- vcf_to_dosage(path, format = "matrix", verbose = FALSE)
    expect_identical(base::attr(geno, "ploidy"), base::as.integer(p))
    expect_equal(base::as.numeric(geno[, 1]), expected[[p]])
  }
})

test_that("the matrix layout matches what allele_freq_poly expects", {
  skip_if_not_installed("vcfR")
  path <- write_test_vcf(tetraploid, samples)
  geno <- vcf_to_dosage(path, format = "matrix", verbose = FALSE)

  expect_true(base::is.matrix(geno))
  expect_identical(base::rownames(geno), samples)          # individuals in rows
  expect_identical(base::ncol(geno), 3L)                   # markers in columns

  freq <- allele_freq_poly(geno, list(A = c("S1", "S2"), B = c("S3", "S4")),
                           ploidy = 4)
  expect_identical(base::rownames(freq), base::colnames(geno))
  expect_true(base::all(freq >= 0 & freq <= 1, na.rm = TRUE))
})

test_that("phased calls and the ID column are handled", {
  skip_if_not_installed("vcfR")
  path <- write_test_vcf(list(c("0|0", "0|1", "1|1", "./.")), samples,
                         ids = "snp1")
  geno <- vcf_to_dosage(path, verbose = FALSE)

  expect_identical(base::names(geno)[2], "snp1")           # marker_names = "auto"
  expect_equal(base::as.numeric(geno$snp1), c(0, 1, 2, NA))

  pos <- vcf_to_dosage(path, marker_names = "position", verbose = FALSE)
  expect_identical(base::names(pos)[2], "chr1_1010")
})

test_that("calls that disagree with the file ploidy are set to missing", {
  skip_if_not_installed("vcfR")
  path <- write_test_vcf(list(c("0/0/0/0", "0/0/0/1", "0/1", "1/1/1/1"),
                              c("0/0/1/1", "0/0/0/0", "0/1/1/1", "1/1/1/1")),
                         samples)
  expect_warning(geno <- vcf_to_dosage(path, verbose = FALSE),
                 "ploidy other than 4")
  expect_true(base::is.na(geno[[2]][3]))
})

test_that("a requested ploidy that contradicts the file warns", {
  skip_if_not_installed("vcfR")
  path <- write_test_vcf(tetraploid, samples)
  expect_warning(vcf_to_dosage(path, ploidy = 2, verbose = FALSE),
                 "does not match")
})

test_that("multiallelic markers are dropped or collapsed", {
  skip_if_not_installed("vcfR")
  path <- write_test_vcf(list(c("0/0/0/0", "0/1/1/2", "0/0/2/2", "1/1/1/1"),
                              c("0/0/0/1", "0/0/1/1", "1/1/1/1", "0/0/0/0")),
                         samples, alts = c("T,G", "T"))

  expect_warning(dropped <- vcf_to_dosage(path, verbose = FALSE),
                 "multiallelic")
  expect_identical(base::ncol(dropped), 2L)                 # id + one marker

  kept <- vcf_to_dosage(path, biallelic_only = FALSE, verbose = FALSE)
  expect_equal(base::as.numeric(kept[[2]]), c(0, 3, 2, 4))  # any non-reference allele
})

test_that("reference dosage and a FORMAT dosage field give consistent results", {
  skip_if_not_installed("vcfR")
  path <- write_test_vcf(tetraploid, samples)
  alt  <- vcf_to_dosage(path, format = "matrix", verbose = FALSE)
  ref  <- vcf_to_dosage(path, format = "matrix", ref_dosage = TRUE, verbose = FALSE)
  expect_true(base::all(alt + ref == 4, na.rm = TRUE))

  ds_path <- write_test_vcf(tetraploid, samples, format = "GT:DS")
  ds <- vcf_to_dosage(ds_path, ploidy = 4, dosage_field = "DS",
                      format = "matrix", verbose = FALSE)
  expect_equal(base::as.vector(ds), base::as.vector(alt))
  expect_error(vcf_to_dosage(ds_path, dosage_field = "DS", verbose = FALSE),
               "ploidy must be supplied")
})

test_that("call-rate filters remove markers and individuals", {
  skip_if_not_installed("vcfR")
  path <- write_test_vcf(list(c("0/0/0/0", "./././.", "./././.", "./././."),
                              c("0/0/0/1", "0/0/1/1", "0/1/1/1", "1/1/1/1"),
                              c("0/0/0/0", "0/0/0/1", "./././.", "1/1/1/1")),
                         samples)

  expect_identical(base::ncol(vcf_to_dosage(path, min_marker_call_rate = 0.5,
                                            verbose = FALSE)), 3L)
  expect_identical(base::nrow(vcf_to_dosage(path, min_sample_call_rate = 0.9,
                                            verbose = FALSE)), 1L)
})

test_that("invalid input is rejected", {
  skip_if_not_installed("vcfR")
  path <- write_test_vcf(tetraploid, samples)

  expect_error(vcf_to_dosage("does_not_exist.vcf"), "not found")
  expect_error(vcf_to_dosage(42), "must be a path")
  expect_error(vcf_to_dosage(path, ploidy = 3.5), "ploidy must be")
  expect_error(vcf_to_dosage(path, min_marker_call_rate = 2), "between 0 and 1")
  expect_error(vcf_to_dosage(path, marker_names = "id", verbose = FALSE), "ID column")
})

test_that("a vcfR object is accepted directly", {
  skip_if_not_installed("vcfR")
  path <- write_test_vcf(tetraploid, samples)
  from_path <- vcf_to_dosage(path, format = "matrix", verbose = FALSE)
  from_obj  <- vcf_to_dosage(vcfR::read.vcfR(path, verbose = FALSE),
                             format = "matrix", verbose = FALSE)
  expect_identical(from_obj, from_path)
})
