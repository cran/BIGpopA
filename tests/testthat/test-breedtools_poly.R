#context("BreedTools - solve_composition_poly")

#  Shared fixtures 
ref_file      <- system.file("test_ref.txt",  package = "BIGpopA")
val_file      <- system.file("test_test.txt", package = "BIGpopA")
ref_ids_file  <- system.file("ref_ids.txt",   package = "BIGpopA")

reference     <- read.table(ref_file,     header = TRUE, row.names = 1, sep = "\t")
validation    <- read.table(val_file,     header = TRUE, row.names = 1, sep = "\t")
reference_ids <- read.table(ref_ids_file, header = TRUE, sep = "\t")

ref_ids <- lapply(as.list(reference_ids), as.character)
freq    <- allele_freq_poly(reference, ref_ids, ploidy = 4)

#  1. Baseline (tetraploid, no ped/groups) 
test_that("baseline tetraploid composition returns expected values", {
  prediction <- as.data.frame(solve_composition_poly(validation, freq, ploidy = 4))
  
  expect_equal(round(mean(as.numeric(freq)),          6), 0.888889, tolerance = 0.01)
  expect_equal(round(mean(as.numeric(prediction$R2)), 6), 0.841454, tolerance = 0.01)
  expect_true(nrow(prediction) == 175)
})

#  2. Default ploidy (diploid) 
test_that("default ploidy = 2 runs without error and returns a data.frame", {
  freq_dip   <- allele_freq_poly(reference, ref_ids, ploidy = 2)
  prediction <- solve_composition_poly(validation, freq_dip)  # ploidy defaults to 2
  
  expect_s3_class(as.data.frame(prediction), "data.frame")
  expect_true(nrow(prediction) > 0)
})

#  3. Ploidy parameter changes output 
test_that("ploidy = 4 scales dosage differently than ploidy = 2", {
  # Use the same freq matrix for both — only ploidy differs
  # Source shows Y is divided by ploidy: Y <- Y / ploidy
  # So with identical input, a higher ploidy should yield different compositions
  pred_2 <- as.data.frame(solve_composition_poly(validation, freq, ploidy = 2))
  pred_4 <- as.data.frame(solve_composition_poly(validation, freq, ploidy = 4))
  
  breed_cols <- intersect(colnames(freq), colnames(pred_2))
  skip_if(length(breed_cols) < 1, "Could not identify breed proportion columns.")
  
  expect_false(isTRUE(all.equal(
    pred_2[, breed_cols, drop = FALSE],
    pred_4[, breed_cols, drop = FALSE]
  )))
})

#  4. Small inline example (tetraploid) 
test_that("inline tetraploid example from documentation runs correctly", {
  allele_freqs_matrix <- matrix(
    c(0.625, 0.500,
      0.500, 0.500,
      0.500, 0.500,
      0.750, 0.500,
      0.625, 0.625),
    nrow = 5, ncol = 2, byrow = TRUE,
    dimnames = list(paste0("SNP", 1:5), c("VarA", "VarB"))
  )
  
  val_geno_matrix <- matrix(
    c(2, 1, 2, 3, 4,
      3, 4, 2, 3, 0),
    nrow = 2, ncol = 5, byrow = TRUE,
    dimnames = list(paste0("Test", 1:2), paste0("SNP", 1:5))
  )
  
  composition <- solve_composition_poly(Y = val_geno_matrix,
                                        X = allele_freqs_matrix,
                                        ploidy = 4)
  
  expect_s3_class(as.data.frame(composition), "data.frame")
  expect_equal(nrow(composition), 2)
})

#  5. Compositions sum to 1 (full dataset) 
test_that("breed composition proportions sum to 1 for each animal", {
  prediction <- as.data.frame(solve_composition_poly(validation, freq, ploidy = 4))
  breed_cols <- colnames(freq)  # reference population names come directly from X
  
  available <- intersect(breed_cols, colnames(prediction))
  skip_if(length(available) < 2, "Could not identify breed proportion columns.")
  
  row_sums <- rowSums(prediction[, available, drop = FALSE])
  expect_true(all(abs(row_sums - 1) < 1e-6))
})

#  6. groups argument returns a named list 
test_that("groups argument returns a named list of data.frames", {
  skip_if_not(existsFunction("QPseparate", where = asNamespace("BIGpopA")),
              "QPseparate not available — skipping groups test.")
  
  all_ids <- rownames(validation)
  half    <- floor(length(all_ids) / 2)
  groups  <- list(GroupA = all_ids[1:half],
                  GroupB = all_ids[(half + 1):length(all_ids)])
  
  result <- solve_composition_poly(validation, freq, groups = groups, ploidy = 4)
  
  expect_type(result, "list")
  expect_named(result, c("GroupA", "GroupB"))
})

#  7. ped argument: basic composition with pedigree 
test_that("ped argument runs and returns a data.frame", {
  ped_file <- system.file("test_ped.txt", package = "BIGpopA")
  skip_if_not(file.exists(ped_file), "Pedigree test file not available.")
  
  ped    <- read.table(ped_file, header = TRUE, sep = "\t")
  result <- solve_composition_poly(validation, freq, ped = ped, ploidy = 4)
  
  expect_s3_class(as.data.frame(result), "data.frame")
  expect_true(nrow(result) > 0)
})

#  8. mia flag returns MIA data.frame 
test_that("mia = TRUE returns maternally inherited allele data.frame", {
  ped_file <- system.file("test_ped.txt", package = "BIGpopA")
  skip_if_not(file.exists(ped_file), "Pedigree test file not available.")
  
  ped    <- read.table(ped_file, header = TRUE, sep = "\t")
  result <- solve_composition_poly(validation, freq, ped = ped,
                                   mia = TRUE, ploidy = 4)
  
  expect_s3_class(as.data.frame(result), "data.frame")
  expect_true(ncol(result) > 0)
})

#  9. sire flag returns sire genotype data.frame 
test_that("sire = TRUE returns sire genotype data.frame", {
  ped_file <- system.file("test_ped.txt", package = "BIGpopA")
  skip_if_not(file.exists(ped_file), "Pedigree test file not available.")
  
  ped    <- read.table(ped_file, header = TRUE, sep = "\t")
  result <- solve_composition_poly(validation, freq, ped = ped,
                                   sire = TRUE, ploidy = 4)
  
  expect_s3_class(as.data.frame(result), "data.frame")
})

#  10. dam flag returns dam genotype data.frame 
test_that("dam = TRUE returns dam genotype data.frame", {
  ped_file <- system.file("test_ped.txt", package = "BIGpopA")
  skip_if_not(file.exists(ped_file), "Pedigree test file not available.")
  
  ped    <- read.table(ped_file, header = TRUE, sep = "\t")
  result <- solve_composition_poly(validation, freq, ped = ped,
                                   dam = TRUE, ploidy = 4)
  
  expect_s3_class(as.data.frame(result), "data.frame")
})

#  11. Extra SNPs in Y not in X are silently dropped 
test_that("extra SNPs in Y not present in X are ignored without error", {
  extra_snps <- matrix(
    rbinom(nrow(validation) * 3, 4, 0.5),
    nrow = nrow(validation),
    dimnames = list(rownames(validation), c("FAKE1", "FAKE2", "FAKE3"))
  )
  validation_extra <- cbind(validation, extra_snps)
  
  expect_no_error(
    solve_composition_poly(validation_extra, freq, ploidy = 4)
  )
})

#  12. Output row count matches number of animals in Y 
test_that("number of rows in output matches number of animals in Y", {
  prediction <- as.data.frame(solve_composition_poly(validation, freq, ploidy = 4))
  expect_equal(nrow(prediction), nrow(validation))
})