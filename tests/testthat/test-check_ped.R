# tests/testthat/test-check_ped.R
library(testthat)

write_ped <- function(df) {
  f <- tempfile(fileext = ".txt")
  utils::write.table(df, f, sep = "\t", row.names = FALSE, quote = FALSE)
  f
}

context("check_ped – Pedigree Quality Checks")


# Return structure


test_that("check_ped returns a named list of length 6", {
  ped <- data.frame(
    id            = c("A", "B", "C"),
    male_parent   = c("0", "A", "A"),
    female_parent = c("0", "0", "0")
  )
  out <- check_ped(write_ped(ped), seed = 1, verbose = FALSE)
  expect_type(out, "list")
  expect_length(out, 6)
  expect_named(out, c(
    "exact_duplicates",
    "conflicting_trios",
    "inconsistent_sex_roles",
    "missing_parents",
    "dependencies",
    "corrected_pedigree"
  ))
})

test_that("check_ped report components are data.frames", {
  ped <- data.frame(
    id            = c("A", "B", "C"),
    male_parent   = c("0", "A", "A"),
    female_parent = c("0", "0", "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE)
  expect_true(is.data.frame(out$exact_duplicates))
  expect_true(is.data.frame(out$conflicting_trios))
  expect_true(is.data.frame(out$inconsistent_sex_roles))
  expect_true(is.data.frame(out$missing_parents))
  expect_true(is.data.frame(out$dependencies))
  expect_true(is.data.frame(out$corrected_pedigree))
})

test_that("corrected_pedigree has lowercase column names and no row_number", {
  ped <- data.frame(
    id            = c("A", "B"),
    male_parent   = c("0", "A"),
    female_parent = c("0", "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE)
  expect_true(all(c("id", "male_parent", "female_parent") %in%
                    names(out$corrected_pedigree)))
  expect_false("row_number" %in% names(out$corrected_pedigree))
})


# Clean pedigree


test_that("clean pedigree produces no issues", {
  ped <- data.frame(
    id            = c("G1", "G2", "P1"),
    male_parent   = c("0",  "0",  "G1"),
    female_parent = c("0",  "0",  "G2")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE)
  expect_equal(nrow(out$exact_duplicates),      0)
  expect_equal(nrow(out$conflicting_trios),      0)
  expect_equal(nrow(out$inconsistent_sex_roles), 0)
  expect_equal(nrow(out$missing_parents),        0)
  expect_equal(nrow(out$dependencies),           0)
  expect_equal(nrow(out$corrected_pedigree),     3)
})


# Check 1: Exact duplicates


test_that("check_ped detects exact duplicates", {
  ped <- data.frame(
    id            = c("A", "A", "B"),
    male_parent   = c("0", "0", "A"),
    female_parent = c("0", "0", "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE)
  expect_equal(nrow(out$exact_duplicates), 2)
  expect_true(all(out$exact_duplicates$id == "A"))
})

test_that("exact duplicates are collapsed to one row in corrected_pedigree", {
  ped <- data.frame(
    id            = c("A", "A", "B"),
    male_parent   = c("0", "0", "A"),
    female_parent = c("0", "0", "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE)
  expect_equal(sum(out$corrected_pedigree$id == "A"), 1)
})

test_that("exact duplicates do not appear in conflicting_trios", {
  ped <- data.frame(
    id            = c("A", "A", "B"),
    male_parent   = c("0", "0", "A"),
    female_parent = c("0", "0", "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE)
  expect_equal(nrow(out$conflicting_trios), 0)
})


# Check 2: Conflicting trios


test_that("check_ped detects conflicting trios", {
  ped <- data.frame(
    id            = c("A", "A", "B"),
    male_parent   = c("X", "Y", "A"),
    female_parent = c("M", "M", "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE)
  expect_equal(nrow(out$conflicting_trios), 2)
  expect_true(all(out$conflicting_trios$id == "A"))
})

test_that("correct_conflicting_trios = TRUE: conflicting field -> '0', consistent kept", {
  ped <- data.frame(
    id            = c("A", "A", "B"),
    male_parent   = c("X", "Y", "A"),
    female_parent = c("M", "M", "0")
  )
  out   <- check_ped(write_ped(ped), verbose = FALSE,
                     correct_conflicting_trios = TRUE)
  a_row <- out$corrected_pedigree[out$corrected_pedigree$id == "A", ]
  expect_equal(nrow(a_row), 1)
  expect_equal(a_row$male_parent,   "0")
  expect_equal(a_row$female_parent, "M")
})

test_that("correct_conflicting_trios = FALSE leaves conflicting rows as-is", {
  ped <- data.frame(
    id            = c("A", "A", "B"),
    male_parent   = c("X", "Y", "A"),
    female_parent = c("M", "M", "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE,
                   correct_conflicting_trios = FALSE)
  expect_equal(sum(out$corrected_pedigree$id == "A"), 2)
})


# Check 3: Missing parents


test_that("check_ped detects missing parents", {
  ped <- data.frame(
    id            = c("A", "B"),
    male_parent   = c("0", "X"),
    female_parent = c("0", "Y")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE)
  expect_equal(nrow(out$missing_parents), 2)
  expect_true("X" %in% out$missing_parents$id)
  expect_true("Y" %in% out$missing_parents$id)
})

test_that("missing parents are added as founder rows in corrected_pedigree", {
  ped <- data.frame(
    id            = c("A", "B"),
    male_parent   = c("0", "X"),
    female_parent = c("0", "Y")
  )
  out   <- check_ped(write_ped(ped), verbose = FALSE)
  x_row <- out$corrected_pedigree[out$corrected_pedigree$id == "X", ]
  y_row <- out$corrected_pedigree[out$corrected_pedigree$id == "Y", ]
  expect_true(nrow(x_row) > 0)
  expect_true(nrow(y_row) > 0)
  expect_equal(x_row$male_parent,   "0")
  expect_equal(x_row$female_parent, "0")
  expect_equal(y_row$male_parent,   "0")
  expect_equal(y_row$female_parent, "0")
})

test_that("a missing parent referenced multiple times is added only once", {
  ped <- data.frame(
    id            = c("B", "C"),
    male_parent   = c("X", "X"),
    female_parent = c("0", "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE)
  expect_equal(sum(out$corrected_pedigree$id == "X"), 1)
})


# Check 4: Inconsistent sex roles


test_that("check_ped detects inconsistent sex roles", {
  ped <- data.frame(
    id            = c("child1", "child2", "P"),
    male_parent   = c("P",      "0",      "0"),
    female_parent = c("0",      "P",      "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE)
  expect_true("inconsistent_sex_roles" %in% names(out))
  expect_gt(nrow(out$inconsistent_sex_roles), 0)
  expect_true(any(out$inconsistent_sex_roles$male_parent   == "P" |
                    out$inconsistent_sex_roles$female_parent == "P"))
})

test_that("correct_inconsistent_sex_roles = TRUE zeros out conflicting parent references", {
  ped <- data.frame(
    id            = c("child1", "child2", "P"),
    male_parent   = c("P",      "0",      "0"),
    female_parent = c("0",      "P",      "0")
  )
  out  <- check_ped(write_ped(ped), verbose = FALSE,
                    correct_inconsistent_sex_roles = TRUE)
  corr <- out$corrected_pedigree
  expect_false(any(corr$male_parent   == "P"))
  expect_false(any(corr$female_parent == "P"))
})

test_that("correct_inconsistent_sex_roles = FALSE leaves conflicting references", {
  ped <- data.frame(
    id            = c("child1", "child2", "P"),
    male_parent   = c("P",      "0",      "0"),
    female_parent = c("0",      "P",      "0")
  )
  out  <- check_ped(write_ped(ped), verbose = FALSE,
                    correct_inconsistent_sex_roles = FALSE)
  corr <- out$corrected_pedigree
  expect_true(any(corr$male_parent == "P" | corr$female_parent == "P"))
})


# Check 5: Dependencies / cycles


test_that("individual that is its own parent is logged as a dependency", {
  ped <- data.frame(
    id            = c("A", "B"),
    male_parent   = c("A", "0"),
    female_parent = c("0", "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE)
  expect_gt(nrow(out$dependencies), 0)
})

test_that("self-parent individual is still present in corrected_pedigree", {
  ped <- data.frame(
    id            = c("A", "B"),
    male_parent   = c("A", "0"),
    female_parent = c("0", "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE)
  expect_true("A" %in% out$corrected_pedigree$id)
})

test_that("check_ped detects a direct two-node cycle", {
  ped <- data.frame(
    id            = c("A", "B"),
    male_parent   = c("B", "A"),
    female_parent = c("0", "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE)
  expect_gt(nrow(out$dependencies), 0)
})

test_that("cycle-involved IDs are still present in corrected_pedigree", {
  ped <- data.frame(
    id            = c("A", "B"),
    male_parent   = c("B", "A"),
    female_parent = c("0", "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE)
  expect_true("A" %in% out$corrected_pedigree$id)
  expect_true("B" %in% out$corrected_pedigree$id)
})

test_that("three-node chain cycle is detected in dependencies", {
  ped <- data.frame(
    id            = c("A", "B", "C"),
    male_parent   = c("C", "A", "B"),
    female_parent = c("0", "0", "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE)
  expect_gt(nrow(out$dependencies), 0)
})

test_that("dependencies data.frame has a 'dependency' column of type character", {
  ped <- data.frame(
    id            = c("A", "B"),
    male_parent   = c("B", "A"),
    female_parent = c("0", "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE)
  expect_true("dependency" %in% names(out$dependencies))
  expect_true(is.character(out$dependencies$dependency))
})


# Input validation


test_that("check_ped errors when required columns are missing", {
  bad_df <- data.frame(
    animal_id = c("a", "b"),
    parent1   = c("0", "a"),
    parent2   = c("0", "0")
  )
  expect_error(
    check_ped(write_ped(bad_df), verbose = FALSE),
    regexp = "missing required column"
  )
})

test_that("non-existent file path raises an error", {
  expect_error(
    check_ped("non_existent_file_xyz.txt", verbose = FALSE)
  )
})

test_that("invalid input type raises a descriptive error for check_ped", {
  expect_error(
    check_ped(list(id = "A"), verbose = FALSE),
    regexp = "file path"
  )
})


# Column name flexibility


test_that("check_ped accepts mixed-case column names (ID, Male_Parent, Female_Parent)", {
  ped <- data.frame(
    ID            = c("A", "B", "C"),
    Male_Parent   = c("0", "A", "A"),
    Female_Parent = c("0", "0", "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE)
  expect_length(out, 6)
  expect_true(all(c("id", "male_parent", "female_parent") %in%
                    names(out$corrected_pedigree)))
})

test_that("check_ped accepts all-uppercase column names", {
  ped <- data.frame(
    ID            = c("A", "B"),
    MALE_PARENT   = c("0", "A"),
    FEMALE_PARENT = c("0", "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE)
  expect_length(out, 6)
})

test_that("check_ped accepts columns in any order", {
  ped <- data.frame(
    female_parent = c("0", "0"),
    male_parent   = c("0", "A"),
    id            = c("A", "B")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE)
  expect_length(out, 6)
  expect_equal(nrow(out$corrected_pedigree), 2)
})

test_that("extra columns beyond the three required are preserved in corrected_pedigree", {
  ped <- data.frame(
    id            = c("A", "B"),
    male_parent   = c("0", "A"),
    female_parent = c("0", "0"),
    cohort        = c(2020L, 2021L)
  )
  out <- check_ped(write_ped(ped), verbose = FALSE)
  expect_true("cohort" %in% names(out$corrected_pedigree))
})


# In-memory input


test_that("check_ped accepts a data.frame directly", {
  ped <- data.frame(
    id            = c("A", "B", "C"),
    male_parent   = c("0", "A", "A"),
    female_parent = c("0", "0", "0")
  )
  out <- check_ped(ped, verbose = FALSE)
  expect_length(out, 6)
  expect_true(all(c("id", "male_parent", "female_parent") %in%
                    names(out$corrected_pedigree)))
})

test_that("check_ped accepts a data.table directly", {
  ped <- data.table::data.table(
    id            = c("A", "B", "C"),
    male_parent   = c("0", "A", "A"),
    female_parent = c("0", "0", "0")
  )
  out <- check_ped(ped, verbose = FALSE)
  expect_length(out, 6)
  expect_true(all(c("id", "male_parent", "female_parent") %in%
                    names(out$corrected_pedigree)))
})

test_that("in-memory and file-path inputs produce identical corrected_pedigree", {
  ped <- data.frame(
    id            = c("A", "B", "C"),
    male_parent   = c("0", "A", "A"),
    female_parent = c("0", "0", "0")
  )
  out_file <- check_ped(write_ped(ped), verbose = FALSE)
  out_mem  <- check_ped(ped,            verbose = FALSE)
  expect_identical(out_file$corrected_pedigree,
                   out_mem$corrected_pedigree)
})


# Verbosity, seed, and side effects


test_that("verbose = FALSE suppresses console output", {
  ped <- data.frame(
    id            = c("A", "B"),
    male_parent   = c("0", "A"),
    female_parent = c("0", "0")
  )
  expect_silent(check_ped(write_ped(ped), verbose = FALSE))
})

test_that("verbose = TRUE produces console output", {
  ped <- data.frame(
    id            = c("A", "B"),
    male_parent   = c("0", "A"),
    female_parent = c("0", "0")
  )
  expect_output(check_ped(write_ped(ped), verbose = TRUE))
})

test_that("check_ped returns invisibly", {
  ped <- data.frame(
    id            = c("A", "B"),
    male_parent   = c("0", "A"),
    female_parent = c("0", "0")
  )
  expect_invisible(check_ped(write_ped(ped), verbose = FALSE))
})

test_that("seed produces reproducible results", {
  ped <- data.frame(
    id            = c("A", "B", "C"),
    male_parent   = c("0", "A", "A"),
    female_parent = c("0", "0", "0")
  )
  f    <- write_ped(ped)
  out1 <- check_ped(f, seed = 42, verbose = FALSE)
  out2 <- check_ped(f, seed = 42, verbose = FALSE)
  expect_identical(out1$corrected_pedigree, out2$corrected_pedigree)
})

test_that("seed = NULL runs without error", {
  ped <- data.frame(
    id            = c("A", "B"),
    male_parent   = c("0", "A"),
    female_parent = c("0", "0")
  )
  expect_no_error(check_ped(write_ped(ped), seed = NULL, verbose = FALSE))
})

test_that("no output files are written to disk", {
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  old_wd  <- getwd()
  setwd(tmp_dir)
  on.exit({ setwd(old_wd); unlink(tmp_dir, recursive = TRUE) }, add = TRUE)
  ped <- data.frame(
    id            = c("A", "B"),
    male_parent   = c("0", "A"),
    female_parent = c("0", "0")
  )
  check_ped(write_ped(ped), verbose = FALSE)
  expect_length(list.files(tmp_dir), 0)
})


# Correction flag interactions


test_that("both correction flags FALSE still returns a corrected_pedigree", {
  ped <- data.frame(
    id            = c("A", "A", "child1", "child2", "P"),
    male_parent   = c("X", "Y", "P",      "0",      "0"),
    female_parent = c("M", "M", "0",      "P",      "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE,
                   correct_conflicting_trios      = FALSE,
                   correct_inconsistent_sex_roles = FALSE)
  expect_true(is.data.frame(out$corrected_pedigree))
  expect_gt(nrow(out$corrected_pedigree), 0)
})


# Integration test


test_that("integration test with bundled fixture file", {
  ped_file <- system.file("check_ped_test.txt", package = "BIGpopA")
  skip_if(ped_file == "", "Bundled fixture file not found; skipping.")
  
  raw        <- utils::read.table(ped_file, header = TRUE)
  names(raw) <- tolower(names(raw))
  
  if ("sire"   %in% names(raw)) names(raw)[names(raw) == "sire"]   <- "male_parent"
  if ("dam"    %in% names(raw)) names(raw)[names(raw) == "dam"]    <- "female_parent"
  if ("animal" %in% names(raw)) names(raw)[names(raw) == "animal"] <- "id"
  
  out <- check_ped(raw, seed = 101919, verbose = FALSE)
  
  expect_length(out, 6)
  expect_gt(nrow(out$inconsistent_sex_roles), 0)
  
  conflicting_ids <- unique(c(
    out$inconsistent_sex_roles$male_parent,
    out$inconsistent_sex_roles$female_parent
  ))
  expect_true(any(c("grandfather2", "grandfather3") %in% conflicting_ids))
  expect_equal(nrow(out$missing_parents), 13)
})