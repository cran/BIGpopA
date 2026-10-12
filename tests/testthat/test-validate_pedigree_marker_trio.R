# tests/testthat/test-validate_pedigree_marker_trio.R

library(testthat)
library(data.table)

mt_data <- function() {
  # one trio (C = P1 x P2), 11 markers, one status per marker
  geno <- data.frame(
    id = c("P1", "P2", "C"),
    M1  = c(0L,  2L,  1L),   # match
    M2  = c(0L,  2L,  0L),   # offspring impossible with the female alone
    M3  = c(0L,  0L,  2L),   # impossible with either parent
    M4  = c(2L,  0L,  0L),   # impossible with the male alone
    M5  = c(NA,  1L,  1L),   # male missing
    M6  = c(0L,  2L,  NA),   # progeny missing
    M7  = c(0L,  NA,  1L),   # female missing
    M8  = c(0L,  0L,  1L),   # each parent fine alone, pair impossible
    M9  = c(NA,  NA,  1L),   # both parents missing
    M10 = c(NA,  1L,  NA),   # progeny and male missing
    M11 = c(NA,  NA,  NA)    # all missing
  )
  ped <- data.frame(id = "C", male_parent = "P1", female_parent = "P2",
                    stringsAsFactors = FALSE)
  list(geno = geno, ped = ped)
}

test_that("marker_trio_table is NULL unless requested", {
  d   <- mt_data()
  out <- validate_pedigree(d$ped, d$geno, verbose = FALSE, plot_results = FALSE)
  expect_null(out$marker_trio_table)
})

test_that("marker_trio_table gives the expected status per marker", {
  d   <- mt_data()
  out <- validate_pedigree(d$ped, d$geno, verbose = FALSE, plot_results = FALSE,
                           marker_trio_table = TRUE)
  tab <- as.data.frame(out$marker_trio_table)
  expect_equal(names(tab), c("marker", "C_P1_P2"))
  expect_equal(tab$marker, paste0("M", 1:11))
  expect_equal(tab$C_P1_P2, c(
    "match", "mismatch_female", "mismatch_both", "mismatch_male",
    "missing_male", "missing_progeny", "missing_female",
    "mismatch_combination", "missing_male_female",
    "missing_progeny_male", "missing_progeny_male_female"))
})

test_that("one column per trio with both parents known and genotyped", {
  d <- mt_data()
  d$ped <- rbind(d$ped,
                 data.frame(id = "P1", male_parent = "0", female_parent = "0"),
                 data.frame(id = "X",  male_parent = "P1", female_parent = "P2"))
  out <- validate_pedigree(d$ped, d$geno, verbose = FALSE, plot_results = FALSE,
                           marker_trio_table = TRUE)
  expect_equal(setdiff(names(out$marker_trio_table), "marker"), "C_P1_P2")
})
