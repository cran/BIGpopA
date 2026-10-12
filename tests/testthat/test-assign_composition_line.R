# tests/testthat/test-assign_composition_line.R

comp <- rbind(
  S1 = c(A = 0.998, B = 0.000, C = 0.002, R2 = 0.99),
  S2 = c(A = 0.529, B = 0.471, C = 0.000, R2 = 0.97),
  S3 = c(A = 0.340, B = 0.330, C = 0.330, R2 = 0.95)
)

test_that("default threshold is 50 and weak assignments are Undetermined", {
  out <- assign_composition_line(comp)
  expect_equal(names(out), c("ID", "A", "B", "C", "Predicted line"))
  expect_equal(out$ID, c("S1", "S2", "S3"))
  expect_equal(out$`Predicted line`, c("A", "A", "Undetermined"))
  expect_equal(out$A, c(99.8, 52.9, 34.0))
})

test_that("threshold 0 always assigns the highest population", {
  out <- assign_composition_line(comp, assign_threshold = 0)
  expect_equal(out$`Predicted line`, c("A", "A", "A"))
})

test_that("a proportion equal to the threshold is assigned", {
  m <- rbind(S1 = c(A = 0.5, B = 0.5))
  expect_equal(assign_composition_line(m, 50)$`Predicted line`, "A")
})

test_that("R2 is dropped and invalid thresholds are rejected", {
  expect_false("R2" %in% names(assign_composition_line(comp)))
  expect_error(assign_composition_line(comp, assign_threshold = 101))
  expect_error(assign_composition_line(comp, assign_threshold = "a"))
})

test_that("works on solve_composition_poly() output", {
  X <- matrix(c(0.625, 0.5, 0.5, 0.5, 0.5, 0.5, 0.75, 0.5, 0.625, 0.625),
              nrow = 5, byrow = TRUE,
              dimnames = list(paste0("SNP", 1:5), c("VarA", "VarB")))
  Y <- matrix(c(2, 1, 2, 3, 4, 3, 4, 2, 3, 0), nrow = 2, byrow = TRUE,
              dimnames = list(paste0("Test", 1:2), paste0("SNP", 1:5)))
  out <- assign_composition_line(solve_composition_poly(Y, X, ploidy = 4))
  expect_equal(nrow(out), 2)
  expect_true(all(out$`Predicted line` %in% c("VarA", "VarB", "Undetermined")))
})
