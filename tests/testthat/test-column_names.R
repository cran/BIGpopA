#### Column names are matched ignoring case ####

base::set.seed(7)
n_snp <- 30
p1 <- base::sample(0:2, n_snp, replace = TRUE)
p2 <- base::sample(0:2, n_snp, replace = TRUE)
gamete <- function(g) base::ifelse(g == 1, base::sample(0:1, length(g), replace = TRUE), g / 2)
geno <- base::data.frame(
  id = c("P1", "P2", "Off1"),
  base::rbind(p1, p2, gamete(p1) + gamete(p2)),
  row.names = NULL
)
base::names(geno)[-1] <- base::paste0("SNP", base::seq_len(n_snp))

ped_lower <- base::data.frame(id = "Off1", male_parent = "P1", female_parent = "P2")
ped_mixed <- stats::setNames(ped_lower, c("ID", "Male_Parent", "FEMALE PARENT"))

test_that(".standardize_names renames only the requested columns", {
  df  <- base::data.frame(ID = 1, Sex = "M", SNP_A = 0)
  out <- BIGpopA:::.standardize_names(df, c("id", "sex"))
  expect_identical(base::names(out), c("id", "sex", "SNP_A"))

  dt  <- data.table::as.data.table(df)
  out <- BIGpopA:::.standardize_names(dt, "id")
  expect_identical(base::names(out)[1], "id")
  expect_identical(base::names(dt)[1],  "ID")    # caller's data.table untouched
})

test_that("validate_pedigree accepts pedigree and genotype columns in any case", {
  geno_upper <- geno
  base::names(geno_upper)[1] <- "Id"
  run <- function(p, g) validate_pedigree(p, g, verbose = FALSE,
                                          plot_results = FALSE)$full_results
  expect_equal(run(ped_mixed, geno_upper), run(ped_lower, geno))
})

test_that("find_parentage accepts parents and progeny columns in any case", {
  parents_lower <- base::data.frame(id = c("P1", "P2"), sex = c("M", "F"))
  parents_mixed <- stats::setNames(parents_lower, c("ID", "Sex"))
  run <- function(par, prog) find_parentage(geno, par, prog, verbose = FALSE,
                                            plot_results = FALSE)$full_results
  expect_equal(run(parents_mixed, base::data.frame(ID = "Off1")),
               run(parents_lower, base::data.frame(id = "Off1")))
})

test_that("check_ped accepts column names in any case", {
  res <- check_ped(ped_mixed, verbose = FALSE)
  expect_true(base::all(c("id", "male_parent", "female_parent") %in%
                          base::names(res$corrected_pedigree)))
})
