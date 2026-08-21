# Tests for arbitrary-ploidy support in find_parentage() and validate_pedigree().
# Even ploidy uses the polysomic gamete-range Mendelian test; odd ploidy falls
# back to the homozygosity-only opposite-homozygote exclusion.

# ---- helper -----------------------------------------------------------------
make_geno <- function(ids, rows, n_markers) {
  mk    <- paste0("S", seq_len(n_markers))
  g     <- data.frame(id = ids, stringsAsFactors = FALSE)
  g[mk] <- do.call(rbind, rows)
  g
}

# ---- even ploidy: tetraploid (ploidy = 4) -----------------------------------
test_that("find_parentage best_pair handles tetraploids (ploidy = 4)", {
  n <- 12
  geno <- make_geno(
    ids  = c("P1", "P2", "OffTrue", "OffFalse"),
    rows = list(rep(0L, n),   # P1 homozygous A
                rep(4L, n),   # P2 homozygous B
                rep(2L, n),   # OffTrue : 0 x 4 can only produce dosage 2
                rep(0L, n)),  # OffFalse: dosage 0 is impossible from 0 x 4
    n_markers = n
  )
  parents <- data.frame(id = c("P1", "P2"), sex = c("M", "F"),
                        stringsAsFactors = FALSE)
  progeny <- data.frame(id = c("OffTrue", "OffFalse"), stringsAsFactors = FALSE)

  res <- find_parentage(geno, parents, progeny, method = "best_pair",
                        ploidy = 4, verbose = FALSE, plot_results = FALSE)
  fr <- as.data.frame(res$full_results)

  expect_equal(fr$status[fr$id == "OffTrue"], "pass")
  expect_equal(as.numeric(fr$mendelian_error_pct[fr$id == "OffTrue"]), 0)
  expect_equal(fr$status[fr$id == "OffFalse"], "high_error")
  expect_equal(as.numeric(fr$mendelian_error_pct[fr$id == "OffFalse"]), 100)
})

test_that("validate_pedigree handles tetraploids (ploidy = 4)", {
  n <- 12
  geno <- make_geno(
    ids  = c("P1", "P2", "OffTrue", "OffFalse"),
    rows = list(rep(0L, n), rep(4L, n), rep(2L, n), rep(0L, n)),
    n_markers = n
  )
  ped <- data.frame(
    id            = c("OffTrue", "OffFalse"),
    male_parent   = c("P1", "P1"),
    female_parent = c("P2", "P2"),
    stringsAsFactors = FALSE
  )
  res <- validate_pedigree(ped, geno, ploidy = 4,
                           verbose = FALSE, plot_results = FALSE)
  fr <- as.data.frame(res$full_results)

  expect_equal(fr$status[fr$id == "OffTrue"], "pass")
  expect_equal(fr$trio_mendelian_error_pct[fr$id == "OffTrue"], 0)
  expect_equal(fr$status[fr$id == "OffFalse"], "fail")
})

# ---- odd ploidy: triploid (ploidy = 3, homozygosity-only) -------------------
test_that("validate_pedigree uses opposite-homozygote exclusion for odd ploidy", {
  n <- 12
  geno <- make_geno(
    ids  = c("Sire", "Dam", "OffOK", "OffBad"),
    rows = list(rep(0L, n),   # Sire  homozygous A
                rep(0L, n),   # Dam   homozygous A
                rep(0L, n),   # OffOK  homozygous A  -> consistent
                rep(3L, n)),  # OffBad homozygous B  -> opposite homozygote
    n_markers = n
  )
  ped <- data.frame(
    id            = c("OffOK", "OffBad"),
    male_parent   = c("Sire", "Sire"),
    female_parent = c("Dam", "Dam"),
    stringsAsFactors = FALSE
  )
  res <- validate_pedigree(ped, geno, ploidy = 3,
                           verbose = FALSE, plot_results = FALSE)
  fr <- as.data.frame(res$full_results)

  expect_equal(fr$status[fr$id == "OffOK"], "pass")
  expect_equal(fr$trio_mendelian_error_pct[fr$id == "OffOK"], 0)
  expect_equal(fr$status[fr$id == "OffBad"], "fail")
  expect_equal(fr$trio_mendelian_error_pct[fr$id == "OffBad"], 100)
})

test_that("odd ploidy counts only homozygous-informative markers", {
  n   <- 12
  hom <- 3                                   # only 3 homozygous-informative markers
  sire <- c(rep(0L, hom), rep(1L, n - hom))  # rest heterozygous (dosage 1 / 2)
  dam  <- c(rep(0L, hom), rep(2L, n - hom))
  off  <- c(rep(0L, hom), rep(1L, n - hom))
  geno <- make_geno(c("Sire", "Dam", "Off"), list(sire, dam, off), n)
  ped  <- data.frame(id = "Off", male_parent = "Sire",
                     female_parent = "Dam", stringsAsFactors = FALSE)

  res <- validate_pedigree(ped, geno, ploidy = 3, min_markers = 10,
                           verbose = FALSE, plot_results = FALSE)
  fr <- as.data.frame(res$full_results)

  expect_equal(fr$trio_markers_tested[fr$id == "Off"], hom)  # heterozygous markers ignored
  expect_equal(fr$status[fr$id == "Off"], "low_markers")     # 3 < min_markers (10)
})

# ---- ploidy validation ------------------------------------------------------
test_that("ploidy must be an integer >= 2", {
  geno <- make_geno(c("A", "B", "C"),
                    list(c(0L, 2L, 0L), c(2L, 0L, 2L), c(1L, 1L, 1L)), 3)
  ped     <- data.frame(id = "C", male_parent = "A", female_parent = "B",
                        stringsAsFactors = FALSE)
  parents <- data.frame(id = c("A", "B"), sex = c("M", "F"),
                        stringsAsFactors = FALSE)
  progeny <- data.frame(id = "C", stringsAsFactors = FALSE)

  expect_error(validate_pedigree(ped, geno, ploidy = 1,
                                 verbose = FALSE, plot_results = FALSE),
               "ploidy must be an integer")
  expect_error(validate_pedigree(ped, geno, ploidy = 2.5,
                                 verbose = FALSE, plot_results = FALSE),
               "ploidy must be an integer")
  expect_error(find_parentage(geno, parents, progeny, ploidy = 0,
                              verbose = FALSE, plot_results = FALSE),
               "ploidy must be an integer")
})

# ---- diploid default unchanged ----------------------------------------------
test_that("default ploidy = 2 reproduces explicit diploid results", {
  n <- 12
  geno <- make_geno(
    ids  = c("P1", "P2", "Off"),
    rows = list(rep(0L, n), rep(2L, n), rep(1L, n)),  # 0 x 2 -> offspring must be 1
    n_markers = n
  )
  ped <- data.frame(id = "Off", male_parent = "P1", female_parent = "P2",
                    stringsAsFactors = FALSE)

  default2  <- validate_pedigree(ped, geno, verbose = FALSE, plot_results = FALSE)
  explicit2 <- validate_pedigree(ped, geno, ploidy = 2,
                                 verbose = FALSE, plot_results = FALSE)

  expect_equal(as.data.frame(default2$full_results),
               as.data.frame(explicit2$full_results))
  expect_equal(as.data.frame(default2$full_results)$status[1], "pass")
})
