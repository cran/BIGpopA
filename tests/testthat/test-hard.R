# tests/testthat/test-hard.R
# Randomised, oracle-based and invariance tests. The oracles below are written
# independently of the package code (explicit gamete enumeration, brute force
# over all candidate pairs), so a pass means the package agrees with the
# definition of the test, not with itself.

library(testthat)
library(data.table)

# ==============================================================================
# Helpers
# ==============================================================================

# gamete allele-B counts a parent with dosage g can transmit at ploidy P,
# found by explicit enumeration of k alleles drawn from g B and (P - g) A
oracle_gametes <- function(g, P) {
  h <- P / 2
  Filter(function(k) k <= g && (h - k) <= (P - g), 0:h)
}

# TRUE if offspring dosage o is possible from parent dosages gm and gf
oracle_valid <- function(o, gm, gf, P) {
  km <- oracle_gametes(gm, P)
  kf <- oracle_gametes(gf, P)
  any(outer(km, kf, "+") == o)
}

# returns c(tested, mismatches) for one trio over all markers (NA = skip)
oracle_trio <- function(o, m, f, P) {
  ok <- !is.na(o) & !is.na(m) & !is.na(f)
  bad <- vapply(which(ok), function(j) !oracle_valid(o[j], m[j], f[j], P), logical(1))
  c(tested = sum(ok), mismatches = sum(bad))
}

# simulate a child dosage from two parents (random meiosis)
sim_child <- function(m, f, P) {
  gam <- function(g) sum(sample(rep(c(1L, 0L), c(g, P - g)), P / 2))
  vapply(seq_along(m), function(j) gam(m[j]) + gam(f[j]), numeric(1))
}

# random panel of parents, valid and random children, with optional NA
make_panel <- function(P, n_markers = 60, n_trios = 24, na_rate = 0.1, seed = 1) {
  set.seed(seed)
  mk <- paste0("S", seq_len(n_markers))
  par_ids <- paste0("P", seq_len(2 * n_trios))
  par <- matrix(sample(0:P, length(par_ids) * n_markers, replace = TRUE),
                nrow = length(par_ids), dimnames = list(par_ids, mk))
  kids <- matrix(NA_real_, nrow = n_trios, ncol = n_markers,
                 dimnames = list(paste0("C", seq_len(n_trios)), mk))
  ped <- data.frame(id = rownames(kids),
                    male_parent = par_ids[seq(1, 2 * n_trios, by = 2)],
                    female_parent = par_ids[seq(2, 2 * n_trios, by = 2)],
                    stringsAsFactors = FALSE)
  for (i in seq_len(n_trios)) {
    m <- par[ped$male_parent[i], ]
    f <- par[ped$female_parent[i], ]
    kids[i, ] <- if (i %% 2 == 0) sim_child(m, f, P)            # valid trio
                 else sample(0:P, n_markers, replace = TRUE)    # random child
  }
  all_g <- rbind(par, kids)
  if (na_rate > 0)
    all_g[sample(length(all_g), round(na_rate * length(all_g)))] <- NA
  geno <- data.frame(id = rownames(all_g), all_g, check.names = FALSE,
                     stringsAsFactors = FALSE, row.names = NULL)
  list(geno = geno, ped = ped, mat = all_g)
}

quiet_validate <- function(ped, geno, ...) {
  validate_pedigree(ped, geno, verbose = FALSE, plot_results = FALSE, ...)
}

# ==============================================================================
# 1. validate_pedigree agrees with brute-force Mendelian oracle
# ==============================================================================

for (P in c(2, 4, 6, 8)) {
  test_that(sprintf("validate_pedigree matches the gamete oracle (ploidy %d, with NA)", P), {
    d   <- make_panel(P, seed = 100 + P)
    out <- quiet_validate(d$ped, d$geno, ploidy = P, min_markers = 3)$full_results
    out <- as.data.frame(out)

    for (i in seq_len(nrow(d$ped))) {
      o  <- oracle_trio(d$mat[d$ped$id[i], ],
                        d$mat[d$ped$male_parent[i], ],
                        d$mat[d$ped$female_parent[i], ], P)
      r  <- out[out$id == d$ped$id[i], ]
      expect_equal(r$trio_markers_tested, unname(o["tested"]))
      expect_equal(r$trio_mendelian_error_pct,
                   round(100 * unname(o["mismatches"]) / unname(o["tested"]), 2),
                   tolerance = 0.011)
      expect_equal(r$status,
                   if (100 * o["mismatches"] / o["tested"] <= 5) "pass" else "fail")
    }
  })
}

test_that("valid simulated trios never fail (no false positives, no NA)", {
  for (P in c(2, 4, 6, 8)) {
    d   <- make_panel(P, na_rate = 0, seed = 7 + P)
    out <- as.data.frame(quiet_validate(d$ped, d$geno, ploidy = P)$full_results)
    valid_ids <- d$ped$id[seq(2, nrow(d$ped), by = 2)]
    expect_true(all(out$trio_mendelian_error_pct[out$id %in% valid_ids] == 0),
                info = paste("ploidy", P))
  }
})

# ==============================================================================
# 2. Odd ploidy (hand-computed)
# ==============================================================================

test_that("triploid trio uses homozygous opposite-allele test only", {
  geno <- data.frame(
    id = c("M", "F", "C"),
    M1 = c(3L, 0L, 0L),   # male 3 vs offspring 0 -> testable, mismatch
    M2 = c(0L, 0L, 0L),   # testable, match
    M3 = c(3L, 1L, 3L),   # male homozygous same as offspring -> match
    M4 = c(1L, 2L, 1L)    # offspring heterozygous -> untestable
  )
  ped <- data.frame(id = "C", male_parent = "M", female_parent = "F")
  r <- as.data.frame(quiet_validate(ped, geno, ploidy = 3, min_markers = 1)$full_results)
  expect_equal(r$trio_markers_tested, 3)
  expect_equal(r$trio_mendelian_error_pct, round(100 / 3, 2), tolerance = 0.011)
})

# ==============================================================================
# 3. Invariance: marker order, trio order and input type do not matter
# ==============================================================================

test_that("results do not depend on marker order or trio order", {
  d    <- make_panel(4, seed = 33)
  base <- as.data.frame(quiet_validate(d$ped, d$geno, ploidy = 4)$full_results)
  set.seed(5)
  geno2 <- d$geno[sample(nrow(d$geno)), c("id", sample(setdiff(names(d$geno), "id")))]
  ped2  <- d$ped[sample(nrow(d$ped)), ]
  perm  <- as.data.frame(quiet_validate(ped2, geno2, ploidy = 4)$full_results)
  cols  <- c("id", "trio_mendelian_error_pct", "trio_markers_tested", "status")
  expect_equal(base[order(base$id), cols], perm[order(perm$id), cols],
               ignore_attr = TRUE)
})

test_that("file input and in-memory input give identical results", {
  d  <- make_panel(4, seed = 34)
  fg <- tempfile(fileext = ".txt"); fp <- tempfile(fileext = ".txt")
  fwrite(d$geno, fg, sep = "\t"); fwrite(d$ped, fp, sep = "\t")
  a <- as.data.frame(quiet_validate(d$ped, d$geno, ploidy = 4)$full_results)
  b <- as.data.frame(quiet_validate(fp, fg, ploidy = 4)$full_results)
  expect_equal(a, b, ignore_attr = TRUE)
})

test_that("founders as file, vector and data.frame give the same result", {
  d <- make_panel(2, seed = 35)
  ped <- rbind(d$ped, data.frame(id = c("P1", "P2"), male_parent = "0", female_parent = "0"))
  ff <- tempfile(fileext = ".txt"); writeLines(c("P1", "P2"), ff)
  a <- quiet_validate(ped, d$geno, founders_file = ff)$full_results
  b <- quiet_validate(ped, d$geno, founders_file = c("P1", "P2"))$full_results
  cc <- quiet_validate(ped, d$geno, founders_file = data.frame(id = c("P1", "P2")))$full_results
  expect_equal(a, b)
  expect_equal(a, cc)
  expect_equal(sum(a$status == "founders"), 2)
})

# ==============================================================================
# 4. marker_trio_table consistent with the oracle
# ==============================================================================

test_that("marker_trio_table mismatch and missing counts agree with the oracle", {
  P <- 4
  d <- make_panel(P, seed = 41)
  out <- quiet_validate(d$ped, d$geno, ploidy = P, marker_trio_table = TRUE)
  tab <- as.data.frame(out$marker_trio_table)
  expect_equal(tab$marker, setdiff(names(d$geno), "id"))

  for (i in seq_len(nrow(d$ped))) {
    col <- paste(d$ped$id[i], d$ped$male_parent[i], d$ped$female_parent[i], sep = "_")
    lab <- tab[[col]]
    o  <- d$mat[d$ped$id[i], ]
    m  <- d$mat[d$ped$male_parent[i], ]
    f  <- d$mat[d$ped$female_parent[i], ]
    any_na <- unname(is.na(o) | is.na(m) | is.na(f))
    expect_equal(grepl("^missing", lab), any_na)
    expect_equal(sum(grepl("^mismatch", lab)),
                 unname(oracle_trio(o, m, f, P)["mismatches"]))
    expect_equal(sum(lab == "match"),
                 unname(oracle_trio(o, m, f, P)["tested"] -
                        oracle_trio(o, m, f, P)["mismatches"]))
  }
})

# ==============================================================================
# 5. marker_summary consistent with the trio results
# ==============================================================================

test_that("marker_summary totals equal the per-trio totals", {
  d <- make_panel(4, seed = 51)
  out <- quiet_validate(d$ped, d$geno, ploidy = 4)
  ms  <- as.data.frame(out$marker_summary)
  fr  <- as.data.frame(out$full_results)
  tested_trios <- fr[fr$trio_markers_tested > 0, ]
  expect_equal(sum(ms$trios_tested), sum(tested_trios$trio_markers_tested))
  expect_equal(sum(ms$trios_mismatch),
               sum(round(tested_trios$trio_markers_tested *
                         tested_trios$trio_mendelian_error_pct / 100)))
  expect_true(all(ms$trios_mismatch <= ms$trios_tested))
  expect_false(is.unsorted(rev(ms$mismatch_pct[!is.na(ms$mismatch_pct)])))
})

# ==============================================================================
# 6. find_parentage best_pair against brute force
# ==============================================================================

test_that("best_pair error equals the brute-force minimum over all pairs", {
  skip_on_cran()
  P <- 4
  set.seed(61)
  n_m <- 50; n_male <- 5; n_fem <- 5
  mk <- paste0("S", seq_len(n_m))
  males <- matrix(sample(0:P, n_male * n_m, TRUE), n_male,
                  dimnames = list(paste0("MA", 1:n_male), mk))
  fems  <- matrix(sample(0:P, n_fem * n_m, TRUE), n_fem,
                  dimnames = list(paste0("FE", 1:n_fem), mk))
  true_pairs <- data.frame(id = paste0("K", 1:6),
                           m = rep(rownames(males), length.out = 6)[c(1, 2, 3, 4, 5, 1)],
                           f = rep(rownames(fems),  length.out = 6)[c(2, 3, 4, 5, 1, 4)],
                           stringsAsFactors = FALSE)
  kids <- t(vapply(seq_len(nrow(true_pairs)), function(i)
    sim_child(males[true_pairs$m[i], ], fems[true_pairs$f[i], ], P), numeric(n_m)))
  dimnames(kids) <- list(true_pairs$id, mk)
  all_g <- rbind(males, fems, kids)
  all_g[sample(length(all_g), 15)] <- NA                       # a few missing calls
  geno <- data.frame(id = rownames(all_g), all_g, check.names = FALSE, row.names = NULL)
  parents <- data.frame(id = c(rownames(males), rownames(fems)),
                        sex = rep(c("M", "F"), c(n_male, n_fem)))
  progeny <- data.frame(id = true_pairs$id)

  res <- as.data.frame(find_parentage(geno, parents, progeny, method = "best_pair",
                                      ploidy = P, min_markers = 3, verbose = FALSE,
                                      plot_results = FALSE)$full_results)

  for (k in true_pairs$id) {
    errs <- outer(rownames(males), rownames(fems), Vectorize(function(a, b) {
      o <- oracle_trio(all_g[k, ], all_g[a, ], all_g[b, ], P)
      100 * o["mismatches"] / o["tested"]
    }))
    expect_equal(as.numeric(res$mendelian_error_pct[res$id == k]), round(min(errs), 2),
                 tolerance = 0.011, info = k)
    expect_equal(res$tied_candidates[res$id == k], sum(abs(errs - min(errs)) < 1e-9),
                 info = k)
  }
  # the simulated pair always has the lowest error unless another pair ties
  unique_best <- res$tied_candidates == 1
  expect_equal(res$male_parent[unique_best],   true_pairs$m[match(res$id[unique_best], true_pairs$id)])
  expect_equal(res$female_parent[unique_best], true_pairs$f[match(res$id[unique_best], true_pairs$id)])
})

test_that("top_n = 3 returns non-decreasing errors per progeny", {
  P <- 4
  d <- make_panel(P, na_rate = 0, seed = 71)
  pid <- grep("^P", d$geno$id, value = TRUE)
  parents <- data.frame(id = pid, sex = rep(c("M", "F"), length.out = length(pid)))
  progeny <- data.frame(id = d$ped$id[1:4])
  r <- as.data.frame(find_parentage(d$geno, parents, progeny, method = "best_pair",
                                    ploidy = P, top_n = 3, verbose = FALSE,
                                    plot_results = FALSE)$full_results)
  e <- r[, c("mendelian_error_pct", "mendelian_error_pct_2", "mendelian_error_pct_3")]
  e[] <- lapply(e, as.numeric)
  expect_true(all(e[, 1] <= e[, 2] & e[, 2] <= e[, 3]))
})

# ==============================================================================
# 7. check_ped invariants on a larger random pedigree
# ==============================================================================

make_big_ped <- function(n = 400, seed = 81) {
  set.seed(seed)
  ids <- sprintf("I%04d", seq_len(n))
  sires <- ids[seq(1, n, by = 2)]            # fixed sex roles
  dams  <- ids[seq(2, n, by = 2)]
  ped <- data.frame(id = ids, male_parent = "0", female_parent = "0",
                    stringsAsFactors = FALSE)
  for (i in 21:n) {                          # first 20 are founders; parents come earlier
    ped$male_parent[i]   <- sample(sires[sires %in% ids[seq_len(i - 1)]], 1)
    ped$female_parent[i] <- sample(dams[dams %in% ids[seq_len(i - 1)]], 1)
  }
  ped
}

test_that("check_ped on a clean random pedigree changes nothing", {
  ped <- make_big_ped()
  out <- check_ped(ped, verbose = FALSE)$corrected_pedigree
  key <- function(x) sort(paste(x$id, x$male_parent, x$female_parent))
  expect_equal(key(out), key(ped))
})

test_that("check_ped removes injected duplicates and is idempotent", {
  ped <- make_big_ped()
  dup <- ped[sample(nrow(ped), 40), ]
  messy <- rbind(ped, dup)[sample(nrow(ped) + 40), ]
  once  <- check_ped(messy, verbose = FALSE)$corrected_pedigree
  twice <- check_ped(once,  verbose = FALSE)$corrected_pedigree
  key <- function(x) sort(paste(x$id, x$male_parent, x$female_parent))
  expect_equal(nrow(once), nrow(ped))
  expect_false(any(duplicated(once)))
  expect_equal(key(once), key(ped))
  expect_equal(key(twice), key(once))
})

test_that("every non-zero parent has its own row after check_ped", {
  ped <- make_big_ped()
  ped$male_parent[ped$id == "I0100"] <- "I0003"   # make sure I0003 is used as a parent
  ped <- ped[ped$id != "I0003", ]                 # then drop its own row
  out <- check_ped(ped, verbose = FALSE)$corrected_pedigree
  parents <- setdiff(c(out$male_parent, out$female_parent), "0")
  expect_true(all(parents %in% out$id))
  expect_true("I0003" %in% out$id)
})

test_that("a cycle inside a large pedigree is reported and nothing is dropped", {
  ped <- make_big_ped()
  ped$male_parent[ped$id == "I0399"] <- "I0001"   # two sires that are each
  ped$male_parent[ped$id == "I0001"] <- "I0399"   # other's sire: a cycle
  res <- check_ped(ped, verbose = FALSE)
  expect_gt(nrow(res$dependencies), 0)
  expect_true(all(ped$id %in% res$corrected_pedigree$id))
})

# ==============================================================================
# 8. Composition: exact recovery of known mixtures
# ==============================================================================

test_that("solve_composition_poly recovers noise-free mixtures", {
  set.seed(91)
  P  <- 4
  nm <- 80
  # lines are fixed for one allele at every marker: line frequencies are 0 or 1
  a  <- sample(0:1, nm, replace = TRUE)
  X  <- cbind(LineA = a, LineB = 1 - a)
  rownames(X) <- paste0("S", seq_len(nm))
  w  <- rbind(c(1, 0), c(0, 1), c(0.5, 0.5), c(0.25, 0.75), c(0.75, 0.25))
  Y  <- t(apply(w, 1, function(wi) P * as.vector(X %*% wi)))   # integer dosages
  dimnames(Y) <- list(paste0("A", 1:5), rownames(X))
  est <- suppressWarnings(solve_composition_poly(Y, X, ploidy = P))  # R2 warns when sd = 0
  expect_equal(unname(est[, 1:2]), w, tolerance = 1e-3)
  expect_equal(unname(rowSums(est[, 1:2])), rep(1, 5), tolerance = 1e-8)
})

test_that("allele_freq_poly matches colMeans / ploidy for random data", {
  set.seed(92)
  P <- 6
  g <- matrix(sample(0:P, 30 * 40, TRUE), 30,
              dimnames = list(paste0("I", 1:30), paste0("S", 1:40)))
  pops <- list(One = paste0("I", 1:10), Two = paste0("I", 11:30))
  f <- allele_freq_poly(g, pops, ploidy = P)
  expect_equal(unname(f[, "One"]), unname(colMeans(g[1:10, ]) / P))
  expect_equal(unname(f[, "Two"]), unname(colMeans(g[11:30, ]) / P))
  expect_true(all(f >= 0 & f <= 1))
})

# ==============================================================================
# 9. Scale
# ==============================================================================

test_that("validate_pedigree handles 300 trios x 3000 markers (ploidy 4) in time", {
  skip_on_cran()
  d <- make_panel(4, n_markers = 3000, n_trios = 300, na_rate = 0.02, seed = 99)
  t0  <- Sys.time()
  out <- quiet_validate(d$ped, d$geno, ploidy = 4)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  expect_equal(nrow(out$full_results), 300)
  expect_equal(nrow(out$marker_summary), 3000)
  expect_lt(secs, 120)
})
