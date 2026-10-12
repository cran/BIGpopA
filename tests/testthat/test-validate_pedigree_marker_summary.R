# tests/testthat/test-validate_pedigree_marker_summary.R
# Independent check of validate_pedigree()$marker_summary against a brute-force
# per-marker oracle (gamete-range rule), on simulated trios with known errors.

library(testthat)
library(data.table)

ms_sim <- function(ploidy, nm = 200, seed = 7) {
  set.seed(seed + ploidy)
  p    <- runif(nm, 0.2, 0.8)
  draw <- function() rbinom(nm, ploidy, p)
  gam  <- function(g) rhyper(length(g), g, ploidy - g, ploidy / 2)
  G <- list()
  for (id in c("M1", "M2", "M3", "F1", "F2", "F3")) G[[id]] <- draw()
  cross <- function(m, f) gam(G[[m]]) + gam(G[[f]])
  G[["C1"]] <- cross("M1", "F1")
  c2 <- cross("M2", "F2"); c2[sample(nm, 6)] <- NA          # missing calls
  G[["C2"]] <- c2
  c3 <- cross("M3", "F3"); idx <- sample(nm, 10); c3[idx] <- draw()[idx]
  G[["C3"]] <- c3                                           # genotyping errors
  G[["C4"]] <- cross("M1", "F2")                            # clean
  G[["C5"]] <- cross("M2", "F3")                            # parent F3 genotyped
  mat <- do.call(rbind, G)
  colnames(mat) <- paste0("m", seq_len(nm))
  geno <- data.table(id = rownames(mat), as.data.table(mat))
  ped  <- data.table(
    id            = c("M1", "M2", "M3", "F1", "F2", "F3", "C1", "C2", "C3", "C4", "C5", "C6", "C7"),
    male_parent   = c(NA, NA, NA, NA, NA, NA, "M1", "M2", "M3", "M1", "M2", "M1", "M1"),
    female_parent = c(NA, NA, NA, NA, NA, NA, "F1", "F2", "F3", "F2", "F3", NA, "FX"))
  # C6: one parent unknown; C7: female parent "FX" has no genotype -> neither is a full trio
  list(geno = geno, ped = ped, GM = mat, ploidy = ploidy)
}

ms_oracle <- function(pop) {
  ped <- as.data.frame(pop$ped); GM <- pop$GM; h <- pop$ploidy / 2
  keep <- !is.na(ped$male_parent) & !is.na(ped$female_parent) &
    ped$id %in% rownames(GM) & ped$male_parent %in% rownames(GM) &
    ped$female_parent %in% rownames(GM)
  ped <- ped[keep, , drop = FALSE]
  tested <- mism <- numeric(ncol(GM))
  for (i in seq_len(nrow(ped))) {
    m <- GM[ped$male_parent[i], ]; f <- GM[ped$female_parent[i], ]; o <- GM[ped$id[i], ]
    ok  <- !is.na(m) & !is.na(f) & !is.na(o)
    lo  <- pmax(0, m - h) + pmax(0, f - h)
    hi  <- pmin(m, h) + pmin(f, h)
    bad <- ok & (o < lo | o > hi)
    tested <- tested + ok
    mism   <- mism + bad
  }
  data.frame(marker = colnames(GM), trios_tested = tested, trios_mismatch = mism,
             stringsAsFactors = FALSE)
}

run_ms <- function(pop)
  suppressWarnings(validate_pedigree(pop$ped, pop$geno, min_markers = 10,
                                     verbose = FALSE, plot_results = FALSE,
                                     ploidy = pop$ploidy))

for (pl in c(2, 4)) {
  test_that(paste0("marker_summary counts match oracle (ploidy ", pl, ")"), {
    pop <- ms_sim(pl)
    res <- run_ms(pop)
    ms  <- as.data.frame(res$marker_summary)
    ex  <- ms_oracle(pop)

    expect_true(all(c("marker", "trios_tested", "trios_mismatch", "mismatch_pct") %in% names(ms)))
    expect_equal(nrow(ms), ncol(pop$GM))
    expect_setequal(ms$marker, ex$marker)

    m <- merge(ms, ex, by = "marker", suffixes = c("", ".ex"))
    expect_equal(as.numeric(m$trios_tested),   as.numeric(m$trios_tested.ex))
    expect_equal(as.numeric(m$trios_mismatch), as.numeric(m$trios_mismatch.ex))

    exp_pct <- ifelse(m$trios_tested.ex > 0,
                      round(m$trios_mismatch.ex / m$trios_tested.ex * 100, 2), NA)
    expect_equal(m$mismatch_pct, exp_pct)

    # mismatches never exceed tested; at most the number of full trios
    expect_true(all(ms$trios_mismatch <= ms$trios_tested))
    expect_true(max(ms$trios_tested) <= 5)   # C1-C5 only
  })

  test_that(paste0("marker_summary is sorted by mismatch_pct desc (ploidy ", pl, ")"), {
    ms <- as.data.frame(run_ms(ms_sim(pl))$marker_summary)
    pct <- ms$mismatch_pct
    pct[is.na(pct)] <- -Inf
    expect_false(is.unsorted(rev(pct)))
    # ties on pct broken by mismatch count desc
    same <- which(diff(pct) == 0)
    expect_true(all(diff(ms$trios_mismatch)[same] <= 0))
  })
}

test_that("a clean population gives zero mismatches everywhere", {
  pop <- ms_sim(2)
  pop$ped <- pop$ped[is.na(male_parent) | id %in% c("C1", "C4")]
  ms <- as.data.frame(run_ms(pop)$marker_summary)
  expect_equal(sum(ms$trios_mismatch), 0)
  expect_equal(sum(ms_oracle(pop)$trios_mismatch), 0)
})

test_that("a known injected error is attributed to the right marker", {
  pop <- ms_sim(2)
  g <- data.table::copy(pop$geno)
  g[id == "M1", m1 := 0L]; g[id == "F1", m1 := 0L]; g[id == "C1", m1 := 2L]
  pop$geno <- g
  pop$GM["M1", "m1"] <- 0; pop$GM["F1", "m1"] <- 0; pop$GM["C1", "m1"] <- 2
  pop$ped <- pop$ped[is.na(male_parent) | id == "C1"]
  ms <- as.data.frame(run_ms(pop)$marker_summary)
  expect_equal(ms$trios_mismatch[ms$marker == "m1"], 1)
  expect_equal(ms$trios_tested[ms$marker == "m1"], 1)
  expect_equal(ms$mismatch_pct[ms$marker == "m1"], 100)
})
