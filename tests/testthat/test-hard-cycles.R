# tests/testthat/test-hard-cycles.R
# Hard tests for the cycle check of check_ped(). The oracle below is written
# independently of the package code (breadth-first search of every
# individual's ancestors), so a pass means check_ped() agrees with the
# definition of a cycle (an individual that is its own ancestor), not with
# itself. Cycles are judged on the pedigree after the corrections, so the
# oracle always reads out$corrected_pedigree.

# ==============================================================================
# Helpers
# ==============================================================================

# Ancestors of every id by breadth-first search over the recorded parents
# (self-parenting is ignored: check_ped() reports it with its own message).
# Groups are ids that are each other's ancestors.
oracle_cycles <- function(ped) {
  ids     <- unique(as.character(ped$id))
  parents <- lapply(ids, function(i) {
    p <- unique(c(ped$male_parent[ped$id == i], ped$female_parent[ped$id == i]))
    setdiff(p[p %in% ids], i)
  })
  names(parents) <- ids

  anc <- lapply(ids, function(i) {
    seen  <- character(0)
    queue <- parents[[i]]
    while (length(queue) > 0) {
      x     <- queue[1]
      queue <- queue[-1]
      if (x %in% seen) next
      seen  <- c(seen, x)
      queue <- c(queue, parents[[x]])
    }
    seen
  })
  names(anc) <- ids

  in_cycle <- ids[vapply(ids, function(i) i %in% anc[[i]], logical(1))]
  groups   <- list()
  left     <- in_cycle
  while (length(left) > 0) {
    i <- left[1]
    g <- left[vapply(left, function(j) j %in% anc[[i]] && i %in% anc[[j]], logical(1))]
    groups[[length(groups) + 1]] <- sort(unique(c(i, g)))
    left <- setdiff(left, g)
  }
  list(parents = parents, groups = groups)
}

# Members, loop path and "more than one loop" flag of each cycle row
parse_cycles <- function(dep) {
  body  <- sub("^Cycle detected involving IDs: ", "",
               dep[startsWith(dep, "Cycle detected involving IDs: ")])
  multi <- grepl(" (more than one loop, e.g. ", body, fixed = TRUE)
  lapply(seq_along(body), function(k) {
    if (multi[k]) {
      members <- strsplit(sub(" \\(more than one loop.*$", "", body[k]), ", ", fixed = TRUE)[[1]]
      loop    <- strsplit(sub("^.*e\\.g\\. (.*)\\)$", "\\1", body[k]), " -> ", fixed = TRUE)[[1]]
    } else {
      loop    <- strsplit(body[k], " -> ", fixed = TRUE)[[1]]
      members <- unique(loop)
    }
    list(members = sort(members), loop = loop, multi = multi[k])
  })
}

member_sets <- function(out) {
  vapply(parse_cycles(out$dependencies$dependency),
         function(f) paste(f$members, collapse = " "), character(1))
}

# A real loop: closed, and each id is a recorded parent of the id before it
is_loop <- function(loop, parents) {
  length(loop) >= 3 && loop[1] == loop[length(loop)] &&
    all(vapply(seq_len(length(loop) - 1),
               function(k) loop[k + 1] %in% parents[[loop[k]]], logical(1)))
}

# check_ped() output against the oracle
expect_cycles_match <- function(out, info = "") {
  orc   <- oracle_cycles(out$corrected_pedigree)
  found <- parse_cycles(out$dependencies$dependency)

  # one row per group, and the same groups
  expect_equal(length(found), length(orc$groups), info = info)
  expect_true(setequal(member_sets(out),
                       vapply(orc$groups, paste, character(1), collapse = " ")), info = info)
  expect_false(any(duplicated(out$dependencies$dependency)), info = info)

  for (f in found) {
    # the path is a real loop inside its group
    expect_true(is_loop(f$loop, orc$parents), info = info)
    expect_true(all(f$loop %in% f$members), info = info)
    # "more than one loop" exactly when the group has more parent links than ids
    n_links <- sum(vapply(f$members, function(m) sum(orc$parents[[m]] %in% f$members), numeric(1)))
    expect_equal(f$multi, n_links > length(f$members), info = info)
    # a single loop visits every member once
    if (!f$multi) {
      expect_true(setequal(f$loop[-length(f$loop)], f$members), info = info)
      expect_false(any(duplicated(f$loop[-length(f$loop)])), info = info)
    }
  }
}

# Random pedigree: parents are mostly earlier individuals (a valid pedigree),
# plus some parents drawn from anyone, which can close loops of any size and
# tie them together. Every parent has its own row and no one is its own parent.
random_ped <- function(n, back_rate, seed) {
  set.seed(seed)
  ids  <- sprintf("I%03d", seq_len(n))
  pick <- function(i) {
    r <- stats::runif(1)
    if (i <= 5 || r < 0.15) return("0")
    if (r < 0.15 + back_rate) return(sample(ids[-i], 1))
    sample(ids[seq_len(i - 1)], 1)
  }
  data.frame(id            = ids,
             male_parent   = vapply(seq_len(n), pick, character(1)),
             female_parent = vapply(seq_len(n), pick, character(1)),
             stringsAsFactors = FALSE)
}

no_corrections <- function(ped) {
  check_ped(ped, verbose = FALSE, correct_conflicting_trios = FALSE,
            correct_inconsistent_sex_roles = FALSE)
}

# ==============================================================================
# 1. Random pedigrees agree with the ancestor oracle
# ==============================================================================

test_that("cycle rows match the ancestor oracle on random pedigrees", {
  for (seed in 1:40) {
    back <- c(0.01, 0.03, 0.08, 0.30)[seed %% 4 + 1]   # from rare loops to one big tangle
    ped  <- random_ped(n = 150, back_rate = back, seed = seed)
    if (seed %% 3 == 0) {
      # repeated ids with other male parents: conflicting trios (and some duplicates)
      extra             <- ped[sample(nrow(ped), 10), ]
      extra$male_parent <- sample(ped$id, 10)
      ped               <- rbind(ped, extra)
    }
    expect_cycles_match(no_corrections(ped), info = paste("seed", seed, "no corrections"))
    expect_cycles_match(check_ped(ped, verbose = FALSE), info = paste("seed", seed, "defaults"))
  }
})

# ==============================================================================
# 2. Invariance
# ==============================================================================

test_that("shuffling the rows finds the same cycles", {
  for (seed in 1:15) {
    ped <- random_ped(n = 120, back_rate = 0.04, seed = seed)
    set.seed(1000 + seed)
    expect_setequal(member_sets(no_corrections(ped[sample(nrow(ped)), ])),
                    member_sets(no_corrections(ped)))
  }
})

test_that("renaming the individuals finds the same cycles", {
  for (seed in 1:15) {
    ped <- random_ped(n = 120, back_rate = 0.04, seed = seed)
    set.seed(2000 + seed)
    new_id      <- stats::setNames(paste0("Z", sample(99999, nrow(ped))), ped$id)
    new_id["0"] <- "0"
    renamed <- data.frame(id            = unname(new_id[ped$id]),
                          male_parent   = unname(new_id[ped$male_parent]),
                          female_parent = unname(new_id[ped$female_parent]),
                          stringsAsFactors = FALSE)
    expected <- vapply(strsplit(member_sets(no_corrections(ped)), " "),
                       function(m) paste(sort(unname(new_id[m])), collapse = " "), character(1))
    expect_setequal(member_sets(no_corrections(renamed)), expected)
  }
})

# ==============================================================================
# 3. One row per cycle, whatever surrounds it
# ==============================================================================

test_that("ancestors and descendants of a cycle add no rows", {
  set.seed(5)
  ring <- paste0("R", 1:5)
  anc  <- paste0("A", 1:10)
  desc <- sprintf("D%02d", 1:50)
  desc_male <- vapply(seq_along(desc), function(i) sample(c(ring, desc[seq_len(i - 1)]), 1),
                      character(1))
  ped <- data.frame(
    id            = c(anc, ring, desc),
    male_parent   = c(rep("0", 10), ring[c(2:5, 1)], desc_male),
    female_parent = c(rep("0", 10), sample(anc, 5), sample(anc, 50, replace = TRUE)),
    stringsAsFactors = FALSE
  )
  out <- check_ped(ped, verbose = FALSE)
  expect_identical(out$dependencies$dependency,
                   "Cycle detected involving IDs: R1 -> R2 -> R3 -> R4 -> R5 -> R1")
})

test_that("loops sharing individuals are one row whose example loop is real", {
  ped <- data.frame(
    id            = c("A", "B", "C", "D", "E"),
    male_parent   = c("B", "C", "A", "A", "D"),
    female_parent = c("E", "0", "0", "0", "0")
  )
  out <- check_ped(ped, verbose = FALSE)
  expect_identical(out$dependencies$dependency,
                   "Cycle detected involving IDs: A, B, C, D, E (more than one loop, e.g. A -> B -> C -> A)")
  expect_cycles_match(out)
})

test_that("a parent used as both sire and dam inside a loop is still one simple loop", {
  ped <- data.frame(id = c("X", "Y"), male_parent = c("Y", "X"), female_parent = c("Y", "0"))
  out <- check_ped(ped, verbose = FALSE, correct_inconsistent_sex_roles = FALSE)
  expect_identical(out$dependencies$dependency, "Cycle detected involving IDs: X -> Y -> X")
})

test_that("self-parenting is reported once, and a loop through that individual is still found", {
  ped <- data.frame(
    id            = c("S", "B", "C", "T"),
    male_parent   = c("S", "C", "B", "0"),
    female_parent = c("0", "B", "0", "T")
  )
  out <- no_corrections(ped)
  expect_setequal(out$dependencies$dependency, c(
    "Dependency: Individual S cannot be its own parent",
    "Dependency: Individual B cannot be its own parent",
    "Dependency: Individual T cannot be its own parent",
    "Cycle detected involving IDs: B -> C -> B"
  ))
})

# ==============================================================================
# 4. Cycles are judged after the corrections
# ==============================================================================

test_that("collapsing a conflicting trio can break a loop", {
  ped <- data.frame(
    id            = c("A",  "A",  "B",  "F1"),
    male_parent   = c("B",  "C",  "A",  "0"),
    female_parent = c("F1", "F1", "F1", "0")
  )
  fixed <- check_ped(ped, verbose = FALSE)
  kept  <- check_ped(ped, verbose = FALSE, correct_conflicting_trios = FALSE)
  expect_false(any(startsWith(fixed$dependencies$dependency, "Cycle")))
  expect_identical(kept$dependencies$dependency, "Cycle detected involving IDs: A -> B -> A")
})

test_that("clearing an inconsistent sex role can break a loop", {
  ped <- data.frame(
    id            = c("A", "B", "C"),
    male_parent   = c("B", "0", "A"),
    female_parent = c("0", "A", "0")
  )
  fixed <- check_ped(ped, verbose = FALSE)
  kept  <- check_ped(ped, verbose = FALSE, correct_inconsistent_sex_roles = FALSE)
  expect_equal(nrow(fixed$dependencies), 0)
  expect_identical(kept$dependencies$dependency, "Cycle detected involving IDs: A -> B -> A")
})

# ==============================================================================
# 5. Size and depth (the previous search was recursive and path by path)
# ==============================================================================

test_that("a 5000-generation lineage and a 5000-member cycle need no recursion", {
  skip_on_cran()
  n    <- 5000
  ids  <- sprintf("G%05d", seq_len(n))
  open <- data.frame(id = ids, male_parent = c("0", ids[-n]), female_parent = "0",
                     stringsAsFactors = FALSE)
  elapsed <- system.time(out_open <- check_ped(open, verbose = FALSE))[["elapsed"]]
  expect_equal(nrow(out_open$dependencies), 0)
  expect_lt(elapsed, 60)

  closed <- open
  closed$male_parent[1] <- ids[n]    # the founder's sire is the last generation
  out_closed <- check_ped(closed, verbose = FALSE)
  expect_identical(out_closed$dependencies$dependency,
                   paste("Cycle detected involving IDs:",
                         paste(c(ids[1], rev(ids[-1]), ids[1]), collapse = " -> ")))
})

test_that("many shared ancestors and many separate cycles stay fast and exact", {
  skip_on_cran()
  # 40 generations of full-sib mating (paths to the founders double every
  # generation) plus 200 separate 3-member loops
  gen  <- 40
  sibs <- data.frame(
    id            = paste0("G", rep(1:gen, each = 2), c("_m", "_f")),
    male_parent   = c("0", "0", rep(paste0("G", 1:(gen - 1), "_m"), each = 2)),
    female_parent = c("0", "0", rep(paste0("G", 1:(gen - 1), "_f"), each = 2)),
    stringsAsFactors = FALSE
  )
  loops <- do.call(rbind, lapply(1:200, function(k) {
    ids <- paste0("L", k, c("a", "b", "c"))
    data.frame(id = ids, male_parent = ids[c(2, 3, 1)], female_parent = "0",
               stringsAsFactors = FALSE)
  }))
  elapsed <- system.time(out <- check_ped(rbind(sibs, loops), verbose = FALSE))[["elapsed"]]
  expect_lt(elapsed, 60)
  expect_equal(nrow(out$dependencies), 200)
  expect_cycles_match(out)
})
