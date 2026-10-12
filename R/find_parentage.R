#' Find Parentage Assignments for Progeny
#'
#' Assigns the most likely parent(s) to each progeny from SNP genotype data
#' using Mendelian error rates or homozygous mismatch rates. Parents or progeny
#' absent from the genotype file are removed with a warning.
#'
#' @param genotypes_file Genotypes as any of: a path to a TSV/CSV/TXT file; a
#'   path to a VCF file (`.vcf` or `.vcf.gz`) or a `vcfR` object, converted
#'   with [vcf_to_dosage()] using `ploidy`; a path to a PLINK `.ped` file
#'   (diploid only), converted with [ped_to_dosage()]; or a data.frame /
#'   data.table with an 'id' column followed by marker columns coded as
#'   allele-B dosage (0, 1, ..., ploidy).
#' @param parents_file Path to a TSV/CSV/TXT file, OR a data.frame /
#'   data.table with an 'id' column and an optional 'sex' column
#'   ('M', 'F', or 'A'). If absent, all parents are treated as ambiguous.
#' @param progeny_file Path to a TSV/CSV/TXT file, OR a data.frame /
#'   data.table with an 'id' column. For `method = "fill_pedigree"` it must
#'   be a pedigree with `id`, `male_parent` and `female_parent` columns (e.g.
#'   the `corrected_pedigree` from [validate_pedigree()]); see Details.
#' @param method Character. One of "best_male_parent", "best_female_parent",
#'   "best_match", "best_pair" (default), or "fill_pedigree" (requires a
#'   pedigree as `progeny_file`; see Details).
#' @param min_markers Integer. Minimum markers required; fewer flags
#'   low_markers (default: 10).
#' @param error_threshold Numeric. Maximum mismatch percentage; exceeded values
#'   flag high_error (default: 5.0). Must be between 0 and 100.
#' @param show_ties Logical. If TRUE, tied best pairs are appended as suffix
#'   columns. Default is TRUE.
#' @param allow_parent_selfing Logical. If FALSE, candidate pairs with identical
#'   male and female parent IDs are excluded. Applies only when method is
#'   "best_pair". Default is FALSE.
#' @param exclude_self_match Logical. If TRUE, each progeny ID is excluded from
#'   its own candidate parent set, preventing self-matches when progeny are also
#'   present in the parents file. Default is TRUE.
#' @param verbose Logical. If TRUE, prints progress and summary. Default is TRUE.
#' @param plot_results Logical. If TRUE, plots the Mendelian error distribution.
#'   Requires ggplot2. Default is TRUE.
#' @param ploidy Integer >= 2. Ploidy level of the species (2 = diploid,
#'   4 = tetraploid, ...). Genotypes must be coded as allele-B dosage
#'   (0, 1, ..., ploidy). Even ploidy uses the polysomic gamete-range Mendelian
#'   test (assumes autopolyploid inheritance; conservative for allopolyploids).
#'   Odd ploidy (e.g. triploid), where balanced gametes are undefined, falls
#'   back to a model-free opposite-homozygote exclusion evaluated on
#'   homozygous-informative markers only (reduced power). Default is 2.
#' @param founders_file Optional. Founder IDs as a path to a one-column file, a
#'   character vector of ids, or a data.frame (first column). Only used when `progeny_file` is a pedigree: founders with
#'   both parents unknown are skipped and not searched. Default is NULL.
#' @param top_n Integer >= 1. Number of best candidates to report per progeny.
#'   With the default of 1, only the best assignment is reported (plus tied
#'   best pairs when `show_ties = TRUE` and `method = "best_pair"`). With
#'   `top_n > 1`, candidates are ranked by error (ties broken by the number of
#'   markers tested) and the best `top_n` are reported: rank 1 in the base
#'   columns and ranks 2 to `top_n` in suffix columns (`male_parent_2`,
#'   `female_parent_2`, `mendelian_error_pct_2`, `markers_tested_2` for
#'   "best_pair"; `best_match_2`, `mendelian_error_pct_2`, `markers_tested_2`
#'   for the other methods). Ties are then ranked like any other candidate, so
#'   `show_ties` has no effect. Progeny with fewer than `top_n` candidates
#'   have `NA` in the remaining ranks. Status refers to rank 1 only.
#'
#' @details **Method "fill_pedigree" ("fill in the blanks").** The
#'   `progeny_file` must have `male_parent` and `female_parent` columns. A parent coded `0`, `NA` or
#'   empty is treated as unknown and each row is handled on its own:
#'   a known male searches only females, a known female searches only males,
#'   and a row with both unknown gets the full pair search (unless listed in
#'   `founders_file`). Rows with both parents known are skipped. The search is
#'   the "best_pair" Mendelian trio test with the known parent held
#'   fixed, and candidates come from `parents_file` (sex M/F/A as before); the
#'   known parent only needs genotypes. Two columns are added to the results:
#'   `search_mode` ("pair", "known_male" or "known_female") and
#'   `known_parent_error_pct` (homozygous mismatch between the known parent and
#'   the progeny; a high value suggests the known parent is wrong). A known
#'   parent absent from the genotype file gives
#'   `status = "known_parent_no_genotype"`. Without `method = "fill_pedigree"`,
#'   any `male_parent` / `female_parent` columns in `progeny_file` are ignored.
#'
#' @return A named list (returned invisibly) with elements:
#' \describe{
#'   \item{pass}{Progeny with a confident parentage assignment.}
#'   \item{high_error}{Progeny whose best assignment exceeds the error threshold.}
#'   \item{low_markers}{Progeny with insufficient markers for a valid assignment.}
#'   \item{full_results}{Complete data.table with all progeny and all output columns.
#'     `tied_candidates` gives the number of candidates (pairs for "best_pair")
#'     that share the lowest error for each progeny; values above 1 mean the
#'     best assignment is ambiguous.}
#'   \item{plot}{ggplot object if plot_results = TRUE, otherwise NULL.}
#' }
#'
#' @examples
#' \donttest{
#' geno_df <- data.frame(
#'   id  = c("P1", "P2", "P3", "Off1", "Off2"),
#'   S1  = c(0L, 2L, 0L, 1L, 0L),
#'   S2  = c(2L, 0L, 2L, 1L, 2L),
#'   S3  = c(0L, 2L, 0L, 1L, 0L),
#'   S4  = c(2L, 0L, 2L, 1L, 2L),
#'   S5  = c(0L, 2L, 0L, 1L, 0L),
#'   S6  = c(2L, 0L, 2L, 1L, 2L),
#'   S7  = c(0L, 2L, 0L, 1L, 0L),
#'   S8  = c(2L, 0L, 2L, 1L, 2L),
#'   S9  = c(0L, 2L, 0L, 1L, 0L),
#'   S10 = c(2L, 0L, 2L, 1L, 2L)
#' )
#'
#' parents_df <- data.frame(
#'   id  = c("P1", "P2", "P3"),
#'   sex = c("M",  "F",  "F"),
#'   stringsAsFactors = FALSE
#' )
#'
#' progeny_df <- data.frame(
#'   id = c("Off1", "Off2"),
#'   stringsAsFactors = FALSE
#' )
#'
#' results <- find_parentage(
#'   genotypes_file = geno_df,
#'   parents_file   = parents_df,
#'   progeny_file   = progeny_df,
#'   method         = "best_pair",
#'   verbose        = FALSE,
#'   plot_results   = FALSE
#' )
#' print(results$full_results)
#' }
#'
#' @author Josue Chinchilla-Vargas
#'
#' @importFrom data.table fread copy CJ rbindlist set data.table as.data.table is.data.table
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline scale_x_continuous scale_y_continuous scale_fill_manual labs theme_classic theme
#' @export
find_parentage <- function(genotypes_file, parents_file, progeny_file,
                           method                = "best_pair",
                           min_markers           = 10,
                           error_threshold       = 5.0,
                           show_ties             = TRUE,
                           allow_parent_selfing   = FALSE,
                           exclude_self_match    = TRUE,
                           verbose               = TRUE,
                           plot_results          = TRUE,
                           ploidy                = 2,
                           top_n                 = 1,
                           founders_file         = NULL) {
  
  ## silence R CMD check NOTEs
  id <- sex <- male_parent <- female_parent <- NULL
  mendelian_error_pct <- plot_status <- status <- NULL
  
  #### Input Validation ####
  allowed_methods <- c("best_male_parent", "best_female_parent",
                       "best_match", "best_pair", "fill_pedigree")
  if (!method %in% allowed_methods)
    stop("Method must be one of: ", paste(allowed_methods, collapse = ", "))
  if (min_markers < 1)
    stop("min_markers must be a positive integer.")
  if (error_threshold < 0 || error_threshold > 100)
    stop("error_threshold must be between 0 and 100.")
  .check_ploidy(ploidy)
  if (!is.numeric(top_n) || length(top_n) != 1 || is.na(top_n) ||
      top_n < 1 || top_n != round(top_n))
    stop("top_n must be a single positive integer.")
  top_n <- as.integer(top_n)

  # Accept file path OR in-memory data.frame / data.table
  read_flex <- function(x, label, ...) {
    if (is.character(x) && length(x) == 1) {
      if (!file.exists(x))
        stop("Error reading input files. Ensure paths are correct and files are TXT/TSV/CSV.")
      data.table::fread(x, sep = "auto", ...)
    } else if (is.data.frame(x) || data.table::is.data.table(x)) {
      data.table::as.data.table(x)
    } else {
      stop(label, " must be a file path (character) or a data.frame / data.table.")
    }
  }
  
  tryCatch({
    genos              <- .read_genotypes(genotypes_file, ploidy = ploidy,
                                          format = "data.frame", verbose = verbose)
    all_parents        <- read_flex(parents_file,   "parents_file")
    progeny_candidates <- read_flex(progeny_file,   "progeny_file")
  }, error = function(e) {
    stop("Error reading input files. Ensure paths are correct and files are TXT/CSV/TSV or VCF. (",
         conditionMessage(e), ")")
  })
  
  # Column names are matched ignoring case (e.g. ID, Sex)
  all_parents        <- .standardize_names(all_parents,        c("id", "sex"))
  progeny_candidates <- .standardize_names(progeny_candidates,
                                           c("id", "male_parent", "female_parent"))

  #### Pedigree ("fill in the blanks") mode ####
  is_unknown <- function(x) is.na(x) | trimws(as.character(x)) %in% c("", "0", "NA")
  method_reported <- method
  pedigree_mode   <- method == "fill_pedigree"
  if (pedigree_mode) {
    if (!all(c("male_parent", "female_parent") %in% names(progeny_candidates)))
      stop("method = 'fill_pedigree' requires a progeny_file with id, male_parent ",
           "and female_parent columns.")
    method <- "best_pair"
    progeny_candidates[, male_parent   := as.character(male_parent)]
    progeny_candidates[, female_parent := as.character(female_parent)]
    founder_ids <- .read_founders(founders_file)
    um         <- is_unknown(progeny_candidates$male_parent)
    uf         <- is_unknown(progeny_candidates$female_parent)
    is_founder <- um & uf & progeny_candidates$id %in% founder_ids
    if (verbose)
      cat(sprintf("Pedigree mode: %d complete trios skipped, %d founders skipped, %d rows with a missing parent to search.\n",
                  sum(!um & !uf), sum(is_founder), sum((um | uf) & !is_founder)))
    progeny_candidates <- progeny_candidates[(um | uf) & !is_founder]
    if (nrow(progeny_candidates) == 0)
      stop("No progeny with a missing parent to fill in.")
  }

  valid_ids       <- genos$id
  removed_parents <- setdiff(all_parents$id, valid_ids)
  if (length(removed_parents) > 0) {
    warning("The following parent IDs were not in the genotype file and will not be analyzed: ",
            paste(removed_parents, collapse = ", "), call. = FALSE)
    all_parents <- all_parents[id %in% valid_ids]
  }
  
  removed_progeny <- setdiff(progeny_candidates$id, valid_ids)
  if (length(removed_progeny) > 0) {
    warning("The following progeny IDs were not in the genotype file and will not be analyzed: ",
            paste(removed_progeny, collapse = ", "), call. = FALSE)
    progeny_candidates <- progeny_candidates[id %in% valid_ids]
  }
  
  if (!"sex" %in% colnames(all_parents)) {
    warning("No 'sex' column in parents file. All parents treated as ambiguous ('A').")
    all_parents[, sex := "A"]
  }
  
  all_parents[, sex := toupper(sex)]
  all_parents <- unique(all_parents, by = c("id", "sex"))
  male_parent_candidates   <- unique(all_parents[sex %in% c("M", "A", "NA"), .SD],
                                           by = "id")
  female_parent_candidates <- unique(all_parents[sex %in% c("F", "A", "NA")],
                                           by = "id")
  
  if (nrow(male_parent_candidates) == 0 &&
      method %in% c("best_male_parent", "best_pair"))
    warning("No valid male parent candidates remain after filtering.", call. = FALSE)
  if (nrow(female_parent_candidates) == 0 &&
      method %in% c("best_female_parent", "best_pair"))
    warning("No valid female parent candidates remain after filtering.", call. = FALSE)
  if (nrow(progeny_candidates) == 0)
    stop("No valid progeny candidates remain after filtering.")
  
  #### Pre-compute genotype matrices once ####
  genos_mat   <- as.matrix(genos, rownames = "id")
  genos_hom   <- data.table::copy(genos)
  marker_cols <- setdiff(names(genos_hom), "id")
  for (col in marker_cols)
    genos_hom[get(col) > 0 & get(col) < ploidy, (col) := NA_integer_]
  genos_hom_mat <- as.matrix(genos_hom, rownames = "id")
  
  #### Status helper ####
  assign_status <- function(markers, error_pct) {
    ifelse(markers < min_markers, "low_markers",
                 ifelse(error_pct > error_threshold, "high_error", "pass"))
  }
  
  #### Logic for Homozygous Matching Methods ####
  if (method %in% c("best_male_parent", "best_female_parent", "best_match")) {
    parent_ids <- switch(method,
                               "best_male_parent"   = male_parent_candidates$id,
                               "best_female_parent" = female_parent_candidates$id,
                               "best_match"         = union(male_parent_candidates$id,
                                                                  female_parent_candidates$id))
    parent_genos  <- genos_hom_mat[rownames(genos_hom_mat) %in% parent_ids,            , drop = FALSE]
    progeny_genos <- genos_hom_mat[rownames(genos_hom_mat) %in% progeny_candidates$id, , drop = FALSE]
    
    n_progeny  <- nrow(progeny_genos)
    results_dt <- data.table::data.table(
      id                  = rownames(progeny_genos),
      best_match          = NA_character_,
      mendelian_error_pct = NA_real_,
      markers_tested      = NA_integer_,
      status              = NA_character_,
      tied_candidates     = NA_integer_
    )
    extra_list <- vector("list", n_progeny)
    
    for (i in seq_len(n_progeny)) {
      prog_id     <- rownames(progeny_genos)[i]
      progeny_vec <- progeny_genos[i, ]
      progeny_mat <- matrix(progeny_vec,
                                  nrow = nrow(parent_genos),
                                  ncol = ncol(parent_genos),
                                  byrow = TRUE)
      mismatches       <- rowSums(parent_genos != progeny_mat, na.rm = TRUE)
      comparisons      <- rowSums(!is.na(parent_genos) & !is.na(progeny_mat))
      percent_mismatch <- (mismatches / comparisons) * 100
      percent_mismatch[is.nan(percent_mismatch)] <- NA
      
      if (exclude_self_match) {
        self_idx <- rownames(parent_genos) == prog_id
        percent_mismatch[self_idx] <- NA_real_
      }
      
      if (top_n > 1) {
        valid_idx  <- which(!is.na(percent_mismatch))
        ranked_idx <- valid_idx[order(percent_mismatch[valid_idx],
                                            -comparisons[valid_idx])]
        top_idx    <- ranked_idx[seq_len(min(top_n, length(ranked_idx)))]
        best_idx   <- top_idx[seq_len(min(1L, length(top_idx)))]
      } else {
        best_idx <- which.min(percent_mismatch)
      }
      if (length(best_idx) == 0) {
        data.table::set(results_dt, i, "markers_tested", 0L)
        data.table::set(results_dt, i, "status",         "low_markers")
        next
      }
      
      best_markers <- comparisons[best_idx]
      best_error   <- round(percent_mismatch[best_idx], 2)
      data.table::set(results_dt, i, "best_match",          rownames(parent_genos)[best_idx])
      data.table::set(results_dt, i, "mendelian_error_pct", best_error)
      data.table::set(results_dt, i, "markers_tested",      as.integer(best_markers))
      data.table::set(results_dt, i, "status",              assign_status(best_markers, best_error))
      data.table::set(results_dt, i, "tied_candidates",
                      as.integer(sum(percent_mismatch == percent_mismatch[best_idx],
                                                 na.rm = TRUE)))
      
      if (top_n > 1 && length(top_idx) > 1) {
        rank_row <- list(id = prog_id)
        for (k in seq(2, length(top_idx))) {
          rank_row[[paste0("best_match_",           k)]] <- rownames(parent_genos)[top_idx[k]]
          rank_row[[paste0("mendelian_error_pct_",  k)]] <- round(percent_mismatch[top_idx[k]], 2)
          rank_row[[paste0("markers_tested_",       k)]] <- as.integer(comparisons[top_idx[k]])
        }
        extra_list[[i]] <- data.table::as.data.table(rank_row)
      }
    }
    
    # Add ranks 2..top_n as suffix columns, keeping the progeny order
    extra_dt <- data.table::rbindlist(Filter(Negate(is.null), extra_list),
                                      fill = TRUE, use.names = TRUE)
    if (top_n > 1) {
      for (k in seq(2, top_n)) {
        cols <- paste0(c("best_match_", "mendelian_error_pct_", "markers_tested_"), k)
        for (col in cols) {
          vals <- if (nrow(extra_dt) > 0 && col %in% names(extra_dt))
            extra_dt[[col]][match(results_dt$id, extra_dt$id)] else
              rep(if (grepl("^best_match", col)) NA_character_ else
                if (grepl("^markers", col)) NA_integer_ else NA_real_, nrow(results_dt))
          data.table::set(results_dt, j = col, value = vals)
        }
      }
    }
    final_df <- results_dt
  }
  
  #### Logic for Best Pair Method ####
  if (method == "best_pair") {
    parent_pairs <- data.table::CJ(male_parent   = male_parent_candidates$id,
                                   female_parent = female_parent_candidates$id)
    if (!allow_parent_selfing) {
      parent_pairs <- parent_pairs[male_parent != female_parent]
      if (verbose) cat("Parent selfing is disallowed. Pairs with identical parents are removed.\n")
    }
    if (nrow(parent_pairs) == 0) stop("No valid parent pairs to test.")
    
    male_parent_genos_mat   <- genos_mat[parent_pairs$male_parent,   , drop = FALSE]
    female_parent_genos_mat <- genos_mat[parent_pairs$female_parent, , drop = FALSE]
    progeny_ids <- progeny_candidates$id
    n_all       <- length(progeny_ids)
    if (pedigree_mode) {
      kmale <- ifelse(is_unknown(progeny_candidates$male_parent),
                            NA_character_, progeny_candidates$male_parent)
      kfem  <- ifelse(is_unknown(progeny_candidates$female_parent),
                            NA_character_, progeny_candidates$female_parent)
    } else {
      kmale <- kfem <- rep(NA_character_, n_all)
    }
    search_mode <- ifelse(!is.na(kmale), "known_male",
                                ifelse(!is.na(kfem), "known_female", "pair"))
    pair_idx    <- which(search_mode == "pair")
    # Full pair matrices are only built for progeny with no known parent
    progeny_mat <- genos_mat[progeny_ids[pair_idx], , drop = FALSE]
    n_progeny   <- nrow(progeny_mat)
    n_pairs     <- nrow(parent_pairs)
    
    mismatch_mat <- matrix(
      vapply(seq_len(n_progeny), function(j) {
        progeny_vec <- progeny_mat[j, ]
        progeny_pair_mat <- matrix(progeny_vec,
                                         nrow = n_pairs,
                                         ncol = ncol(progeny_mat),
                                         byrow = TRUE)
        rowSums(
          .mend_mismatch(male_parent_genos_mat, female_parent_genos_mat,
                         progeny_pair_mat, ploidy),
          na.rm = TRUE
        )
      }, numeric(n_pairs)),
      nrow = n_pairs, ncol = n_progeny
    )
    
    comparison_mat <- matrix(
      vapply(seq_len(n_progeny), function(j) {
        progeny_vec <- progeny_mat[j, ]
        progeny_pair_mat <- matrix(progeny_vec,
                                         nrow = n_pairs,
                                         ncol = ncol(progeny_mat),
                                         byrow = TRUE)
        rowSums(
          .mend_testable(male_parent_genos_mat, female_parent_genos_mat,
                         progeny_pair_mat, ploidy)
        )
      }, numeric(n_pairs)),
      nrow = n_pairs, ncol = n_progeny
    )
    
    pct_mismatch_mat <- (mismatch_mat / comparison_mat) * 100
    pct_mismatch_mat[is.nan(pct_mismatch_mat)] <- NA
    
    results_dt <- data.table::data.table(
      id                  = progeny_ids,
      male_parent         = NA_character_,
      female_parent       = NA_character_,
      mendelian_error_pct = NA_character_,
      markers_tested      = NA_integer_,
      status              = NA_character_,
      tied_candidates     = NA_integer_
    )
    
    # Known-parent bookkeeping (pedigree mode only)
    known_id_vec <- ifelse(search_mode == "known_male", kmale,
                                 ifelse(search_mode == "known_female", kfem,
                                              NA_character_))
    if (pedigree_mode) {
      known_err <- rep(NA_real_, n_all)
      for (j in which(search_mode != "pair")) {
        kid <- known_id_vec[j]
        if (kid %in% rownames(genos_hom_mat)) {
          kh <- genos_hom_mat[kid, ]
          ph <- genos_hom_mat[progeny_ids[j], ]
          n_cmp <- sum(!is.na(kh) & !is.na(ph))
          known_err[j] <- if (n_cmp == 0) NA_real_ else
            round(sum(kh != ph, na.rm = TRUE) / n_cmp * 100, 2)
        }
      }
      data.table::set(results_dt, j = "search_mode",            value = search_mode)
      data.table::set(results_dt, j = "known_parent_error_pct", value = known_err)
      data.table::set(results_dt, which(search_mode == "known_male"),
                      "male_parent",   known_id_vec[search_mode == "known_male"])
      data.table::set(results_dt, which(search_mode == "known_female"),
                      "female_parent", known_id_vec[search_mode == "known_female"])
    }
    
    results_list <- vector("list", n_all)
    for (j in seq_len(n_all)) {
      prog_id <- progeny_ids[j]
      mode_j  <- search_mode[j]
      
      if (mode_j == "pair") {
        pairs_j          <- parent_pairs
        pj               <- match(j, pair_idx)
        percent_mismatch <- pct_mismatch_mat[, pj]
        comparisons      <- comparison_mat[,  pj]
      } else {
        known_id <- known_id_vec[j]
        if (!known_id %in% rownames(genos_mat)) {
          warning("Known parent '", known_id, "' of progeny '", prog_id,
                  "' is not in the genotype file; this progeny cannot be tested.",
                  call. = FALSE)
          data.table::set(results_dt, j, "markers_tested", 0L)
          data.table::set(results_dt, j, "status",         "known_parent_no_genotype")
          next
        }
        cand <- if (mode_j == "known_male") female_parent_candidates$id else
          male_parent_candidates$id
        if (!allow_parent_selfing) cand <- cand[cand != known_id]
        if (exclude_self_match)    cand <- cand[cand != prog_id]
        cand <- cand[cand %in% rownames(genos_mat)]
        if (length(cand) == 0) {
          data.table::set(results_dt, j, "markers_tested", 0L)
          data.table::set(results_dt, j, "status",         "low_markers")
          next
        }
        pairs_j <- if (mode_j == "known_male")
          data.table::data.table(male_parent = known_id, female_parent = cand) else
            data.table::data.table(male_parent = cand,   female_parent = known_id)
        pm <- genos_mat[pairs_j$male_parent,   , drop = FALSE]
        pf <- genos_mat[pairs_j$female_parent, , drop = FALSE]
        pg <- matrix(genos_mat[prog_id, ], nrow = nrow(pairs_j),
                           ncol = ncol(genos_mat), byrow = TRUE)
        mm          <- rowSums(.mend_mismatch(pm, pf, pg, ploidy), na.rm = TRUE)
        comparisons <- rowSums(.mend_testable(pm, pf, pg, ploidy))
        percent_mismatch <- (mm / comparisons) * 100
        percent_mismatch[is.nan(percent_mismatch)] <- NA
      }
      
      if (exclude_self_match) {
        self_pair_idx <- pairs_j$male_parent == prog_id |
          pairs_j$female_parent == prog_id
        percent_mismatch[self_pair_idx] <- NA_real_
      }
      
      min_mismatch_val <- min(percent_mismatch, na.rm = TRUE)
      
      if (is.infinite(min_mismatch_val)) {
        data.table::set(results_dt, j, "markers_tested", 0L)
        data.table::set(results_dt, j, "status",         "low_markers")
        next
      }
      
      if (top_n > 1) {
        # Rank all candidate pairs by error, then by markers tested; keep top_n
        valid_idx    <- which(!is.na(percent_mismatch))
        ranked_idx   <- valid_idx[order(percent_mismatch[valid_idx],
                                              -comparisons[valid_idx])]
        best_indices <- ranked_idx[seq_len(min(top_n, length(ranked_idx)))]
      } else {
        best_indices <- which(percent_mismatch == min_mismatch_val)
        if (length(best_indices) > 1) {
          best_markers_per_pair <- comparisons[best_indices]
          max_markers           <- max(best_markers_per_pair)
          best_indices          <- best_indices[best_markers_per_pair == max_markers]
        }
      }
      
      best_pairs   <- pairs_j[best_indices]
      best_markers <- comparisons[best_indices[1]]
      best_error   <- round(min_mismatch_val, 2)
      a_status     <- assign_status(best_markers, best_error)
      
      if (top_n == 1 && !show_ties && nrow(best_pairs) > 1) {
        warning("Progeny '", prog_id, "' has ", nrow(best_pairs),
                " tied best pairs. Only one is reported as show_ties=FALSE.",
                call. = FALSE)
      }
      
      report_all    <- show_ties || top_n > 1
      num_to_report <- min(nrow(best_pairs),
                                 if (report_all) nrow(best_pairs) else 1)
      
      data.table::set(results_dt, j, "male_parent",         best_pairs$male_parent[1])
      data.table::set(results_dt, j, "female_parent",       best_pairs$female_parent[1])
      data.table::set(results_dt, j, "mendelian_error_pct", sprintf("%.2f", min_mismatch_val))
      data.table::set(results_dt, j, "markers_tested",      as.integer(best_markers))
      data.table::set(results_dt, j, "status",              a_status)
      data.table::set(results_dt, j, "tied_candidates",
                      as.integer(sum(percent_mismatch == min_mismatch_val,
                                                 na.rm = TRUE)))
      
      if (report_all && num_to_report > 1) {
        tie_row <- list(id = prog_id)
        for (k in seq(2, num_to_report)) {
          tie_row[[paste0("male_parent_",         k)]] <- best_pairs$male_parent[k]
          tie_row[[paste0("female_parent_",       k)]] <- best_pairs$female_parent[k]
          tie_row[[paste0("mendelian_error_pct_", k)]] <-
            if (top_n > 1) round(percent_mismatch[best_indices[k]], 2) else min_mismatch_val
          tie_row[[paste0("markers_tested_",      k)]] <- comparisons[best_indices[k]]
        }
        results_list[[j]] <- data.table::as.data.table(tie_row)
      }
    }
    
    tie_rows <- data.table::rbindlist(
      Filter(Negate(is.null), results_list),
      fill      = TRUE,
      use.names = TRUE
    )
    if (nrow(tie_rows) > 0) {
      final_df <- merge(results_dt, tie_rows, by = "id", all.x = TRUE)
      for (col in names(final_df))
        data.table::set(final_df, which(final_df[[col]] == ""), col, NA_character_)
    } else {
      final_df <- results_dt
    }
  }
  
  #### Compile named list ####
  output_list <- list(
    pass         = final_df[status == "pass"],
    high_error   = final_df[status == "high_error"],
    low_markers  = final_df[status == "low_markers"],
    full_results = final_df,
    plot         = NULL
  )
  
  #### Verbose output ####
  if (verbose) {
    total_progeny <- nrow(final_df)
    cat("\n=== Parentage Assignment Report ===\n")
    cat(sprintf("\nTotal progeny evaluated: %d\n", total_progeny))
    cat(sprintf("Method: %s  |  Error threshold: %s%%  |  Minimum markers: %s\n",
                method_reported, error_threshold, min_markers))
    
    n_pass <- nrow(output_list$pass)
    if (n_pass > 0) {
      cat(sprintf("\n%d progeny passed (%.1f%%).\n",
                              n_pass, (n_pass / total_progeny) * 100))
    } else {
      cat("\nNo progeny passed.\n")
    }
    
    n_high <- nrow(output_list$high_error)
    if (n_high > 0) {
      cat(sprintf("\n%d progeny flagged high_error (%.1f%%):\n",
                              n_high, (n_high / total_progeny) * 100))
      print(output_list$high_error)
    } else {
      cat("\nNo progeny flagged for high error.\n")
    }
    
    n_low <- nrow(output_list$low_markers)
    if (n_low > 0) {
      cat(sprintf("\n%d progeny flagged low_markers (%.1f%%):\n",
                              n_low, (n_low / total_progeny) * 100))
      print(output_list$low_markers)
    } else {
      cat("\nNo progeny flagged for low marker count.\n")
    }
    
    cat("\nFull results are included in the returned list as $full_results.\n")
  }
  
  #### Plot Results ####
  if (plot_results) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("ggplot2 is required for plot_results = TRUE. Please install it.",
              call. = FALSE)
    } else {
      plot_df <- final_df[!is.na(final_df$mendelian_error_pct)]
      plot_df$mendelian_error_pct <- as.numeric(plot_df$mendelian_error_pct)
      plot_df$plot_status <- ifelse(
        plot_df$status == "pass",        "pass",
        ifelse(
          plot_df$status == "high_error",  "high_error",
          ifelse(
            plot_df$status == "low_markers", "low_markers", "other")))
      
      n_total <- nrow(plot_df)
      n_pass  <- sum(plot_df$status == "pass",        na.rm = TRUE)
      n_high  <- sum(plot_df$status == "high_error",  na.rm = TRUE)
      n_low   <- sum(plot_df$status == "low_markers", na.rm = TRUE)
      
      threshold_label <- paste0(
        "Error Threshold: ", error_threshold, "%  |  ",
        "Pass: ", n_pass, "  |  ",
        "High Error: ", n_high, "  |  ",
        "Low Markers: ", n_low
      )
      
      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = mendelian_error_pct, fill = plot_status)
      ) +
        ggplot2::geom_histogram(binwidth = 1, color = "white", alpha = 0.9) +
        ggplot2::geom_vline(xintercept = error_threshold,
                            linetype = "dashed", color = "black", linewidth = 1) +
        ggplot2::scale_x_continuous(breaks = seq(0, 100, by = 5)) +
        ggplot2::scale_y_continuous() +
        ggplot2::scale_fill_manual(
          values = c("pass"        = "#339900",
                     "high_error"  = "#cc3333",
                     "low_markers" = "#F1C40F",
                     "other"       = "#BDC3C7"),
          labels = c("pass"        = "Pass",
                     "high_error"  = "High Error",
                     "low_markers" = "Low Markers",
                     "other"       = "Other")
        ) +
        ggplot2::labs(
          title    = "Parentage Mendelian Error Distribution",
          subtitle = paste0("Progeny Tested: ", n_total,
                                  "\n \n", threshold_label),
          x        = "Mendelian Error (%)",
          y        = "Number of Progeny",
          fill     = "Status"
        ) +
        ggplot2::theme_classic(base_size = 13) +
        ggplot2::theme(legend.position = "top")
      
      print(p)
      output_list$plot <- p
    }
  }
  
  return(invisible(output_list))
}
