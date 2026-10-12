#' Validate Pedigree Trios Using Mendelian Error Analysis
#'
#' Validates parent-offspring trios against SNP genotype data using Mendelian
#' error rates. Identifies incorrect parentage assignments and outputs a
#' corrected pedigree. Founder trios
#' (both parents coded as 0) are preserved unchanged if a founders file is
#' supplied. Trios absent from the genotype file are retained as
#' no_genotype_data.
#'
#' @param pedigree_file Path to the pedigree file (TSV/CSV/TXT), OR a
#'   data.frame / data.table with columns: id, male_parent, female_parent.
#' @param genotypes_file Genotypes as any of: a path to a TSV/CSV/TXT file; a
#'   path to a VCF file (`.vcf` or `.vcf.gz`) or a `vcfR` object, converted
#'   with [vcf_to_dosage()] using `ploidy`; a path to a PLINK `.ped` file
#'   (diploid only), converted with [ped_to_dosage()]; or a data.frame /
#'   data.table with an id column followed by marker columns coded as
#'   allele-B dosage (0, 1, ..., ploidy).
#' @param founders_file Optional. Founder IDs as a path to a one-column file, a
#'   character vector of ids, or a data.frame (first column). Founders with both parents coded as 0 are left unchanged.
#'   Defaults to NULL.
#' @param trio_error_threshold Numeric. Maximum Mendelian error percentage to
#'   classify a trio as pass (default: 5.0). Must be between 0 and 100.
#' @param min_markers Integer. Minimum non-missing markers required to evaluate
#'   a trio (default: 10).
#' @param single_parent_error_threshold Numeric. Maximum homozygous-marker
#'   mismatch percentage for a parent to be considered acceptable
#'   (default: 2.0). Must be between 0 and 100.
#' @param verbose Logical. If TRUE, prints progress, summary, and results to
#'   the console (default: TRUE).
#' @param plot_results Logical. If TRUE, prints a histogram of trio Mendelian
#'   error percentages with a threshold line (default: TRUE).
#' @param ploidy Integer >= 2. Ploidy level of the species (2 = diploid,
#'   4 = tetraploid, ...). Genotypes must be coded as allele-B dosage
#'   (0, 1, ..., ploidy). Even ploidy uses the polysomic gamete-range Mendelian
#'   test (assumes autopolyploid inheritance; conservative for allopolyploids).
#'   Odd ploidy (e.g. triploid), where balanced gametes are undefined, falls
#'   back to a model-free opposite-homozygote exclusion evaluated on
#'   homozygous-informative markers only (reduced power). Default is 2.
#' @param marker_trio_table Logical. If TRUE, the returned list includes
#'   `marker_trio_table`, a marker-by-trio table of per-marker statuses (see
#'   Details). Default is FALSE.
#'
#' @return An invisible named list with the following elements:
#' \describe{
#'   \item{pass}{Trios that passed the Mendelian error threshold.}
#'   \item{fail}{Trios that failed the Mendelian error threshold.}
#'   \item{low_markers}{Trios with insufficient markers for evaluation.}
#'   \item{no_genotype_data}{Trios absent from the genotype file.}
#'   \item{founders}{Trios identified as founders.}
#'   \item{missing_parents}{Trios with one or both parents coded as 0 (non-founders).}
#'   \item{full_results}{Complete data.table with all trios and all output columns.}
#'   \item{corrected_pedigree}{Pedigree table after applying recommended corrections.}
#'   \item{marker_summary}{Per-marker Mendelian mismatch counts over every trio
#'     with both parents known and genotyped: `marker`, `trios_tested`,
#'     `trios_mismatch` and `mismatch_pct`, sorted worst first. Markers with a
#'     much higher rate than the rest (e.g. paralogous or miscalled markers)
#'     are candidates for removal before rerunning.}
#'   \item{marker_trio_table}{If `marker_trio_table = TRUE`, a data.table with
#'     one row per marker (`marker`) and one column per trio with both parents
#'     known and genotyped, named `id_male_parent_female_parent`; otherwise NULL.
#'     See Details for the status values.}
#'   \item{plot}{ggplot object if plot_results = TRUE, otherwise NULL.}
#' }
#'
#' @details The cells of `marker_trio_table` take one of these values: `match`;
#'   `mismatch_male` or `mismatch_female` (the offspring is impossible with that
#'   parent alone, whatever the other parent contributes); `mismatch_both` (the
#'   offspring is impossible with either parent alone); `mismatch_combination`
#'   (each parent is compatible alone but the pair is not, for example a 0 x 0
#'   cross with a heterozygous offspring); `untestable` (odd ploidy, marker not
#'   informative); `missing_*` (no call, named after every individual without a
#'   call: `missing_progeny`, `missing_male`, `missing_female`, `missing_male_female`,
#'   `missing_progeny_male`, `missing_progeny_female` or `missing_progeny_male_female`).
#'
#'   The `recommended_correction` column of `full_results` takes these
#'   values:
#' \describe{
#'   \item{none}{Trio passes, or a parent is missing / no genotype data.}
#'   \item{remove_male_parent, remove_female_parent}{The trio fails and only one
#'     parent fails the single-parent homozygous check; that parent is set to 0
#'     in `corrected_pedigree`.}
#'   \item{remove_both}{The trio fails and both parents fail the homozygous
#'     check; both are set to 0.}
#'   \item{unresolved}{The trio fails but neither parent fails the homozygous
#'     check (for example a recorded self whose offspring carries alleles the
#'     parent lacks). The pedigree is likely wrong but the faulty parent cannot
#'     be identified, so `corrected_pedigree` is left unchanged for manual review.}
#'   \item{keep_both}{Only with low markers: the trio is within the threshold
#'     and both parents pass.}
#' }
#'   Low-marker trios get the same values prefixed with `low_markers_`. In the
#'   plot, every trio above the error threshold is coloured as a failure;
#'   unresolved trios are shown in grey.
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
#' ped_df <- data.frame(
#'   id            = c("Off1", "Off2"),
#'   male_parent   = c("P1",   "P1"),
#'   female_parent = c("P2",   "P3"),
#'   stringsAsFactors = FALSE
#' )
#'
#' results <- validate_pedigree(
#'   pedigree_file  = ped_df,
#'   genotypes_file = geno_df,
#'   verbose        = FALSE,
#'   plot_results   = FALSE
#' )
#' print(results$full_results)
#' }
#'
#' @author Josue Chinchilla-Vargas
#'
#' @importFrom data.table fread copy data.table set rbindlist as.data.table is.data.table
#' @importFrom dplyr case_when
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline scale_x_continuous scale_y_continuous scale_fill_manual labs theme_classic theme
#' @export
validate_pedigree <- function(pedigree_file, genotypes_file,
                              founders_file                 = NULL,
                              trio_error_threshold          = 5.0,
                              min_markers                   = 10,
                              single_parent_error_threshold = 2.0,
                              verbose                       = TRUE,
                              plot_results                  = TRUE,
                              ploidy                        = 2,
                              marker_trio_table             = FALSE) {

  ## silence R CMD check NOTEs
  id <- male_parent <- female_parent <- status <- trio_mendelian_error_pct <- NULL
  plot_status <- recommended_correction <- NULL
  marker <- trios_tested <- trios_mismatch <- mismatch_pct <- NULL

  #### Input validation ####
  if (trio_error_threshold < 0 || trio_error_threshold > 100)
    stop("trio_error_threshold must be between 0 and 100")
  if (single_parent_error_threshold < 0 || single_parent_error_threshold > 100)
    stop("single_parent_error_threshold must be between 0 and 100")
  .check_ploidy(ploidy)

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
    pedigree <- read_flex(pedigree_file, "pedigree_file", colClasses = "character")
    genos    <- .read_genotypes(genotypes_file, ploidy = ploidy,
                                format = "data.frame", verbose = verbose)
  }, error = function(e) {
    stop("Error reading input files. Ensure paths are correct and files are TXT/TSV/CSV or VCF. (",
         conditionMessage(e), ")")
  })

  #### Check required columns ####
  # Column names are matched ignoring case (e.g. ID, Male_Parent, FEMALE_PARENT)
  pedigree <- .standardize_names(pedigree, c("id", "male_parent", "female_parent"))
  required_ped_cols <- c("id", "male_parent", "female_parent")
  missing_cols <- setdiff(required_ped_cols, names(pedigree))
  if (length(missing_cols) > 0)
    stop("Pedigree file missing required columns: ",
         paste(missing_cols, collapse = ", "))
  if (!"id" %in% names(genos))
    stop("Genotypes file must have an 'id' column")

  pedigree[, male_parent   := as.character(male_parent)]
  pedigree[, female_parent := as.character(female_parent)]
  original_pedigree <- data.table::copy(pedigree)

  #### Read founders list ####
  founder_ids <- .read_founders(founders_file)

  #### Build genotype matrices ####
  genos_mat   <- as.matrix(genos, rownames = "id")
  genos_hom   <- data.table::copy(genos)
  marker_cols <- setdiff(names(genos_hom), "id")
  for (col in marker_cols)
    genos_hom[get(col) > 0 & get(col) < ploidy, (col) := NA_integer_]
  genos_hom_mat <- as.matrix(genos_hom, rownames = "id")

  #### Identify trios missing from the genotype file ####
  valid_ids <- as.character(genos$id)
  has_geno <- pedigree[id %in% valid_ids &
                         (male_parent %in% valid_ids   | male_parent == "0") &
                         (female_parent %in% valid_ids | female_parent == "0")]
  no_geno_rows <- pedigree[!(id %in% valid_ids) |
                             (!(male_parent %in% valid_ids)   & male_parent != "0") |
                             (!(female_parent %in% valid_ids) & female_parent != "0")]
  if (nrow(no_geno_rows) > 0 && verbose)
    cat(sprintf("Found %d trios with missing genotype data; flagged as no_genotype_data.\n",
                nrow(no_geno_rows)))
  pedigree <- has_geno
  if (nrow(pedigree) == 0)
    stop("No valid trios remain after filtering for genotype availability.")

  #### Per-marker mismatch counters (filled in the trio loop) ####
  marker_names      <- colnames(genos_mat)
  marker_tested_n   <- numeric(length(marker_names))
  marker_mismatch_n <- numeric(length(marker_names))

  #### Main trio evaluation loop ####
  results_list <- lapply(seq_len(nrow(pedigree)), function(i) {
    prog_id          <- pedigree$id[i]
    male_parent_id   <- pedigree$male_parent[i]
    female_parent_id <- pedigree$female_parent[i]
    correction_decision     <- "none"
    error_pct               <- NA_real_
    status                  <- "no_data"
    markers_tested          <- 0L
    male_parent_error_pct   <- NA_real_
    female_parent_error_pct <- NA_real_

    if (male_parent_id == "0" && female_parent_id == "0" &&
        prog_id %in% founder_ids) {
      status              <- "founders"
      correction_decision <- "none"
    } else {
      if (male_parent_id == "0" && female_parent_id == "0") {
        status              <- "missing_both_parents"
        correction_decision <- "none"
      } else if (male_parent_id == "0" && female_parent_id != "0") {
        status              <- "missing_male_parent"
        correction_decision <- "none"
      } else if (male_parent_id != "0" && female_parent_id == "0") {
        status              <- "missing_female_parent"
        correction_decision <- "none"
      } else {
        progeny_vec       <- genos_mat[prog_id, ]
        male_parent_vec   <- genos_mat[male_parent_id, ]
        female_parent_vec <- genos_mat[female_parent_id, ]
        mm_vec <- .mend_mismatch(male_parent_vec, female_parent_vec,
                                 progeny_vec, ploidy)
        tt_vec <- .mend_testable(male_parent_vec, female_parent_vec,
                                 progeny_vec, ploidy)
        mismatches     <- sum(mm_vec, na.rm = TRUE)
        markers_tested <- sum(tt_vec)
        marker_tested_n   <<- marker_tested_n   + as.numeric(tt_vec)
        marker_mismatch_n <<- marker_mismatch_n +
          as.numeric((mm_vec %in% TRUE) & tt_vec)
        if (markers_tested == 0) {
          status              <- "no_data"
          correction_decision <- "none"
        } else {
          error_pct <- (mismatches / markers_tested) * 100
          if (markers_tested < min_markers) {
            status <- "low_markers"
          } else if (error_pct <= trio_error_threshold) {
            status              <- "pass"
            correction_decision <- "none"
          } else {
            status <- "fail"
          }
          if (status %in% c("fail", "low_markers")) {
            progeny_hom       <- genos_hom_mat[prog_id, ]
            male_parent_hom   <- genos_hom_mat[male_parent_id, ]
            female_parent_hom <- genos_hom_mat[female_parent_id, ]
            male_comparisons <- sum(!is.na(male_parent_hom) &
                                            !is.na(progeny_hom))
            male_parent_error_pct <- if (male_comparisons == 0) NA_real_ else
              round((sum(male_parent_hom != progeny_hom, na.rm = TRUE) /
                             male_comparisons) * 100, 2)
            female_comparisons <- sum(!is.na(female_parent_hom) &
                                              !is.na(progeny_hom))
            female_parent_error_pct <- if (female_comparisons == 0) NA_real_ else
              round((sum(female_parent_hom != progeny_hom, na.rm = TRUE) /
                             female_comparisons) * 100, 2)
            male_acceptable   <- !is.na(male_parent_error_pct) &&
              male_parent_error_pct <= single_parent_error_threshold
            female_acceptable <- !is.na(female_parent_error_pct) &&
              female_parent_error_pct <= single_parent_error_threshold
            if (male_acceptable && female_acceptable) {
              # Neither parent is excluded by the homozygous check. If the trio
              # still fails (e.g. a recorded self whose offspring carries alleles
              # absent from the parent), the pedigree is wrong but the faulty
              # parent cannot be identified, so the trio is flagged and the parents are not kept.
              correction_decision <- if (error_pct > trio_error_threshold) "unresolved" else "keep_both"
            } else if (male_acceptable && !female_acceptable) {
              correction_decision    <- "remove_female_parent"
            } else if (!male_acceptable && female_acceptable) {
              correction_decision  <- "remove_male_parent"
            } else {
              correction_decision    <- "remove_both"
            }
            if (status == "low_markers")
              correction_decision <- paste0("low_markers_", correction_decision)
          }
        }
      }
    }

    data.table::data.table(
      id                              = prog_id,
      orig_male_parent                = male_parent_id,
      orig_female_parent              = female_parent_id,
      trio_mendelian_error_pct        = round(error_pct, 2),
      trio_markers_tested             = markers_tested,
      status                          = status,
      recommended_correction          = correction_decision,
      male_parent_hom_error_pct       = male_parent_error_pct,
      female_parent_hom_error_pct     = female_parent_error_pct
    )
  })

  final_df <- data.table::rbindlist(results_list)

  #### Append no_genotype_data rows ####
  if (nrow(no_geno_rows) > 0) {
    no_geno_df <- data.table::data.table(
      id                              = no_geno_rows$id,
      orig_male_parent                = no_geno_rows$male_parent,
      orig_female_parent              = no_geno_rows$female_parent,
      trio_mendelian_error_pct        = NA_real_,
      trio_markers_tested             = 0L,
      status                          = "no_genotype_data",
      recommended_correction          = "none",
      male_parent_hom_error_pct       = NA_real_,
      female_parent_hom_error_pct     = NA_real_
    )
    final_df <- data.table::rbindlist(list(final_df, no_geno_df))
  }

  #### Build corrected pedigree ####
  corrected_pedigree <- data.table::copy(original_pedigree)
  for (i in seq_len(nrow(final_df))) {
    prog_id  <- final_df$id[i]
    decision <- final_df$recommended_correction[i]
    row_idx  <- which(corrected_pedigree$id == prog_id)
    # low_markers_* decisions apply the same removal as the plain decision
    base_decision <- sub("^low_markers_", "", decision)
    if (base_decision %in% c("remove_male_parent", "remove_both"))
      data.table::set(corrected_pedigree, row_idx, "male_parent",   "0")
    if (base_decision %in% c("remove_female_parent", "remove_both"))
      data.table::set(corrected_pedigree, row_idx, "female_parent", "0")
  }

  #### Summary output ####
  if (verbose) {
    total_trios   <- nrow(final_df)
    status_counts <- table(final_df$status)
    cat("\n--- Trio Validation Summary ---\n")
    cat(sprintf("Total trios in pedigree: %d\n", total_trios))
    for (s in names(status_counts))
      cat(sprintf("%-28s: %d (%.1f%%)\n", s, status_counts[s],
                  (status_counts[s] / total_trios) * 100))
    cat(sprintf("Error threshold: %s%%\n", trio_error_threshold))
    cat(sprintf("Homozygous threshold: %s%%\n", single_parent_error_threshold))
    cat(sprintf("Minimum markers required: %s\n\n", min_markers))
    corrections <- table(final_df$recommended_correction)
    cat("Correction summary:\n")
    for (decision in names(corrections))
      if (decision != "none")
        cat(sprintf("  %s: %d\n", decision, corrections[decision]))
    cat("\n")
    print(final_df)
  }

  #### Plot results ####
  p <- NULL
  if (plot_results) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("ggplot2 is required for plot_results = TRUE. Please install it.", call. = FALSE)
    } else {
      # Colour by the trio error first, so no trio above the threshold is shown as pass
      plot_df <- final_df %>%
        dplyr::filter(!is.na(trio_mendelian_error_pct)) %>%
        dplyr::mutate(plot_status = dplyr::case_when(
          trio_mendelian_error_pct <= trio_error_threshold ~ "pass",
          recommended_correction %in% c("remove_male_parent",
                                        "remove_female_parent",
                                        "low_markers_remove_male_parent",
                                        "low_markers_remove_female_parent") ~ "fail_one_parent",
          recommended_correction %in% c("remove_both",
                                        "low_markers_remove_both") ~ "fail_both_parents",
          recommended_correction %in% c("unresolved",
                                        "low_markers_unresolved") ~ "fail_unresolved",
          TRUE ~ "other"
        ))
      n_total <- nrow(plot_df)
      n_fail  <- sum(plot_df$trio_mendelian_error_pct > trio_error_threshold)
      n_pass  <- sum(plot_df$trio_mendelian_error_pct <= trio_error_threshold)
      threshold_label <- paste0(
        "Mendelian Error Threshold: ", trio_error_threshold, "%  |  ",
        "Lost: ", n_fail, " trios  |  ",
        "Kept: ", n_pass, " trios"
      )
      p <- ggplot2::ggplot(plot_df,
                           ggplot2::aes(x = trio_mendelian_error_pct,
                                        fill = plot_status)) +
        ggplot2::geom_histogram(binwidth = 1, color = "white", alpha = 0.9) +
        ggplot2::geom_vline(xintercept = trio_error_threshold,
                            linetype = "dashed", color = "black", linewidth = 1) +
        ggplot2::scale_x_continuous(breaks = seq(0, 100, by = 5)) +
        ggplot2::scale_y_continuous(breaks = seq(0, 100, by = 5)) +
        ggplot2::scale_fill_manual(
          values = c("pass"              = "#339900",
                     "fail_one_parent"   = "#F1C40F",
                     "fail_both_parents" = "#cc3333",
                     "fail_unresolved"   = "#808080",
                     "other"             = "#BDC3C7"),
          labels = c("pass"              = "Pass",
                     "fail_one_parent"   = "Fail - One Parent",
                     "fail_both_parents" = "Fail - Both Parents",
                     "fail_unresolved"   = "Fail - Parents Unresolved",
                     "other"             = "Other")
        ) +
        ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2)) +   # keeps every label visible
        ggplot2::labs(
          title    = "Trio Mendelian Error Distribution",
          subtitle = paste0("Trios with Genotype Data Tested: ", n_total,
                            "\n \n", threshold_label),
          x        = "Mendelian Error (%)",
          y        = "Number of Trios",
          fill     = "Status"
        ) +
        ggplot2::theme_classic(base_size = 13) +
        ggplot2::theme(legend.position = "top")
      print(p)
    }
  }

  #### Optional marker x trio status table ####
  marker_trio_out <- NULL
  if (marker_trio_table) {
    full_trios <- pedigree[male_parent != "0" & female_parent != "0"]
    status_list <- lapply(seq_len(nrow(full_trios)), function(i)
      .marker_trio_status(genos_mat[full_trios$male_parent[i],   ],
                          genos_mat[full_trios$female_parent[i], ],
                          genos_mat[full_trios$id[i],            ],
                          ploidy))
    names(status_list) <- make.unique(
      paste(full_trios$id, full_trios$male_parent, full_trios$female_parent, sep = "_"))
    marker_trio_out <- data.table::as.data.table(
      c(list(marker = marker_names), status_list))
  }

  #### Per-marker mismatch summary ####
  marker_summary <- data.table::data.table(
    marker         = marker_names,
    trios_tested   = as.integer(marker_tested_n),
    trios_mismatch = as.integer(marker_mismatch_n)
  ) %>%
    dplyr::mutate(mismatch_pct = ifelse(trios_tested > 0,
                                        round(trios_mismatch / trios_tested * 100, 2),
                                        NA_real_)) %>%
    dplyr::arrange(dplyr::desc(mismatch_pct), dplyr::desc(trios_mismatch), marker)

  #### Compile and return named list ####
  output_list <- list(
    pass               = final_df[status == "pass"],
    fail               = final_df[status == "fail"],
    low_markers        = final_df[status == "low_markers"],
    no_genotype_data   = final_df[status == "no_genotype_data"],
    founders           = final_df[status == "founders"],
    missing_parents    = final_df[status %in% c("missing_both_parents",
                                                "missing_male_parent",
                                                "missing_female_parent")],
    full_results       = final_df,
    corrected_pedigree = corrected_pedigree,
    marker_summary     = marker_summary,
    marker_trio_table  = marker_trio_out,
    plot               = p
  )
  return(invisible(output_list))
}
