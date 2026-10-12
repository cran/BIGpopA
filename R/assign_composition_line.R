#' Assign Each Sample to a Reference Line
#'
#' Builds the composition table shown in the Familia BreedTools tab from the
#' output of [solve_composition_poly()]: one row per sample with the estimated
#' percentage of each reference population and the predicted line. A sample is
#' assigned to the population with its highest proportion only if that
#' proportion is at least `assign_threshold`; otherwise it is labelled
#' `"Undetermined"`.
#'
#' @param composition Matrix or data.frame returned by
#'   [solve_composition_poly()]: samples in rows, reference populations in
#'   columns (proportions summing to 1). An `R2` column, if present, is dropped.
#' @param assign_threshold Numeric between 0 and 100. Minimum percentage the
#'   highest population must reach for the sample to be assigned to it
#'   (default: 50). Use 0 to always assign the highest population.
#'
#' @return A data.frame with an `ID` column, one column per reference
#'   population holding the estimated percentage (0 to 100, rounded to one
#'   decimal) and a `Predicted line` column. Ties between populations go to the
#'   first column.
#'
#' @examples
#' comp <- rbind(
#'   Sample01 = c(Duroc = 0.998, Hampshire = 0.000, Landrace = 0.002, R2 = 0.991),
#'   Sample02 = c(Duroc = 0.529, Hampshire = 0.471, Landrace = 0.000, R2 = 0.973),
#'   Sample03 = c(Duroc = 0.340, Hampshire = 0.330, Landrace = 0.330, R2 = 0.952)
#' )
#' assign_composition_line(comp)
#' assign_composition_line(comp, assign_threshold = 0)
#'
#' @export
assign_composition_line <- function(composition, assign_threshold = 50) {
  if (!is.numeric(assign_threshold) || length(assign_threshold) != 1 ||
      is.na(assign_threshold) || assign_threshold < 0 || assign_threshold > 100)
    stop("assign_threshold must be a single number between 0 and 100.")

  comp <- as.data.frame(composition, check.names = FALSE)
  comp <- comp[, !colnames(comp) %in% "R2", drop = FALSE]
  if (ncol(comp) == 0) stop("composition has no population columns.")
  comp[] <- lapply(comp, as.numeric)
  mat <- as.matrix(comp)

  top_value <- apply(mat, 1, max, na.rm = TRUE)
  top_line  <- colnames(mat)[max.col(mat, ties.method = "first")]
  predicted <- ifelse(top_value * 100 >= assign_threshold, top_line, "Undetermined")

  out <- data.frame(ID = rownames(mat), round(mat * 100, 1),
                    check.names = FALSE, stringsAsFactors = FALSE,
                    row.names = NULL)
  out[["Predicted line"]] <- predicted
  out
}
