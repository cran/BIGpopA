#' Check and Correct Common Pedigree Errors
#'
#' Reads a 3-column pedigree file (id, male_parent, female_parent) and performs
#' quality checks, optionally correcting detected errors. Exact duplicates and
#' missing parents are always corrected. Conflicting trios and inconsistent sex
#' roles are corrected when their respective arguments are TRUE. Cycles are
#' reported only and must be resolved manually.
#'
#' @param ped.file Path to the pedigree text file (TSV/CSV/TXT), OR a
#'   data.frame / data.table with columns: id, male_parent, female_parent.
#' @param seed Optional integer seed for reproducibility. Pass NULL (default)
#'   to skip setting a seed.
#' @param verbose Logical. If TRUE (default), prints the report to the console.
#' @param correct_conflicting_trios Logical. If TRUE (default), sets conflicting
#'   male_parent and female_parent to 0 and collapses to one row per ID.
#' @param correct_inconsistent_sex_roles Logical. If TRUE (default), sets
#'   male_parent and female_parent to 0 for rows involving IDs found as both,
#'   then removes any resulting exact duplicates.
#'
#' @return An invisible named list of data frames:
#' \describe{
#'   \item{exact_duplicates}{Exact duplicate rows found in the input.}
#'   \item{conflicting_trios}{IDs with conflicting male_parent or female_parent assignments.}
#'   \item{inconsistent_sex_roles}{Rows where a conflicting ID appears as male_parent or female_parent.}
#'   \item{missing_parents}{Parent IDs absent from id, added as founders.}
#'   \item{dependencies}{Individuals recorded as their own parent, and cycles
#'     (individuals that are their own ancestors). Each cycle is listed once,
#'     as a path in which every id is a recorded parent of the one before it,
#'     back to the first id (e.g. "A -> B -> A"). Must be resolved manually.}
#'   \item{corrected_pedigree}{Corrected pedigree table.}
#' }
#'
#' @examples
#' # Self-contained example using a data.frame
#' ped_df <- data.frame(
#'   id            = c("A", "B", "C", "C", "D"),
#'   male_parent   = c("0", "0", "A", "A", "B"),
#'   female_parent = c("0", "0", "B", "B", "C"),
#'   stringsAsFactors = FALSE
#' )
#' ped_errors <- check_ped(ped.file = ped_df, seed = 101919, verbose = FALSE)
#' names(ped_errors)
#' head(ped_errors$corrected_pedigree)
#'
#' \donttest{
#' library(data.table)
#' ped_dt <- data.table(id = c("A", "B", "C"),
#'                      male_parent   = c("0", "0", "A"),
#'                      female_parent = c("0", "0", "B"))
#' ped_errors <- check_ped(ped.file = ped_dt, verbose = FALSE)
#' }
#'
#' @author Josue Chinchilla-Vargas
#'
#' @importFrom dplyr %>% mutate filter group_by ungroup summarize distinct bind_rows select first n n_distinct if_else row_number
#' @importFrom stats setNames
#' @importFrom utils read.table
#' @importFrom janitor clean_names
#' @export
check_ped <- function(ped.file,
                      seed                           = NULL,
                      verbose                        = TRUE,
                      correct_conflicting_trios      = TRUE,
                      correct_inconsistent_sex_roles = TRUE) {
  
  #### setup ####
  if (!is.null(seed)) set.seed(seed)
  
  # Accept file path OR in-memory data.frame / data.table
  if (is.character(ped.file) && length(ped.file) == 1 && file.exists(ped.file)) {
    data <- utils::read.table(ped.file, header = TRUE)
    data <- janitor::clean_names(data)
  } else if (is.data.frame(ped.file) || data.table::is.data.table(ped.file)) {
    data <- as.data.frame(ped.file)
    data <- janitor::clean_names(data)
  } else {
    stop("ped.file must be a valid file path (character) or a data.frame / data.table.")
  }
  
  # clean_names() already lower-cases; this also covers names it would not map
  data <- .standardize_names(data, c("id", "male_parent", "female_parent"))

  required_cols <- c("id", "male_parent", "female_parent")
  missing_cols  <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(
      "Input is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      ".\nExpected columns: id, male_parent, female_parent."
    )
  }
  
  extra_cols <- setdiff(names(data), required_cols)
  data       <- data[, c(required_cols, extra_cols)]
  
  data <- data %>%
    dplyr::mutate(
      id            = as.character(id),
      male_parent   = as.character(male_parent),
      female_parent = as.character(female_parent)
    )
  
  data <- data %>% dplyr::mutate(row_number = dplyr::row_number(), .before = id)
  
  errors          <- list()
  missing_parents <- data.frame(
    row_number    = integer(),
    id            = character(),
    male_parent   = character(),
    female_parent = character(),
    stringsAsFactors = FALSE
  )
  
  #### check 1: exact duplicates (always fixed) ####
  exact_duplicates <- data[
    duplicated(data %>% dplyr::select(-row_number)) |
      duplicated(data %>% dplyr::select(-row_number), fromLast = TRUE),
  ]
  if (nrow(exact_duplicates) > 0) {
    data <- data %>%
      dplyr::select(-row_number) %>%
      dplyr::distinct() %>%
      dplyr::mutate(row_number = dplyr::row_number(), .before = id)
  }
  
  #### check 2: conflicting trios ####
  repeated_ids <- data %>%
    dplyr::group_by(id) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::ungroup()
  
  conflicting_trios_ids <- repeated_ids %>%
    dplyr::group_by(id) %>%
    dplyr::filter(dplyr::n_distinct(male_parent) > 1 |
                    dplyr::n_distinct(female_parent) > 1) %>%
    dplyr::ungroup()
  
  if (correct_conflicting_trios && nrow(conflicting_trios_ids) > 0) {
    data <- data %>%
      dplyr::group_by(id) %>%
      dplyr::summarize(
        row_number    = dplyr::first(row_number),
        male_parent   = if (dplyr::n_distinct(male_parent)   > 1) "0" else dplyr::first(male_parent),
        female_parent = if (dplyr::n_distinct(female_parent) > 1) "0" else dplyr::first(female_parent),
        .groups = "drop"
      ) %>%
      dplyr::select(row_number, id, male_parent, female_parent)
  }
  
  conflicting_trios <- conflicting_trios_ids
  
  #### check 3: missing parents (always fixed) ####
  for (i in seq_len(nrow(data))) {
    id            <- data$id[i]
    male_parent   <- data$male_parent[i]
    female_parent <- data$female_parent[i]
    
    # Include already-queued missing parents to avoid adding the same ID twice
    all_ids <- c(data$id, missing_parents$id)
    
    if (male_parent != "0" && male_parent != id && !male_parent %in% all_ids) {
      missing_parents <- rbind(
        missing_parents,
        data.frame(row_number = data$row_number[i], id = male_parent,
                   male_parent = "0", female_parent = "0",
                   stringsAsFactors = FALSE)
      )
    }
    if (female_parent != "0" && female_parent != id && !female_parent %in% all_ids) {
      missing_parents <- rbind(
        missing_parents,
        data.frame(row_number = data$row_number[i], id = female_parent,
                   male_parent = "0", female_parent = "0",
                   stringsAsFactors = FALSE)
      )
    }
    
    if (male_parent == id || female_parent == id) {
      errors <- append(errors, paste("Dependency: Individual", id,
                                     "cannot be its own parent"))
    }
  }
  
  missing_parents <- dplyr::distinct(missing_parents)
  if (nrow(missing_parents) > 0) {
    data <- dplyr::bind_rows(data, missing_parents)
  }
  
  #### check 4: inconsistent sex roles ####
  male_ids            <- unique(data$male_parent[data$male_parent   != "0"])
  female_ids          <- unique(data$female_parent[data$female_parent != "0"])
  conflicting_sex_ids <- intersect(male_ids, female_ids)
  
  inconsistent_sex_roles <- data %>%
    dplyr::filter(male_parent %in% conflicting_sex_ids |
                    female_parent %in% conflicting_sex_ids)
  
  if (correct_inconsistent_sex_roles && length(conflicting_sex_ids) > 0) {
    data <- data %>%
      dplyr::mutate(
        male_parent   = dplyr::if_else(male_parent   %in% conflicting_sex_ids, "0", male_parent),
        female_parent = dplyr::if_else(female_parent %in% conflicting_sex_ids, "0", female_parent)
      ) %>%
      dplyr::distinct(id, male_parent, female_parent, .keep_all = TRUE)
  }
  
  #### check 5: dependencies (cycles) -- reported only ####
  # Individuals that are their own ancestors form a group tied by cycles (a
  # strongly connected component of the id -> parent graph). Tarjan's
  # algorithm finds each group once, in linear time and without recursion, so
  # a cycle is reported once however many individuals it has and however many
  # descendants they have. Self-parenting is reported above and skipped here.
  cyc_ids   <- unique(data$id)
  n_ids     <- length(cyc_ids)
  e_child   <- match(c(data$id, data$id), cyc_ids)
  e_parent  <- match(c(data$male_parent, data$female_parent), cyc_ids)
  keep_edge <- !is.na(e_parent) & e_child != e_parent   # parents with a row, no self-parenting
  e_child   <- e_child[keep_edge]
  e_parent  <- e_parent[keep_edge]
  up        <- split(e_parent, factor(e_child, levels = seq_len(n_ids)))   # parents of each id

  index    <- rep(NA_integer_, n_ids)   # visit order
  low      <- integer(n_ids)            # lowest visit order reachable
  on_stack <- logical(n_ids)
  stk_pos  <- integer(n_ids)            # position of each id in the stack
  stk      <- integer(n_ids)            # ids of the groups being built
  path     <- integer(n_ids)            # current search path (no recursion)
  next_up  <- integer(n_ids)            # next parent to visit, per id
  top      <- 0L
  depth    <- 0L
  counter  <- 0L
  groups   <- list()

  for (s in seq_len(n_ids)) {
    if (!is.na(index[s])) next
    todo <- s                                   # id to visit next
    while (!is.na(todo) || depth > 0L) {
      if (!is.na(todo)) {
        # visit: number it and put it on both stacks
        counter        <- counter + 1L
        index[todo]    <- counter
        low[todo]      <- counter
        top            <- top + 1L
        stk[top]       <- todo
        stk_pos[todo]  <- top
        on_stack[todo] <- TRUE
        depth          <- depth + 1L
        path[depth]    <- todo
        next_up[todo]  <- 1L
        todo           <- NA_integer_
      }
      v <- path[depth]
      if (next_up[v] <= length(up[[v]])) {
        w          <- up[[v]][next_up[v]]
        next_up[v] <- next_up[v] + 1L
        if (is.na(index[w])) {
          todo <- w
        } else if (on_stack[w]) {
          low[v] <- min(low[v], index[w])
        }
      } else {
        # all parents of v done: close v, and its group if v is the root
        depth <- depth - 1L
        if (depth > 0L) low[path[depth]] <- min(low[path[depth]], low[v])
        if (low[v] == index[v]) {
          members           <- stk[stk_pos[v]:top]
          on_stack[members] <- FALSE
          top               <- stk_pos[v] - 1L
          if (length(members) > 1) groups[[length(groups) + 1]] <- sort(members)
        }
      }
    }
  }

  # One line per group, in pedigree order: the loop as a path in which each id
  # is a recorded parent of the id before it, back to the first id
  groups <- groups[order(vapply(groups, min, integer(1)))]
  for (members in groups) {
    in_group <- e_child %in% members & e_parent %in% members
    g_child  <- e_child[in_group]
    g_parent <- e_parent[in_group]
    loop     <- members[1]
    v        <- g_parent[match(members[1], g_child)]
    while (!v %in% loop) {
      loop <- c(loop, v)
      v    <- g_parent[match(v, g_child)]
    }
    loop <- c(loop[match(v, loop):length(loop)], v)
    # a single loop when every id has one parent inside the group
    single <- !any(duplicated(g_child[!duplicated(cbind(g_child, g_parent))]))
    cycle  <- if (single) {
      paste(cyc_ids[loop], collapse = " -> ")
    } else {
      paste0(paste(cyc_ids[members], collapse = ", "), " (more than one loop, e.g. ",
             paste(cyc_ids[loop], collapse = " -> "), ")")
    }
    errors <- append(errors, paste("Cycle detected involving IDs:", cycle))
  }
  
  #### compile findings ####
  input_ped_report <- list(
    exact_duplicates       = exact_duplicates,
    conflicting_trios      = conflicting_trios,
    inconsistent_sex_roles = inconsistent_sex_roles,
    missing_parents        = missing_parents,
    dependencies           = data.frame(dependency = as.character(unique(unlist(errors))),
                                        stringsAsFactors = FALSE),   # keeps the column when empty
    corrected_pedigree     = data %>% dplyr::select(-row_number)
  )
  
  #### output ####
  if (verbose) {
    cat("\n=== Pedigree Quality Check Report ===\n")
    
    if (nrow(exact_duplicates) > 0) {
      cat("\nExact duplicate trios detected (removed in corrected pedigree):\n")
      print(exact_duplicates)
    } else cat("\nNo exact duplicate trios found.\n")
    
    if (nrow(conflicting_trios) > 0) {
      cat("\nConflicting trios detected:\n")
      print(conflicting_trios)
      if (correct_conflicting_trios) {
        cat("  -> parents set to 0 and collapsed to one row in corrected pedigree.\n")
      } else {
        cat("  -> correct_conflicting_trios = FALSE: left as-is in corrected pedigree.\n")
      }
    } else cat("\nNo conflicting trios found.\n")
    
    if (nrow(missing_parents) > 0) {
      cat("\nParents missing as IDs (added as founders in corrected pedigree):\n")
      print(missing_parents)
    } else cat("\nNo missing parents found.\n")
    
    if (nrow(inconsistent_sex_roles) > 0) {
      cat("\nIDs found as both male_parent and female_parent:\n")
      print(inconsistent_sex_roles)
      if (correct_inconsistent_sex_roles) {
        cat("  -> parent fields set to 0 for conflicting IDs in corrected pedigree.\n")
      } else {
        cat("  -> correct_inconsistent_sex_roles = FALSE: left as-is.\n")
      }
    } else cat("\nNo IDs found as both male_parent and female_parent.\n")
    
    if (nrow(input_ped_report$dependencies) > 0) {
      cat("\nDependencies detected (must be resolved manually):\n")
      print(input_ped_report$dependencies)
    } else cat("\nNo dependencies detected.\n")
    
    cat("\nThe corrected pedigree is included in the returned list as corrected_pedigree.\n")
  }
  
  invisible(input_ped_report)
}