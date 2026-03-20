#' Combine nTARP BestCluster solutions into one table and assign stable final cluster IDs
#'
#' Consolidates a set of nTARP "best cluster" solutions (branching runs) that contain
#' `LabeledClusters` into a single ID-level data frame. For each run, the function adds a
#' per-run cluster assignment column, constructs a concatenated `ClusterPath`, and assigns
#' a stable numeric `FinalClusterID` based on unique `ClusterPath` values.
#'
#' @param contextual_variables_df Optional `data.frame` of contextual variables to merge in.
#' If supplied, it must contain a column named `id_name` (default `"mcid"`). The merge is a
#' left join on the provided `ids` vector, so every ID in `ids` appears in the output.
#' @param nTARP_best_clusters A list of run objects. Each run object must either:
#' \itemize{
#'   \item contain `LabeledClusters` at the top level, or
#'   \item contain a top-level list element that contains `LabeledClusters`.
#' }
#' `LabeledClusters` must be a named list with `"Cluster 1"` and `"Cluster 2"` entries containing IDs.
#' The list should be named (recommended). Names are used as run identifiers (column suffixes).
#' If unnamed, runs are labeled sequentially (`run1`, `run2`, ...).
#' @param ids A vector of IDs (character or coercible) that defines the universe of rows in the output.
#'
#' @return A `data.frame` with one row per ID in `ids`, optionally merged with contextual variables,
#' plus per-run cluster columns, `ClusterPath`, and `FinalClusterID`.
#' @examples
#' data <- data.frame(X1 = c(0.5, -0.2, 0.1, 0.3, -0.1, 0.2, 5.2, 4.8, 5.1, 5.0,
#'                          -4.5, -5.2, -4.8, -5.1, -4.9, -5.3, 0.0, 0.2, 5.3, -5.0),
#'   X2 = c(0.3, -0.1, 0.2, 0.1, 0.0, 0.2, 5.0, 4.9, 5.3, 5.1,
#'          5.0, 5.2, 4.7, 4.9, 5.1, 4.8, -0.2, 0.0, 5.2, -4.9),
#'   X3 = c(0.4, 0.0, 0.1, -0.1, 0.2, 0.0, 5.1, 4.7, 5.2, 5.0,
#'          -5.0, -4.8, -5.3, -5.1, -4.9, -5.2, 0.1, 0.3, 5.0, -5.1)
#' )
#' nTARP_result <- nTARP_bisecting(data = data,number_of_projections = 100,withinss_threshold = 0.36)
#' result <- build_solution_from_labeled_clusters(nTARP_best_clusters = nTARP_result$BestClusters,
#' ids = 1:10, contextual_variables_df = data)
#' str(result)
#' @export

build_solution_from_labeled_clusters <- function(nTARP_best_clusters,
                                                ids = NULL,
                                                contextual_variables_df = NULL
                                                )
{
  # ---------- helpers ----------
  .get_run_names <- function(x) {
    nms <- names(x)
    if (is.null(nms) || any(nms == "")) {
      nms <- paste0("run", seq_along(x))
    }
    nms
  }

  .infer_branch_paths_from_names <- function(run_names) {
    # Tries to support the legacy "0", "01", "0112" pattern when run_names are numeric-like.
    # If names aren't numeric-like, returns all NA (treat as roots).
    numeric_like <- grepl("^[0-9]+$", run_names)
    if (!all(numeric_like)) return(rep(NA_character_, length(run_names)))

    out <- rep(NA_character_, length(run_names))
    for (i in seq_along(run_names)) {
      code <- run_names[i]
      if (identical(code, "0")) out[i] <- NA_character_ else out[i] <- sub("^0", "", code)
    }
    out
  }

  .extract_labeledclusters_df <- function(run_obj) {
    labeled <- NULL

    if (is.list(run_obj) && "LabeledClusters" %in% names(run_obj)) {
      labeled <- run_obj$LabeledClusters
    } else if (is.list(run_obj)) {
      for (obj in run_obj) {
        if (is.list(obj) && "LabeledClusters" %in% names(obj)) {
          labeled <- obj$LabeledClusters
          break
        }
      }
    }

    if (is.null(labeled)) {
      stop("No `LabeledClusters` found in one of the provided `nTARP_best_clusters` elements.")
    }
    if (!is.list(labeled) || !all(c("Cluster 1", "Cluster 2") %in% names(labeled))) {
      stop("`LabeledClusters` must be a named list with `Cluster 1` and `Cluster 2`.")
    }

    cl1 <- as.character(unlist(labeled[["Cluster 1"]]))
    cl2 <- as.character(unlist(labeled[["Cluster 2"]]))

    data.frame(
      id = c(cl1, cl2),
      cluster = c(rep(1L, length(cl1)), rep(2L, length(cl2))),
      stringsAsFactors = FALSE
    )
  }

  .order_runs_tree <- function(run_names, branch_paths) {
    # Deterministic: roots first, then increasing depth, then lexicographic
    depth <- ifelse(is.na(branch_paths), 0L, nchar(branch_paths))
    order(depth, branch_paths, run_names)
  }

  # ---------- validate inputs ----------

  if (!is.list(nTARP_best_clusters) || length(nTARP_best_clusters) == 0) {
    stop("`nTARP_best_clusters` must be a non-empty list of run objects.")
  }
  if (is.null(ids) || length(ids) == 0) {
    stop("`ids` must be a non-empty vector of IDs.")
  }
  ids <- as.character(ids)

  if (!is.null(contextual_variables_df)) {
    if (!is.data.frame(contextual_variables_df)) {
      stop("`contextual_variables_df` must be a data.frame or NULL.")
    }
  }

  run_names <- .get_run_names(nTARP_best_clusters)
  branch_paths <- .infer_branch_paths_from_names(run_names)
  branch_paths <- as.character(branch_paths)
  ord <- .order_runs_tree(run_names, branch_paths)
  nTARP_best_clusters <- nTARP_best_clusters[ord]
  run_names <- run_names[ord]
  branch_paths <- branch_paths[ord]

  # ---------- initialize output ---------------------------------------

  result <- data.frame(ID = ids)

  result$ClusterPath <- NA_character_

  # ---------- core logic: merge per run + update ClusterPath ----------

  for (i in seq_along(nTARP_best_clusters)) {
    run_id <- run_names[i]
    parent_path <- branch_paths[i]
    colname <- run_id
    df <- .extract_labeledclusters_df(nTARP_best_clusters[[i]])
    names(df) <- c("ID", colname)
    df[["ID"]] <- as.character(df[["ID"]])
    # left-join WITHOUT merge (preserves order)
    m <- match(result[["ID"]], df[["ID"]])
    result[[colname]] <- df[[colname]][m]
    # Update ClusterPath
    if (is.na(parent_path)) {
      # root run: set where ClusterPath is empty
      idx <- !is.na(result[[colname]]) & is.na(result$ClusterPath)
    } else {
      # child run: append only when current path matches parent branch
      idx <- !is.na(result[[colname]]) & result$ClusterPath == parent_path
    }
  }

  # ---------- stable numeric IDs based on path ----------
  result$ClusterPath <- apply(result[, 3:(ncol(result))], 1, function(x) paste0(x, collapse = ""))
  result$ClusterPath <- gsub("NA","",result$ClusterPath)
  cluster_names <- apply(result[, 3:(ncol(result))], 1, function(x) paste0(x, collapse = ""))
  cluster_names <- as.numeric(factor(cluster_names))
  result$FinalClusterID <- cluster_names

  if (!is.null(contextual_variables_df)) {
    #result <- merge(result, contextual_variables_df, by = ID, all.x = TRUE, sort = FALSE)
    result <- cbind(contextual_variables_df,result)
  }
  return(result)
}
