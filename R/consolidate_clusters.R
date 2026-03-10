#' Merge clusters as post-processing
#'
#' This function combines user-specified clusters after `nTARP` has run.
#' The input to the function is the output of the
#' `build_solution_from_labeled_clusters` function, along with the two
#' cluster labels the user would like to merge. The output is returned
#' in the same format as `build_solution_from_labeled_clusters`, with
#' the final solution labels renamed to reflect the merging.
#'
#' @param cluster_path_matrix Data frame — output of
#'   `build_solution_from_labeled_clusters` showing which branch each
#'   observation belongs to from the `nTARP` clustering
#' @param first_cluster_to_combine Numeric — label of the first cluster to merge
#' @param second_cluster_to_combine Numeric — label of the second cluster to merge
#'
#' @return A data frame in the same format as the first argument, with the
#'   final column showing the cluster IDs relabeled based on the chosen merge.
#' @export

consolidate_clusters <- function(cluster_path_matrix,first_cluster_to_combine,second_cluster_to_combine)
{
  clusters_to_combine <- c(first_cluster_to_combine,second_cluster_to_combine)
  final_clusters <- cluster_path_matrix$FinalClusterID
  relevant_indices <- which(final_clusters %in% clusters_to_combine)
  final_clusters[relevant_indices] <- "x"
  final_clusters <- as.numeric(factor(final_clusters))
  cluster_path_matrix$FinalClusterID <- final_clusters
  return(cluster_path_matrix)
}
