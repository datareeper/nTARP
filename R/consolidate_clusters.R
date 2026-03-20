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
#' @examples
#' data <- data.frame(X1 = c(0.5, -0.2, 0.1, 0.3, -0.1, 0.2, 5.2, 4.8, 5.1, 5.0,
#'                          -4.5, -5.2, -4.8, -5.1, -4.9, -5.3, 0.0, 0.2, 5.3, -5.0),
#'   X2 = c(0.3, -0.1, 0.2, 0.1, 0.0, 0.2, 5.0, 4.9, 5.3, 5.1,
#'          5.0, 5.2, 4.7, 4.9, 5.1, 4.8, -0.2, 0.0, 5.2, -4.9),
#'   X3 = c(0.4, 0.0, 0.1, -0.1, 0.2, 0.0, 5.1, 4.7, 5.2, 5.0,
#'          -5.0, -4.8, -5.3, -5.1, -4.9, -5.2, 0.1, 0.3, 5.0, -5.1)
#' )
#' nTARP_result <- nTARP_bisecting(data = data,number_of_projections = 100,
#' withinss_threshold = 0.36, minimum_cluster_size_percent = 30)
#' result <- build_solution_from_labeled_clusters(nTARP_best_clusters = nTARP_result$BestClusters,
#' ids = 1:20, contextual_variables_df = data)
#' str(result)
#' result <- consolidate_clusters(result,first_cluster_to_combine = 1,second_cluster_to_combine = 2)
#' str(result)
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
