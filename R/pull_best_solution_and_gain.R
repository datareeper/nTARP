#' Get IDs from clusterings
#'
#' #' @keywords  internal
#'
#' Pulls the IDs used to identify observations in cluster solutions so that
#' they can be used for other purposes, including tracking observations
#' between solutions.
#'
#' @param gains_list Numeric vector — gain values from the
#'   `calculate_gains_for_solution` output (`Gains`)
#' @param full_gains_list List — gain values from the
#'   `calculate_gains_for_solution` output (`FullGainInformation`)
#' @param clusterings List — a list of cluster solutions from the `nTARP`
#'   function output (`Clusterings`)
#' @param ids Numeric or character vector — identifying labels for individuals in the clusters
#' @param contextual_variable Vector of integers or characters — variable used as the basis
#'   for comparing clusters
#'
#' @return A list with two entries:
#'   - One containing the IDs in cluster 1
#'   - One containing the IDs in cluster 2


pull_best_solution_and_gain <- function(gains_list,full_gains_list,clusterings,ids,contextual_variable)
{
  clusters_to_get <- which(gains_list == max(gains_list))
  #There might be more than one solution with the maximum gain, we'll select the first one (it doesn't matter much which one is used).
  if (length(clusters_to_get) > 1)
  {
    clusters_to_get <- clusters_to_get[[1]]
  }
  best_solution_distribution <- full_gains_list[clusters_to_get]
  clusters_to_get <- clusterings[c(clusters_to_get,clusters_to_get+1)]
  cluster_1_size <- length(clusters_to_get[[1]])
  cluster_2_size <- length(clusters_to_get[[2]])
  number_of_observations <- cluster_1_size+cluster_2_size
  clusters <- rep(NA,number_of_observations)
  clusters[clusters_to_get[[1]]] <-  1
  clusters[clusters_to_get[[2]]] <-  2

  labeled_clusters <- get_ids_from_clusterings(ids,clusters)
  quality_metrics <- evaluate_cluster_solution(clusters,contextual_variable)
  output <- list(best_solution_distribution,
                 cluster_1_size,
                 cluster_2_size,
                 clusters,
                 labeled_clusters
  )
  names(output) <- c("BestSolutionDistribution",
                     "Cluster1Size",
                     "Cluster2Size",
                     "ClustersForSolution",
                     "LabeledClusters"
  )
  return(output)
}
