#' Calculate gain for solutions
#'
#' @keywords  internal
#'
#' Calculates the gain for each cluster solution from the `nTARP`
#' output based on a user-specified contextual variable, in order to determine
#' which clusters are most helpful for explaining the clustering structure.
#' The gains are derived from the Gini index of the original cluster and the
#' two resulting clusters obtained by bisecting it using `nTARP`.
#'
#' The gain for a split is defined as:
#'
#'   Gain = Gini(parent) -
#'          [ (n_1 / n) * Gini(cluster_1) +
#'            (n_2 / n) * Gini(cluster_2) ]
#'
#' where n is the size of the original (parent) cluster, and n_1 and n_2
#' are the sizes of the two resulting clusters.
#'
#' @param clusterings List — a list of cluster solutions from the `nTARP`
#'   function output (`"Clusterings"`)
#' @param contextual_variable Vector of integers or characters — variable to use as the basis
#'   for comparing clusters
#'
#' @return A list containing:
#'   (1) Gains for each pair of clusters, and
#'   (2) Output from the `evaluate_cluster_solution` function containing full
#'       information on each solution's Gini indices and distribution of groups.

calculate_gain_for_solutions <- function(clusterings,contextual_variable)
{
  number_of_clusterings <- length(clusterings)
  number_of_observations <- length(clusterings[[1]]) + length(clusterings[[2]])
  gain <- rep(-Inf,number_of_clusterings)
  index <- 1
  recording_index <- 1
  gain_list <- list()
  while (index < (number_of_clusterings))
  {
    cluster_solution <- rep(NA,number_of_observations)
    cluster_solution[clusterings[[index]]] <- 1
    cluster_solution[clusterings[[index + 1]]] <- 2
    result <- evaluate_cluster_solution(cluster_solution,contextual_variable)
    gain[recording_index ] <- result$Gain
    gain_list[[recording_index]] <- result
    index <- index + 2
    recording_index <- recording_index + 2
  }
  output <- list(gain,gain_list)
  names(output) <- c("Gains","FullGainInformation")
  return(output)
}
