#' Run nTARP repeatedly in a bisecting fashion
#'
#' Repeatedly applies `nTARP` to iteratively bisect a dataset until a minimum
#' cluster size threshold is reached.
#'
#' This function supports two strategies for selecting the optimal split
#' at each step:
#'
#' (1) Within-Cluster Compactness Criterion:
#'     The optimal solution is selected based on the normalized within-cluster
#'     sum of squares (WSS). The split that minimizes normalized WSS is retained.
#'
#' (2) Contextual Purity Criterion:
#'     The optimal solution is selected using a contextual variable. Inspired
#'     by decision tree learning, the algorithm evaluates candidate splits
#'     based on improvements in class purity (i.e., Gini reduction) with respect
#'     to the contextual variable. The split that maximizes purity gain is retained.
#'
#' The process continues recursively (bisecting the largest eligible cluster)
#' until no resulting cluster meets the user-defined minimum size threshold.
#'
#' @param data Numeric matrix — dataset to be clustered using `nTARP`
#' @param number_of_projections Numeric — number of random projections for `nTARP` to try for each run
#' @param withinss_threshold Numeric — maximum value defining what a "quality cluster" is,
#'   based on the solution's normalized within-cluster sum of squares (typically 0.36)
#' @param ids Numeric or character vector — identifying labels for individuals in the clusters
#' @param contextual_variable Vector of integers or characters — variable to use as the basis
#'   for comparing clusters. This is `NULL` by default, which analytically corresponds
#'   to option (1).
#' @param minimum_cluster_size_percent Numeric — minimum size allowable for a cluster
#'   (expressed as a percentage)
#' @examples
#' # 20-point example dataset
#' data <- data.frame(
#'   X1 = c(0.5, -0.2, 0.1, 0.3, -0.1, 0.2, 5.2, 4.8, 5.1, 5.0,
#'          -4.5, -5.2, -4.8, -5.1, -4.9, -5.3, 0.0, 0.2, 5.3, -5.0),
#'   X2 = c(0.3, -0.1, 0.2, 0.1, 0.0, 0.2, 5.0, 4.9, 5.3, 5.1,
#'          5.0, 5.2, 4.7, 4.9, 5.1, 4.8, -0.2, 0.0, 5.2, -4.9),
#'   X3 = c(0.4, 0.0, 0.1, -0.1, 0.2, 0.0, 5.1, 4.7, 5.2, 5.0,
#'          -5.0, -4.8, -5.3, -5.1, -4.9, -5.2, 0.1, 0.3, 5.0, -5.1)
#' )
#'
#' # Run nTARP without contextual variable
#' result1 <- nTARP_bisecting(
#'   data = data,
#'   number_of_projections = 10,
#'   withinss_threshold = 0.36
#' )
#' str(result1)
#'
#' # Add a latent group as contextual variable
#' latent_group <- c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2,
#'                   3, 3, 3, 3, 3, 3, 1, 1, 2, 3)
#'
#' # Run nTARP with contextual variable
#' result2 <- nTARP_bisecting(
#'   data = data,
#'   number_of_projections = 10,
#'   withinss_threshold = 0.36,
#'   contextual_variable = latent_group
#' )
#' str(result2)
#'
#'
#' @return A list containing:
#'   (1) Complete solutions (i.e., outputs from the `nTARP` function),
#'   (2) Clusters with the best gains identified using the
#'       `pull_best_solution_and_gain` function,
#'   (3) Within-cluster sum of squares for each solution,
#'   (4) Gains for each solution (if a contextual variable is used).
#' @export

nTARP_bisecting <- function(data,
                            number_of_projections,
                            withinss_threshold,
                            ids = NULL,
                            minimum_cluster_size_percent = 20,
                            contextual_variable = NULL
                            )
{
  if (is.null(ids))
  {
    ids <- 1:nrow(data)
  }

  if (is.null(contextual_variable))
  {
    output <- nTARP_complete_solution_no_contextual_variable(data,
                                                                 number_of_projections,
                                                                 withinss_threshold,
                                                                 ids,
                                                                 minimum_cluster_size_percent
    )
  }
  else
  {
    output <- nTARP_complete_solution_with_contextual_variable(data,
                                                             number_of_projections,
                                                             withinss_threshold,
                                                             ids,
                                                             contextual_variable,
                                                             minimum_cluster_size_percent
    )
  }
  return(output)
}
