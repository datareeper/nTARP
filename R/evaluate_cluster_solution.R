#' Evaluate cluster solution
#'
#'#' @keywords  internal
#'
#' Calculates the gain for a cluster solution from `nTARP` based on a user-specified
#' contextual variable to determine which clusters are most informative in explaining
#' the data structure. The gain is calculated using the Gini index. A Gini index of 0
#' indicates a perfectly pure clustering where each cluster fully separates the groups
#' defined by the contextual variable.
#'
#' The Gini index for a cluster is defined as:
#'
#'   Gini = 1 - sum(p_i^2)
#'
#' where p_i is the proportion of observations in the cluster belonging to group i
#' of the contextual variable.
#'
#' The gain for a solution is calculated as the weighted reduction in Gini index:
#'
#'   Gain = Gini(parent) -
#'          [ (n_1 / n) * Gini(cluster_1) +
#'            (n_2 / n) * Gini(cluster_2) ]
#'
#' where n is the size of the parent cluster, and n_1 and n_2 are the sizes of the
#' resulting clusters.
#'
#' @param cluster_member_vector Numeric vector — cluster assignments from `nTARP`
#' @param contextual_variable Vector of integers or characters — variable used as the basis
#'   for comparing clusters
#'
#' @return A list containing:
#'   (1) Gain for the solution,
#'   (2) Distribution of the contextual variable in cluster 1 (frequency),
#'   (3) Distribution of the contextual variable in cluster 2 (frequency),
#'   (4) Distribution of the contextual variable in cluster 1 (percentage),
#'   (5) Distribution of the contextual variable in cluster 2 (percentage),
#'   (6) Gini index for the parent cluster,
#'   (7) Gini index for cluster 1,
#'   (8) Gini index for cluster 2.

evaluate_cluster_solution <- function(cluster_member_vector,contextual_variable)
{
  cluster_1_data <- contextual_variable[which(cluster_member_vector == 1)]
  cluster_2_data <- contextual_variable[which(cluster_member_vector == 2)]
  cluster_1_size <- length(cluster_1_data)
  cluster_2_size <- length(cluster_2_data)
  total_sample <- length(cluster_member_vector)
  distribution_of_classes_1 <- table(cluster_1_data)
  distribution_of_classes_2 <- table(cluster_2_data)
  distribution_of_classes_1_count <- distribution_of_classes_1
  distribution_of_classes_2_count <-  distribution_of_classes_2
  distribution_of_classes_1 <- distribution_of_classes_1/sum(distribution_of_classes_1)
  distribution_of_classes_2 <- distribution_of_classes_2/sum(distribution_of_classes_2)
  overall_distribution <- table(contextual_variable)
  distribution_of_classes_overall <- overall_distribution / sum(overall_distribution)
  Gini_index_1 <- 1 - sum(distribution_of_classes_1^2)
  Gini_index_2 <- 1 - sum(distribution_of_classes_2^2)
  Gini_index_parent <- 1 - sum(distribution_of_classes_overall^2)
  Gain <- Gini_index_parent - (cluster_1_size/total_sample)*Gini_index_1 - (cluster_2_size/total_sample)*Gini_index_2
  output <- list(Gain,
                 distribution_of_classes_1_count,
                 distribution_of_classes_2_count,
                 distribution_of_classes_1,
                 distribution_of_classes_2,
                 Gini_index_parent,
                 Gini_index_1,
                 Gini_index_2
  )
  names(output) <- c("Gain",
                     "Cluster1DistributionCount",
                     "Cluster2DistributionCount",
                     "Cluster1DistributionPercent",
                     "Cluster2DistributionPercent",
                     "GiniParent",
                     "GiniCluster1",
                     "GiniCluster2"
  )
  return(output)
}
