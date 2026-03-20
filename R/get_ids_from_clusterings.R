#' Get IDs from clusterings
#'
#' Pulls the IDs used to identify observations in cluster solutions so they
#' can be used for other purposes, including tracking observations between solutions.
#'
#' @param ids Numeric or character vector - identifying labels for individuals in the clusters
#' @param cluster_member_vector Numeric vector - cluster assignments from `nTARP`
#'
#' @return A list with two entries:
#'   - One containing the IDs in cluster 1
#'   - One containing the IDs in cluster 2
#'
#' @keywords internal
#' @noRd

get_ids_from_clusterings <- function(ids, cluster_member_vector)
{
  labeled_clusters <- cluster_member_vector
  for (cluster_index in 1:length(cluster_member_vector))
  {
    labeled_clusters[cluster_index] <- ids[cluster_index]
  }
  first_cluster <- labeled_clusters[which(cluster_member_vector == 1)]
  second_cluster <- labeled_clusters[which(cluster_member_vector == 2)]
  first_cluster <- unname(as.integer(first_cluster))
  second_cluster <- unname(as.integer(second_cluster))
  labeled_clusters <- list(first_cluster,second_cluster)
  names(labeled_clusters) <- c("Cluster 1","Cluster 2")
  return(labeled_clusters)
}
