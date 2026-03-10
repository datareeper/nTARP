#' Run nTARP repeatedly in a bisecting fashion (using normalized within sum of squares)
#'
#' #' @keywords  internal
#'
#' Repeatedly applies `nTARP` to iteratively bisect a dataset until a minimum
#' cluster size threshold is reached, using within-cluster compactness to select
#' optimal splits.
#'
#' At each step, the algorithm evaluates candidate splits based on the
#' normalized within-cluster sum of squares (WSS). The split that minimizes
#' normalized WSS is retained.
#'
#' The process continues recursively, bisecting the largest eligible cluster,
#' until no resulting cluster meets the user-defined minimum size threshold.
#'
#' @param data Numeric matrix — dataset to be clustered using `nTARP`
#' @param number_of_projections Numeric — number of random projections for `nTARP` to try for each run (usually 1000 to start)
#' @param withinss_threshold Numeric — maximum value defining what a "quality cluster" is,
#'   based on the solution's normalized within-cluster sum of squares (typically 0.36)
#' @param ids Numeric or character vector — identifying labels for individuals in the clusters
#' @param minimum_cluster_size_percent Numeric — minimum size allowable for a cluster to be
#'   further bisected (as a percentage)
#'
#' @return A list containing:
#'   (1) Complete solutions (i.e., outputs from the `nTARP` function),
#'   (2) Clusters with the best gains identified using the
#'       `pull_best_solution_and_gain` function,
#'   (3) Within-cluster sum of squares for each solution

nTARP_complete_solution_no_contextual_variable <- function(data,
                                                            number_of_projections,
                                                            withinss_threshold,
                                                            ids,
                                                            minimum_cluster_size_percent
)
{
  #Initialize the outputs, which will be lists containing the information for
  #each cluster solution
  complete_solution_list <- list()
  best_cluster_solution_list <- list()
  withinss_list <- list()

  #Next, we'll initialize the other variables we need to track. These variables
  #will be used to name the clusters and record which solutions we've found.
  sample_size <- nrow(data)
  traveled_clusters <- "0"
  initial_cluster <- "0"
  assign(initial_cluster,NULL)
  next_cluster <- initial_cluster

  #We'll print a message that tells the user a solution is attempting to
  #be found
  message(paste("Finding solution for cluster:",next_cluster, sep = " "))

  #We'll start by breaking down the first cluster
  cluster_solution <- nTARP(data,number_of_projections,withinss_threshold,ids)
  assign(initial_cluster,cluster_solution)
  number_of_observations <- nrow(data)
  labeled_clusters <- get_ids_from_clusterings(ids,
                                                    cluster_solution$OptimalSolution$cluster
                                                    )

  best_solution <- list(BestSolutionWSS = cluster_solution$OptimalWithinss,
                     Cluster1Size = cluster_solution$OptimalSolution$size[1],
                     Cluster2Size = cluster_solution$OptimalSolution$size[2],
                     ClustersForSolution = unname(as.integer(cluster_solution$OptimalSolution$cluster)),
                     LabeledClusters = labeled_clusters
  )

  best_cluster_output_name <- paste(next_cluster,"BestCluster",sep="_")
  assign(best_cluster_output_name,best_solution)

  #We'll record the solution in the lists we established earlier.
  complete_solution_list <- append(complete_solution_list,list(cluster_solution))
  names(complete_solution_list)[length(complete_solution_list)] <- next_cluster
  best_cluster_solution_list <- append(best_cluster_solution_list,list(best_solution))
  names(best_cluster_solution_list)[length(best_cluster_solution_list)] <- next_cluster
  withinss_list <- append(withinss_list,list(cluster_solution$AllWithinss))
  names(withinss_list)[length(withinss_list)] <- next_cluster

  cluster1size <- best_solution$Cluster1Size
  cluster2size <- best_solution$Cluster2Size
  smallest_cluster <- min(cluster1size,cluster2size)

  #Next, we'll break down the first branch of the clusters until we get
  #to a cluster smaller than a user-defined size.
  while(smallest_cluster >= minimum_cluster_size_percent*sample_size/100)
  {
    if (cluster1size >= cluster2size) {
        cluster_to_split <- 1
      } else {
        cluster_to_split <- 2
      }
    next_cluster <- paste(next_cluster,cluster_to_split,sep="")
    assign(next_cluster,NULL)
    message(paste("Finding solution for cluster:",next_cluster, sep = " "))

    data_to_use <- cluster_solution$OriginalData
    data_to_use <- as.data.frame(data_to_use)
    ID <- data_to_use$ids

    #First we'll pull the data from the appropriate cluster
    relevant_observations_in_cluster <- unname(which(best_solution$ClustersForSolution == cluster_to_split))
    relevant_data_in_cluster_to_split <- data_to_use[relevant_observations_in_cluster,]
    relevant_data_in_cluster_to_split <- relevant_data_in_cluster_to_split[,-1]
    relevant_data_in_cluster_to_split[] <- lapply(relevant_data_in_cluster_to_split, as.numeric)
    ID <- ID[relevant_observations_in_cluster]

        cluster_solution <- nTARP(relevant_data_in_cluster_to_split,
                              number_of_projections,
                              withinss_threshold,
                              ID)
    assign(next_cluster,cluster_solution)
    number_of_observations <- nrow(relevant_data_in_cluster_to_split)

    labeled_clusters <- get_ids_from_clusterings(ID,
                                                      cluster_solution$OptimalSolution$cluster
    )

    best_solution <- list(BestSolutionWSS = cluster_solution$OptimalWithinss,
                       Cluster1Size = cluster_solution$OptimalSolution$size[1],
                       Cluster2Size = cluster_solution$OptimalSolution$size[2],
                       ClustersForSolution = unname(as.integer(cluster_solution$OptimalSolution$cluster)),
                       LabeledClusters = labeled_clusters
    )

    best_cluster_output_name <- paste(next_cluster,"BestCluster",sep="_")
    assign(best_cluster_output_name ,best_solution)

    if (!(next_cluster %in% names(complete_solution_list))) {
      complete_solution_list[[next_cluster]] <- cluster_solution
    }
    if (!(next_cluster %in% names(best_cluster_solution_list))) {
      best_cluster_solution_list[[next_cluster]] <- best_solution
    }
    if (!(next_cluster %in% names(withinss_list))) {
      withinss_list[[next_cluster]] <- cluster_solution$AllWithinss
    }

    traveled_clusters <- c(traveled_clusters,next_cluster)
    cluster1size <- best_solution$Cluster1Size
    cluster2size <- best_solution$Cluster2Size
    smallest_cluster <- min(cluster1size,cluster2size)
  }
  traveled_clusters <- traveled_clusters[-length(traveled_clusters)]


  #Now we'll break down the clusters iteratively
  while (length(traveled_clusters) > 0)
  {
    #We'll keep a list of the clusters we visit during the breakdown so we can
    #return to a solution that could be broken down further.
    current_travel <- traveled_clusters
    traveled_clusters <- character(0)  # empty the list

    #For each of the cluster solutions, we'll try breaking it down further by returning
    #to previous solutions and running the same process as before.
    for (solution_name in rev(current_travel)) {
      relevant_solution_name_overall <- solution_name
      relevant_solution_name_best <- paste(relevant_solution_name_overall,"BestCluster",sep="_")
      cluster_solution <- get(relevant_solution_name_overall)
      best_solution <- get(relevant_solution_name_best)
      cluster1size <- best_solution$Cluster1Size
      cluster2size <- best_solution$Cluster2Size
      smallest_cluster <- min(cluster1size,cluster2size)

      #We need to reverse the logic this time to make sure the other
      #cluster is broken down.
      if (cluster1size < cluster2size) {
        cluster_to_split <- 1
      } else {
        cluster_to_split <- 2
      }
      next_cluster <- paste(solution_name,cluster_to_split,sep="")

      #This section repeats the same process as in the first breakdown of the
      #initial cluster solution.
      while(smallest_cluster >= minimum_cluster_size_percent*sample_size/100)
      {
        data_to_use <- cluster_solution$OriginalData
        data_to_use <- as.data.frame(data_to_use)
        ID <- data_to_use$ids

        #First we'll pull the data from the appropriate cluster
        relevant_observations_in_cluster <- unname(which(best_solution$ClustersForSolution == cluster_to_split))
        relevant_data_in_cluster_to_split <- data_to_use[relevant_observations_in_cluster,]
        relevant_data_in_cluster_to_split <- relevant_data_in_cluster_to_split[,-1]
        relevant_data_in_cluster_to_split[] <- lapply(relevant_data_in_cluster_to_split, as.numeric)
        ID <- ID[relevant_observations_in_cluster]


        cluster_solution <- nTARP(relevant_data_in_cluster_to_split,
                                  number_of_projections,
                                  withinss_threshold,
                                  ID)
        assign(next_cluster,cluster_solution)
        number_of_observations <- nrow(relevant_data_in_cluster_to_split)

        labeled_clusters <- get_ids_from_clusterings(ID,
                                                          cluster_solution$OptimalSolution$cluster
        )

        best_solution <- list(BestSolutionWSS = cluster_solution$OptimalWithinss,
                           Cluster1Size = cluster_solution$OptimalSolution$size[1],
                           Cluster2Size = cluster_solution$OptimalSolution$size[2],
                           ClustersForSolution = unname(as.integer(cluster_solution$OptimalSolution$cluster)),
                           LabeledClusters = labeled_clusters
        )

        best_cluster_output_name <- paste(next_cluster,"BestCluster",sep="_")
        assign(best_cluster_output_name ,best_solution)

        if (cluster1size >= cluster2size) {
          cluster_to_split <- 1
        } else {
          cluster_to_split <- 2
        }

        if (!(next_cluster %in% names(complete_solution_list))) {
          complete_solution_list[[next_cluster]] <- cluster_solution
          message(paste("Finding solution for cluster:",next_cluster, sep = " "))
        }
        if (!(next_cluster %in% names(best_cluster_solution_list))) {
          best_cluster_solution_list[[next_cluster]] <- best_solution
        }
        if (!(next_cluster %in% names(withinss_list))) {
          withinss_list[[next_cluster]] <- cluster_solution$AllWithinss
        }

        traveled_clusters <- c(traveled_clusters,next_cluster)
        next_cluster <- paste(next_cluster,cluster_to_split,sep="")

        assign(next_cluster,NULL)
        cluster1size <- best_solution$Cluster1Size
        cluster2size <- best_solution$Cluster2Size
        smallest_cluster <- min(cluster1size,cluster2size)
      }
    }
    traveled_clusters <- traveled_clusters[-length(traveled_clusters)]
  }

  #Rename the clusters to ensure they are fetchable using typical R functions
  complete_names <- names(complete_solution_list)
  complete_names <- paste("Cluster",complete_names,sep="_")
  names(complete_solution_list) <- complete_names
  names(best_cluster_solution_list) <- complete_names

  #Form the output as a list
  output <- list(complete_solution_list,
                 best_cluster_solution_list,
                 withinss_list
                 )
  names(output) <- c("SolutionsList",
                     "BestClusters",
                     "Withinss"
                     )

  return(output)
}
