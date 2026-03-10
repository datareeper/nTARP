#' Thresholding After Random Projections (n-TARP) Clustering
#'
#' Implements the n-TARP clustering technique by projecting the data into
#' a one-dimensional space and performing k-means. The data can be either
#' unlabeled or labeled. The only required parameters are the number of
#' projections and the within-cluster sum of squares threshold.
#' Suggested starting values: `number_of_projections = 1000` and
#' `withinss_threshold = 0.36`.
#'
#' @param data Numeric matrix — dataset to be clustered using `nTARP`
#' @param number_of_projections Numeric — number of random projections for `nTARP` to try for each run
#' @param withinss_threshold Numeric — maximum value defining what a "quality cluster" is,
#'   based on the solution's normalized within-cluster sum of squares (typically 0.36)
#' @param ids Numeric or character vector — identifying labels for individuals in the clusters
#'
#' @return A list containing results and supporting data from the k-means clustering analysis:
#'   (1) `OptimalSolution`: the optimal clustering solution, including cluster assignments and centroids,
#'   (2) `OptimalProjection`: the projection vector associated with the optimal solution,
#'   (3) `Threshold`: the threshold used for determining cluster membership or filtering,
#'   (4) `Direction`: indicates where a new data point should be placed if using the result as a classifier,
#'   (5) `OptimalWithinss`: the within-cluster sum of squares for the optimal solution,
#'   (6) `AllWithinss`: the within-cluster sum of squares for all candidate solutions,
#'   (7) `Clusterings`: all clustering solutions generated during analysis,
#'   (8) `OriginalData`: the original dataset used for clustering,
#'   (9) `OriginalIDs`: the identifiers of the original observations.
#'
#' @references Tarun, Y.; Boutin, M. (2018). n-TARP Binary Clustering Code. Purdue University Research Repository. doi:10.4231/R74B2ZJV
#' @importFrom stats runif kmeans var setNames
#' @export

nTARP <- function(data,number_of_projections,withinss_threshold,ids=NULL)
{
  #Prepare the data
  data <- as.matrix(data)
  m <- nrow(data)
  n <- ncol(data)

  #Preparing the withinss parameters
  withinss <- matrix(0,1,number_of_projections)
  optimal_withinss <- Inf

  #Preparing the random vectors
  random_vectors <- matrix(stats::runif(number_of_projections * n, min = -1, max = 1),
                           nrow = number_of_projections,
                           ncol = n)
  projections <- NULL

  #Preparing the list of valid clusterings
  clusterings <- list()
  index <- 1

  #Try random projections for the number of attempts specified
  for (i in 1:number_of_projections)
  {
    random_vector_basis <- do.call(rbind, replicate(m, random_vectors[i,], simplify=FALSE))
    projection <- rowSums(data*random_vector_basis)
    k.means.results <- stats::kmeans(projection,2,iter.max = 10, nstart = 10,"Hartigan-Wong")
    withinss[i] <- sum(k.means.results$withinss)/(var(projection)*m)
    if (withinss[i] < withinss_threshold)
    {
      cluster.members.1 <- which(k.means.results$cluster == 1)
      cluster.members.2 <- which(k.means.results$cluster == 2)
      clusterings[[index]] <- cluster.members.1
      clusterings[[index+1]] <- cluster.members.2
      index <- index + 2
    }
    if (withinss[i] < optimal_withinss)
    {
      optimal_withinss <- withinss[i]
      k_means_results_optimal <- k.means.results
      optimal_projection <- projection
      optimal_projection_vector <- random_vectors[i,]
    }
  }

  if (optimal_withinss > withinss_threshold)
  {
    #If the threshold is violated, the data may not have a viable clustering
    warning("No clusterings found were below the established threshold")
  }

  #Find the values of the projections in the new 1D space and their means
  optimal_projection_1 <- optimal_projection[which(k_means_results_optimal$cluster == 1)]
  optimal_projection_2 <- optimal_projection[which(k_means_results_optimal$cluster == 2)]
  mean.1 <- mean(optimal_projection_1)
  mean.2 <- mean(optimal_projection_2)

  #The threshold separating the clusters is defined as follows, the direction indicates where a new data point should be placed
  if (mean.1 > mean.2)
  {
    threshold <- (min(optimal_projection_1) - max(optimal_projection_2 ))/2 +   max(optimal_projection_2)
    direction <- 1
  }
  if (mean.1 < mean.2)
  {
    threshold <- (max(optimal_projection_2) - max(optimal_projection_1))/2 + max(optimal_projection_1)
    direction <- 2
  }

  if(!is.null(ids))
  {
    data <- cbind(ids,data)
    data <- as.data.frame(data)
  }

  #Write the output into a list with the necessary values for analysis
  output <- list("OptimalSolution" = k_means_results_optimal,
                 "OptimalProjection" = optimal_projection_vector,
                 "Threshold" = threshold,
                 "Direction" = direction,
                 "OptimalWithinss" = optimal_withinss,
                 "AllWithinss" = withinss,
                 "Clusterings" = clusterings,
                 "OriginalData" = data,
                 "OriginalIDs" = ids
  )
  names(output) <- c("OptimalSolution",
                     "OptimalProjection",
                     "Threshold",
                     "Direction",
                     "OptimalWithinss",
                     "AllWithinss",
                     "Clusterings",
                     "OriginalData",
                     "OriginalIDs"
  )
  return(output)
}
