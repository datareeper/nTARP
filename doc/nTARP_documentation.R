## -----------------------------------------------------------------------------
install.packages("nTARP")
library(nTARP)

## -----------------------------------------------------------------------------
data(iris)
head(iris)

## -----------------------------------------------------------------------------
table(iris$Species)

## -----------------------------------------------------------------------------
set.seed(123532) #random seed for reproducibility
cluster_solution <- nTARP(data = iris[,1:4],
                          number_of_projections = 1000,
                          withinss_threshold = 0.36
                          )

## -----------------------------------------------------------------------------
cluster_solution$OptimalSolution
paste("The optimal within sum of squares was", round(cluster_solution$OptimalWithinss,2))

## -----------------------------------------------------------------------------
print("This is the random vector associated with the best cluster.")
cluster_solution$OptimalProjection
print("This is the boundary where the clusters transition.")
cluster_solution$Threshold
print("The *direction* tells us which cluster has the larger mean, so in this case, points with values larger than the threshold should be assigned to cluster 2.")
cluster_solution$Direction

## -----------------------------------------------------------------------------
observation_to_classify <- iris[1,1:4]
observation_after_projection <- sum(observation_to_classify * cluster_solution$OptimalProjection)
print("The resulting 1-D projection is...")
observation_after_projection 
paste("Our direction, ", cluster_solution$Direction,", tells us that cluster 2 has a mean bigger than the threshold, so since the projected value ", round(observation_after_projection,2)," is larger than the threshold ",round(cluster_solution$Threshold,2)," we should assign it to cluster 2.", sep="")

## -----------------------------------------------------------------------------
hist(cluster_solution$AllWithinss)

## -----------------------------------------------------------------------------
one_feasible_solution <- list(cluster_solution$Clusterings[[1]],cluster_solution$Clusterings[[2]])
one_feasible_solution 

## -----------------------------------------------------------------------------
length(cluster_solution$Clusterings)/2

## -----------------------------------------------------------------------------
set.seed(123532) #random seed for reproducibility
cluster_solution_bisecting <- nTARP_bisecting(data = iris[,1:4],
                                    number_of_projections = 1000,
                                    withinss_threshold = 0.36,
                                    minimum_cluster_size_percent = 20 #default value
                                    )

## -----------------------------------------------------------------------------
cluster_solution_bisecting$BestClusters

## -----------------------------------------------------------------------------
labeled_clusters <- build_solution_from_labeled_clusters(nTARP_best_clusters= cluster_solution_bisecting$BestClusters, ids = 1:nrow(iris))
labeled_clusters[sample(1:150,10,replace = FALSE),] #We'll take a sampling of the observations

## -----------------------------------------------------------------------------
table(labeled_clusters$FinalClusterID,iris$Species)

## -----------------------------------------------------------------------------
labeled_clusters <- consolidate_clusters(cluster_path_matrix = labeled_clusters,
                                         first_cluster_to_combine = 1,
                                         second_cluster_to_combine = 2)
table(labeled_clusters$FinalClusterID,iris$Species)

## -----------------------------------------------------------------------------
#Merge clusters two more times
labeled_clusters <- consolidate_clusters(cluster_path_matrix = labeled_clusters,
                                         first_cluster_to_combine = 1,
                                         second_cluster_to_combine = 2)
labeled_clusters <- consolidate_clusters(cluster_path_matrix = labeled_clusters,
                                         first_cluster_to_combine = 1,
                                         second_cluster_to_combine = 2)
table(labeled_clusters$FinalClusterID,iris$Species)

## -----------------------------------------------------------------------------
set.seed(123532) #random seed for reproducibility
cluster_solution_bisecting_contextual <- nTARP_bisecting(data = iris[,1:4],
                                                         number_of_projections = 1000,
                                                         withinss_threshold = 0.36,
                                                         minimum_cluster_size_percent = 20,
                                                         contextual_variable = iris$Species
                                                         )
labeled_clusters <- build_solution_from_labeled_clusters(nTARP_best_clusters= cluster_solution_bisecting_contextual$BestClusters, ids = 1:nrow(iris))

## -----------------------------------------------------------------------------
labeled_clusters_with_context <- build_solution_from_labeled_clusters(nTARP_best_clusters= cluster_solution_bisecting_contextual$BestClusters, ids = 1:nrow(iris), contextual_variables_df = iris)
head(labeled_clusters_with_context)

## -----------------------------------------------------------------------------
table(labeled_clusters$FinalClusterID,iris$Species)

## -----------------------------------------------------------------------------
labeled_clusters <- consolidate_clusters(cluster_path_matrix = labeled_clusters,
                                         first_cluster_to_combine = 1,
                                         second_cluster_to_combine = 2)
table(labeled_clusters$FinalClusterID,iris$Species)

## -----------------------------------------------------------------------------
if (!requireNamespace("HDclassif", quietly = TRUE)) {
  install.packages("HDclassif")
}
library(HDclassif)
data(wine)
head(wine)
wine[,2:14] <- scale(wine[,2:14])

## -----------------------------------------------------------------------------
table(wine$class)

## -----------------------------------------------------------------------------
set.seed(123532) #random seed for reproducibility
cluster_solution_bisecting_wine <- nTARP_bisecting(data = wine[,2:14],
                                    number_of_projections = 100^2,
                                    withinss_threshold = 0.36,
                                    minimum_cluster_size_percent = 20 #default value
                                    )
labeled_clusters_wine <- build_solution_from_labeled_clusters(nTARP_best_clusters= cluster_solution_bisecting_wine$BestClusters, ids = 1:nrow(wine))
table(labeled_clusters_wine$FinalClusterID,wine$class)

## -----------------------------------------------------------------------------
#Merge clusters two more times
print("This is the result of combining clusters 1 and 2, during which 3-6 were relabled to be 1-4 with the merged cluster being labeled 5.")
labeled_clusters_wine <- consolidate_clusters(cluster_path_matrix = labeled_clusters_wine,
                                         first_cluster_to_combine = 1,
                                         second_cluster_to_combine = 2)
table(labeled_clusters_wine$FinalClusterID,wine$class)
#Merge 1 and 2 again
print("This is the result of combining clusters 3 and 4 (that were relabeled 1 and 2 in the last step). Now there are 4 clusters left.")
labeled_clusters_wine <- consolidate_clusters(cluster_path_matrix = labeled_clusters_wine,
                                         first_cluster_to_combine = 1,
                                         second_cluster_to_combine = 2)
table(labeled_clusters_wine$FinalClusterID,wine$class)
#Merge 1 and 2 again

print("We combine that last two clusters, giving us 3 clusters at the end.")
labeled_clusters_wine <- consolidate_clusters(cluster_path_matrix = labeled_clusters_wine,
                                         first_cluster_to_combine = 1,
                                         second_cluster_to_combine = 2)
table(labeled_clusters_wine$FinalClusterID,wine$class)

## -----------------------------------------------------------------------------
set.seed(123532) #random seed for reproducibility
cluster_solution_bisecting_wine <- nTARP_bisecting(data = wine[,2:14],
                                    number_of_projections = 100^2,
                                    withinss_threshold = 0.36,
                                    minimum_cluster_size_percent = 20,
                                    contextual_variable = wine$class
                                    )
labeled_clusters_wine <- build_solution_from_labeled_clusters(nTARP_best_clusters= cluster_solution_bisecting_wine$BestClusters, ids = 1:nrow(wine))
table(labeled_clusters_wine$FinalClusterID,wine$class)

## -----------------------------------------------------------------------------
#Merge clusters two more times
print("This is the result of combining clusters 1 and 2, during which 3-5 were relabled to be 1-3 with the merged cluster being labeled 4.")
labeled_clusters_wine <- consolidate_clusters(cluster_path_matrix = labeled_clusters_wine,
                                         first_cluster_to_combine = 1,
                                         second_cluster_to_combine = 2)
table(labeled_clusters_wine$FinalClusterID,wine$class)
#Merge 1 and 2 again
print("We combine that last two clusters, giving us 3 clusters at the end.")
labeled_clusters_wine <- consolidate_clusters(cluster_path_matrix = labeled_clusters_wine,
                                         first_cluster_to_combine = 1,
                                         second_cluster_to_combine = 2)
table(labeled_clusters_wine$FinalClusterID,wine$class)

