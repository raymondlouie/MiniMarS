#' Find performance measure of predicted markers
#'
#' @param markers_sel Selected markers for classification
#' @param input_matrix_train Feature training matrix with cells as columns, and features as rows.
#' @param input_matrix_test Feature test matrix with cells as columns, and features as rows.
#' @param unique_clusters_sample Unique clusters
#' @param clusters_num_train Cluster annotation for training set (numerical)
#' @param clusters_num_test Cluster annotation for testing set (numerical)
#' @param clusters_train Cluster annotation for training set
#' @param clusters_test Cluster annotation for test set
#' @param method List of methods to find cluster markers.
#' \itemize{
#'   \item \code{xgBoost}
#'   \item \code{geneBasis}
#'   \item \code{all}: Use all methods
#' }
#' @param nrounds Number of rounds used in xgBoost algorithm.
#' @param nthread Number of threads used in xgBoost algorithm.
#'
#' @return A list containing performance of testing set, for each cluster
#' @export
performanceAllMarkers <- function (list_markers,
                                   final_out,
                                   method="all",
                                   nrounds=1500,
                                   nthread=6,
                                   verbose=FALSE,
                                   ...) {

    input_matrix_train = final_out$training_matrix
    input_matrix_test = final_out$test_matrix
    unique_clusters = final_out$unique_clusters_sample
    clusters_num_train = final_out$training_clusters_num
    clusters_num_test = final_out$test_clusters_num
    clusters_train = final_out$training_clusters
    clusters_test = final_out$test_clusters

    list_performance = c()
    for (i in 1:length(list_markers)){
        markers_sel = list_markers[[i]]
        markers_sel= markers_sel[!is.na(markers_sel)]
        # print("ite list_markers")
        if (verbose){
            print(names(list_markers)[[i]])
            print(table(clusters_num_train))
            print(table(clusters_train))
            print(table(clusters_num_test))
            print(table(clusters_test))
            print(dim(input_matrix_train))
            print(dim(input_matrix_test))
        }


        list_performance[[names(list_markers)[[i]]]]=performanceMarkers(markers_sel,
                                                                        t(as.matrix(input_matrix_train)),
                                                                        t(as.matrix(input_matrix_test)),
                                                                        unique_clusters,
                                                                        clusters_num_train,
                                                                        clusters_num_test,
                                                                        clusters_train,
                                                                        clusters_test,
                                                                        method=method,
                                                                        nrounds=1500,
                                                                        nthread=6,
                                                                        verbose=verbose)

    }

    return(list_performance)


}
