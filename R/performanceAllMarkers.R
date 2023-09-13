#' Find performance measure of predicted markers
#'
#' @param list_markers Output of `findClusterMarkers` function.
#' @param final_out Output of `processSubsampling` function.
#' @param method List of methods to find the performance of the cluster markers.
#' \itemize{
#'   \item \code{xgBoost}
#'   \item \code{geneBasis}
#'   \item \code{all}: Use all methods
#' }
#' @param nrounds Number of rounds used in xgBoost algorithm.
#' @param nthread Number of threads used in xgBoost algorithm.
#'
#' @return A list containing the performance of testing set, for each marker identification method.
#' @export
performanceAllMarkers <- function (list_markers,
                                   final_out,
                                   method="all",
                                   nrounds=1500,
                                   nthread=6,
                                   testSet = "test"
                                   verbose=FALSE,
                                   ...) {
    
    input_matrix_train = final_out$training_matrix
    unique_clusters = final_out$unique_clusters_sample
    clusters_num_train = final_out$training_clusters_num
    clusters_train = final_out$training_clusters
    
    
    if (testSet %in% "test"){
        input_matrix_test = final_out$test_matrix
        clusters_num_test = final_out$test_clusters_num
        clusters_test = final_out$test_clusters
    } else{
        input_matrix_test = final_out$valid_matrix
        clusters_num_test = final_out$valid_clusters_num
        clusters_test = final_out$valid_clusters
    }
    
    list_performance = c()
    for (i in 1:length(list_markers)){
        markers_sel = list_markers[[i]]
        markers_sel= markers_sel[!is.na(markers_sel)]
        method_name = names(list_markers)[[i]]
        
        # message(method_name)
        if (method_name != "runtime_secs"){
            
            list_performance[[names(list_markers)[[i]]]]=performanceMarkers(markers_sel,
                                                                            t(as.matrix(input_matrix_train)),
                                                                            t(as.matrix(input_matrix_test)),
                                                                            unique_clusters,
                                                                            clusters_num_train,
                                                                            clusters_num_test,
                                                                            clusters_train,
                                                                            clusters_test,
                                                                            method=method,
                                                                            nrounds=nrounds,
                                                                            nthread=nthread,
                                                                            method_marker_name=method_name,
                                                                            verbose=verbose)
            
        }
    }
    
    return(list_performance)
    
    
}
