#' Find performance measure of the user's own markers
#'
#' @param user_markers Markers input from the users.
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
#' @return A list containing the performance of testing set, for user's own marker panel.
#' @export
performanceOwnMarkers <- function (user_markers,
                                   final_out,
                                   method="all",
                                   nrounds=1500,
                                   nthread=6,
                                   testSet = "test",
                                   verbose=FALSE,
                                   ...) {
  
    
    input_matrix_train = final_out$training_matrix
    unique_clusters = final_out$unique_clusters_sample
    clusters_num_train = final_out$training_clusters_num
    clusters_train = final_out$training_clusters

    # Check if all the customised markers are included in the data
    if(!is.character(user_markers)){
      stop("Please input the correct format (character) of the marker list")
    }

    if(!all(user_markers %in% colnames(input_matrix_train))){
      message("Not all input markers are included in the data.")
      message("Only the overlapped marker(s) are processed: ")
      print(user_markers[user_markers %in% colnames(input_matrix_train)])
      message("The following marker(s) are excluded: ")
      print(user_markers[!user_markers %in% colnames(input_matrix_train)])
      # Update the marker list
      user_markers = user_markers[user_markers %in% colnames(input_matrix_train)]
      }

  
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
    list_performance = performanceMarkers(user_markers,
                                          t(as.matrix(input_matrix_train)),
                                          t(as.matrix(input_matrix_test)),
                                          unique_clusters,
                                          clusters_num_train,
                                          clusters_num_test,
                                          clusters_train,
                                          clusters_test,
                                          method = method,
                                          nrounds = nrounds,
                                          nthread = nthread,
                                          method_marker_name = "customised",
                                          verbose = verbose)
        
    
    return(list_performance)
    
    
}
