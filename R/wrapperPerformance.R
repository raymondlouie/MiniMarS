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
performanceMarkers <- function (markers_sel,
                                input_matrix_train,
                                input_matrix_test,
                                unique_clusters_sample,
                                clusters_num_train,
                                clusters_num_test,
                                clusters_train,
                                clusters_test,
                                method="all",
                                nrounds=1500,
                                nthread=6,
                                method_marker_name="xgBoost",
                                verbose=FALSE,
                                ...) {

    all_methods = c("xgBoost","geneBasis")

    if (method == "all"){
        method = all_methods
    }

    diff_methods = setdiff(method,all_methods)

    if (length(diff_methods) >  0) {
        warning(paste0(method, " not found. Using remaining or all methods."))
        method = intersect(method,all_methods)
        if (length(method) == 0) {
            method = all_methods
        }
    }

    if (length(clusters_num_test) != dim(input_matrix_test)[2]) {
        stop("Number of cells in cluster_num_test does not match the number of cells in the test matrix.")
    }

    if (length(clusters_test) != dim(input_matrix_test)[2]) {
        stop("Number of cells in cluster_test does not match the number of cells in the test matrix.")
    }

    if (length(clusters_num_train) != dim(input_matrix_train)[2]) {
        stop("Number of cells in cluster_num_test does not match the number of cells in the training matrix.")
    }

    if (length(clusters_train) != dim(input_matrix_train)[2]) {
        stop("Number of cells in cluster_train does not match the number of cells in the training matrix.")
    }

    list_measures = list()

    for (i in 1:length(method)) {
        curr_method = method[[i]]

        if (verbose){
            message("Calculating performance of marker selection method:",method_marker_name,  "using performance method:",curr_method,
                    ".\n")
        }

        if (curr_method == "geneBasis") {
            performance_df=geneBasisPerformance(markers_sel,
                                                t(input_matrix_test),
                                                clusters_test,
                                                unique_clusters_sample,
                                                ...)
        }

        if (curr_method == "xgBoost") {
            performance_df=xgboostPerformance(markers_sel,
                                              t(input_matrix_train),
                                              t(input_matrix_test),
                                              clusters_num_train,
                                              clusters_num_test,
                                              clusters_train,
                                              clusters_test,
                                              unique_clusters_sample,
                                              nrounds=nrounds,
                                              nthread=nthread,
                                              ...)
        }


        list_measures[[curr_method]] = performance_df


    }


    return(list_measures)


}
