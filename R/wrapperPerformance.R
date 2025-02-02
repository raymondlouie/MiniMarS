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

    # Remove NA markers, which may occur due to geneBasis
    indexNA  = which(is.na(markers_sel))
    if (length(indexNA)>0){
        markers_sel = markers_sel[-indexNA]
    }

    all_methods = c("xgBoost","geneBasis")

    if (method == "all"){
        method = all_methods
    }

    diff_methods = setdiff(method,all_methods)

    if (length(diff_methods) >  0) {
        warning(paste0(method, " not found. Using remaining or all methods.\n"))
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
    runtime_secs = list()
    for (i in 1:length(method)) {
        curr_method = method[[i]]

        if (verbose){
            message("Calculating performance of marker selection method ", paste(method_marker_name,collapse=", "),
                    " using performance method ",curr_method,
                    ".\n")
        }

        if (curr_method == "geneBasis") {
            start_time <- Sys.time()

            test_prediction=geneBasisPredictCell(markers_sel,
                                                 t(input_matrix_test),
                                                 clusters_test,
                                                 unique_clusters_sample)

            performance_df <- calculatePerformance (test_prediction,
                                                    unique_clusters_sample)
            performance_cluster_df = performance_df[[1]]
            performance_all_df = performance_df[[2]]

            end_time <- Sys.time()
            runtime_secs[i] <- as.numeric(end_time-start_time, units="secs")
            names(runtime_secs)[i] <- "geneBasis"
        }

        if (curr_method == "xgBoost") {
            start_time <- Sys.time()

            test_prediction=xgboostPredictCell(markers_sel,
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
            performance_df <- calculatePerformance (test_prediction,
                                                    unique_clusters_sample)
            performance_cluster_df = performance_df[[1]]
            performance_all_df = performance_df[[2]]

            end_time <- Sys.time()
            runtime_secs[i] <- as.numeric(end_time-start_time, units="secs")
            names(runtime_secs)[i] <- "xgBoost"
        }


        list_measures[[paste0(curr_method,"_performance_cluster")]] = performance_cluster_df
        list_measures[[paste0(curr_method,"_performance_all")]] = performance_all_df

        list_measures[['runtime']]=runtime_secs

    }


    return(list_measures)


}



#' Wrapper function for xgboostPerformance
#'
#' @param markers_sel Single cell experiment object
#' @param input_matrix_train Feature training matrix with cells as columns, and features as rows.
#' @param input_matrix_test Feature test matrix with cells as columns, and features as rows.
#' @param clusters_num_train Cluster annotation for training set (numerical)
#' @param clusters_num_test Cluster annotation for test set (numerical)
#' @param clusters_train Cluster annotation for training set
#' @param clusters_test Cluster annotation for test set
#' @param unique_clusters_sample Unique clusters
#'
#' @return The performance of the input markers, as determined by a xgBoost algorithm
#' @export
xgboostPredictCell <- function (markers_sel,
                                input_matrix_train,
                                input_matrix_test,
                                clusters_num_train,
                                clusters_num_test,
                                clusters_train,
                                clusters_test,
                                unique_clusters_sample,
                                ...){

    clusters_num_train = as.numeric(clusters_num_train)
    clusters_num_test = as.numeric(clusters_num_test)

    xgboost_train = xgboost::xgb.DMatrix(data=as.matrix(input_matrix_train[, markers_sel]),
                                         label=clusters_num_train)

    xgboost_test = xgboost::xgb.DMatrix(data=as.matrix(input_matrix_test[, markers_sel]),
                                        label=clusters_num_test)

    # train a model using our training data
    numberOfClasses <- length(unique(clusters_num_train))

    xgb_params <- list("objective" = "multi:softprob",
                       "eval_metric" = "mlogloss",
                       "num_class" = numberOfClasses)

    built.model <- xgboost::xgb.train(params = xgb_params,
                                      data = xgboost_train,
                                      ...)

    importance_train.mt <- xgboost::xgb.importance(colnames(xgboost_train),
                                                   model = built.model)


    pred_test <- predict(built.model,
                         newdata = xgboost_test)

    test_prediction <- matrix(pred_test, nrow = numberOfClasses,
                              ncol=length(pred_test)/numberOfClasses) %>%
        t() %>%
        data.frame() %>%
        mutate(true = clusters_num_test ,
               predict_num = max.col(., "last")-1) ##predicted

    test_prediction$celltype <- clusters_test
    test_prediction$mapped_celltype <- unlist(lapply(test_prediction$predict_num,
                                                     function (x) unique_clusters_sample[which(names(unique_clusters_sample) %in% x)]))


    return(test_prediction)

}



#' Wrapper function for geneBasisPerformance
#'
#' @param markers_sel Selected markers
#' @param input_matrix_test Feature test matrix with cells as columns, and features as rows.
#' @param clusters_test Cluster annotation for test set
#' @param unique_clusters_sample Unique clusters
#'
#' @return True positive rate calculated using get_celltype_mapping from geneBasisR
#' @export
geneBasisPredictCell <- function (markers_sel,
                                  input_matrix_test,
                                  clusters_test,
                                  unique_clusters_sample,
                                  ...){

    sce_test  <- SingleCellExperiment::SingleCellExperiment(list(counts=t(input_matrix_test)),
                                                            colData=data.frame(cell_type=clusters_test))
    # logcounts(sce_test) <- log2(t(input_matrix_test) + 1)
    SingleCellExperiment::logcounts(sce_test) <- t(input_matrix_test)

    cluster_map = geneBasisR::get_celltype_mapping(sce_test ,
                                                   genes.selection = markers_sel,
                                                   celltype.id = "cell_type",
                                                   return.stat = T)

    test_prediction=cluster_map$mapping


    return(test_prediction)

}


#' Wrapper function for geneBasisPerformance
#'
#' @param test_prediction Selected markers
#' @param unique_clusters_sample Unique clusters
#'
#' @return List of two elements: performance of each cluster, and performance over all clusters. Metrics include precision_weighted, precision_macro, recall_weighted, recall_macro, F1_macro, F1_weighted, precision_micro
#' @export
calculatePerformance <- function (test_prediction,
                                  unique_clusters_sample,
                                  ...){

    list_TP = list()
    list_FP = list()
    list_FN = list()
    list_precision = list()
    list_recall = list()
    list_F1 = list()
    list_weight= list()

    for (i in 1:length(unique_clusters_sample)){
        curr_cluster = unique_clusters_sample[[i]]
        curr_test_predict = test_prediction[which(test_prediction$celltype %in% curr_cluster),]
        TP = length(which(curr_test_predict$celltype == curr_test_predict$mapped_celltype))
        FP= length(which(curr_test_predict$celltype !=curr_test_predict$mapped_celltype))

        curr_test_predict2 = test_prediction[which(test_prediction$mapped_celltype %in% curr_cluster),]
        FN = length(which(curr_test_predict2$celltype != curr_test_predict2$mapped_celltype))

        curr_F1 = 2*TP/(2*TP + FP + FN)
        curr_precision = TP/(TP+FP)
        curr_recall = TP/(TP+FN)

        curr_df = data.frame(cluster = curr_cluster,
                             precision = curr_precision,
                             recall = curr_recall,
                             F1 = curr_F1)

        list_TP[[i]] = TP
        list_FP[[i]] = FP
        list_FN[[i]] = FN
        list_precision[[i]] =curr_precision
        list_recall[[i]] =curr_recall
        list_F1[[i]] = curr_F1
        list_weight[[i]] = dim(curr_test_predict2)[1]

        if (i==1){
            performance_cluster_df = curr_df
        } else{
            performance_cluster_df = rbind(performance_cluster_df,curr_df)

        }
    }


    performance_all_df = data.frame(
        precision_weighted = weighted.mean(unlist(list_precision),unlist(list_weight)),
        precision_macro = mean(unlist(list_precision)),
        recall_weighted  = weighted.mean(unlist(list_recall),unlist(list_weight)),
        recall_macro = mean(unlist(list_recall)),
        F1_macro = mean(unlist(list_F1)),
        F1_weighted  = weighted.mean(unlist(list_F1),unlist(list_weight)),
        precision_micro = sum(unlist(list_TP))/(sum(unlist(list_TP)+sum(unlist(list_FP))))
        # recall_micro = sum(unlist(list_TP))/(sum(unlist(list_TP)+sum(unlist(list_FN))))
    )

    return(list(performance_cluster_df,performance_all_df))
}
