#' Find the most informative cluster markers
#'
#' @param matrix_all Feature matrix with cells as rows, and features as columns.
#' @param clusters Cluster annotation for each cell.
#' @param clusters_sel Selected clusters
#' @param num_markers Number of markers to output.
#' @param subsample_num Number of cells to sub-sample
#' @param sub.seed Seed number of sub-sample
#' @param train_test_ratio Ratio of training to test data
#' @param method_cluster List of methods to find cluster markers.
#' \itemize{
#'   \item \code{citeFUSE}
#'   \item \code{sc2marker}
#'   \item \code{geneBasis}
#'   \item \code{xgBoost}
#'   \item \code{all}: Use all methods
#' }
#' @param method_performance List of methods to evaluate performance.
#' \itemize{
#'   \item \code{xgBoost}
#'   \item \code{geneBasis}
#'   \item \code{all}: Use all methods
#' }
#'
#' @param cluster_proportional Sub-sampling equally or proportional
#' \itemize{
#'   \item \code{equal}
#'   \item \code{proportional}
#' }
#' @return A list containing
#' \itemize{
#'   \item \code{markers}: Informative markers for each method
#'   \item \code{performance}: Performance of each method
#' }
#' @export
findClusterPerformance <- function (matrix_all,
                                    clusters,
                                    clusters_sel = "ALLCLUSTER",
                                    num_markers,
                                    subsample_num,
                                    sub.seed = 42,
                                    train_test_ratio = 0.9,
                                    method_cluster = "all",
                                    method_performance = "all",
                                    cluster_proportion = "proportional",
                                    verbose = FALSE,
                                    ...) {

    # print("findClusterPerformance")
    # print(dim(matrix_all))

    if (length(clusters) != dim(matrix_all)[1]){
        stop("Number of cluster annotation cells is not equal to the number of cells in feature matrix.")
    }

    if (length(clusters_sel) == 1 && clusters_sel == "ALLCLUSTER"){
        clusters_sel = unique(clusters)
    }

    if (length(setdiff(clusters_sel,
                       unique(clusters)
    ))>0){
        warning(paste0("Invalid cluster name. Using all clusters"))
        clusters_sel = unique(clusters)
    }  else{
        index_keep = which(clusters %in% clusters_sel)
        clusters = clusters[index_keep]
        matrix_all = matrix_all[index_keep,]
    }

    if (subsample_num > nrows(matrix_all)){
        warning(paste0("Number of sub-samples more than number of cells. Using all cells."))


        subsample_num=nrow(matrix_all)


    }
    set.seed(sub.seed)


    unique_clusters = unique(clusters)



    if (cluster_proportion == "proportional"){
        temp_table = table(clusters)
        temp_table = temp_table[match(unique_clusters,
                                      names(temp_table))]
        cluster_ratio = as.numeric(temp_table / sum(temp_table))
    } else{
        cluster_ratio = rep(1 / length(unique_clusters),
                            length(unique_clusters))
    }


    list_index_training = list()
    list_index_test= list()
    icount=1
    for (i in 1:length(unique_clusters)){
        curr_cluster = unique_clusters[[i]]
        index_temp = which(clusters %in% curr_cluster)

        # index of cells for training and testing
        # index_training_temp = index_temp[1:round(train_test_ratio*length(index_temp))]
        # index_test_temp = setdiff(index_temp,
        #                           index_training_temp)

        # index of cells for training and testing
        index_test_temp = index_temp[1:ceiling((1 - train_test_ratio) * length(index_temp))]
        index_training_temp = setdiff(index_temp,
                                      index_test_temp)

        # subsample
        index_training_sample_temp = index_training_temp[1:min(max(ceiling(train_test_ratio*subsample_num*cluster_ratio[[i]]),
                                                                   4),
                                                               length(index_training_temp))]
        index_test_sample_temp = index_test_temp[1:min(max(ceiling((1-train_test_ratio)*subsample_num*cluster_ratio[[i]]),4),
                                                       length(index_test_temp))]

        if (length(index_test_sample_temp) > 3 & length(index_training_sample_temp) > 3){

            list_index_training[[icount]] = index_training_sample_temp
            list_index_test[[icount]] = index_test_sample_temp
            icount = icount + 1
        }


    }

    index_train = unlist(list_index_training)
    index_test = unlist(list_index_test)

    sample_index = c(index_train,index_test)
    matrix_sample = matrix_all[sample_index,]
    clusters_sample = clusters[sample_index]

    # create numeric form of clusters, used in XgBoost
    unique_clusters = unique(clusters_sample)
    num_clust = length(unique_clusters)
    label <- 0:(num_clust - 1)
    names(unique_clusters) = label



    input_matrix_train = matrix_all[index_train,]
    clusters_train = clusters[index_train]
    clusters_num_train = unlist(lapply(clusters_train,
                                       function (x) as.numeric(names(unique_clusters)[which(as.character(unique_clusters) %in% x)])))


    input_matrix_test = matrix_all[index_test,]
    clusters_test = clusters[index_test]
    clusters_num_test= unlist(lapply(clusters_test,
                                     function (x) as.numeric(names(unique_clusters)[which(as.character(unique_clusters) %in% x)])))


    if (length(unique(clusters_test))<2){
        stop("Not enough data in test set. Increase sample size")

    }
    if (length(unique(clusters_train))<2){
        stop("Not enough data in training set. Increase sample size")

    }


    if (verbose){
        print(table(clusters_test))
        print(table(clusters_train))


    }

    list_markers = findClusterMarkers(t(as.matrix(input_matrix_train)),
                                      clusters_train,
                                      num_markers,
                                      method=method_cluster,
                                      verbose=verbose)


    # calculate most occurring markers
    for (i in 1:length(list_markers)){
        curr_df = data.frame(markers = list_markers[[i]],
                             method = names(list_markers)[[i]])
        if (i == 1){
            total_df = curr_df
        } else{
            total_df = rbind(total_df,curr_df)
        }

    }

    table_compare = data.frame(table(total_df$markers))
    table_compare = table_compare[order(table_compare$Freq,decreasing = TRUE),]
    list_markers[["consensus"]] = as.character(table_compare$Var1[1:num_markers])


    if (verbose){
        print("Find performance metrics")
    }

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


        list_performance[[names(list_markers)[[i]]]] = performanceMarkers(markers_sel,
                                                                        t(as.matrix(input_matrix_train)),
                                                                        t(as.matrix(input_matrix_test)),
                                                                        unique_clusters,
                                                                        clusters_num_train,
                                                                        clusters_num_test,
                                                                        clusters_train,
                                                                        clusters_test,
                                                                        method = method_performance,
                                                                        nrounds = 1500,
                                                                        nthread = 6,
                                                                        verbose = verbose)

    }

    output = list()
    output$markers = list_markers
    output$performance = list_performance

    return(output)

}
