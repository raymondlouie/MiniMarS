#' Find the most informative cluster markers
#'
#' @param matrix_all Feature matrix with cells as columns, and features as rows.
#' @param clusters_all Cluster annotation for each cell.
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
#' @return A list containing
#' \itemize{
#'   \item \code{markers}: Informative markers for each method
#'   \item \code{performance}: Performance of each method
#' }
#' @export
findClusterPerformance <- function (matrix_all,
                                    clusters_all,
                                    num_markers,
                                    subsample_num,
                                    sub.seed=42,
                                    train_test_ratio,
                                    method_cluster="all",
                                    method_performance="all",
                                    verbose=FALSE,
                                    ...) {

    # print("findClusterPerformance")
    # print(dim(matrix_all))

    if (subsample_num > dim(matrix_all)[1]){
        warning(paste0("Number of sub-samples more than number of cells. Using all cells"))
        subsample_num=dim(matrix_all)[1]

    }

    # subsample
    set.seed(sub.seed)
    sample_index = sample(x = 1:dim(matrix_all)[1],size = subsample_num,replace = FALSE)
    matrix_sample= matrix_all[sample_index,]
    clusters_sample = clusters_all[sample_index]


    # create numeric form of clusters, used in XgBoost
    unique_clusters_sample = unique(clusters_sample)
    num_clust= length(unique_clusters_sample)
    label <- 0:(num_clust-1)
    names(unique_clusters_sample) = label
    clusters_newlabel_sample = unlist(lapply(clusters_sample,
                                             function (x) as.numeric(names(unique_clusters_sample)[which(as.character(unique_clusters_sample) %in% x)])))


    # divide data into training and test sets
    index_train = sample(x = 1:dim(matrix_sample)[1],size = round(train_test_ratio*dim(matrix_sample)[1]),replace = FALSE)
    index_test = setdiff(1:dim(matrix_sample)[1] , index_train)

    input_matrix_train = matrix_all[index_train,]
    clusters_train = clusters_sample[index_train]
    clusters_num_train = clusters_newlabel_sample[index_train]

    clusters_train_df = data.frame(table(clusters_train))
    remove_cells = which(clusters_train %in% clusters_train_df$clusters_train[which(clusters_train_df$Freq <6)])
    if (length(remove_cells)>0){
        input_matrix_train = input_matrix_train[-remove_cells,]
        clusters_train = clusters_train[-remove_cells]
        clusters_num_train = clusters_num_train[-remove_cells]
    }

    input_matrix_test = matrix_all[index_test,]
    clusters_test = clusters_sample[index_test]
    clusters_num_test = clusters_newlabel_sample[index_test]

    clusters_test_df = data.frame(table(clusters_test))
    remove_cells = which(clusters_test %in% clusters_test_df$clusters_test[which(clusters_test_df$Freq <6)])
    if (length(remove_cells)>0){
        input_matrix_test = input_matrix_test[-remove_cells,]
        clusters_test = clusters_test[-remove_cells]
        clusters_num_test = clusters_num_test[-remove_cells]
    }

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


    # calculate most occuring markers
    for (i in 1:length(list_markers)){
        curr_df = data.frame(markers = list_markers[[i]],
                             method = names(list_markers)[[i]])
        if (i==1){
            total_df = curr_df
        } else{
            total_df = rbind(total_df,curr_df)
        }

    }

    table_compare = data.frame(table(total_df$markers))
    table_compare = table_compare[order(table_compare$Freq,decreasing=TRUE),]
    list_markers[["consensus"]] = as.character(table_compare$Var1[1:num_markers])


    if (verbose){
        print("Find performance metrics")
    }

    list_performance = c()
    for (i in 1:length(list_markers)){
        markers_sel = list_markers[[i]]
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
                                                                        unique_clusters_sample,
                                                                        clusters_num_train,
                                                                        clusters_num_test,
                                                                        clusters_train,
                                                                        clusters_test,
                                                                        method=method_performance,
                                                                        nrounds=1500,
                                                                        nthread=6,
                                                                        verbose=verbose)

    }

    output = list()
    output$markers = list_markers
    output$performance = list_performance

    return(output)

}
