#' Process input format (feature matrix, Seurat object, SCE object) to a format usable in other functions
#'
#' @param sc_object Input matrix and cell type object (required). Can be
#' \itemize{
#'   \item \code{Matrix:} Feature matrix (cells as rows, features as columns). cluster annotation has to be provided in `clusters_all` input.
#'   \item \code{SCE object}: SCE object, with the cluster annotation specified by the `sce_cluster` column name in the `colData` slot.
#'   \item \code{Seurat object}: Seurat object, with the cluster annotation in the `active.idents` slot.
#' }
#' @param clusters_all Cluster annotation, used if the input `sc_object` is a feature matrix.
#' @param sce_cluster Cluster annotation column name in the 'colData' slot, to be used if the input `sc_object` is a SCE object.
#' @param seurat_assay Seurat assay slot, if the input is a Seurat object. If not specified, default assay is used.
#' @param seurat_data Seurat data slot, if the input is a Seurat object. If not specified, the `data` slot will be used.
#'
#' @return List containing
#' \itemize{
#'   \item \code{matrix}: Feature matrix (cells as rows, features as columns)
#'   \item \code{clusters}: Cluster annotation vector
#' }

#' @export
processInputFormat =function(sc_object,
                             clusters_all=NULL,
                             sce_cluster = NULL,
                             seurat_assay=NULL,
                             seurat_slot=NULL,
                             verbose=TRUE,
                             ...){

    if (inherits(x = sc_object, what = c("matrix", "Matrix", "dgCMatrix"))) {
        if (verbose){
            message("Matrix input.\n")
        }
        sc_matrix <- sc_object
        if (is.null(clusters_all)){
            stop("Please provide a cluster annotation.")
        }
        sc_clusters = clusters_all

        if (length(clusters_all) != dim(sc_matrix)[1]){

            if (length(clusters_all) == dim(sc_matrix)[2]){
                warning(cat("The feature matrix will be transposed to have cells as rows, features as columns.\n"))
                sc_matrix= t(as.matrix(sc_matrix))
            } else{
                stop("Number of annotated cells is not equal to the number of cells in feature matrix. Please check the length of the cluster vector or the dimension of the feature matrix (cells as rows, features as columns).")
            }
        }
    }
    # Seurat object as input
    if (is(sc_object,"Seurat")) {
        if (verbose){
            cat("Seurat input.\n")
        }
        if (!requireNamespace("Seurat", quietly = TRUE)) {
            stop("Please install Seurat.")
        }

        if (is.null(seurat_assay)) {
            seurat_assay = DefaultAssay(sc_object)
            cat(paste0("The default Seurat assay used will be `", seurat_assay, "`. Please change the default assay if this is not correct."),'\n')
        } else{
            seurat_assay = DefaultAssay(sc_object)
        }

        if (is.null(seurat_slot)) {
            seurat_slot= "data"
            sc_matrix <- t(as.matrix(Seurat::GetAssayData(sc_object, assay = seurat_assay, slot = seurat_slot)))
            cat(paste0("The default Seurat slot used will be `", seurat_slot, "`. Please provide an alternative slot input if this is not correct."),'\n')
        } else{
            seurat_slot= "data"
            sc_matrix <- t(as.matrix(Seurat::GetAssayData(sc_object, assay = seurat_assay, slot = seurat_slot)))
        }

        cat("The active idents will be used for the cluster annotation. If this is not correct, please change the active idents to the correct annotation.")

        sc_clusters=Seurat::Idents(sc_object)
    }
    # SingleCellExperiment object as input
    if (is(sc_object,"SingleCellExperiment")) {
        if (verbose){
            cat("SCE input.\n")
        }
        if ("counts" %in% SummarizedExperiment::assayNames(sc_object)) {
            cat("The `counts` assay is used.",'\n')
            sc_matrix <- t(as.matrix(SingleCellExperiment::counts(sc_object)))
        } else {
            stop("SingleCellExperiment object must contain an assay named `counts`.")
        }
        if (is.null(sce_cluster)) {
            stop("Please specify the column name in the `colData` assay in the SingleCellExperiment object corresponding to the cluster annotation",'\n')
        }
        sc_clusters = sc_object@colData[,sce_cluster]
    }

    if (verbose){
        cat("Feature matrix dimension:", dim(sc_object)[1],"x",dim(sc_object)[2], ". Cluster annotation vector length:", length(sc_clusters),".")
    }

    sc_out = list(matrix=sc_matrix,
                  clusters=sc_clusters)

    return(sc_out)

}

#' Select clusters to compare with (based on user input)
#'
#' @param sc_out Output of function `processInputFormat`.
#' @param clusters_sel Subset of clusters to identify. Default is to use all clusters.
#'
#' @return List containing
#' \itemize{
#'   \item \code{matrix}: Feature matrix (cells as rows, features as columns)
#'   \item \code{clusters}: Cluster annotation vector
#' }
#' @export
processClusterSelection =function(sc_out,
                                  clusters_sel="all_clusters",
                                  verbose=TRUE,
                                  ...){

    sc_clusters = sc_out$clusters
    sc_matrix = sc_out$matrix

    if (length(clusters_sel)==1 && clusters_sel == "all_clusters"){
        if (verbose){
            cat("Using all clusters.\n")
        }
        clusters_sel = unique(sc_clusters)
    }
    if (length(clusters_sel)==1){
        stop("Please select more than one cluster.")
    }
    if (length(setdiff(clusters_sel,
                       unique(sc_clusters)))>0){
        warning(cat("Invalid cluster name. Using all clusters."))
        clusters_sel = unique(sc_clusters)
    }  else{
        index_keep = which(sc_clusters %in% clusters_sel)
        sc_clusters = sc_clusters[index_keep]
        sc_matrix = sc_matrix[index_keep,]
    }


    if (verbose){
        cat("Using clusters:",paste(unique(sc_clusters),collapse=", "),"\n")
    }

    cluster_selection_out = list(matrix=sc_matrix,
                                 clusters=sc_clusters)

    return(cluster_selection_out)

}


#' Sub-sample and split data into training and test set
#'
#' @param Output of function `processInputFormat`.
#' @param subsample_num Number of cells after sub-sammpling.
#' @param train_test_ratio Training to test data ratio.
#' @param cluster_proportion
#' \itemize{
#'   \item \code{proportional:} (default) Same proportion of cells in each cluster, for the training and data sets, compared to the original cluster proportion.
#'   \item \code{equal}: Equal proportion of cells in each cluster, for the training and data sets, compared to the original cluster proportion.
#' }
#'
#' @return List containing
#' \itemize{
#'   \item \code{training_matrix}:
#'   \item \code{training_clusters}:
#'   \item \code{training_clusters_num}:
#'   \item \code{test_matrix}:
#'   \item \code{test_clusters}:
#'   \item \code{test_clusters_num}:
#' }
#' @export
processSubsampling =function(cluster_selection_out,
                             subsample_num=1000,
                             train_test_ratio = 0.9,
                             cluster_proportion= "proportional",
                             verbose=TRUE,
                             ...){

    clusters_all = cluster_selection_out$clusters
    matrix_all = cluster_selection_out$matrix

    if (train_test_ratio <= 0 | train_test_ratio>=1){
        stop("Please choose a training:test ratio between 0 and 1.")
    }

    numCells = dim(matrix_all)[1]
    if (subsample_num > numCells){
        warning(paste0("Number of sub-samples more than number of cells. Using all cells."))
        subsample_num=numCells
    }

    unique_clusters_all = unique(clusters_all)
    if (cluster_proportion == "proportional"){
        temp_table = table(clusters_all)
        temp_table = temp_table[match(unique_clusters_all,
                                      names(temp_table))]
        cluster_ratio = as.numeric(temp_table/sum(temp_table))
    } else{
        cluster_ratio = rep(1/length(unique_clusters_all),
                            length(unique_clusters_all))
    }

    list_index_training = list()
    list_index_test= list()
    icount=1

    # Sub-sample for each cluster. Ensure each training and testing data has at least four cells in the cluster, otherwise errors occur downstream.
    for (i in 1:length(unique_clusters_all)){
        curr_cluster = unique_clusters_all[[i]]
        index_temp= which(clusters_all %in% curr_cluster)

        # index of cells for training and testing
        index_test_temp = index_temp[1:ceiling((1-train_test_ratio)*length(index_temp))]
        index_training_temp = setdiff(index_temp,
                                      index_test_temp)

        # Sub-sample to obtain training and test data
        index_training_sample_temp = index_training_temp[1:min(max(ceiling(train_test_ratio*subsample_num*cluster_ratio[[i]]),
                                                                   4),
                                                               length(index_training_temp))]
        index_test_sample_temp = index_test_temp[1:min(max(ceiling((1-train_test_ratio)*subsample_num*cluster_ratio[[i]]),4),
                                                       length(index_test_temp))]

        if (length(index_test_sample_temp)>3 & length(index_training_sample_temp)>3){
            list_index_training[[icount]]=index_training_sample_temp
            list_index_test[[icount]]=index_test_sample_temp
            icount=icount+1
        }
    }

    index_train = unlist(list_index_training)
    index_test = unlist(list_index_test)

    sample_index = c(index_train,index_test)
    clusters_sample = clusters_all[sample_index]

    if (length(index_train)<10){
        stop("Not enough cells in dataset.")
    }

    if (length(index_test)<10){
        stop("Not enough cells in dataset.")
    }

    # Create numeric form of clusters, used in XgBoost
    unique_clusters_sample = unique(clusters_sample)
    num_clust= length(unique_clusters_sample)
    label <- 0:(num_clust-1)
    names(unique_clusters_sample) = label

    # Form training cluster and matrix
    input_matrix_train = matrix_all[index_train,]
    clusters_train = clusters_all[index_train]
    clusters_num_train = unlist(lapply(clusters_train,
                                       function (x) as.numeric(names(unique_clusters_sample)[which(as.character(unique_clusters_sample) %in% x)])))

    # Form test cluster and matrix
    input_matrix_test = matrix_all[index_test,]
    clusters_test = clusters_all[index_test]
    clusters_num_test= unlist(lapply(clusters_test,
                                     function (x) as.numeric(names(unique_clusters_sample)[which(as.character(unique_clusters_sample) %in% x)])))

    if (length(unique(clusters_test))<2){
        stop("Less than two clusters after sub-sampling in the testing data. Please increase sub-sample size.")

    }
    if (length(unique(clusters_train))<2){
        stop("Less than two clusters after sub-sampling in the training data. Please increase sub-sample size.")
    }

    if (verbose){
        print(table(clusters_test))
        print(table(clusters_train))
    }

    final_out = list(training_matrix=input_matrix_train,
                     training_clusters=clusters_train,
                     training_clusters_num = clusters_num_train,
                     test_matrix=input_matrix_test,
                     test_clusters=clusters_test,
                     test_clusters_num = clusters_num_test,
                     unique_clusters_sample=unique_clusters_sample)

    return(final_out)
}
