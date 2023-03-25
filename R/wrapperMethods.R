#' Wrapper function for CiteFuse
#'
#' @param sce Single cell experiment object
#' @param num_markers Number of markers to output
#'
#' @return The most informative markers determined by CiteFuse
#' @export
citeFuseWrapper <- function (sce,
                             num_markers=15,
                             ...){

    sce_alt <- SummarizedExperiment(list(raw=sce@assays@data$counts))
    altExp(sce, "protein") <- sce_alt

    set.seed(2020)
    sce <- CiteFuse::importanceADT(sce,
                                   group = factor(sce$cell_type),
                                   altExp_name ="protein",
                                   exprs_value = "raw")

    importance_scores = sort(sce@metadata$importanceADT,decreasing=TRUE)
    return(names(importance_scores[1:num_markers]))

}

#' Wrapper function for sc2marker
#'
#' @param input_matrix Marker matrix (cell vs markers)
#' @param clusters Cell type annotation
#' @param num_markers Number of markers to output
#'
#' @return The most informative markers determined by sc2marker
#' @export
sc2markerWrapper <- function (input_matrix,
                              clusters,
                              num_markers=15,
                              ...){


    seurat_object = Seurat::CreateSeuratObject(input_matrix,
                                               meta.data =data.frame(cell_type=clusters) )
    Seurat::Idents(object = seurat_object)=clusters
    all.markers <- sc2marker::Detect_single_marker_all(seurat_object, ...)

    unique_clusters = names(all.markers)
    num_clusters = length(unique_clusters)

    # old code
    # list_markers= list()
    # for (i in 1:num_markers){
    #     ii = (i-1) %% num_clusters
    #     curr_df = all.markers[[ii+1]]
    #     curr_df = curr_df[which(curr_df$direction %in% "+"),]
    #     index_remove = which(curr_df$gene %in% unlist(list_markers))
    #     if (length(index_remove)>0){
    #         curr_df = curr_df[-index_remove,]
    #     }
    #     list_markers[[i]] = curr_df$gene[[1]]
    # }
    message(num_markers)
    icount=1
    list_markers= list()
    for (i in 1:num_markers){
        ii = (i-1) %% num_clusters
        curr_df = all.markers[[ii+1]]
        curr_df = curr_df[which(curr_df$direction %in% "+"),]
        index_remove = which(curr_df$gene %in% unlist(list_markers))
        if (length(index_remove)>0){
            curr_df = curr_df[-index_remove,]
        }
    
        if (dim(curr_df)[1]>0){
            list_markers[[icount]] = curr_df$gene[[1]]
            icount=icount+1
        }
    }
    return(unlist(list_markers))

}

#' Wrapper function for geneBasis
#'
#' @param sce Single cell experiment object
#' @param clusters Cell type annotation
#' @param num_markers Number of markers to output
#'
#' @return The most informative markers determined by geneBasis
#' @export
geneBasisWrapper <- function (sce,
                              clusters,
                              num_markers=15,
                              ...){

    sce = geneBasisR::retain_informative_genes(sce,
                                               ...)
    geneBasis_num_markers = dim(sce)[1]
    if (geneBasis_num_markers<num_markers){
        warning("\n Number of markers from geneBasis is less than the number of input markers. Reducing number of markers. \n")
        num_markers = geneBasis_num_markers-1
    }
    marker_output = geneBasisR::gene_search(sce, n_genes_total = num_markers,...)
    return(marker_output$gene)

}


#' Wrapper function for xgBoost
#'
#' @param input_matrix Marker matrix with cells as rows, and features as columns.
#' @param clusters Cell type annotation
#' @param num_markers Number of markers to output
#'
#' @return The most informative markers determined by xgBoost
#' @export
xgBoostWrapper <- function (input_matrix, clusters,num_markers, nrounds=1500,nthread=6, ...){

    unique_clusters = unique(clusters)
    num_clust= length(unique_clusters)
    label <- 0:(num_clust-1)
    names(unique_clusters) = label
    clusters_newlabel = unlist(lapply(clusters,
                                      function (x) as.numeric(names(unique_clusters)[which(as.character(unique_clusters) %in% x)])))

    # convert features to numbers, because xgb.importance seems to have trouble with greek letters
    marker_num = 1:dim(input_matrix)[2]
    names(marker_num) = colnames(input_matrix)
    colnames(input_matrix) = marker_num

    fstat=apply(input_matrix,2,function (x) na.omit(anova(aov(x~as.factor(clusters)))$"F value"))
    fstat <- fstat[order(fstat, decreasing = T)]
    markers_fstat <- names(fstat)[1:min(num_markers*3,dim(input_matrix)[2])]

    xgboost_train = xgboost::xgb.DMatrix(data=input_matrix[, markers_fstat],
                                         label=clusters_newlabel)

    xgb_params <- list("objective" = "multi:softprob",
                       "eval_metric" = "mlogloss",
                       "num_class" = num_clust)

    built.model <- xgboost::xgb.train(params = xgb_params,
                                      data = xgboost_train,
                                      nrounds = nrounds,
                                      nthread = nthread)
    importance.mt <- xgboost::xgb.importance(colnames(xgboost_train), model = built.model)
    markers_xgboost = importance.mt$Feature


    if (length(markers_xgboost)<num_markers){
        # warning(paste0("XgBoost produced less markers than selected: ", paste(markers_xgboost,collapse=", "),
        #                ". Adding additional markers from the fstat algorithm."))
        message("XgBoost produced less markers than selected. Adding additional markers from the fstat algorithm.")
        additional_markers = setdiff(markers_fstat,markers_xgboost)
        markers_xgboost = c(markers_xgboost,additional_markers)
        markers_xgboost = markers_xgboost[1:num_markers]
    }

    markers_xgboost_label = unlist(lapply(markers_xgboost,
                                          function (x) names(marker_num)[which(marker_num %in% x)]))

    return(markers_xgboost_label[1:num_markers])

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
xgboostPerformance <- function (markers_sel,
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

    pred_test <- predict(built.model,
                         newdata = xgboost_test)

    test_prediction <- matrix(pred_test, nrow = numberOfClasses,
                              ncol=length(pred_test)/numberOfClasses) %>%
        t() %>%
        data.frame() %>%
        mutate(true = clusters_num_test ,
               predict_num = max.col(., "last")-1) ##predicted

    test_prediction$true_lab <- clusters_test
    test_prediction$predict_lab <- unlist(lapply(test_prediction$predict_num,
                                                 function (x) unique_clusters_sample[which(names(unique_clusters_sample) %in% x)]))

    for (i in 1:length(unique_clusters_sample)){
        curr_cluster = unique_clusters_sample[[i]]
        curr_test_predict = test_prediction[which(test_prediction$true_lab %in% curr_cluster),]
        number_correct= which(curr_test_predict$true_lab==curr_test_predict$predict_lab)
        curr_TP = length(number_correct)/dim(curr_test_predict)[1]
        curr_df = data.frame(cluster = curr_cluster,
                             TP = curr_TP)
        if (i==1){
            performance_xgBoost_df = curr_df
        } else{
            performance_xgBoost_df = rbind(performance_xgBoost_df,curr_df)

        }
    }

    return(performance_xgBoost_df)

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
geneBasisPerformance <- function (markers_sel,
                                  input_matrix_test,
                                  clusters_test,
                                  unique_clusters_sample,
                                  ...){

    sce_test  <- SingleCellExperiment::SingleCellExperiment(list(counts=t(input_matrix_test)),
                                                            colData=data.frame(cell_type=clusters_test))
    # logcounts(sce_test) <- log2(t(input_matrix_test) + 1)
    logcounts(sce_test) <- t(input_matrix_test)

    cluster_map = geneBasisR::get_celltype_mapping(sce_test ,
                                                   genes.selection = markers_sel,
                                                   celltype.id = "cell_type",
                                                   return.stat = T)

    test_prediction=cluster_map$mapping
    test_stat = cluster_map$stat

    for (i in 1:length(unique_clusters_sample)){
        curr_cluster = unique_clusters_sample[[i]]
        curr_test_predict = test_prediction[which(test_prediction$celltype %in% curr_cluster),]
        number_correct= which(curr_test_predict$celltype==curr_test_predict$mapped_celltype)
        curr_TP = length(number_correct)/dim(curr_test_predict)[1]
        curr_df = data.frame(cluster = curr_cluster,
                             TP = curr_TP)
        if (i==1){
            performance_genebasis_df = curr_df
        } else{
            performance_genebasis_df = rbind(performance_genebasis_df,curr_df)

        }
    }

    return(performance_genebasis_df)

}
