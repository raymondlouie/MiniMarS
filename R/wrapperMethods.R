#' Wrapper function for CiteFuse
#'
#' @param sce Single cell experiment object
#' @param num_markers Number of markers to output
#'
#' @return The most informative markers determined by CiteFuse
#' @export
citeFuseWrapper <- function (sce,
                             num_markers=15,
                             min_freq=3,
                             ...){

    # Remove cells with very low library size, which causes issues in CiteFuse
    totalCount = colSums(sce@assays@data$counts)
    # index_remove = which(totalCount < quantile(totalCount,0.1))
    index_remove = which(totalCount < (-60))

    if (length(index_remove)>0){
        message(paste0(length(index_remove), " cell(s) with low library size have been removed.\n"))
        sce = sce[,-index_remove]
    }

    table_df = data.frame(table(sce$cell_type))

    cell_remove_type = table_df$Var1[which(table_df$Freq < min_freq)]
    index_remove= which(sce$cell_type %in% cell_remove_type)

    if (length(index_remove)>0){
        sce= sce[,-index_remove]
    }
    sce$cell_type  = droplevels(factor(sce$cell_type))
    sce_alt <- SummarizedExperiment::SummarizedExperiment(list(raw=sce@assays@data$counts))
    SingleCellExperiment::altExp(sce, "protein") <- sce_alt

    # set.seed(2020)
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
#' @param method Differential algorithm
#'
#' @return The most informative markers determined by sc2marker
#' @export
seuratWrapper <- function (input_matrix,
                           clusters,
                           num_markers=15,
                           method = "bimod",
                           ...){


    seurat_object = Seurat::CreateSeuratObject(input_matrix)
    Seurat::Idents(object = seurat_object)=clusters

    markers_df = Seurat::FindAllMarkers(seurat_object,
                                        slot="counts",
                                        test.use = method,
                                        only.pos=TRUE)

    num_markers_each = floor(num_markers/length(unique(clusters)))

    # Calculate initial markers
    markers_df1=markers_df %>%
        group_by(cluster) %>%
        slice_max(n = num_markers_each, order_by = avg_log2FC)

    # Calculate markers to fill up, based on the top FC

    markers_df2=markers_df %>%
        group_by(cluster) %>%
        slice_max(n = num_markers_each+1, order_by = avg_log2FC)

    markers_df3 = markers_df2[which(markers_df2$gene %in%
                                        setdiff(markers_df2$gene,markers_df1$gene)),]

    first_list = unique(markers_df1$gene)
    second_list = unique(markers_df3$gene)

    final_list = c(first_list,second_list[1:(num_markers-length(first_list))])


    return(final_list)

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


  # message("here")
  print(dim(input_matrix))
  print(table(clusters))

  seurat_object = Seurat::CreateSeuratObject(input_matrix,
                                             meta.data =data.frame(cell_type=clusters) )
  # print("here2")

  Seurat::Idents(object = seurat_object)=clusters
  # seurat_object@assays$RNA@counts = input_matrix
  # seurat_object@assays$RNA@data = input_matrix

  # message("here2")

  all.markers <- sc2marker::Detect_single_marker_all(seurat_object, ...)
  
  # message("here3")
  # all.markers = "test"
  # print("here")

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
      message(i)
      message(paste(curr_df$gene,collapse=", "))
      # message(paste0(capture.output(curr_df), collapse = "\n"))
      # print(curr_df)
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
geneBasisWrapper <- function (sce, clusters, num_markers=15, ...) {
    # Save the current warning setting
    original_warn_setting <- getOption("warn")

    # Set the warning level to 1
    options(warn = 1)
  
    num_markers_original = num_markers

    sce = geneBasisR::retain_informative_genes(sce, ...)
    cat("\n")
    geneBasis_num_markers = dim(sce)[1]
    if (geneBasis_num_markers <= num_markers) {
        num_markers = geneBasis_num_markers - 1
        warning("\nNumber of markers from geneBasis is no larger than the number of input markers.\nReducing the number of markers.\n", paste0(num_markers, " markers used now."))
            }
    marker_output = geneBasisR::gene_search(sce, n_genes_total = num_markers, ...)
    
    # Restore the original warning setting
    options(warn = original_warn_setting)
    
    return(c(marker_output$gene, rep(NA, num_markers_original - num_markers)))
}




#' Wrapper function for xgBoost
#'
#' @param input_matrix Marker matrix with cells as rows, and features as columns.
#' @param clusters Cell type annotation
#' @param num_markers Number of markers to output
#'
#' @return The most informative markers determined by xgBoost
#' @export
fstatWrapper <- function (input_matrix, clusters,num_markers,  ...){

    # unique_clusters = unique(clusters)
    # num_clust= length(unique_clusters)
    # label <- 0:(num_clust-1)
    # names(unique_clusters) = label
    # clusters_newlabel = unlist(lapply(clusters,
    #                                   function (x) as.numeric(names(unique_clusters)[which(as.character(unique_clusters) %in% x)])))

    # convert features to numbers, because xgb.importance seems to have trouble with greek letters
    # marker_num = 1:dim(input_matrix)[2]
    # names(marker_num) = colnames(input_matrix)
    # colnames(input_matrix) = marker_num

    fstat=apply(input_matrix,2,function (x) na.omit(anova(aov(x~as.factor(clusters)))$"F value"))
    fstat <- fstat[order(unlist(fstat), decreasing = T)]
    # markers_fstat <- names(fstat)[1:min(num_markers*3,dim(input_matrix)[2])]


    return(names(fstat)[1:num_markers])

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
    fstat <- fstat[order(unlist(fstat), decreasing = T)]
    markers_fstat <- names(fstat)[1:min(num_markers*3,length(fstat))]

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







#' Wrapper function for Consensus function

calculateConsensus <- function (list_markers,
                                input_matrix_train,
                                clusters_train,
                                num_markers=15,
                                method = "other",
                                list_weight_num = list_weight_num,
                                verbose=TRUE,
                                ...){


    for (i in 1:length(list_markers)){
        curr_df = data.frame(markers = list_markers[[i]],
                             method = names(list_markers)[[i]])
        if (method =="weighted"){
            curr_df$weight = list_weight_num[which(names(list_weight_num) %in% names(list_markers)[[i]])]
        } else{
            curr_df$weight=1
        }
        if (i==1){
            total_df = curr_df
        } else{
            total_df = rbind(total_df,curr_df)
        }

    }

    table_compare = data.frame(table(total_df$markers))
    table_compare = table_compare[order(table_compare$Freq,decreasing=TRUE),]
    table_compare$finalAdd = table_compare$Var1
    # rownames(table_compare) = NULL
    table_temp = table_compare

    if (method =="fstat"){
        fstat=apply(input_matrix_train,2,
                    function (x) na.omit(anova(aov(x~as.factor(clusters_train)))$"F value"))
        temp_gain <- fstat[order(unlist(fstat), decreasing = T)]

        tempValue = rep(0,dim(table_compare)[1])
        common_names = intersect(table_compare$Var1,names(temp_gain))
        tempValue[match(common_names,table_compare$Var1)] = temp_gain[match(common_names,
                                                                            names(temp_gain))]
        table_compare$finalAdd =table_compare$Freq + tempValue

        if (verbose){
            message("Calculating consensus using majority and fstat to resolve ties.")
        }


    } else if (method=="xgBoost"){

        unique_clusters = unique(clusters_train)
        num_clust= length(unique_clusters)
        label <- 0:(num_clust-1)
        names(unique_clusters) = label
        clusters_num_train = unlist(lapply(clusters_train,
                                           function (x) as.numeric(names(unique_clusters)[which(as.character(unique_clusters) %in% x)])))


        xgboost_train = xgboost::xgb.DMatrix(data=as.matrix(input_matrix_train[, table_compare$Var1]),
                                             label=clusters_num_train)

        # train a model using our training data
        numberOfClasses <- length(unique(clusters_num_train))

        xgb_params <- list("objective" = "multi:softprob",
                           "eval_metric" = "mlogloss",
                           "num_class" = numberOfClasses)

        built.model <- xgboost::xgb.train(params = xgb_params,
                                          data = xgboost_train,
                                          nrounds=2)

        importance_train.mt <- xgboost::xgb.importance(colnames(xgboost_train),
                                                       model = built.model)
        temp_gain = importance_train.mt$Gain
        names(temp_gain) = importance_train.mt$Feature

        tempValue = rep(0,dim(table_compare)[1])
        common_names = intersect(table_compare$Var1,names(temp_gain))
        tempValue[match(common_names,table_compare$Var1)] = temp_gain[match(common_names,
                                                                            names(temp_gain))]

        table_compare$finalAdd =table_compare$Freq + tempValue

        if (verbose){
            message("Calculating consensus using majority and xgBoost to resolve ties.")
        }

    } else if (method=="weighted"){

        # table_compare= table_weighted$wt[which(table_weighted$var %in% table_compare$Var1)]
        table_compare=aggregate(x = list("wt" = total_df$weight),
                                by = list("var" = total_df$markers),
                                FUN = sum)
        colnames(table_compare) = c("Var1","finalAdd")
        print(table_compare)

        if (verbose){
            message("Calculating consensus using weighted.")
            print(table_compare)
        }


    } else{
        if (verbose){
            message("Calculating consensus using simple majority with random ties.")
        }

    }
    table_compare = table_compare[order(table_compare$finalAdd,
                                        decreasing=TRUE),]

    consensus_markers = as.character(table_compare$Var1[1:num_markers])

    return(consensus_markers)

}
