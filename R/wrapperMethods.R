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


#' Wrapper function for Seurat
#'
#' @param input_matrix Marker matrix (cell vs markers)
#' @param clusters Cell type annotation
#' @param num_markers Number of markers to output
#' @param method Differential algorithm used in Seruat's FindAllMarkers
#'
#' @return The most informative markers determined by Seurat, based on FC
#' @export
seuratWrapper <- function (input_matrix,
                           clusters,
                           num_markers=15,
                           method = "bimod",
                           ...){
    
    
    seurat_object = Seurat::CreateSeuratObject(input_matrix)
    Seurat::Idents(object = seurat_object)=clusters
    
    # Convert if assay is V5
    if(class(seurat_object@assays$RNA)[1]=="Assay5"){
        # seurat_object[["RNA"]] <- as(object = seurat_object[["RNA"]], Class = "Assay")
        seurat_object[["RNA"]] <- SeuratObject::CreateAssayObject(counts = seurat_object[["RNA"]]$counts)
    }
    
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
    
    
    print(dim(input_matrix))
    print(table(clusters))
    
    seurat_object = Seurat::CreateSeuratObject(input_matrix,
                                               meta.data =data.frame(cell_type=clusters) )
    
    Seurat::Idents(object = seurat_object)=clusters
    
    # Convert if assay is V5
    if(class(seurat_object@assays$RNA)[1]=="Assay5"){
        # seurat_object[["RNA"]] <- as(object = seurat_object[["RNA"]], Class = "Assay")
        seurat_object[["RNA"]] <- SeuratObject::CreateAssayObject(counts = seurat_object[["RNA"]]$counts)
    }
    
    all.markers <- sc2marker::Detect_single_marker_all(seurat_object, ...)
    
    unique_clusters = names(all.markers)
    num_clusters = length(unique_clusters)
    
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
#' 
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



#' Wrapper function for fstat
#'
#' @param input_matrix Marker matrix with cells as rows, and features as columns.
#' @param clusters Cell type annotation
#' @param num_markers Number of markers to output
#'
#' @return The most informative markers determined by fstat
#' @export
fstatWrapper <- function (input_matrix, clusters,num_markers,  ...){
    
    fstat=apply(input_matrix,2,function (x) na.omit(anova(aov(x~as.factor(clusters)))$"F value"))
    fstat <- fstat[order(unlist(fstat), decreasing = T)]
    
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


#' Wrapper function to calculate consensus
#' 
#' @param list_markers List of string vector, with each element of list corresponding to a marker method
#' @param input_matrix_train Training feature matrix (cells as rows, features as columns)
#' @param clusters_train Training cluster annotation vector
#' @param num_markers Number of markers 
#' @param method Consensus method, based on the methods in `list_markers`
#' \itemize{
#'   \item \code{other} (default) Consensus calculated based on majority rule.
#'   \item \code{weighted}: The contribution of each marker is weighted according `list_weight_num`
#'   \item \code{fstat}: The contribution of each marker is weighted according to the fstat algorithm.
#'   \item \code{xgBoost}: The contribution of each marker is weighted according xgBoost.
#' }
#' @param list_weight_num Weight of each marker. Only used if `method` is weighted.
#'
#' @return Consensus markers
#' @export

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
    table_compare$finalAdd = table_compare$Freq
    # rownames(table_compare) = NULL
    table_temp = table_compare
    
    if (method =="fstat"){
        if (verbose){
            message("Calculating consensus using majority and fstat to resolve ties.")
        }
        
        fstat=apply(input_matrix_train,2,
                    function (x) na.omit(anova(aov(x~as.factor(clusters_train)))$"F value"))
        temp_gain <- fstat[order(unlist(fstat), decreasing = T)]
        
        tempValue = rep(0,dim(table_compare)[1])
        common_names = intersect(table_compare$Var1,names(temp_gain))
        tempValue[match(common_names,table_compare$Var1)] = as.numeric(temp_gain[match(common_names,
                                                                                       names(temp_gain))])
        print(table_compare$Freq)
        print(tempValue)
        table_compare$finalAdd =table_compare$Freq + tempValue
        
        
        
    } else if (method=="xgBoost"){
        
        if (verbose){
            message("Calculating consensus using majority and xgBoost to resolve ties.")
        }
        
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
        
        
        
    } else if (method=="weighted"){
        
        if (verbose){
            message("Calculating consensus using weighted.")
            print(table_compare)
        }
        
        # table_compare= table_weighted$wt[which(table_weighted$var %in% table_compare$Var1)]
        table_compare=aggregate(x = list("wt" = total_df$weight),
                                by = list("var" = total_df$markers),
                                FUN = sum)
        colnames(table_compare) = c("Var1","finalAdd")
        print(table_compare)
        
        
        
        
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

#' Wrapper function to calculate consensus whilst selecting certain methods
#' 
#' @param list_markers_temp List of string vector, with each element of list corresponding to a marker method
#' @param final_out List of matrix and cluster information produced by function `processSubsampling`
#' @param num_markers Number of markers 
#' @param chosen_measure The performance measure used to choose the methods used in the consensus. Options are precision_weighted, precision_macro, recall_weighted, recall_macro, F1_macro, F1_weighted, precision_micro
#' @param list_performance_valid Performance of the chosen markers using the validation matrix. Each element in list corresponds to a method.
#' @param metric_thres Threshold above which a method will be considered in the consensus calculation, based on `chosen_measure`
#' @param metric_topnum Number of methods to consider in consensus (based on sorted `chosen_measure`).
#' @return Consensus markers
#' @export



calculateConsensus_wrap <- function(list_markers_temp,
                                    final_out,
                                    num_markers=num_markers,
                                    chosen_measure= "F1_macro",
                                    list_performance_valid=list(),
                                    metric_thres = 0,
                                    metric_topnum = 1,
                                    verbose=TRUE,
                                    ...){
    
    input_matrix = t(as.matrix(final_out$training_matrix))
    
    
    if (length(list_performance_valid)==0){
        list_performance_valid = performanceAllMarkers(list_markers_temp,
                                                       final_out=final_out,
                                                       # method="all",
                                                       num_markers=num_markers,
                                                       nrounds=1500,
                                                       nthread=6,
                                                       testSet = "valid",
                                                       verbose=FALSE)
        message("Validation list manually created")
        
    } else{
        message("Validation list provided")
        
    }
    clusters = final_out$training_clusters
    
    if (length(clusters) != dim(input_matrix)[2]) {
        stop("Number of clusters do not match the dimension of the input matrix.")
    }
    
    # print(colnames(input_matrix))
    if (is.null(colnames(input_matrix))){
        warning("No cell names in matrix. Manually assigning names. Please manually check the input matrix matches up with cluster input.\n")
        colnames(input_matrix) = 1:dim(input_matrix)[2]
        
    }
    names(clusters) = colnames(input_matrix)
    
    
    index_remove = grep("consensus",names(list_performance_valid))
    if (length(index_remove)>0){
        list_performance_valid = list_performance_valid[-index_remove]
    }
    
    list_weight = list()
    for (jj in 1:length(list_performance_valid)){
        curr_list_measures = list_performance_valid[[jj]]$xgBoost_performance_all
        list_weight[[names(list_performance_valid)[[jj]]]] = as.numeric(curr_list_measures[which(names(curr_list_measures) %in% chosen_measure)])
    }
    list_weight_num = as.numeric(list_weight)
    list_weight_num= list_weight_num/sum(list_weight_num)
    names(list_weight_num) = names(list_performance_valid)
    
    names(list_weight) <- names(list_performance_valid)
    v_weight <- do.call(c, list_weight)
    num_total_methods = length(v_weight)
    topMethods <- names(v_weight[order(v_weight, decreasing = T)][1:min(metric_topnum,length(v_weight))])
    thresMethods = names(v_weight)[which(v_weight>metric_thres)]
    keepMethods = intersect(topMethods,thresMethods)
    
    if (length(keepMethods)==0){
        keepMethods = names(v_weight[order(v_weight, decreasing = T)][1:min(1,length(v_weight))])
        # message(paste0("No methods passed threshold. Choosing top method: ",keepMethods))
        message(paste0("No methods passed threshold. Choosing top method"))
        
        # } else{
        # message(paste("Methods used to calculate consensus:",keepMethods,collapse=", "))
        
    }
    message(paste0("Methods used to calculate consensus using measure ",
                   chosen_measure, ", threshold: ",
                   metric_thres, 
                   " and keeping the top ",
                   min(metric_topnum,num_total_methods), 
                   " methods: ",
                   paste0(keepMethods,collapse=", ")))
    
    
    ##only use marker sets from keep_methods
    list_markers_temp <- list_markers_temp[names(list_markers_temp) %in% keepMethods]
    
    #the weight vector needs to be adjusted accordingly
    list_weight_num <- list_weight_num[names(list_weight_num) %in% keepMethods]
    
    
    # if (verbose){
    #     message(paste0("Weighted list is ", paste0(list_weight_num, collapse=", ")))
    # }
    
    i = 1
    runtime_secs = list()
    start_time <- Sys.time()
    list_markers= list()
    
    if (length(keepMethods)>1){
        
        list_markers[["consensus_weighted"]] = calculateConsensus(list_markers_temp,
                                                                  t(input_matrix),
                                                                  clusters,
                                                                  num_markers=num_markers,
                                                                  method = "weighted",
                                                                  list_weight_num = list_weight_num,
                                                                  verbose=TRUE)
        
        end_time <- Sys.time()
        runtime_secs[i] <- as.numeric(end_time-start_time, units="secs")
        names(runtime_secs)[i] <- "consensus_weighted"
        
        i = length(runtime_secs)+1
        start_time <- Sys.time()
        list_markers[["consensus_naive"]] = calculateConsensus(list_markers_temp,
                                                               t(input_matrix),
                                                               clusters,
                                                               num_markers=num_markers,
                                                               verbose=TRUE)
        end_time <- Sys.time()
        runtime_secs[i] <- as.numeric(end_time-start_time, units="secs")
        names(runtime_secs)[i] <- "consensus_naive"
        
        i = length(runtime_secs)+1
        start_time <- Sys.time()
        list_markers[["consensus_fstat"]] =  calculateConsensus(list_markers_temp,
                                                                t(input_matrix),
                                                                clusters,
                                                                num_markers=num_markers,
                                                                method = "fstat",
                                                                verbose=TRUE)
        end_time <- Sys.time()
        runtime_secs[i] <- as.numeric(end_time-start_time, units="secs")
        names(runtime_secs)[i] <- "consensus_fstat"
        
        i = length(runtime_secs)+1
        start_time <- Sys.time()
        list_markers[["consensus_xgboost"]] = calculateConsensus(list_markers_temp,
                                                                 t(input_matrix),
                                                                 clusters,
                                                                 num_markers=num_markers,
                                                                 method = "xgBoost",
                                                                 verbose=TRUE)
        end_time <- Sys.time()
        runtime_secs[i] <- as.numeric(end_time-start_time, units="secs")
        names(runtime_secs)[i] <- "consensus_xgboost"
        
    } else{
        message("Only one method chosen.")
        names(list_markers_temp) = paste0("consensusTop","_",names(list_markers_temp))
        list_markers = list_markers_temp
        runtime_secs = c("consensusTop" = 0)
        
    }
    
    fstat=apply(t(input_matrix),2,function (x) na.omit(anova(aov(x~as.factor(clusters)))$"F value"))
    fstat <- fstat[order(unlist(fstat), decreasing = T)]
    
    for (i in 1:length(list_markers)){
        # Order markers according to fstat
        list_markers[[i]] = intersect(names(fstat),list_markers[[i]])
        
    }
    
    list_markers[["runtime_secs"]] <- runtime_secs
    
    return(list_markers)
    
}


#' Find the minimum number of markers to satisfy performance threshold
#' 
#' @param final_out List of matrix and cluster information produced by function `processSubsampling`
#' @param list_markers_test List of markers to test 
#' @param chosen_measure The performance measure used to choose the methods used in the consensus. Options are precision_weighted, precision_macro, recall_weighted, recall_macro, F1_macro, F1_weighted, precision_micro
#' @param threshold Minimum threshold for `chosen_measure`
#' @return List containing 
#' - Markers corresponding to best method
#' - Performance of these markers
#' - Markers corresponding to all methods
#' - Performance of all methods
#' - Time of all methods
#' @export

minMarker <- function (final_out,
                       list_markers_test=c(5,10,15,20,25,30,40),
                       chosen_measure = "F1_macro",
                       threshold  = 0.8,
                       seed=44,
                       ...){
    
    list_all = list()
    
    for (i in 1:length(list_markers_test)){
        
        numMarkers = list_markers_test[[i]]
        
        list_markers_time = findClusterMarkers(final_out,
                                               num_markers = numMarkers,
                                               method = "all",
                                               verbose = TRUE)
        
        list_time = list_markers_time$runtime_secs
        names(list_time) = names(list_markers_time)[which(!(names(list_markers_time) %in% c("consensus",
                                                                                            "runtime_secs")))]
        list_markers = list_markers_time[which(!(names(list_markers_time) %in% c("runtime_secs")))]
        
        
        list_markers_time_consensus= calculateConsensus_wrap(list_markers,
                                                             final_out,
                                                             num_markers=numMarkers)
        
        list_time_all = c(list_time,list_markers_time_consensus$runtime_secs)
        
        list_markers_all = c(list_markers,list_markers_time_consensus)
        list_markers_all = list_markers_all[which(!(names(list_markers_all) %in% c("runtime_secs")))]
        
        
        
        list_performance_all = performanceAllMarkers(list_markers_all,
                                                     final_out = final_out,
                                                     method = "xgBoost",
                                                     nrounds = 1500,
                                                     nthread = 6,
                                                     verbose = TRUE)
        
        curr_performance = list_performance_all[[grep("Top",names(list_performance_all))]]
        curr_performance_metric = curr_performance$xgBoost_performance_all[[chosen_measure]]
        list_all[[i]]= list(markersTop = list_markers_all[[grep("Top",names(list_markers_all))]],
                            performanceTop = curr_performance,
                            markersAll=list_markers_all,
                            performanceAll=list_performance_all,
                            timeAll = list_time_all)
        if (curr_performance_metric>threshold){
            break;
        }
    }
    names(list_all) = paste0(list_markers_test[1:length(list_all)], " markers")
    
    if (curr_performance_metric>threshold){
        
        message(paste0("Threshold reached with ", numMarkers, " markers and ", chosen_measure,
                       " score of ", curr_performance_metric," (user threshold=", threshold,")."))
        
    } else{
        list_all= NULL
        
        message("Threshold not reached.")
        
    }
    return(list_all)
}

#' Find the minimum number of markers to satisfy performance threshold across all clusters
#' 
#' @param final_out List of matrix and cluster information produced by function `processSubsampling`
#' @param list_markers_test List of markers to test 
#' @param chosen_measure The performance measure used to choose the methods used in the consensus. Options are precision_weighted, precision_macro, recall_weighted, recall_macro, F1_macro, F1_weighted, precision_micro
#' @param threshold Minimum threshold over all clusters for `chosen_measure`
#' @param clusters_sel Clusters to consider.
#' @return List containing 
#' - Markers corresponding to best method
#' - Performance of these markers
#' - Markers corresponding to all methods
#' - Performance of all methods
#' - Time of all methods
#' 
#' @export

minMarker_clusters <- function (final_out,
                                list_markers_test=c(5,10,15,20,25,30,40),
                                chosen_measure = "F1",
                                threshold  = 0.8,
                                clusters_sel="all",
                                seed=44,
                                ...){
    
    list_all = list()
    for (i in 1:length(list_markers_test)){
        
        numMarkers = list_markers_test[[i]]
        
        list_markers_time = findClusterMarkers(final_out,
                                               num_markers = numMarkers,
                                               method = "all",
                                               verbose = TRUE)
        
        list_time = list_markers_time$runtime_secs
        names(list_time) = names(list_markers_time)[which(!(names(list_markers_time) %in% c("consensus",
                                                                                            "runtime_secs")))]
        list_markers = list_markers_time[which(!(names(list_markers_time) %in% c("runtime_secs")))]
        
        
        list_markers_time_consensus= calculateConsensus_wrap(list_markers,
                                                             final_out,
                                                             num_markers=numMarkers)
        
        list_time_all = c(list_time,list_markers_time_consensus$runtime_secs)
        
        list_markers_all = c(list_markers,list_markers_time_consensus)
        list_markers_all = list_markers_all[which(!(names(list_markers_all) %in% c("runtime_secs")))]
        
        
        
        list_performance_all = performanceAllMarkers(list_markers_all,
                                                     final_out = final_out,
                                                     method = "xgBoost",
                                                     nrounds = 1500,
                                                     nthread = 6,
                                                     verbose = TRUE)
        
        curr_performance = list_performance_all[[grep("Top",names(list_performance_all))]]
        curr_performance_metric = curr_performance$xgBoost_performance_cluster[,chosen_measure]
        names(curr_performance_metric) = curr_performance$xgBoost_performance_cluster$cluster
        
        if (clusters_sel == "all"){
            clusters_sel2 = curr_performance$xgBoost_performance_cluster$cluster
        } else{
            clusters_sel2 = intersect(clusters_sel,curr_performance$xgBoost_performance_cluster$cluster)
        }
        if (length(clusters_sel2)==0){
            warning("Please choose valid clusters.")
        }
        curr_performance_metric= curr_performance_metric[which(names(curr_performance_metric) %in% clusters_sel2)]
        message(paste0("Clusters selected: ", paste0(clusters_sel2,collapse=", ")))
        
        list_all[[i]]= list(markersTop = list_markers_all[[grep("Top",names(list_markers_all))]],
                            performanceTop = curr_performance,
                            markersAll=list_markers_all,
                            performanceAll=list_performance_all,
                            timeAll = list_time_all)
        
        if (min(curr_performance_metric)>threshold){
            break;
        }
        
    }
    names(list_all) = paste0(list_markers_test[1:length(list_all)], " markers")
    
    
    if (min(curr_performance_metric)>threshold){
        message(paste0("Threshold reached with ", numMarkers, " markers and ", chosen_measure,
                       " score (minimum over all clusters) of ", min(curr_performance_metric)," ( user threshold=", threshold,")."))
        
    } else{
        message("Threshold not reached.")
        list_all= NULL
    }
    
}    
