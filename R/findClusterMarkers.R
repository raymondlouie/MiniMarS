#' Find the most informative cluster markers
#'
#' @param final_out List of matrix and cluster information produced by function `processSubsampling`
#' @param num_markers Number of markers 
#' @param method List of methods to find cluster markers.
#' \itemize{
#'   \item \code{citeFuse}
#'   \item \code{sc2marker}
#'   \item \code{geneBasis}
#'   \item \code{xgBoost}
#'   \item \code{seurat}
#'   \item \code{all}: Use all methods
#' }
#' 
#' @param metric_thres Methods to consider in consensus, which are above this threshold based on F1_macro
#' @param metric_topnum Number of methods to consider in consensus, based on sorted F1_macro.
#'
#' @return A list containing
#' \itemize{
#'   \item \code{markers}: Informative markers for each method
#'   \item \code{score_min}: Minimum of cell type scores as predicted by geneBasis
#'   \item \code{score_median}: Median of cell type scores as predicted by geneBasis
#'   \item \code{celltype_stat}: Cell type stats as predicted by geneBasis
#'   \item \code{runtime}: Computation time for each method
#' }
#' @export
findClusterMarkers <- function (final_out,
                                num_markers = 15,
                                method = "all",
                                metric_thres = 0,
                                metric_topnum = 2,
                                verbose = FALSE,
                                ...) {
    
    input_matrix = t(as.matrix(final_out$training_matrix))
    clusters = final_out$training_clusters
    
    num_markers_original = num_markers
    
    all_methods = c("citeFuse","sc2marker","geneBasis","xgBoost","fstat",
                    "seurat_wilcox","seurat_bimod","seurat_roc","seurat_t","seurat_LR")

    
    if ((length(method)==1 && method == "all") | (identical(sort(method), sort(all_methods)))){
        message("Using all methods.")
        method = all_methods
    } else{
        method = all_methods[which(all_methods %in% method)]
        
        # Recognize the input methods
        if (length(method) == 0) {
            warning("No method or invalid method selected. Using all methods.\n")
            method = all_methods
        } else {
            message(paste0("\nUsing the following method(s): ", paste(method, collapse = ", ")))
        }
    }
    

    # Print out the methods used after checks
    message(paste0("Methods used in this analysis: ", paste(method, collapse = ", "),"\n"))
    
    if (length(clusters) != dim(input_matrix)[2]) {
        stop("Number of clusters do not match the dimension of the input matrix.")
    }
    
    
    # print(colnames(input_matrix))
    if (is.null(colnames(input_matrix))){
        warning("No cell names in matrix. Manually assigning names. Please manually check the input matrix matches up with cluster input.\n")
        colnames(input_matrix) = 1:dim(input_matrix)[2]
        
    }
    names(clusters) = colnames(input_matrix)
    
    sce  <- SingleCellExperiment::SingleCellExperiment(list(counts=input_matrix),
                                                       colData=data.frame(cell_type=clusters))
    
    # logcounts(sce) <- log2(input_matrix + 1)
    SingleCellExperiment::logcounts(sce) <- input_matrix
    
    list_markers = list()
    runtime_secs <- c()
    
    for (i in 1:length(method)) {
        curr_method = method[[i]]
        # curr_method_old = method_old[[i]]
        
        if (verbose){
            message(paste0("\nCaclulating markers using ", curr_method,".\n"))
        }
        
        
        if (curr_method == "citeFuse") {
            start_time <- Sys.time()
            curr_markers = citeFuseWrapper(sce,
                                           num_markers,
                                           subsample = TRUE)
            end_time <- Sys.time()
            runtime_secs[i] <- as.numeric(end_time-start_time, units="secs")
            names(runtime_secs)[i] <- "citeFuse"
        }
        
        if (curr_method == "sc2marker") {
            start_time <- Sys.time()
            if (packageVersion("Seurat")>="5"){
                message("sc2marker doesn't work with Seurat v5. Please try Seurat v4.")
            } else{
                curr_markers = sc2markerWrapper(input_matrix,
                                                clusters,
                                                num_markers)
                end_time <- Sys.time()
                runtime_secs[i] <- as.numeric(end_time-start_time, units="secs")
                names(runtime_secs)[i] <- "sc2marker"
            }
        }
        
        if (curr_method == "geneBasis") {
            start_time <- Sys.time()
            SingleCellExperiment::logcounts(sce) <- input_matrix-min(input_matrix)
            
            curr_markers = geneBasisWrapper(sce,
                                            clusters,
                                            num_markers)
            end_time <- Sys.time()
            runtime_secs[i] <- as.numeric(end_time-start_time, units="secs")
            names(runtime_secs)[i] <- "geneBasis"
        }
        
        if (curr_method == "xgBoost") {
            start_time <- Sys.time()
            curr_markers = xgBoostWrapper(t(input_matrix),
                                          clusters,
                                          num_markers)
            end_time <- Sys.time()
            runtime_secs[i] <- as.numeric(end_time-start_time, units="secs")
            names(runtime_secs)[i] <- "xgBoost"
        }
        
        if (curr_method == "fstat") {
            start_time <- Sys.time()
            curr_markers = fstatWrapper(t(input_matrix),
                                        clusters,
                                        num_markers)
            end_time <- Sys.time()
            runtime_secs[i] <- as.numeric(end_time-start_time, units="secs")
            names(runtime_secs)[i] <- "fstat"
        }
        
        if (length(grep("seurat",curr_method))>0){
            start_time <- Sys.time()
            curr_method2 = gsub("seurat_","",curr_method)
            curr_markers = seuratWrapper(input_matrix,
                                         clusters,
                                         num_markers,
                                         curr_method2)
            end_time <- Sys.time()
            runtime_secs[i] <- as.numeric(end_time-start_time, units="secs")
            names(runtime_secs)[i] <- curr_method
            
        }
        
        
        
        list_markers[[curr_method]] = curr_markers
        
    }
    
    list_markers_temp = list_markers
    
    # calculate consensus if more than one method chosen
    if (length(method)>1){
        
        i = length(runtime_secs)+1
        start_time <- Sys.time()
        
        list_performance_valid = performanceAllMarkers(list_markers_temp,
                                                       final_out=final_out,
                                                       method="all",
                                                       nrounds=1500,
                                                       nthread=6,
                                                       testSet = "valid",
                                                       verbose=FALSE)
        
        chosen_measure = "F1_macro"
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
        topMethods <- names(v_weight[order(v_weight, decreasing = T)][1:min(metric_topnum,length(v_weight))])
        thresMethods = names(v_weight)[which(v_weight>metric_thres)]
        keepMethods = intersect(topMethods,thresMethods)
        
        ##only use marker sets from top two methods (with the highest chosen_measure)
        list_markers_temp <- list_markers_temp[names(list_markers_temp) %in% keepMethods]
        
        #the weight vector needs to be adjusted accordingly
        list_weight_num <- list_weight_num[names(list_weight_num) %in% keepMethods]
        
        
        if (verbose){
            message(paste0("Weighted list is ", list_weight_num))
        }
        
        
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
