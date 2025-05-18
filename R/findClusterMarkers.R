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
    # all_methods = c("citeFuse","sc2marker","xgBoost","fstat",
    #                 "seurat_wilcox","seurat_bimod","seurat_roc","seurat_t","seurat_LR")
    
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
    list_errors = list()
    for (i in 1:length(method)) {
        curr_method = method[[i]]
        # curr_method_old = method_old[[i]]
        
        if (verbose){
            message(paste0("\nCaclulating markers using ", curr_method,".\n"))
        }
        
        flag=0
        if (curr_method == "citeFuse") {
            start_time <- Sys.time()
            
            # curr_markers=NULL
            tryCatch({
                curr_markers <- citeFuseWrapper(sce,
                                                num_markers,
                                                subsample = TRUE)
            }, error = function(e) {
                error_txt <<- paste0(curr_method, " has error: ", e$message)
                curr_markers <<- NULL  # Optional: to avoid using undefined variable later
                flag <<- 1
            })
            
            
        }
        
        if (curr_method == "sc2marker") {
            start_time <- Sys.time()
            
            
            tryCatch({
                curr_markers = sc2markerWrapper(input_matrix,
                                                clusters,
                                                num_markers)
            }, error = function(e) {
                error_txt <<- paste0(curr_method, " has error: ", e$message)
                curr_markers <<- NULL  # Optional: to avoid using undefined variable later
                flag <<- 1
            })
            
            
            end_time <- Sys.time()
            runtime_secs[i] <- as.numeric(end_time-start_time, units="secs")
            names(runtime_secs)[i] <- "sc2marker"
        }
        
        if (curr_method == "geneBasis") {
            start_time <- Sys.time()
            SingleCellExperiment::logcounts(sce) <- input_matrix-min(input_matrix)
            
            
            tryCatch({
                curr_markers = geneBasisWrapper(sce,
                                                clusters,
                                                num_markers)
            }, error = function(e) {
                error_txt <<- paste0(curr_method, " has error: ", e$message)
                curr_markers <<- NULL  # Optional: to avoid using undefined variable later
                flag <<- 1
            })
            
            
            end_time <- Sys.time()
            runtime_secs[i] <- as.numeric(end_time-start_time, units="secs")
            names(runtime_secs)[i] <- "geneBasis"
        }
        
        if (curr_method == "xgBoost") {
            start_time <- Sys.time()
            
            tryCatch({
                curr_markers =     curr_markers = xgBoostWrapper(t(input_matrix),
                                                                 clusters,
                                                                 num_markers)
            }, error = function(e) {
                error_txt <<- paste0(curr_method, " has error: ", e$message)
                curr_markers <<- NULL  # Optional: to avoid using undefined variable later
                flag <<- 1
            })
            
            end_time <- Sys.time()
            runtime_secs[i] <- as.numeric(end_time-start_time, units="secs")
            names(runtime_secs)[i] <- "xgBoost"
        }
        
        if (curr_method == "fstat") {
            start_time <- Sys.time()
            
            
            tryCatch({
                curr_markers = fstatWrapper(t(input_matrix),
                                            clusters,
                                            num_markers)
            }, error = function(e) {
                error_txt <<- paste0(curr_method, " has error: ", e$message)
                curr_markers <<- NULL  # Optional: to avoid using undefined variable later
                flag <<- 1
            })
            
            end_time <- Sys.time()
            runtime_secs[i] <- as.numeric(end_time-start_time, units="secs")
            names(runtime_secs)[i] <- "fstat"
        }
        
        if (length(grep("seurat",curr_method))>0){
            start_time <- Sys.time()
            curr_method2 = gsub("seurat_","",curr_method)
            
            tryCatch({
                curr_markers = seuratWrapper(input_matrix,
                                             clusters,
                                             num_markers,
                                             curr_method2)
            }, error = function(e) {
                error_txt <<- paste0(curr_method, " has error: ", e$message)
                curr_markers <<- NULL  # Optional: to avoid using undefined variable later
                flag <<- 1
            })
            
            end_time <- Sys.time()
            runtime_secs[i] <- as.numeric(end_time-start_time, units="secs")
            names(runtime_secs)[i] <- curr_method
            
        }
        
        
        if (flag == 0){
            list_markers[[curr_method]] = curr_markers
        } else{
            list_errors[[curr_method]] = error_txt
            message(error_txt)
        }
    }
    
    # if (verbose){
        message(cat("Method with no errors used: ",paste0(names(list_markers),collapse=", ")))
        message(cat("Methods with errors removed: ",paste0(names(list_errors),collapse=", ")))
    # }
    
    list_markers_temp = list_markers
    
    
    fstat=apply(t(input_matrix),2,function (x) na.omit(anova(aov(x~as.factor(clusters)))$"F value"))
    fstat <- fstat[order(unlist(fstat), decreasing = T)]
    
    for (i in 1:length(list_markers)){
        # Order markers according to fstat
        list_markers[[i]] = intersect(names(fstat),list_markers[[i]])
        
    }
    
    list_markers[["runtime_secs"]] <- runtime_secs
    
    return(list_markers)
    
}
