#' Find the most informative cluster markers
#'
#' @param input_matrix Feature matrix with cells as columns, and features as rows.
#' @param clusters Cluster annotation for each cell.
#' @param num_markers Number of markers to output.
#' @param method List of methods to find cluster markers.
#' \itemize{
#'   \item \code{citeFuse}
#'   \item \code{sc2marker}
#'   \item \code{geneBasis}
#'   \item \code{xgBoost}
#'   \item \code{all}: Use all methods
#' }
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
findClusterMarkers <- function (input_matrix,
                                clusters,
                                num_markers = 15,
                                method = "all",
                                verbose = FALSE,
                                ...) {

    input_matrix= t(as.matrix(input_matrix))

    num_markers_original = num_markers

    all_methods = c("citeFuse","sc2marker","geneBasis","xgBoost","fstat",
                    "seurat_wilcox","seurat_bimod","seurat_roc","seurat_t","seurat_LR")

    # # put to lower in case of user typo
    # method_old = all_methods
    #     all_methods=unlist(lapply(all_methods, tolower))
    # method=unlist(lapply(method, tolower))

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

    # # Recognize the input methods
    # if (length(method) == 0) {
    #   stop("No method or invalid method selected.")
    #
    # # } else if((length(method) == 1 & all(method=="all") | all(method %in% all_methods))) {
    # } else if (method=="all") {
    #   method = all_methods
    # } else {
    #   diff_methods = setdiff(method, all_methods)
    #   if(length(diff_methods) >  0 & length(diff_methods) < length(method)){
    #     warning(paste0(paste(diff_methods, collapse = ", "), " not found. Using the remaining method(s)."))
    #     method = intersect(method,all_methods)
    #   } else if(length(diff_methods) != 0){
    #     stop(paste0("No available method selected.\nPlease select at least one method from the following: ",
    #                 paste(all_methods, collapse = ", ")))
    #   }
    # }

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
            curr_markers = sc2markerWrapper(input_matrix,
                                            clusters,
                                            num_markers)
            end_time <- Sys.time()
            runtime_secs[i] <- as.numeric(end_time-start_time, units="secs")
            names(runtime_secs)[i] <- "sc2marker"
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



    # calculate consensus if more than one method chosen
    if (length(method)>1){

        i = length(runtime_secs)+1
        start_time <- Sys.time()
        list_markers[["consensus_naive"]] = calculateConsensus(list_markers,
                                                               t(input_matrix),
                                                               clusters,
                                                               num_markers=num_markers,
                                                               verbose=TRUE)
        end_time <- Sys.time()
        runtime_secs[i] <- as.numeric(end_time-start_time, units="secs")
        names(runtime_secs)[i] <- "consensus_naive"

        i = length(runtime_secs)+1
        start_time <- Sys.time()
        list_markers[["consensus_fstat"]] =  calculateConsensus(list_markers,
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
        list_markers[["consensus_xgboost"]] = calculateConsensus(list_markers,
                                                                 t(input_matrix),
                                                                 clusters,
                                                                 num_markers=num_markers,
                                                                 method = "xgboost",
                                                                 verbose=TRUE)
        end_time <- Sys.time()
        runtime_secs[i] <- as.numeric(end_time-start_time, units="secs")
        names(runtime_secs)[i] <- "consensus_xgboost"

    }

    list_markers[["runtime_secs"]] <- runtime_secs

    return(list_markers)

}
