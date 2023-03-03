#' Find the most informative cluster markers
#'
#' @param input_matrix Feature matrix with cells as columns, and features as rows.
#' @param clusters Cluster annotation for each cell.
#' @param num_markers Number of markers to output.
#' @param method List of methods to find cluster markers.
#' \itemize{
#'   \item \code{citeFUSE}
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

    all_methods = c("citeFuse","sc2marker","geneBasis","xgBoost")

    # Recognize the input methods
    if (length(method) == 0) {
      stop("No method selected.")

    } else if((length(method) == 1 & all(method=="all") | all(method %in% all_methods))) {
      method = all_methods
    } else {
      diff_methods = setdiff(method, all_methods)
      if(length(diff_methods) >  0 & length(diff_methods) < length(method)){
        warning(paste0(paste(diff_methods, collapse = ", "), " not found. Using the remaining method(s)."))
        method = intersect(method,all_methods)
      } else if(length(diff_methods) != 0){
        stop(paste0("No available method selected.\nPlease select at least one method from the following: ",
                    paste(all_methods, collapse = ", ")))
      }
    }

    # Print out the methods used after checks
    message(paste0("Methods used in this analysis: ", paste(method, collapse = ", ")))

    if (length(clusters) != dim(input_matrix)[2]) {
        stop("Number of clusters do not match the dimension of the input matrix.")
    }


    sce  <- SingleCellExperiment::SingleCellExperiment(list(counts=input_matrix),
                                                       colData=data.frame(cell_type=clusters))

    logcounts(sce) <- log2(input_matrix + 1)

    list_markers = list()

    for (i in 1:length(method)) {
        curr_method = method[[i]]

        if (verbose){
            message(paste0("Caclulating markers using ", curr_method,".\n"))
        }


        if (curr_method == "citeFuse") {
            curr_markers = citeFuseWrapper(sce,
                                           num_markers,
                                           subsample = TRUE)
        }

        if (curr_method == "sc2marker") {
            curr_markers = sc2markerWrapper(input_matrix,
                                            clusters,
                                            num_markers)
        }

        if (curr_method == "geneBasis") {
            curr_markers = geneBasisWrapper(sce,
                                            clusters,
                                            num_markers)
        }

        if (curr_method == "xgBoost") {
            curr_markers = xgBoostWrapper(t(input_matrix),
                                          clusters,
                                          num_markers)
        }

        list_markers[[curr_method]] = curr_markers

    }

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

    return(list_markers)

}
