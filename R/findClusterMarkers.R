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
                                method="all") {

    num_markers_original=num_markers
    # num_markers = min(2*num_markers,
    # dim(input_matrix)[1])

    all_methods = c("citeFuse","sc2marker","geneBasis")

    if (method == "all"){
        method = all_methods
    }

    diff_methods = setdiff(method,all_methods)

    if (length(diff_methods) >  0) {
        warning(paste0(method, " not found. Using remaining or all methods."))
        method = intersect(method,all_methods)
        if (length(method) == 0) {
            method = all_methods
        }
    }

    if (length(clusters) != dim(input_matrix)[2]) {
        stop("Number of clusters do not match the dimension of the input matrix.")
    }

    sce  <- SingleCellExperiment(list(counts=input_matrix),
                                 colData=data.frame(cell_type=clusters))
    logcounts(sce) <- log2(input_matrix + 1)

    list_markers = list()
    list_celltype_pred = list()
    list_celltype_min = list()
    list_celltype_median = list()
    list_celltype_stats = list()
    for (i in 1:length(method)) {
        curr_method = method[[i]]

        if (curr_method == "citeFuse") {
            curr_markers = citeFuseWrapper(sce,
                                           num_markers,
                                           subsample=TRUE)
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

        list_markers[[curr_method]] = curr_markers

        cluster_map = get_celltype_mapping(sce ,
                                           genes.selection = curr_markers[1:num_markers_original],
                                           celltype.id = "cell_type",
                                           return.stat = T)

        list_celltype_pred[[curr_method]]=cluster_map$stat
        list_celltype_min[[curr_method]]=min(cluster_map$stat$frac_correctly_mapped)
        list_celltype_median[[curr_method]]=median(cluster_map$stat$frac_correctly_mapped)
        list_celltype_stats[[curr_method]]=cluster_map$stat

    }

    output_list = list(markers = list_markers,
                       score_min =list_celltype_min,
                       score_median = list_celltype_median,
                       celltype_stat = list_celltype_stats)

    return(output_list)


}
