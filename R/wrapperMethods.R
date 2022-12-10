#' Wrapper function for citeFUSE
#'
#' @param sce Single cell experiment object
#' @param num_markers Number of markers to output
#'
#' @return The most informative markers determined by citeFUSE
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
#' @param sce Single cell experiment object
#' @param num_markers Number of markers to output
#'
#' @return The most informative markers determined by sc2marker
#' @export
sc2markerWrapper <- function (input_matrix,
                              clusters,
                              num_markers=15, ...){

    seurat_object = Seurat::CreateSeuratObject(input_matrix,
                                               meta.data =data.frame(cell_type=clusters) )
    names(clusters) = colnames(input_matrix)
    Seurat::Idents(object = seurat_object)=clusters
    all.markers <- sc2marker::Detect_single_marker_all(seurat_object, ...)

    unique_clusters = names(all.markers)
    num_clusters = length(unique_clusters)

    list_markers= list()
    for (i in 1:num_markers){
        ii = (i-1) %% num_clusters
        print(ii)
        curr_df = all.markers[[ii+1]]
        curr_df = curr_df[which(curr_df$direction %in% "+"),]
        index_remove = which(curr_df$gene %in% unlist(list_markers))
        if (length(index_remove)>0){
            curr_df = curr_df[-index_remove,]
        }
        list_markers[[i]] = curr_df$gene[[1]]
    }
    return(unlist(list_markers))

}

#' Wrapper function for geneBasis
#'
#' @param sce Single cell experiment object
#' @param num_markers Number of markers to output
#'
#' @return The most informative markers determined by geneBasis
#' @export
geneBasisWrapper <- function (sce,
                              clusters,
                              num_markers=15,
                              ...){

    sce = geneBasisR::retain_informative_genes(sce,...)
    marker_output = geneBasisR::gene_search(sce, n_genes_total = num_markers,...)
    return(marker_output$gene)

}
