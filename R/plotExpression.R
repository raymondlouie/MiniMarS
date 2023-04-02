#' Visualize gene expression data for a list of markers using UMAP or violin plot
#'
#' @param list_markers A list of gene markers to be plotted (required).
#' @param input_obj An input object containing gene expression data and cluster information (required).
#' @param plot_type The type of plot to generate, either "umap" or "violin" (default is "umap").
#' @param markers_plot A vector of markers to be plotted (default is NA, which will use all markers in the `list_markers` object).
#' @param order_method A method for ordering the markers (default is alphabetical).
#' @param cluster_levels A vector of cluster levels (default will use the unique clusters in the input object).
#' @param marker_name_plot A logical value indicating whether to display marker names on the plot (default is TRUE).
#' @param return_ggplot A logical value indicating whether to return a ggplot object (default is FALSE).
#' @param title_size The size of the title text (default is 30).
#' @param text_size The size of the text (default is 30).
#' @param point_size The size of the points in the plot (default is 3).
#' @param ... Additional arguments passed to the underlying plotting functions.
#'
#' @return If return_ggplot is TRUE, a ggplot object with the generated plot.
#'
#' @examples
#' # Given an input object with gene expression data and cluster information, and a list of markers:
#' plotExpression(list_markers, input_obj, plot_type = "umap")
#'
#' @export

plotExpression =function(list_markers,
                         input_obj,
                         plot_type = "umap",
                         markers_plot = NA,
                         order_method=NA,
                         cluster_levels=NA,
                         marker_name_plot=TRUE,
                         return_ggplot=FALSE,
                         title_size=30,
                         text_size=30,
                         point_size=3,
                         ...){

    if (plot_type == "umap"){
        matrix_log = log10(1+input_obj$matrix)

        umap<-data.frame(uwot::umap(as.matrix(matrix_log),
                                    n_neighbors=30,
                                    n_epochs=200,
                                    init = "spectral",
                                    n_threads = 8,
                                    min_dist = 0.1,
                                    verbose = TRUE))

        colnames(umap)= c("UMAP_1","UMAP_2")
        umap$clusters = factor(input_obj$clusters)


        if (is.na(cluster_levels)){
            cluster_levels = sort(unique(input_obj$clusters))
        }


        list_markers_df = reshape2::melt(list_markers)
        total_table_df = data.frame(table(unlist(list_markers)))
        total_table_df = total_table_df[order(total_table_df$Freq,decreasing=TRUE),]
        total_table_df$methods = unlist(lapply(total_table_df$Var1,
                                               function (x) paste0(unique(list_markers_df$L1[which(list_markers_df$value %in% x)]),collapse=", ")))


        if (is.na(markers_plot)){
            markers_plot = as.character(total_table_df$Var1[1:length(list_markers[[1]])])

        }
        total_table_df= total_table_df[which(total_table_df$Var1 %in% markers_plot),]

        plot_mat = as.data.frame(matrix_log[,as.character(total_table_df$Var1)],
                                 , check.names = FALSE)
        colnames(plot_mat)=as.character(total_table_df$Var1)

        if (marker_name_plot){
            tempName = unlist(lapply(colnames(plot_mat),
                                     function (x) paste0(x," \n (",
                                                         total_table_df$methods[which(as.character(total_table_df$Var1) %in% x)],
                                                         ")")))
            colnames(plot_mat) =tempName
        }

        plot_mat$clusters = input_obj$clusters
        plot_mat$UMAP_1 = umap$UMAP_1
        plot_mat$UMAP_2 = umap$UMAP_2


        plot_df_melt = reshape2::melt(plot_mat,
                            measure.vars = setdiff(colnames(plot_mat),c("UMAP_1","UMAP_2","clusters")))

        getPalette = colorRampPalette(brewer.pal(9, "Set1"))
        p_umap_cluster = ggplot(umap,aes_string(x="UMAP_1",y="UMAP_2",col="clusters"))+
            geom_point(show.legend = TRUE,
                       size=point_size*1.5
            ) +
            guides(colour = guide_legend(override.aes = list(size=7))) +
            scale_colour_manual(values=getPalette(length(unique(input_obj$clusters))))+
            theme_classic() +
            theme(legend.key.size= unit(1, "cm"),
                  axis.text=element_text(size=text_size*1.5),
                  axis.title=element_text(size=text_size*1.5),
                  legend.text = element_text(size = text_size),
                  legend.title = element_text(size = text_size))+
            labs(colour="Clusters")

        p1 = Seurat::LabelClusters(p_umap_cluster,
                                   id = "clusters",
                                   size=text_size/2)

        print(p1)

        p_umap = ggplot(plot_df_melt,aes_string(x="UMAP_1",y="UMAP_2",col="value"))+
            geom_point(show.legend = TRUE,
                       size=point_size
            ) +
            scale_colour_gradientn(name="", colours=c("blue","cyan","green","yellow","red"))+

            theme_classic() +
            theme(legend.key.size= unit(1, "cm"),
                  axis.text=element_text(size=30),
                  axis.title=element_text(size=30,face="bold"),
                  strip.text = element_text(size = title_size),
                  legend.text = element_text(size = text_size),
                  legend.title = element_text(size = text_size))+
            labs(colour="Expression (Log10)")+
            facet_wrap(~variable,ncol=2)

        print(p_umap)

        if (return_ggplot){
            return(p_umap)
        }

    }

    if (plot_type == "violin"){

        if (is.na(cluster_levels)){
            cluster_levels = sort(unique(input_obj$clusters))
        }


        list_markers_df = reshape2::melt(list_markers)
        total_table_df = data.frame(table(unlist(list_markers)))
        total_table_df = total_table_df[order(total_table_df$Freq,decreasing=TRUE),]
        total_table_df$methods = unlist(lapply(total_table_df$Var1,
                                               function (x) paste0(unique(list_markers_df$L1[which(list_markers_df$value %in% x)]),collapse=", ")))


        if (is.na(markers_plot)){
            markers_plot = as.character(total_table_df$Var1[1:length(list_markers[[1]])])

        }
        total_table_df= total_table_df[which(total_table_df$Var1 %in% markers_plot),]

        plot_mat = as.data.frame(log10(1+input_obj$matrix[,as.character(total_table_df$Var1)]),
                                 , check.names = FALSE)
        colnames(plot_mat)=as.character(total_table_df$Var1)

        # plot_df_truncate =plot_mat[,which(colnames(plot_mat) %in% as.character(total_table_df$Var1))]

        if (marker_name_plot){
            tempName = unlist(lapply(colnames(plot_mat),
                                     function (x) paste0(x," \n (",
                                                         total_table_df$methods[which(as.character(total_table_df$Var1) %in% x)],
                                                         ")")))
            colnames(plot_mat) =tempName
        }

        plot_mat$clusters = input_obj$clusters


        plot_df_melt = reshape2::melt(plot_mat,
                            measure.vars = setdiff(colnames(plot_mat),"clusters"))

        plot_df_melt$clusters = factor(plot_df_melt$clusters,
                                       levels = cluster_levels)

        p_violin= ggplot(plot_df_melt, aes_string(y = "value", x = "clusters", fill = "clusters")) +
            geom_violin(position = position_dodge(width = 4),scale = "width",draw_quantiles = c(0.5)) +
            geom_jitter(shape=16, position=position_jitter(0.2),size=2) +
            theme_bw()+
            theme(
                axis.title=element_text(size = text_size+20),
                axis.text.x=element_text(size=text_size,
                                         angle=85,
                                         hjust=1),
                legend.text = element_text(size = text_size),
                legend.title = element_text(size = text_size),
                axis.text.y=element_text(size=text_size),
                strip.text = element_text(size = title_size))+
            facet_wrap(~variable,ncol=2) +
            ylab("Expression (log10)")+
            labs(fill="Clusters") + xlab("Clusters")
        print(p_violin)

        if (return_ggplot){
            return(p_violin)
        }

    }
}
