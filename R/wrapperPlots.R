#' Plot markers identified by each method on a heatmap
#'
#' @param list_markers List of selected markers from each method, from output of findClusterMarkers function.
#' @param text_size Text size
#' @param color_tile Color of each tile
#' @param fill_tile Fill of each tile
#' @param coord_ratio Heatmap tile coordinate ratio
#'
#' @return A ggplot output
#' @export
plotMarkers =function(list_markers,
                      text_size = 20,
                      color_tile = "blue",
                      fill_tile = "blue",
                      coord_ratio  = 0.2,
                      ...){
    
    # Impute the component that has less number of markers than the set number
    max_len <- max(lengths(list_markers))
    list_markers <- lapply(list_markers, `length<-`, max_len)
  
    # create a rectangular data frame 
    markers_df = reshape2::melt(data.frame(list_markers),
                                measure.vars = names(list_markers))
    table_markers_df = data.frame(table(markers_df$value))
    table_markers_df = table_markers_df[order(table_markers_df$Freq,decreasing=TRUE),]
    
    # Order markers according to most frequent.
    markers_df$value = factor(markers_df$value,
                              levels = rev(table_markers_df$Var1))
    markers_df = markers_df[!is.na(markers_df$value),]
    
    p1=ggplot(markers_df, aes(x = variable, y = value),...) +
        geom_tile(color = color_tile, fill = fill_tile) + theme_classic()+
        theme(axis.text = element_text(size = text_size),
              axis.title = element_text(size = text_size),
              axis.text.x = element_text(angle = 45, hjust = 1, size = text_size))+
        xlab("Method") + ylab("Feature markers") +
        coord_fixed(ratio=coord_ratio)
    
    return(p1)
}

#' Plot performance of each marker
#'
#' @param list_performance List of performance metrics, from output of performanceAllMarkers function.
#' @param tile_text_size Tile text size
#' @param text_size Text size
#'
#' @return A list of ggplot outputs
#' @export
plotPerformance <- function(list_performance,
                            metric="F1", #other options are "precision" and "recall"
                            tile_text_size = 3,
                            text_size = 8,
                            ...){
    
   
    #extracting performance_cluster
    use.ls <- list()
    for(i in 1:length(list_performance)){
        ls.here <- list_performance[[i]]
        ls.here[[1]]$marker_method <- names(list_performance)[i]
        use.ls[[i]] <- ls.here[[1]]
    }
    names(use.ls) <- names(list_performance)
    
    curr_plot <- do.call(rbind, use.ls)
    names(curr_plot)[1] <- c("Clusters")
    curr_plot[,2:4] <- apply(curr_plot[,2:4], 2, function(x){round(x, digits=3)})
    
    if(metric == "F1"){
        colnames(curr_plot)[which(colnames(curr_plot) == "F1")] <- "Metric"
    }else if(metric == "precision"){
        colnames(curr_plot)[which(colnames(curr_plot) == "precision")] <- "Metric"
    }else if(metric == "recall"){
        colnames(curr_plot)[which(colnames(curr_plot) == "recall")] <- "Metric"
    }
    
    
    p1 <- ggplot(curr_plot, aes(x = marker_method, y = Clusters, fill=Metric)) +
        geom_tile(color="black") + theme_bw()+
        scale_fill_gradient(low = "white", high = "red") +
        
        theme(axis.text = element_text(size = text_size),
              axis.title = element_text(size = text_size),
              axis.text.x = element_text(angle = 45, hjust = 1, size = text_size),
              plot.title = element_text(size = text_size),
              
              legend.position="none")+
        xlab("Method") + ylab("Clusters") +
        coord_fixed(ratio=0.2)+
        geom_text(aes(marker_method, Clusters, label = Metric), colour = "black", check_overlap = FALSE,
                  size = tile_text_size)  +
        ggtitle(paste0("performance by cluster (xgBoost): ", metric)) 
    # print(p1)
    # list_p1[[i]] = p1
    
    # }
    
    return(p1)
    # return(list_p1)
    
}
