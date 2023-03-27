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
                      text_size=40,
                      color_tile = "blue",
                      fill_tile = "blue",
                      coord_ratio  = 0.2,
                      ...){
    
    
    markers_df = reshape2::melt(data.frame(list_markers),
                                measure.vars = names(list_markers))
    table_markers_df = data.frame(table(markers_df$value))
    table_markers_df = table_markers_df[order(table_markers_df$Freq,decreasing=TRUE),]
    
    # Order markers according to most frequent.
    markers_df$value = factor(markers_df$value,
                              levels=rev(table_markers_df$Var1))
    markers_df = markers_df[!is.na(markers_df$value),]
    
    p1=ggplot(markers_df, aes(x = variable, y = value),...) +
        geom_tile(color=color_tile,fill=fill_tile) + theme_classic()+
        theme(axis.text = element_text(size=text_size),
              axis.title = element_text(size=text_size),
              axis.text.x = element_text(angle=45,hjust=1,size=text_size))+
        xlab("Method") + ylab("Feature markers") +
        coord_fixed(ratio=coord_ratio)
    print(p1)
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
plotPerformance =function(list_performance,
                          tile_text_size=7,
                          text_size=40,
                          ...){
    
    # Create dataframe suitable for plotting
    performance_df = data.frame(list_performance)
    
    for (i in seq(1,dim(performance_df)[2],2)){
        curr_df = performance_df[,c(i,i+1)]
        tempSplit = unlist(lapply(as.character(colnames(curr_df)),
                                  function (x) strsplit(x,split="[.]")[[1]]))
        colnames(curr_df) = c("Clusters","TP")
        curr_df$marker_method=tempSplit[[1]]
        curr_df$performance_method=gsub("_performance",
                                        "",
                                        tempSplit[[2]])
        curr_df$TP
        if (i==1){
            performance_plot_df = curr_df
        } else{
            performance_plot_df = rbind(performance_plot_df,curr_df)
        }
        
    }
    
    
    # performance_melt_df = reshape2::melt(performance_df,
    #                                      id.vars = grep("performance.cluster",colnames(performance_df)))
    # performance_melt_df2 = reshape2::melt(performance_melt_df,
    #                                       id.vars = c("variable","value"))
    # colnames(performance_melt_df2) = c("name","TP","name2","Clusters")
    # performance_plot_df = performance_melt_df2[,c("name", "TP","Clusters")]
    # tempSplit = lapply(as.character(performance_plot_df$name),
    #                    function (x) strsplit(x,split="[.]")[[1]])
    # performance_plot_df$marker_method = unlist(lapply(tempSplit,
    #                                                   function (x) x[[1]]))
    # 
    # performance_plot_df$performance_method = unlist(lapply(tempSplit,
    #                                                        function (x) x[[2]]))
    # performance_plot_df$performance_method = gsub("_performance","",performance_plot_df$performance_method)
    performance_plot_df$TPround = round(performance_plot_df$TP,digits=3)
    
    
    unique_performance = unique(performance_plot_df$performance_method)
    
    performance_plot_df$marker_method = factor(performance_plot_df$marker_method,
                                               levels = names(list_performance))
    
    list_p1 = list()
    for (i in 1:length(unique_performance)){
        
        curr_performance = unique_performance[[i]]
        curr_plot = performance_plot_df[which(performance_plot_df$performance_method %in% curr_performance),]
        
        p1=ggplot(curr_plot, aes(x = marker_method, y = Clusters,fill=TP)) +
            geom_tile(color="black") + theme_bw()+
            scale_fill_gradient(low = "white", high = "red") +
            theme(axis.text = element_text(size=text_size),
                  axis.title = element_text(size=text_size),
                  axis.text.x = element_text(angle=45,hjust=1,text_size),
                  plot.title = element_text(size=text_size),
                  legend.position="none")+
            xlab("Method") + ylab("Clusters") +
            coord_fixed(ratio=0.2)+
            geom_text(aes(marker_method, Clusters, label=TPround), colour = "black", check_overlap = FALSE,
                      size=tile_text_size)  +
            ggtitle(curr_performance)
        list_p1[[i]] = p1
        print(p1)
    }
    
    return(list_p1)
    
}
