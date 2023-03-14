---
title: "Script to test methods"
author: "Raymond Louie"
date: "`r format(Sys.time(), '%d %B, %Y')`"
output:
  html_document:
    number_sections: yes
    self_contained: yes
    toc: yes
    toc_float: yes
  pdf_document:
    toc: yes   
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(cache=FALSE)
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(results = 'asis')
knitr::opts_chunk$set(fig.width = 25)
knitr::opts_chunk$set(fig.height = 10)
knitr::opts_chunk$set(fig.pos = '!h')
knitr::opts_knit$set(root.dir = getwd())
knitr::knit_hooks$set(timeit = local({
    now = NULL
    function(before, options) {
        if (before) {
            now <<- Sys.time()
        } else {
            res = difftime(Sys.time(), now)
            now <<- NULL
            # use options$label if you want the chunk label as well
            paste('Time for this code chunk:', as.character(res))
        }
    }})
)
```
```{r clear, echo=FALSE}
rm(list = ls())
graphics.off()
```
<!-- First we load the libraries: -->

```{r,warning=FALSE,message=FALSE,echo=FALSE}
# library(Matrix)
# library(dplyr)
# library(Seurat)
# library(ggplot2)
# library(grid)
# require(gridExtra)
# library(data.table)
# library(knitr)
# library(kableExtra)
# library(CiteFuse)
# library(scater)
# library(SingleCellExperiment)
# library(geneBasisR)
# library(xgboost)
# library(dplyr)  
# library(caret)
# library(Matrix)
# library(ClusterMarkers)
# source("performanceMarkers.R")
# 
# source("findClusterMarkers.R")
# source("wrapperMethods.R")
# source("findClusterPerformance.R")
```
<!-- install.packages("ClusterMarkers",repos=NULL,type="source") -->


Load libraries and example data.
```{r}
packages_required = c("CiteFuse","sc2marker","geneBasisR","xgboost","dplyr","ClusterMarkers")
packages_required_not_installed=setdiff(packages_required, rownames(installed.packages()))
if (length(packages_required_not_installed)>0){
    stop(paste0("Please install packages",packages_required_not_installed))
}

library(ClusterMarkers)
library(dplyr)

# Check to see if packages are installed



library(SingleCellExperiment)
data(sce)
input_matrix = t(sce@assays@data$counts)
clusters = sce$cell_type
```

The input data can  either be a i) feature matrix (with cluster vectors), ii) Seurat object or SCE object. 

We will first convert the input to the desired format required for downstream analysis, showing all three input data examples:
```{r}
sce_in = processInputFormat(sc_object=sce,
                            sce_cluster="cell_type",
                            verbose=TRUE)

manual_in = processInputFormat(sc_object=input_matrix,
                               clusters_all=clusters,
                               verbose=TRUE)                           

library(Seurat)
sc_object = CreateSeuratObject(input_matrix)
Idents(object = sc_object) <- clusters
seurat_in = processInputFormat(sc_object=sc_object,
                               verbose=TRUE)
unique(clusters)
```

We can now select a subset of clusters to identify markers for, via the `clusters_sel` input.
```{r}
clusters_sel = c("CD4-positive, alpha-beta memory T cell",
                 "naive thymus-derived CD8-positive, alpha-beta T cell")
sc_in = sce_in # as an example, select the SCE input
cluster_selection_out= processClusterSelection(sc_in,
                                               clusters_sel=clusters_sel,
                                               verbose=TRUE)
```                                               

In the next step, we i) Sub-sample  the data, and ii) Divide the data into a training and test set.
```{r}
final_out = processSubsampling(cluster_selection_out,
                               clusters_sel="all_clusters",
                               subsample_num=1000,
                               train_test_ratio = 0.9,
                               cluster_proportion= "proportional",
                               verbose=TRUE)
```

We now find the markers to identify the cluster. There are four methods implemented to identify the clusters using the `method` argument:  "citeFUSE", "sc2marker", "geneBasis" and "xgBoost". The default option is to use "all" methods. 
```{r}
list_markers = findClusterMarkers(final_out$training_matrix,
                                  final_out$training_clusters,
                                  num_markers=15,
                                  method="all",
                                  verbose=TRUE)
```

Finally, we  evaulate the performance of the markers using the test data. There are two methods implemented to test the performance using the `method` argument:  "xgBoost" and "geneBasis". The default option is to use "all" methods. 
```{r}
list_performance = performanceAllMarkers(list_markers,
                                         final_out=final_out,
                                         method="all",
                                         nrounds=1500,
                                         nthread=6,
                                         verbose=TRUE)
```

```{r}
sessionInfo()
```
