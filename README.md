# ClusterMarkers
 
ClusterMarkers finds the markers which best define a cluster, using a number of different pre-exisitng methods.

### Installation

```

## Install development version
devtools::install_github("raymondlouie/ClusterMarkers", ref="devel") 
```

### Example work flow
An example of the `Milo` work flow to get started:

```{r}
data(sce)
input_matrix = sce@assays@data$counts
clusters = sce$cell_type

output_results = findClusterMarkers(input_matrix,
                                    clusters,
                                    num_markers = 15,
                                    method="all")
```


