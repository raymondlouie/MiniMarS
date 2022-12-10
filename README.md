# ClusterMarkers
 
ClusterMarkers finds the markers which best define a cluster, using a number of different pre-exisitng methods.

## Installation

### Packages required to be pre-installed

Please install the following packages first.

```
citeFuse
sc2Marker
geneBasis
Seurat
SingleCellExperiment
```

### Installation

```
devtools::install_github("raymondlouie/ClusterMarkers") 
```

## Example work flow
An example of the `ClusterMarkers` work flow to get started:

```{r}
library(SingleCellExperiment)
data(sce)
input_matrix = sce@assays@data$counts
clusters = sce$cell_type

output_results = findClusterMarkers(input_matrix,
                                    clusters,
                                    num_markers = 15,
                                    method="all")
```


