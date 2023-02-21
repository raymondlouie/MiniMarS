# ClusterMarkers
 
ClusterMarkers finds the markers which best define a cluster, using a number of different pre-exisitng methods.

## Installation

### Packages required to be pre-installed

Please install the following packages first .

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
or download the package [here](https://www.dropbox.com/s/5dz53xqcp5u4sf4/ClusterMarkers_0.1.0.tar.gz?dl=0) and install it using
```
install.packages("~/Downloads/ClusterMarkers_0.1.0.tar.gz", type = "source", repo = NULL)
```

## Example work flow
An example of the `ClusterMarkers` work flow to get started:

```{r}
library(SingleCellExperiment)
data(sce)
input_matrix = sce@assays@data$counts
clusters = sce$cell_type

output_df = findClusterPerformance(input_matrix,
                                   clusters,
                                   num_markers=15,
                                   subsample_num=100,
					     sub.seed=42,
                                   train_test_ratio=0.9,
                                   method_cluster="all",
                                   method_performance="all",
                                   cluster_proportion = "proportional",
                                   verbose=FALSE 
```


