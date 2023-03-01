# ClusterMarkers
 
ClusterMarkers finds the markers which best define a cluster, using a number of different pre-exisitng methods.

## Installation

### Packages required to be pre-installed

Please install the following packages first .

```
CiteFuse
sc2marker
geneBasisR
Seurat
SingleCellExperiment
dplyr
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

Load libraries and example data.
```{r}
library(dplyr)
library(SingleCellExperiment)
data(sce)
# input_matrix = t(sce@assays@data$counts)
# clusters = sce$cell_type
```


The input data can be a feature matrix (with cluster vectors), Seurat object or SCE object. We will first convert the input to the desired format required for downstream analysis.
```{r}
sc_out = processInputFormat(sc_object=sce,
                            sce_cluster="cell_type",
                            verbose=TRUE)

```

The user can now select a subset of clusters to find markers for, via the `clusters_sel` input.
```{r}
cluster_selection_out= processClusterSelection(sc_out,
                                               clusters_sel="all_clusters",
                                               verbose=TRUE)
```                                               

The next step is the sub-sampling of the data, and dividing the data into a training and test set.
```{r}
final_out = processSubsampling(cluster_selection_out,
                               clusters_sel="all_clusters",
                               subsample_num=1000,
                               train_test_ratio = 0.9,
                               cluster_proportion= "proportional",
                               verbose=TRUE)
```

We now find the markers to distinguish the clusters
```{r}

list_markers = findClusterMarkers(final_out$training_matrix,
                                  final_out$training_clusters,
                                  num_markers=15,
                                  method="all",
                                  verbose=TRUE)
```

Finally, we will evaulate the performance of the markers using the test data.
```{r}
list_performance = performanceAllMarkers(list_markers,
                                         final_out=final_out,
                                         method="all",
                                         nrounds=1500,
                                         nthread=6,
                                         verbose=TRUE)
```

<br>

## Public datasets
We provide below public datasets for users to try the `ClusterMarkers` package. The datasets all contain protein features as well as either cell annotation or cluster labels from the corresponding papers or websites.

### Dataset 1
Blood and bone marrow samples of human from the BD Rhapsody assay.<br>
Number of cells in this dataset: around 100K<br>
[Click here to access the data](https://cellxgene.cziscience.com/collections/93eebe82-d8c3-41bc-a906-63b5b5f24a9d)<br>

### Dataset 2
Mucosa-associated lymphoid tissue samples of human from the 10X technology.<br>
Number of cells in this dataset: around 10K<br>
[Click here to access the data](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/malt_10k_protein_v3)<br>

### Dataset 3
Spleen and lymph nodes samples of mouse from the 10X technology.<br>
Number of cells in this dataset: around 40K<br>
[Click here to access the raw data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150599)<br>
[Click here to access the processed data](https://github.com/YosefLab/totalVI_reproducibility/)<br>

### Dataset 4
PBMC (peripheral blood mononuclear cells) samples of human from the 10X technology.<br>
Number of cells in this dataset: around 10K<br>
[Click here to access the data](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_protein_v3)<br>


