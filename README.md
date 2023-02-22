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
 					     clusters_sel=c("CD4-positive, alpha-beta memory T cell","CD8-positive, alpha-beta memory T cell"),
                                   num_markers=15,
                                   subsample_num=100,
					     sub.seed=42,
                                   train_test_ratio=0.9,
                                   method_cluster="all",
                                   method_performance="all",
                                   cluster_proportion = "proportional",
                                   verbose=FALSE 
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
[Click here to access the raw data](https://github.com/YosefLab/totalVI_reproducibility/)<br>
[Click here to access the processed data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150599)<br>

### Dataset 4
PBMC (peripheral blood mononuclear cells) samples of human from the 10X technology.<br>
Number of cells in this dataset: around 10K<br>
[Click here to access the data](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_protein_v3)<br>


