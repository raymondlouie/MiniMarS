# ClusterMarkers
 
ClusterMarkers finds the markers which best define a cluster, using a number of different pre-exisitng methods.

## Installation

### Packages required to be pre-installed

Please install the following packages first. 

```{r}
# CiteFuse
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("CiteFuse")
devtools::install_github("tpq/propr") # propr package required for CiteFuse to run

# sc2marker
if (!require("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("https://github.com/CostaLab/sc2marker", build_vignettes = TRUE)

# geneBasisR
devtools::install_github("MarioniLab/geneBasisR") 

# Seurat
install.packages('Seurat')

# xgboost
install.packages("xgboost", repos="http://dmlc.ml/drat/", type = "source")

# SingleCellExperiment
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")

# dplyr
install.packages("dplyr")
```

### Installation of ClusterMarkers

Please run the following to install the `ClusterMarkers` package:
```
devtools::install_github("raymondlouie/ClusterMarkers") 
```
or download the package [here](https://www.dropbox.com/s/5dz53xqcp5u4sf4/ClusterMarkers_0.1.0.tar.gz?dl=0) and install it using
```
install.packages("~/Downloads/ClusterMarkers_0.1.0.tar.gz", type = "source", repo = NULL)
```

## Example work flow

Here is an example of the `ClusterMarkers` work flow to get started:

Load libraries and example data.
```{r}
# Check to see if packages are installed.
packages_required = c("CiteFuse","sc2marker","geneBasisR","xgboost","dplyr","ClusterMarkers")
packages_required_not_installed=setdiff(packages_required, rownames(installed.packages()))
if (length(packages_required_not_installed)>0){
    stop(paste0("Please install packages ",packages_required_not_installed,"\n"))
}

library(ClusterMarkers)
library(dplyr)
library(SingleCellExperiment)
data(sce)
input_matrix = t(sce@assays@data$counts)
clusters = sce$cell_type
```

The input data can  either be a i) feature matrix (with cluster vectors), ii) Seurat object or SCE object. 

We will first convert the input to the desired format required for downstream analysis, showing all three input data examples:
```{r}
# SCE input example. 
sce_in = processInputFormat(sc_object=sce,
                            sce_cluster="cell_type",
                            verbose=TRUE)

# Feature matrix with cluster vector example.                            
manual_in = processInputFormat(sc_object=input_matrix,
                            clusters_all=clusters,
                            verbose=TRUE)                           

# Seurat input example.
library(Seurat)
sc_object = CreateSeuratObject(input_matrix)
Idents(object = sc_object) <- clusters
seurat_in = processInputFormat(sc_object=sc_object,
                            verbose=TRUE)
```

We now select a subset of clusters (`clusters_sel`) to identify markers for. Default is using all clusters.
```{r}
clusters_sel = c("CD4-positive, alpha-beta memory T cell",
                 "naive thymus-derived CD8-positive, alpha-beta T cell")
sc_in = sce_in # As an example, select the SCE input

cluster_selection_out= processClusterSelection(sc_in,
                                               clusters_sel=clusters_sel,
                                               verbose=TRUE)
```   

In the next step, we i) sub-sample  the data, and ii) divide the data into a training and test set.
```{r}
final_out = processSubsampling(cluster_selection_out,
                               clusters_sel="all_clusters",
                               subsample_num=1000,
                               train_test_ratio = 0.9,
                               cluster_proportion= "proportional",
                               verbose=TRUE)
```

We now find the markers to identify the clusters. There are four methods implemented to identify the clusters using the `method` argument:  "citeFUSE", "sc2marker", "geneBasis" and "xgBoost". The default option is to use "all" methods. 
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

```
 [1] sp_1.4-6                    SeuratObject_4.1.0          Seurat_4.1.1.9002           SingleCellExperiment_1.16.0
 [5] SummarizedExperiment_1.24.0 Biobase_2.54.0              GenomicRanges_1.46.1        GenomeInfoDb_1.30.0        
 [9] IRanges_2.28.0              S4Vectors_0.32.3            BiocGenerics_0.40.0         MatrixGenerics_1.6.0       
[13] matrixStats_0.61.0          dplyr_1.0.7                 ClusterMarkers_0.1.0 
```

<br>

## Public datasets
We provide below public datasets for users to try the `ClusterMarkers` package. The datasets all contain protein features as well as either cell annotation or cluster labels from the corresponding papers or websites.

### Dataset 1
Blood and bone marrow samples of human from the BD Rhapsody assay.<br>
Number of cells in this dataset: around 100,000 (after QC: 49,100 from healthy controls and 31,600 from leukemia patients)<br>
Number of protein features: 97<br>
43 cell types from the data provider<br>
[Click here to access the data](https://cellxgene.cziscience.com/collections/93eebe82-d8c3-41bc-a906-63b5b5f24a9d)<br>

### Dataset 2
Mucosa-associated lymphoid tissue samples of human from the 10X technology.<br>
Number of cells in this dataset: around 10,000 (after QC: 8,500)<br>
Number of protein features: 17<br>
11 cell types from the data provider<br>
[Click here to access the data](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/malt_10k_protein_v3)<br>

### Dataset 3
Spleen and lymph nodes samples of mouse from the 10X technology.<br>
Number of cells in this dataset: around 40,000 (after QC: 16,800 with 111 protein features, 15,800 with 206 protein features)<br>
Number of protein features: 111 and 206 (two sets of data) - need to find the name of each protein feature (16 Mar 2023)<br>
Both sets, 35 cell types from the data provider<br>
[Click here to access the raw data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150599)<br>
[Click here to access the processed data](https://github.com/YosefLab/totalVI_reproducibility/)<br>

### Dataset 4
PBMC (peripheral blood mononuclear cells) samples of human from the 10X technology.<br>
Number of cells in this dataset: around 10,000 (after QC: 7,800)<br>
Number of protein features: 17<br>
11 cell types from the data provider<br>
[Click here to access the data](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_protein_v3)<br>


