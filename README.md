# MiniMarS
 
MiniMarS finds the protein markers that best define biological clusters, using a number of different pre-existing methods for marker selection.

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
install.packages("xgboost")

# SingleCellExperiment
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")

# dplyr
install.packages("dplyr")
```

### Install `MiniMarS`

Please run the following to install the `MiniMarS` package (typical install time: 5-10 minutes):
```
devtools::install_github("https://github.com/raymondlouie/MiniMarS",force=TRUE)

```

or download the package [here](https://www.dropbox.com/scl/fi/2vk9wy1j1tw7j2vutzp2q/MiniMarS_0.3.1.tar.gz?rlkey=kt5dj2oqyozwn9u82v6d9mhhz&dl=0) and install it using the following command
```
install.packages("~/Downloads/MiniMarS_0.3.1.tar.gz", type = "source", repos = NULL)
```

## Standard `MiniMarS` workflow demo

Expected runtime: 1 minute

### Load the libraries
```{r}
# Check to see if all the packages are installed.
packages_required = c("CiteFuse","sc2marker","geneBasisR","xgboost","dplyr","MiniMarS")
packages_required_not_installed=setdiff(packages_required, rownames(installed.packages()))
if (length(packages_required_not_installed)>0){
    stop(paste0("Please install packages ",packages_required_not_installed,"\n"))
}

library(MiniMarS)
library(dplyr)
library(SingleCellExperiment)
```

### Convert the input data to the desired format. 
The input data can be i) an SCE object, ii) a matrix of features and a vector of cell type annotations, or iii) a Seurat object. 

```{r}
# Load the data
```

i) If you have a SCE object `sce`:
```{r}
sc_in = processInputFormat(sc_object = sce,
                           sce_cluster = "cell_type", #pre-defined in the SCE object
                           verbose = TRUE)
```
                               
ii) If you have a feature matrix (feature x cell) `input_matrix` and a vector of cell type annotations for each cell `cluster` (the length of this vector should be the same as the number of columns of the `input_matrix`):
```{r}
sc_in = processInputFormat(sc_object = input_matrix,
                               clusters_all = clusters,
                               verbose = TRUE)
```
                               
iii) If you have a Seurat object:
```{r}
library(Seurat)

sc_in = processInputFormat(sc_object = seurat_object,
                           verbose=TRUE)
```

### Select the clusters you want to identify the markers for. 
```{r}
# If 'clusters_sel' is not defined, then the default is to use all clusters.
clusters_sel = c("CD4-positive, alpha-beta memory T cell",
                 "naive thymus-derived CD8-positive, alpha-beta T cell")
cluster_selection_out= processClusterSelection(sc_in,
                                               clusters_sel = clusters_sel,
                                               verbose = TRUE)
```   

### Sub-sample and divide the dataset into training, validation, and testing sets.
```{r}
final_out = processSubsampling(cluster_selection_out,
                               subsample_num = 100,
                               verbose = TRUE,
                               seed = 8)
```

### Find the minimum number of markers to identify the clusters. 
Finds the minimum number of markers using our wrapper function, which selects the top method by default, out of the following methods: "citeFUSE", "sc2marker", "geneBasis", "xgBoost", "fstat", "seurat_wilcox", "seurat_bimod", "seurat_roc", "seurat_t", "seurat_LR", "consensus_weighted", "consensus_naive", "consensus_fstat", and "consensus_xgboost". 
```{r}
minMarker_result <- minMarker(final_out,
                              list_markers_test=c(5,10,15,20,25,30,40),
                              # numTimes = 50,
                              metric = "F1_macro",
                              cluster = "all",
                              threshold  = 0.7)
```


### Visualise the identified markers and their performance.
```{r}
library(ggplot2)

list_markers_all = minMarker_result[[length(minMarker_result)]]$markersAll
list_performance_all= minMarker_result[[length(minMarker_result)]]$performanceAll

plotMarkers(list_markers_all)
plotPerformance(list_performance_all)

## Visualisation
library(RColorBrewer)
plotExpression(list_markers_all,
                   sc_in,
                   plot_type="violin")

plotExpression(list_markers_all,
                   sc_in,
                   plot_type="umap")

```

```{r}
sessionInfo()
```

```
 [1] sp_1.4-6                    SeuratObject_4.1.0          Seurat_4.1.1.9002           SingleCellExperiment_1.16.0
 [5] SummarizedExperiment_1.24.0 Biobase_2.54.0              GenomicRanges_1.46.1        GenomeInfoDb_1.30.0        
 [9] IRanges_2.28.0              S4Vectors_0.32.3            BiocGenerics_0.40.0         MatrixGenerics_1.6.0       
[13] matrixStats_0.61.0          dplyr_1.0.7                 MiniMarS_0.1.0 
```

<br>

## Using the user's own marker panel

Under some scenarios, the users may want to test their own customised marker panel instead of the predicted ones. We recommend using the following codes to evaluate the performance of the user's customised marker input.
```{r}
user_markers = c("CD80", "CD86", "CD274", "CD273", "CD275", "CD11b", "CD137L", "CD70", "unidentified_marker1","unidentified_marker2")
own_list_performance = performanceOwnMarkers(user_markers,
                                         final_out = final_out,
                                         method = "all",
                                         nrounds = 1500,
                                         nthread = 6,
                                         verbose = TRUE)
print(own_list_performance)
```

<br>

## Public datasets
We provide below public datasets for users to try the `MiniMarS` package. The datasets all contain protein features as well as either cell annotation or cluster labels from the corresponding papers or websites.

### Dataset 1
Human bone marrow samples from the BD Rhapsody assay.<br>
Number of cells in this dataset: around 100,000 (after QC: 49,100 from healthy controls and 31,600 from leukemia patients)<br>
Number of protein features: 97<br>
43 cell types/clusters from the data provider<br>
[Click here to access the data](https://cellxgene.cziscience.com/collections/93eebe82-d8c3-41bc-a906-63b5b5f24a9d)<br>

### Dataset 2
Human mucosa-associated lymphoid tissue samples from the 10X technology.<br>
Number of cells in this dataset: around 10,000 (after QC: 8,500)<br>
Number of protein features: 17<br>
11 cell types/clusters from the data provider<br>
[Click here to access the data](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/malt_10k_protein_v3)<br>

### Dataset 3
Mouse spleen and lymph nodes samples from the 10X technology.<br>
Number of cells in this dataset: around 40,000 (after QC: 16,800 with 111 protein features, 15,800 with 206 protein features)<br>
Number of protein features: 111 and 206 (two sets of data) - need to find the name of each protein feature (16 Mar 2023)<br>
Both sets, 35 cell types/clusters from the data provider<br>
[Click here to access the raw data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150599)<br>
[Click here to access the processed data](https://github.com/YosefLab/totalVI_reproducibility/)<br>

### Dataset 4
Human PBMC (peripheral blood mononuclear cells) samples from the 10X technology.<br>
Number of cells in this dataset: around 10,000 (after QC: 7,800)<br>
Number of protein features: 17<br>
11 cell types/clusters from the data provider<br>
[Click here to access the data](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_protein_v3)<br>

## Citation and Runtime

Please see our [paper](https://www.biorxiv.org/cgi/content/short/2025.07.14.664825v1) for more details, including runtime. 
