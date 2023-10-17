# MiniMarS
 
MiniMarS finds the markers that best define a cluster, using a number of different pre-existing methods.

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

### Installation of MiniMarS

Please run the following to install the `MiniMarS` package from the Development branch:
```
devtools::install_github("raymondlouie/MiniMarS", ref = "Dev")
```

or download the package [here](https://www.dropbox.com/scl/fi/vqpz8h45t2tf5uxrzrrwo/MiniMarS_0.2.0.tar.gz?rlkey=u0y0o03j0xqnz1usp8gjmpz3m&dl=0) and install it using the following command
```
install.packages("~/Downloads/MiniMarS_0.2.0.tar.gz", type = "source", repo = NULL)
```

## Example workflow

Here is an example of the `MiniMarS` workflow to get started:

Load libraries and example data.
```{r}
# Check to see if packages are installed.
packages_required = c("CiteFuse","sc2marker","geneBasisR","xgboost","dplyr","MiniMarS")
packages_required_not_installed=setdiff(packages_required, rownames(installed.packages()))
if (length(packages_required_not_installed)>0){
    stop(paste0("Please install packages ",packages_required_not_installed,"\n"))
}

library(MiniMarS)
library(dplyr)
library(SingleCellExperiment)
```

The input data can  either be i) a feature matrix (with cluster vectors), ii) a Seurat object or an SCE object. 

Here we use the SCE object included in the package as an example
```{r}
# Load the dataset
data(sce)
```

First, we convert the input to the desired format required for downstream analysis, showing all three input data examples:
```{r}
# SCE input example. 
sc_in = processInputFormat(sc_object = sce,
                           sce_cluster = "cell_type",
                           verbose = TRUE)
                               
# Feature matrix with cluster vector example.
# The 'input_matrix' should be formatted as feature x cell matrix
input_matrix = sce@assays@data$counts
# The 'clusters' should be a vector of cell cluster annotations corresponding to each cell (i.e., the row of the input_matrix)
clusters = sce$cell_type
sc_in = processInputFormat(sc_object = input_matrix,
                               clusters_all = clusters,
                               verbose = TRUE)

sc_in_all=sc_in                               
# Seurat input example.
library(Seurat)
# Create a Seurat object or read in the user's own object
sc_object = CreateSeuratObject(input_matrix, assay = "Protein")
Idents(object = sc_object) <- clusters
sc_in = processInputFormat(sc_object = sc_object,
                           verbose=TRUE)
```

Second, we select a subset of clusters (`clusters_sel`) to identify markers. The default is to use all clusters.
```{r}
clusters_sel = c("CD4-positive, alpha-beta memory T cell",
                 "naive thymus-derived CD8-positive, alpha-beta T cell")

cluster_selection_out= processClusterSelection(sc_in,
                                               clusters_sel = clusters_sel,
                                               verbose = TRUE)
```   

Third, we i) sub-sample  the data, and ii) divide the data into a training and test set.
```{r}
final_out = processSubsampling(cluster_selection_out,
                               subsample_num = 1000,
                               train_test_ratio = 0.9,
                               cluster_proportion = "proportional",
                               verbose = TRUE,
                               seed = 8)
```

Fourth, we now find the markers to identify the clusters. There are four methods implemented to identify the clusters using the `method` argument:  "citeFUSE", "sc2marker", "geneBasis" and "xgBoost". The default option is to use "all" methods. 
```{r}
list_markers_time = findClusterMarkers(final_out$training_matrix,
                                  final_out$training_clusters,
                                  num_markers = 15,
                                  method = "all",
                                  verbose = TRUE)

list_time = list_markers_time$runtime_secs
names(list_time) = names(list_markers_time)[which(!(names(list_markers_time) %in% c("consensus",
                                                                                    "runtime_secs")))]
list_markers = list_markers_time[which(!(names(list_markers_time) %in% c("runtime_secs")))]
```

Finally, we  evaluate the performance of the markers using the test data. There are two methods implemented to test the performance using the `method` argument:  "xgBoost" and "geneBasis". The default option is to use "all" methods. 
```{r}
list_performance = performanceAllMarkers(list_markers,
                                         final_out = final_out,
                                         method = "all",
                                         nrounds = 1500,
                                         nthread = 6,
                                         verbose = TRUE)
```

We can print out the identified markers and their performance:
```{r}
library(ggplot2)

plotMarkers(list_markers)
plotPerformance(list_performance)

library(RColorBrewer)
plotExpression(list_markers,
                   sc_in_all,
                   plot_type="violin")

plotExpression(list_markers,
                   sc_in_all,
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


