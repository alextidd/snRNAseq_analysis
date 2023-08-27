# snRNAseq_workflow

snRNAseq_workflow image built as a modification of 'almurphy-scfdev-dev.img'.

Build a sandbox of the image, enter the sandbox and open R:

```
singularity build --sandbox snRNAseq_workflow/ docker://almurphy/scfdev:dev
singularity shell -f --writable snRNAseq_workflow/
INFO:    User not listed in /etc/subuid, trying root-mapped namespace
INFO:    Using fakeroot command combined with root-mapped namespace
Apptainer> R
```

Then in R inside the container, run install.packages() or BiocManager::install() to download the following packages in this order:

* dplyr
* Seurat
* SeuratWrapper
* inferCNV
* SingleR
* scDblFinder
* celldex
* Matrix
* stringr
* reshape2
* RColorBrewer
* ggplot2
* pheatmap
* purrr
* rmarkdown
* xfun
* Monocle3
* rjson
* clustree
* dittoSeq
* ggforce
* ggrepel
* magrittr
* readr
* SummarizedExperiment
* janitor
* gridExtra
* intrinsicDimension
* knitr
* lemon
* patchwork
* rlang
* glmGamPoi (from github)

Install SCPA and its dependencies:

```
devtools::install_version("crossmatch", version = "1.3.1", repos = "http://cran.us.r-project.org")
devtools::install_version("multicross", version = "2.1.0", repos = "http://cran.us.r-project.org")
BiocManager::install("ComplexHeatmap")
BiocManager::install("GSEABase")
BiocManager::install("GSVA")
BiocManager::install("singscore")
install.packages("clustermole")
devtools::install_github("jackbibby1/SCPA")
install.packages("msigdbr")
```

Install a previous version of EWCE from GitHub:

```
devtools::install_github("neurogenomics/EWCE")
```

Then exit the container (`exit`) and convert it to an image:

```
singularity build snRNAseq_workflow.img snRNAseq_workflow/
```

# inferCNV

```
singularity build infercnv.latest.simg docker://trinityctat/infercnv:latest
singularity build --sandbox infercnv.latest/ infercnv.latest.simg
singularity shell -f --writable infercnv.latest/
Apptainer> R
> install.packages("rjson")
> install.packages("readr")
> q()
exit
singularity build infercnv.latest.img infercnv.latest/
```
