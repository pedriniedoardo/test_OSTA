# -------------------------------------------------------------------------
# 3.2.1 Calling Python
# from the command line
# conda create -n scvi python=3.12
# conda activate scvi
# python -m pip install scvi-tools

# python version ----------------------------------------------------------
Sys.setenv(RETICULATE_PYTHON = "/home/edo/micromamba/envs/env_scvi/bin/python")
library(reticulate)
reticulate::use_condaenv(condaenv = "/home/edo/micromamba/envs/env_scvi")

py_config()

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(sceasy)
library(reticulate)
library(zellkonverter)

# sample vigette ----------------------------------------------------------

scvi <- import("scvi")
# dir.create(td <- tempfile())
# this fiils in reticulate. rownload the file from the call in jupyther
# (ad <- scvi$data$cortex(save_path=td))

h5ad_path <- "../../data/cortex.h5ad"
ad <- scvi$data$read_h5ad(h5ad_path)

# 3.2.2 Continuing in R
# We can access any of the variables above in R. For basic outputs, this works out of the box:
unique(ad$obs$cell_type)

sce <- AnnData2SCE(ad)

# 3.2.3 Back to Python
# We can also do the reverse, i.e. go from R’s SingleCellExperiment to Python’s AnnData:
ad <- SCE2AnnData(sce, X_name="X")

# convert the file --------------------------------------------------------
# sample conversion from anndata to as seurat object
cortex.seurat <- sceasy::convertFormat(h5ad_path, from="anndata", to="seurat",outFile='../../out/object/cortex.rds')

# convert the seurat object to SingleCellExperiment
cortex.sce <- as.SingleCellExperiment(cortex.seurat)

# convert the SingleCellExperiment to AnnData
cortex.ad <- SCE2AnnData(cortex.sce)
unique(cortex.ad$obs$cell_type)
