# -------------------------------------------------------------------------
# 3.2.1 Calling Python
# from the command line
# micromamba create -n env_scvi -c conda-forge conda python=3.12 scanpy mkl-include r-base r-essentials r-reticulate
# micromamba activate scvi
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
library(SpatialExperiment)
library(ggspavis)

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

# try also the inverse conversion
cortex.adata <- sceasy::convertFormat(cortex.seurat, from="seurat", to="anndata",outFile='../../out/object/cortex_sceasy.h5ad')
# confirm the creation
dir("../../out/object/")

# convert the seurat object to SingleCellExperiment
cortex.sce <- as.SingleCellExperiment(cortex.seurat)

# convert the SingleCellExperiment to AnnData
cortex.ad <- SCE2AnnData(cortex.sce)
unique(cortex.ad$obs$cell_type)

# -------------------------------------------------------------------------
# sample function to conver a Seurat spatial object into a spatial experiment object

seurat_to_spe <- function(seu, sample_id, img_id) {
  ## Convert to SCE
  sce <- Seurat::as.SingleCellExperiment(seu)
  
  ## Extract spatial coordinates
  spatialCoords <- GetTissueCoordinates(seu)[, 1:2] %>%
    as.matrix(.)
  # spatialCoords <- as.matrix(
  #   seu@images[[img_id]]@coordinates[, c("imagecol", "imagerow")])
  
  ## Extract and process image data
  img <- SpatialExperiment::SpatialImage(
    x = as.raster(seu@images[[img_id]]@image))
  
  imgData <- DataFrame(
    sample_id = sample_id,
    image_id = img_id,
    data = I(list(img)),
    scaleFactor = seu@images[[img_id]]@scale.factors$lowres)
  
  # Convert to SpatialExperiment
  spe <- SpatialExperiment(
    assays = assays(sce),
    rowData = rowData(sce),
    colData = colData(sce),
    metadata = metadata(sce),
    reducedDims = reducedDims(sce),
    altExps = altExps(sce),
    sample_id = sample_id,
    spatialCoords = spatialCoords,
    imgData = imgData
  )
  # indicate all spots are on the tissue
  spe$in_tissue <- 1
  spe$sample_id <- sample_id
  # Return Spatial Experiment object
  spe
}

# test
# SeuratData::InstallData(ds = "stxBrain")
se1 <- SeuratData::LoadData(ds = "stxBrain", type = "posterior1")
# sample plot
SpatialFeaturePlot(se1,features = "nCount_Spatial")

# conver the object
spe <- seurat_to_spe(seu = se1, sample_id = "posterior1", img_id = "posterior1")

# sample plot
plotSpots(spe,point_size = 0.05)
plotSpots(spe,point_size = 0.05,x_coord = "y",y_coord = "x")

plotSpots(spe,point_size = 0.05,annotate = "nCount_Spatial",x_coord = "y",y_coord = "x") + 
  scale_color_gradientn(colors=pals::jet())

