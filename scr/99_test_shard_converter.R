# libraries ---------------------------------------------------------------
library(schard)
library(tidyverse)
library(Seurat)
library(loupeR)

# download the sample files -----------------------------------------------
# download.file('https://datasets.cellxgene.cziscience.com/c5ac5c36-f60c-4680-8018-2d6cb65c0a37.h5ad','../../data/vis.heart.h5ad')
# download.file('https://datasets.cellxgene.cziscience.com/8cc521c8-c4ff-4cba-a07b-cae67a9dcba9.h5ad','../../data/sn.heart.h5ad')
# download.file('https://covid19.cog.sanger.ac.uk/baron16.processed.h5ad','../../data/ba16.h5ad')

# testing vignette --------------------------------------------------------
# load h5ad as Single Cell Experiment
ba16.sce <- schard::h5ad2sce('../../data/ba16.h5ad')

# load h5ad as Seurat
snhx <- schard::h5ad2seurat('../../data/sn.heart.h5ad')

# load all visium samples as single Seurat object
visx01 <- schard::h5ad2seurat_spatial('../../data/vis.heart.h5ad')
visx02 <- schard::h5ad2seurat_spatial('../../data/vis.heart.h5ad',use.raw = T)

# see the difference in loading with or without raw parameter
visx01[["Spatial"]]$counts[1:10,1:10]
visx02[["Spatial"]]$counts[1:10,1:10]

# or load as list of Seurat objects (per slide)
visl <- schard::h5ad2seurat_spatial('../../data/vis.heart.h5ad',simplify = FALSE)

# or load raw counts
snhr <- schard::h5ad2seurat('../../data/sn.heart.h5ad',use.raw = TRUE)

# raw counts for visium
visr <- schard::h5ad2seurat_spatial('../../data/vis.heart.h5ad',use.raw = TRUE)

# check that it works
Seurat::SpatialPlot(visx01,features = 'total_counts')
Seurat::SpatialPlot(visx01,features = 'total_counts',images = 'HCAHeartST11702009')
Seurat::SpatialPlot(visl$HCAHeartST11702010,features = 'total_counts')

# raw counts are different from normolized ones
plot(colSums(visx@assays$Spatial),colSums(visr@assays$Spatial),pch=16)

# the name of reduction is 'Xumap_' (autotranslated from scanpy to Seurat), somehow DimPlot manages to find it, but probably safier to specify it manually with reduction = 'Xumap_'
Seurat::DimPlot(snhx,group.by = 'cell_state')

# testing misc ------------------------------------------------------------
test01 <- schard::h5ad2seurat_spatial('../../data/stardist/Sample_HDP031_C1_nuclei_grouped.h5ad')
Seurat::SpatialPlot(test01,features = 'total_counts')

test02 <- schard::h5ad2seurat_spatial('../../data/stardist/Sample_HDP031_C1_expanded_nuclei.h5ad')
Seurat::SpatialPlot(test02,features = 'total_counts')

GetTissueCoordinates(test01) %>% head()
GetTissueCoordinates(test02) %>% head()

# try to save the cloupe from the object
# import the library
# library("loupeR")

# convert the SeuratObject named `seurat_obj` to a Loupe file
create_loupe_from_seurat(test01,output_dir = "../../out/object/",output_name = "Sample_HDP031_C1_nuclei_grouped")
