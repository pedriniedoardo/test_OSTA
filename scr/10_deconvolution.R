# AIM ---------------------------------------------------------------------
# run sample deconvolution on Visium and VisiumHD dataset

# libraries ---------------------------------------------------------------
# library(OSTA.data)
library(VisiumIO)
library(SpatialExperiment)
library(ggspavis)
library(spacexr)
library(tidyverse)
library(patchwork)
library(DropletUtils)
library(BiocParallel)
library(scran)
library(scater)
library(pheatmap)
library(ComplexHeatmap)

# -------------------------------------------------------------------------
# In this example of Visium breast cancer data (Janesick et al. 2023), we perform Visium cell type deconvolution without Chromium reference and compare the concordance with the provided Visium annotation by 10x Genomics.

# retrieve dataset from OSF repository
# id <- "Visium_HumanBreast_Janesick"
# pa <- OSTA.data_load(id)
# dir.create(td <- tempfile())
# unzip(pa, exdir=td)

# data have been retrieved from a different session because of the current issue of OSTA.data
# location of the data
loc_sample_visium <- "../../data/Visium_HumanBreast_Janesick/"

# read into 'SpatialExperiment'
vis <- TENxVisium(
  spacerangerOut=loc_sample_visium, 
  processing = "filtered", 
  format="h5", 
  images="lowres") |> 
  import()

vis

# retrieve spot annotations & add as metadata
df <- read.csv(file.path(loc_sample_visium, "annotation.csv"))
head(df)
cs <- match(colnames(vis), df$Barcode)
vis$anno <- factor(df$Annotation[cs])

# set gene symbols as feature names
head(rownames(vis))
head(rowData(vis)$Symbol)

rownames(vis) <- make.unique(rowData(vis)$Symbol)
vis

# -------------------------------------------------------------------------
# plot spatial coordinates (spots)
plotSpots(vis)
plotVisium(vis,spots = F)|plotVisium(vis,zoom = T)
plotVisium(vis,)


# Deconvolution comes after proper quality control, as detailed in Chapter 8, and is usually performed on unnormalized and untransformed (i.e., raw) counts. Here, we quickly check if any spot has a count of 0 and remove it.
sub <- list(mt=grep("^MT-", rownames(vis)))
vis <- addPerCellQCMetrics(vis, subsets=sub)
vis$log_sum <- log(vis$sum)
vis@colData

plotSpots(vis,annotate = "log_sum")+scale_color_gradientn(colors=pals::jet())|
  plotSpots(vis,annotate = "subsets_mt_percent")+scale_color_gradientn(colors=pals::jet())

# All spots have non-zero library size, so we can proceed without subsetting. We first visualize the spot-level cell type annotation provided by 10x Genomics.
plotSpots(vis, 
          annotate = "anno", point_size = 1, 
          pal = unname(pals::trubetskoy())) +
  theme(legend.key.size = unit(0, "lines"))

# -------------------------------------------------------------------------
# Now, we load the Chromium reference data for the Visium dataset. To streamline the demonstration, we consolidate some of the cell type annotations provided by 10x Genomics (i.e. Annotation) into more generalized categories (i.e. Annogrp).

# retrieve dataset from OSF repo
# id <- "Chromium_HumanBreast_Janesick"
# pa <- OSTA.data_load(id)
# dir.create(td <- tempfile())
# unzip(pa, exdir=td)
loc_sample_chromium <- "../../data/Chromium_HumanBreast_Janesick/"

# read into 'SingleCellExperiment'
fs <- list.files(loc_sample_chromium, full.names=TRUE)
fs
h5 <- grep("h5$", fs, value=TRUE)
h5

sce <- read10xCounts(h5, col.names=TRUE)
sce

# use gene symbols as feature names
head(rownames(sce))
head(rowData(sce)$Symbol)

rownames(sce) <- make.unique(rowData(sce)$Symbol)

# retrieve cell type labels
csv <- grep("csv$", fs, value=TRUE)
csv
cd <- read.csv(csv, row.names=1)
head(cd)

# simplify annotations
pat <- c("B Cell"="B",
         "T Cell"="T",
         "Mac"="macro",
         "Mast"="mast",
         "DCs"="dendritic",
         "Peri"="perivas",
         "End"="endo",
         "Str"="stromal",
         "Inv"="tumor",
         "Myo"="myoepi",
         "Hyb" = "hybrid")

# -------------------------------------------------------------------------
# my implementation
df_meta <- cd %>%
  rownames_to_column("barcode")

df_pat <- data.frame(Annogrp = pat) %>%
  rownames_to_column("pattern")

df_meta_full <- pmap(list(df_pat$pattern,df_pat$Annogrp), function(pat,ann){
  df_meta %>%
    filter(str_detect(Annotation,pattern = pat)) %>%
    mutate(Annotation2 = ann)
}) %>%
  bind_rows() %>%
  # collapse the hybrid annotations
  group_by(barcode) %>%  
  summarize(size = n(), 
            Annogrp = paste(Annotation2,collapse="|")) %>%
  # column_to_rownames("barcode") %>%
  # if hybrid make it as NA
  mutate(Annogrp = case_when(str_detect(Annogrp, "hybrid")~"hybrid",
                             T~Annogrp)) %>%
  # join the with the full metadata
  left_join(df_meta,y = .,by = c("barcode")) %>%
  column_to_rownames("barcode") %>%
  mutate(Annogrp_final = case_when(is.na(Annogrp)~Annotation,
                                   Annogrp == "hybrid" ~ NA,
                                   T~Annogrp))

# count the annotations
table(df_meta_full$Annogrp_final)

# check meta before update
colData(sce)

sce$Annogrp <- df_meta_full[colnames(sce),"Annogrp_final"]
# check meta after update
colData(sce)
colData(sce)$Annogrp %>% is.na() %>% sum()

# -------------------------------------------------------------------------
# vignette implementation
# lab <- lab2 <- cd$Annotation
# for (. in names(pat))
#   lab2[grep(., lab)] <- pat[.]
# lab2[grepl("Hyb", lab)] <- NA
# lab2 <- gsub("\\s", "", lab2)
# 
# # add as cell metadata
# table(cd$Annogrp <- lab2)
# 
# colData(sce)
# colData(sce)[names(cd)] <- cd[colnames(sce), ]
# colData(sce)
# 
# data.frame(vignette = cd[colnames(sce), ],
#            custom = df_meta_full[colnames(sce),"Annogrp_final"])
# 
# colData(sce)$Annogrp %>% is.na() %>% sum()
# 
# 
# cd[colnames(sce), ] %>%
#   filter(str_detect(Annotation, "Hyb")) %>%
#   head()

# -------------------------------------------------------------------------

# We only keep the Chromium data with an annotation and are not labeled as “Hybrid”, as these correspond to mixed subpopulations.
sce2 <- sce[, !is.na(sce$Annogrp)]
dim(sce2)

# 10.3 RCTD ---------------------------------------------------------------
# Next, we perform deconvolution with spacexr (f.k.a. RCTD)(Cable et al. 2022). By default, run.RCTD()’s doublet_mode="doublet", specifies at most two subpopulations coexist in a data unit, such as a “spot”; here, we set doublet_mode="full" in order to allow for an arbitrary number of subpopulations to be fit instead.

# prep spatial data (Visium)
xy <- data.frame(spatialCoords(vis))
.vis <- SpatialRNA(xy, counts(vis))

# prep reference data (Chromium)
ids <- factor(sce2$Annogrp)
names(ids) <- colnames(sce2)
head(ids)

ref <- Reference(counts(sce2), ids)

# run 'RCDT' deconvolution
res <- create.RCTD(.vis, ref, max_cores=15)
res <- run.RCTD(res, doublet_mode="full")
str(res)

# Weights inferred by RCTD should be normalized such that proportions of cell types sum to 1 in each spot:
# normalize weights to obtain proportions
ws <- normalize_weights(res@results$weights)
# add proportion estimates as spot metadata
ws <- data.frame(as.matrix(ws))
colData(vis)[names(ws)] <- ws[colnames(vis), ]
colData(vis)

# Let’s visualize deconvolution weights in space, i.e., coloring by the proportion of a given cell type estimated to fall within a given spot:
# lapply(names(ws), \(.) 
#        plotSpots(vis, annotate=., point_size = 0.05)) |>
#   wrap_plots(nrow=3) & theme(
#     legend.key.width=unit(0.5, "lines"),
#     legend.key.height=unit(1, "lines")) &
#   scale_color_gradientn(colors=pals::jet())

lapply(names(ws),function(x){
  plotSpots(vis, annotate=x, point_size = 0.05)
}) %>%
  wrap_plots(nrow=3) &
  theme(
    legend.key.width=unit(0.5, "lines"),
    legend.key.height=unit(1, "lines")) &
  scale_color_gradientn(colors=pals::jet())

# The deconvolution result can also be viewed as a heatmap, where rows = cells and columns = clusters:
# pheatmap(ws, show_rownames=FALSE, show_colnames=TRUE,
#          cellwidth=12, treeheight_row=5, treeheight_col=5)

Heatmap(ws,show_row_names = F,
        col = viridis::viridis(option = "turbo",n = 10))

# We see that more than half of the spots are estimated to have a stromal proportion of more than 50%. Few spots have an intense and distinct signal for cancerous subpopulations, DCIS1 and DCIS2.

# For comparison with spot annotations provided by 10x Genomics, we include majority voted cell type from deconvolution as RCTD. Note that, because stromal cells show broad signals across the entire tissue, to better investigate immune cell signals, we remove stromal from the majority vote calculation for an alternative label: RCTD_no_stroma.
# derive majority vote label
ids <- names(ws)[apply(ws, 1, which.max)]
names(ids) <- rownames(ws)
head(ids)

vis$RCTD <- factor(ids[colnames(vis)])

# derive majority vote nouding stromal cells
ws_no_stroma <- ws[, colnames(ws) != "stromal"]
ids_no_stroma <- names(ws_no_stroma)[apply(ws_no_stroma, 1, which.max)]
names(ids_no_stroma) <- rownames(ws)
vis$RCTD_no_stroma <- factor(ids_no_stroma[colnames(vis)])

# We can visualize these three annotations spatially:
lapply(c("anno", "RCTD", "RCTD_no_stroma"), function(x){
  plotSpots(vis, annotate=x)
}) %>%
  wrap_plots(nrow=1) &
  theme(legend.key.size=unit(0, "lines")) &
  scale_color_manual(values=unname(pals::trubetskoy()))

# Note the strong stromal signals and macrophage being the second common cell type for stromal cells. To help characterize subpopulations from deconvolution, we can view their distribution against the provided annotation:
cd <- data.frame(colData(vis))
head(cd)

df <- as.data.frame(with(cd, table(RCTD, anno)))
head(df)
fd <- as.data.frame(with(cd, table(RCTD_no_stroma, anno)))
head(fd)

ggplot(df, aes(Freq, RCTD, fill=anno)) + ggtitle("RCTD") +
  ggplot(fd, aes(Freq, RCTD_no_stroma, fill=anno)) + ggtitle("RCTD_no_stroma") +
  plot_layout(nrow=1, guides="collect") &
  labs(x="Proportion", y=NULL) &
  coord_cartesian(expand=FALSE) &
  geom_col(width=1, col="white", position="fill") &
  scale_fill_manual(values=unname(pals::trubetskoy())) &
  theme_minimal() & theme(aspect.ratio=1,
                          legend.key.size=unit(2/3, "lines"),
                          plot.title=element_text(hjust=0.5))

# Next, we can investigate the agreement between provided annotation against the two deconvolution majority vote labels.
hm <- \(.) pheatmap(., show_rownames=TRUE, show_colnames=TRUE,
                    cellwidth=10, cellheight=10, treeheight_row=5, treeheight_col=5)

hm(prop.table(table(vis$anno, vis$RCTD), 2))
hm(prop.table(table(vis$anno, vis$RCTD_no_stroma), 2))

# Overall, we observe agreement between the provided spot-labels and RCTD deconvolution derived annotations. Before cleaning up stromal, some immune cell types, such as dendritic and mast, never had a chance to have the highest cell type proportion. On the left panel, among among all the spots annotated by RCTD as T cells, nearly all of them are in “immune” type in the provided annotation. Strong agreements are also observed for spots with cell type of “DCIS1”, “DCIS2”, and “Invasive tumor”.

# -------------------------------------------------------------------------
# PC regression
# Next, we prepare the principal components (PCs) needed to perform PC regression:

vis <- logNormCounts(vis)
# Feature selection 
dec <- modelGeneVar(vis)
head(dec)
hvg <- getTopHVGs(dec, prop = 0.1)
head(hvg)

# Dimension reduction 
set.seed(1234)
vis <- runPCA(vis, ncomponents = 20, subset_row = hvg)

# We fit the deconvolution result of each cell type against the first 10 PCs to obtain 10 regressions.
idx <- rownames(ws)
ids <- colnames(ws)
pcs <- reducedDim(vis, "PCA")
pcs <- pcs[idx, seq_len(10)]
pcr <- lapply(ids, \(id) {
  fit <- summary(lm(pcs ~ vis[, idx][[id]]))
  r2 <- sapply(fit, \(.) .$adj.r.squared)
  data.frame(id, pc=seq_along(r2), r2)
}) |> do.call(what=rbind)

# Here we plot the coefficient of determination of the first 10 PCs for each cell type.
pcr$id <- factor(pcr$id, ids)
pal <- pals::trubetskoy()
ggplot(pcr, aes(pc, r2, col=id)) +
  geom_line(show.legend=FALSE) + geom_point() +
  scale_color_manual("predictor", values=unname(pal)) +
  scale_x_continuous(breaks=c(1, seq(5, 20, 5))) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
  labs(x="principal component", y="coeff. of determination") +
  guides(col=guide_legend(override.aes=list(size=2))) +
  coord_cartesian(xlim=c(1, 10)) +
  theme_minimal() + theme(
    panel.grid.minor=element_blank(),
    legend.key.size=unit(0, "lines"))

# Let’s inspect the key drivers of (expression) variability in terms of PCs. Considering deconvolution results from above, we can see that
# PC1 distinguishes stromal, tumor, macrophage from the rest of the tissue;
# PC3 clearly separates DCIS1;
# PC4 captures T cells region;
# PC8 separates DCIS2.

pcs <- reducedDim(vis, "PCA")
colData(vis)[colnames(pcs)] <- pcs

lapply(c("DCIS1", "T", "DCIS2", colnames(pcs)[c(3,4,8)]),function(x){
  plotSpots(vis, annotate=x) +
    scale_color_gradientn(x, colors=pals::jet())
}) %>%
  wrap_plots(nrow=2) &
  theme(
    plot.title=element_blank(),
    legend.key.width=unit(0.5, "lines"),
    legend.key.height=unit(1, "lines"))

# Note that the direction of each PC is irrelevant from how much variation it explains.
# In conclusion, deconvolution-based cell type proportion estimates are able to recapitulate PCs and, in turn, expression variability.
