rm(list=ls())
library(ggspavis)
library(SpatialExperiment)
library(scater)
`%notin%` <- Negate(`%in%`)
setwd("./Data/mouse_cerebellum/")
load("SlideseqV2_ROI.rds")
sce = SingleCellExperiment(assays=list(counts=sp_count), 
                           colData=location)

rd <- S4Vectors::DataFrame(
  symbol = rowData(sce)$Symbol)

(spe <- SpatialExperiment(
  assays = list(counts = assay(sce)),
  colData = colData(sce)
))

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rownames(spe))
# calculate per-spot QC metrics and store in colData
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
colData(spe)$sum <- colSums(counts(spe))
spatialCoords(spe) <- as.matrix(colData(sce)[,1:2])
spatialCoordsNames(spe) <- c("x_coord","y_coord")
spatialData(spe) <- cbind(spatialData(spe), 
                          in_tissue = rep(1,dim(spe)[2]))
colnames(colData(spe))[1:2] <- c("row","col")
spe1 = spe
      
## Library size
# select QC threshold for library size
qc_lib_size <- colData(spe1)$sum < 100
colData(spe1)$qc_lib_size <- qc_lib_size

## Number of expressed features
# select QC threshold for expressed genes
qc_detected <- colData(spe1)$detected <= 60
colData(spe1)$qc_detected <- qc_detected

## Proportion of mitochondrial read
# select QC threshold for mitochondrial read proportion
qc_mito <- colData(spe1)$subsets_mito_percent > 50
colData(spe1)$qc_mito <- qc_mito

## Remove low-quality spots
# combined set of discarded spots
discard <- qc_lib_size | qc_detected | qc_mito

# store in object
# remove combined set of low-quality spots
expr = assays(spe1)$counts
# extract metadata from gobject
metadata = colData(spe1)
sce = SingleCellExperiment(assays=list(counts=expr), 
                           colData=metadata)
libsizes <- colSums(expr)
size.factors <- libsizes/mean(libsizes)
logcounts(sce) <- log2(t(t(expr)/size.factors) + 1)
all(colnames(sce) == names(discard))
colData(sce)$discard <- discard
sce <- sce[, !colData(sce)$discard]

## Gene level
### basic filter: remove undetected genes
# calculate logcounts (log-transformed normalized counts) and store in object
spe1 = sce
spe1 <- logNormCounts(spe1)
qc_undetected_gene <- rowSums(assays(spe1)$logcounts > 0) > 0

### basic filter: remove lowly expressed genes: at least xx non-zero cells:
qc_low_gene <- rowSums(assays(spe1)$counts > 0) >= 100
## remove combined set of low-quality spots

spe1 <- spe1[qc_low_gene & qc_undetected_gene,]
dim(spe1)

## store the filtered dataset
sce_one = spe1
path = './Data/mouse_cerebellum/'
save(sce_one, file = paste0(path,"cerebellum_filtered.rda"))
