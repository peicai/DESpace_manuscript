## Spit sce into count.csv and meta.csv 
## As input for methods in python

############################ mouse cerebellum ####################
load("./DESpace_data/Data/mouse_cerebellum/cerebellum_filtered.rda")
library(SingleCellExperiment)
sce = sce_one
sce <- scuttle::logNormCounts(sce)
save_dir = "./DESpace_data/Data/mouse_cerebellum/"
object = sce
expr = assays(object)$counts
metadata = colData(object)
locs = as.data.frame(cbind(metadata$row,metadata$col))
rownames(locs) <- colnames(expr)
colnames(locs) <- c('row','col')
write.csv(as.matrix(expr), quote = FALSE, row.names = TRUE,
          file = paste0(save_dir, '100_count.csv'))
write.csv(metadata, quote = FALSE, row.names = TRUE,
          file = paste0(save_dir, '100_meta.csv'))
write.csv(locs, quote = FALSE, row.names = TRUE,
          file = paste0(save_dir, '100_location.csv'))