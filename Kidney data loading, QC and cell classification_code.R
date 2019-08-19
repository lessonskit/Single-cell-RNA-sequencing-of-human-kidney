library(Seurat)
library(magrittr)
library(harmony)
library(dplyr)

#Kidney data loading
K1.data <- Read10X(data.dir = "/home/yuzhenyuan/single cell data/kid/kidney1/GRCh38")
K1 <- CreateSeuratObject(raw.data = K1.data, min.cells = 8, min.genes = 200, project = "kidney1")
K21.data <- Read10X(data.dir = "/home/yuzhenyuan/single cell data/kid/kidney2/GRCh38")
K21 <- CreateSeuratObject(raw.data = K21.data, min.cells = 6, min.genes = 200, project = "kidney2")
K22.data <- Read10X(data.dir = "/home/yuzhenyuan/single cell data/kid/kidney3/GRCh38")
K22 <- CreateSeuratObject(raw.data = K22.data, min.cells = 10, min.genes = 200, project = "kidney3")
K2<- MergeSeurat(object1 = K21,object2 = K22,add.cell.id1 = "mer1", add.cell.id2 = "mer2",project = "kidm1")
kid<- MergeSeurat(object1 = K1,object2 = K2,add.cell.id1 = "batch1", add.cell.id2 = "batch2",project = "kidney.merge")

#quality control
mito.genes <- grep(pattern = "^MT-", x = rownames(x = kid@raw.data), value = TRUE)
percent.mito <- Matrix::colSums(kid@raw.data[mito.genes, ]) / Matrix::colSums(kid@raw.data)
kid <- AddMetaData(object = kid,metadata = percent.mito,col.name = "percent.mito")
VlnPlot(object = kid, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
kid <- FilterCells(object = kid, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2000, 0.3))
kid <- NormalizeData(object = kid, normalization.method = "LogNormalize", scale.factor = 10000)
kid <- FindVariableGenes(object = kid, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
s.genes <-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
kid<- CellCycleScoring(object = kid, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = TRUE)
kid <- ScaleData(object = kid, vars.to.regress = c("nUMI", "percent.mito","S.Score", "G2M.Score"))

#Eliminate batch effects with harmony and cell classification
kid@var.genes <- split(row.names(kid@meta.data), kid@meta.data$orig.ident) %>% lapply(function(cells_use) {
    temp_data <- SubsetData(kid, cells_use)
    temp_data %<>% FindVariableGenes(do.plot = FALSE)
    head(rownames(temp_data@hvg.info), 1000)
}) %>% unlist %>% unique
kid %<>% ScaleData()
kid %<>% RunPCA(pc.genes = kid@var.genes, pcs.compute = 20, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
options(repr.plot.height = 3, repr.plot.width = 6)
system.time(kid %<>% RunHarmony("orig.ident", theta = 2, plot_convergence = TRUE))
options(repr.plot.height = 6, repr.plot.width = 12)
p1 <- DimPlot(object = kid, reduction.use = "harmony", pt.size = .1, group.by = "orig.ident", do.return = T)
p2 <- VlnPlot(object = kid, features.plot = "Harmony1", group.by = "orig.ident", do.return = TRUE)
plot_grid(p1,p2)
PrintDim(object = kid, reduction.type = "harmony", dims.print = 1:5, genes.print = 10)
DimHeatmap(object = kid, reduction.type = "harmony", cells.use = 500, dim.use = 1:6, do.balanced = TRUE)
system.time({
kid %<>% RunTSNE(reduction.use = "harmony", dims.use = 1:20, do.fast = T)
kid %<>% FindClusters(reduction.type = "harmony", resolution = 0.6, dims.use = 1:20, force.recalc = TRUE, print.output = FALSE)
})
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7,8,9,10,11,12)
new.cluster.ids <- c(1, 2, 3, 4, 5, 6, 7,8,9,10,11,12,13)
kid@ident <- plyr::mapvalues(x = kid@ident, from = current.cluster.ids, to = new.cluster.ids)
p3 <- TSNEPlot(kid, do.return = T, pt.size = 0.5, group.by = "orig.ident")
p4 <- TSNEPlot(kid, do.label = T, do.return = T, pt.size = 1)
plot_grid(p3, p4)

#Calculating differentially expressed genes (DEGs) and Save rds file
kidney.markers <- FindAllMarkers(object = kid, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(kidney.markers,sep="\t",file="/home/yuzhenyuan/single cell data/kid/har/0.6_20.xls")
saveRDS(kid,file="/home/yuzhenyuan/single cell data/kid/har/0.6_20.rds")

#Some visual figure generation
DoHeatmap(object = kid, genes.use = c("SLC13A3","SLC34A1","GPX3","DCXR","SLC17A3","SLC22A8","SLC22A7","GNLY","NKG7","CD3D","CD3E","LYZ","CD14","KRT8","KRT18","CD24","VCAM1","UMOD","DEFB1","CD79A","CD79B","CLDN8","AQP2","ATP6V1G3","ATP6V0D2","TMEM213"), slim.col.label = TRUE, remove.key = TRUE)
TSNEPlot(kid, do.return = T, pt.size = 1, group.by = "orig.ident")
TSNEPlot(kid, do.label = T, do.return = T, pt.size = 1)
TSNEPlot(kid, do.return = T, pt.size = 1, group.by = "Phase")
VlnPlot(object = kid, point.size.use =0, ident.include= c(1,2,3,4,5,6), features.plot = c("GPX3", "DCXR","SLC13A3","SLC34A1","SLC22A8","SLC22A7"))
VlnPlot(object = kid, point.size.use =0, ident.include= c(12,13), features.plot = c("AQP2", "ATP6V1B1","ATP6V0D2","ATP6V1G3"))
FeaturePlot(object = kid, features.plot = c("NKG7","GNLY","CD3D","IL7R"), cols.use = c("grey", "blue"), reduction.use = "tsne")

#Select a subset of PT cells
PT <- SubsetData(kid, ident.use = c(1,2,3,4,5,6), subset.raw = T)
saveRDS(kid,file="/home/yuzhenyuan/single cell data/kid/har/PT.rds")

