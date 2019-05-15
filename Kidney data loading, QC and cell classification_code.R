library(Seurat)
library(dplyr)

#Kidney data loading
K1.data <- Read10X(data.dir = "/home/yuzhenyuan/single cell data/kid/kidney1/GRCh38")
K1 <- CreateSeuratObject(raw.data = K1.data, min.cells = 8, min.genes = 200, project = "kidney1")
K21.data <- Read10X(data.dir = "/home/yuzhenyuan/single cell data/kid/kidney2/GRCh38")
K21 <- CreateSeuratObject(raw.data = K21.data, min.cells = 6, min.genes = 200, project = "kidney2")
K22.data <- Read10X(data.dir = "/home/yuzhenyuan/single cell data/kid/kidney3/GRCh38")
K22 <- CreateSeuratObject(raw.data = K22.data, min.cells = 10, min.genes = 200, project = "kidney3")
K2<- MergeSeurat(object1 = K21,object2 = K22,add.cell.id1 = "mer1", add.cell.id2 = "mer2",project = "kidm1")
kidney<- MergeSeurat(object1 = K1,object2 = K2,add.cell.id1 = "batch1", add.cell.id2 = "batch2",project = "kidney.merge")

#quality control
mito.genes <- grep(pattern = "^MT-", x = rownames(x = kidney@raw.data), value = TRUE)
percent.mito <- Matrix::colSums(kidney@raw.data[mito.genes, ]) / Matrix::colSums(kidney@raw.data)
kidney <- AddMetaData(object = kidney,metadata = percent.mito,col.name = "percent.mito")
VlnPlot(object = kidney, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
kidney <- FilterCells(object = kidney, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2000, 0.3))
kidney <- NormalizeData(object = kidney, normalization.method = "LogNormalize", scale.factor = 10000)
kidney <- FindVariableGenes(object = kidney, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
s.genes <-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
kidney<- CellCycleScoring(object = kidney, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = TRUE)
kidney <- ScaleData(object = kidney, vars.to.regress = c("nUMI", "percent.mito","S.Score", "G2M.Score"))

#cell classification
kidney <- RunPCA(object = kidney, pc.genes = kidney@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
kidney<- FindClusters(object = kidney, reduction.type = "pca", dims.use = 1:20, resolution = 0.6, print.output = 0, save.SNN = TRUE)
kidney <- RunTSNE(object = kidney, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = kidney,do.label = TRUE)

#Save rds file
saveRDS(kidney,file="/home/yuzhenyuan/single cell data/kid/2000_0.6_0.3.rds")

#read rds file
kidney=readRDS(file="/home/yuzhenyuan/single cell data/kid/2000_0.6_0.3.rds")
kidney.markers <- FindAllMarkers(object = kidney, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(kidney.markers,sep="\t",file="/home/yuzhenyuan/single cell data/kid/2000_0.6_0.3.xls")
DoHeatmap(object = kidney, genes.use = c("SLC13A3","SLC34A1","GPX3","DCXR","SLC17A3","SLC22A8","SLC7A13","SLC16A9","SLC22A7","KRT8","KRT18","CD24","VCAM1","LYZ","CD14","GNLY","NKG7","CD3D","CD3E","CD79A","CD79B","UMOD","DEFB1","CLDN8","AQP2","SLC4A11","ATP6V1G3","ATP6V0D2","TMEM213"), slim.col.label = TRUE, remove.key = TRUE)
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6,7,8,9,10,11,12)
new.cluster.ids <- c("Proximal convoluted tubule 1", "Glomerular parietal epithelial cells", "Proximal convoluted tubule 2", "Proximal tubule", "Proximal straight tubule 1", "NK-T cells", "Monocytes", "Proximal convoluted tubule 3", "Proximal straight tubule 2", "Distal tubule", "B cells", "collecting duct principal cells", "Collecting duct intercalated cells")
kidney@ident <- plyr::mapvalues(x = kidney@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = kidney,do.label = TRUE)

#Select a subset of PT cells
PT <- SubsetData(kidney, ident.use = c(0,2,3,4,7,8), subset.raw = T)
saveRDS(PT,file="/home/yuzhenyuan/single cell data/kid/PT.rds")
