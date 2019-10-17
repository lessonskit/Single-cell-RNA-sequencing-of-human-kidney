library(Seurat)
library(magrittr)
library(harmony)
library(dplyr)

#Kidney data loading
K1.data <- Read10X(data.dir = "/staging/yuzhenyuan/kidney1/kidney1_result/outs/filtered_feature_bc_matrix")
K1 <- CreateSeuratObject(counts = K1.data, project = "kidney1", min.cells = 8, min.features = 200)
K2.data <- Read10X(data.dir = "/staging/yuzhenyuan/kidney2/combine/kidney2_result/outs/filtered_feature_bc_matrix")
K2 <- CreateSeuratObject(counts = K2.data, project = "kidney2", min.cells = 6, min.features = 200)
K3.data <- Read10X(data.dir = "/staging/yuzhenyuan/kidney3/kidney3_result/outs/filtered_feature_bc_matrix")
K3 <- CreateSeuratObject(counts = K3.data, project = "kidney3", min.cells = 10, min.features = 200)
kid <- merge(x = K1, y = list(K2, K3))

#quality control
kid[["percent.mt"]] <- PercentageFeatureSet(kid, pattern = "^MT-")
VlnPlot(kid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(kid, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(kid, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
kid <- subset(kid, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 30)
kid <- NormalizeData(kid, normalization.method = "LogNormalize", scale.factor = 10000)
kid <- NormalizeData(kid)
kid <- FindVariableFeatures(kid, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(kid), 10)
plot1 <- VariableFeaturePlot(kid)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
s.genes <-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
kid <- CellCycleScoring(kid, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
all.genes <- rownames(kid)
kid <- ScaleData(kid, vars.to.regress = c("S.Score", "G2M.Score"), features = all.genes)

#Eliminate batch effects with harmony and cell classification
kid <- RunPCA(kid, pc.genes = kid@var.genes, npcs = 20, verbose = FALSE)
options(repr.plot.height = 2.5, repr.plot.width = 6)
kid <- kid %>% 
    RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(kid, 'harmony')
harmony_embeddings[1:5, 1:5]
kid <- kid %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.25) %>% 
    identity()
new.cluster.ids <- c(1, 2, 3, 4, 5, 6, 7,8,9,10)
names(new.cluster.ids) <- levels(kid)
kid <- RenameIdents(kid, new.cluster.ids)

#Calculating differentially expressed genes (DEGs) and Save rds file
kid.markers <- FindAllMarkers(kid, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(kid.markers,sep="\t",file="/home/yuzhenyuan/Seurat/0.2_20.xls")
saveRDS(kid,file="/home/yuzhenyuan/kid/har/0.25_20.rds")

#Some visual figure generation
DimPlot(kid, reduction = "umap", group.by = "orig.ident", pt.size = .1, split.by = 'orig.ident')
DimPlot(kid, reduction = "umap", group.by = "Phase", pt.size = .1)
DimPlot(kid, reduction = "umap", label = TRUE, pt.size = .1)
DoHeatmap(kid, features = c("SLC13A3","SLC34A1","GPX3","DCXR","SLC17A3","SLC22A8","SLC22A7","GNLY","NKG7","CD3D","CD3E","LYZ","CD14","KRT8","KRT18","CD24","VCAM1","UMOD","DEFB1","CLDN8","AQP2","CD79A","CD79B","ATP6V1G3","ATP6V0D2","TMEM213"))
VlnPlot(kid, pt.size =0, idents= c(1,2,3), features = c("GPX3", "DCXR","SLC13A3","SLC34A1","SLC22A8","SLC22A7"))
VlnPlot(kid, idents= c(8,10), features = c("AQP2", "ATP6V1B1","ATP6V0D2","ATP6V1G3"))

##tSNE Plot
kid <-RunTSNE(kid, reduction = "harmony", dims = 1:20)
TSNEPlot(kid, do.label = T, label = TRUE, do.return = T, pt.size = 1)
TSNEPlot(kid, do.return = T, pt.size = 1, group.by = "orig.ident", split.by = 'orig.ident')
TSNEPlot(kid, do.return = T, pt.size = 1, group.by = "Phase")

#Select a subset of PT cells
PT <- SubsetData(kid, ident.use = c(0,1,2), subset.raw = T)
saveRDS(PT,file="/staging/yuzhenyuan/PT.rds")
