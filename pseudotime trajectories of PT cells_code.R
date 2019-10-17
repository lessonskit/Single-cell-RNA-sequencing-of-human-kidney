library(Seurat)

#PT data loading
PT=readRDS(file="/staging/yuzhenyuan/PT.rds")

#Seurat data convert to monocle data
library(monocle)
data <- as(as.matrix(PT@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = PT@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
my_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
my_cds <- estimateSizeFactors(my_cds)
my_cds <- estimateDispersions(my_cds)
my_cds <- detectGenes(my_cds, min_expr = 0.1)
head(fData(my_cds))
head(pData(my_cds))
pData(my_cds)$UMI <- Matrix::colSums(exprs(my_cds))
disp_table <- dispersionTable(my_cds)
head(disp_table)
table(disp_table$mean_expression>=0.1)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
my_cds <- setOrderingFilter(my_cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(my_cds)
expressed_genes <- row.names(subset(fData(my_cds), num_cells_expressed >= 10))
my_cds_subset <- my_cds[expressed_genes, ]
my_cds_subset
head(pData(my_cds_subset))
my_cds_subset <- detectGenes(my_cds_subset, min_expr = 0.1)
fData(my_cds_subset)$use_for_ordering <- fData(my_cds_subset)$num_cells_expressed > 0.05 * ncol(my_cds_subset)
table(fData(my_cds_subset)$use_for_ordering)
plot_pc_variance_explained(my_cds_subset, return_all = FALSE)
my_cds_subset <- reduceDimension(my_cds_subset,max_components = 2,norm_method = 'log',num_dim = 10,reduction_method = 'tSNE',verbose = TRUE)
my_cds_subset <- clusterCells(my_cds_subset, verbose = FALSE)
plot_rho_delta(my_cds_subset, rho_threshold = 2, delta_threshold = 10)
my_cds_subset <- clusterCells(my_cds_subset,rho_threshold = 2,delta_threshold = 10,skip_rho_sigma = T,verbose = FALSE)
table(pData(my_cds_subset)$Cluster)
plot_cell_clusters(my_cds_subset)
head(pData(my_cds_subset))
clustering_DEG_genes <- differentialGeneTest(my_cds_subset,fullModelFormulaStr = '~Cluster',cores = 22)
dim(clustering_DEG_genes)
library(dplyr)
clustering_DEG_genes %>% arrange(qval) %>% head()
my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
my_cds_subset <- setOrderingFilter(my_cds_subset, ordering_genes = my_ordering_genes)
my_cds_subset <- reduceDimension(my_cds_subset, method = 'DDRTree')
my_cds_subset <- orderCells(my_cds_subset)

#pseudotime trajectories calculated
plot_cell_trajectory(my_cds_subset, color_by = "State")
plot_cell_trajectory(my_cds_subset, color_by = "res.0.6")
plot_cell_trajectory(my_cds_subset, color_by = "orig.ident")
head(pData(my_cds_subset))
my_pseudotime_de <- differentialGeneTest(my_cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 22)
my_pseudotime_de %>% arrange(qval) %>% head()
my_pseudotime_de %>% arrange(qval) %>% head() %>% select(gene_short_name) -> my_pseudotime_gene
plot_cell_trajectory(my_cds_subset, color_by = "Pseudotime")

#"A" stand for top 6 genes of affecting the fate decisions
A=c("AKR1A1","PDZK1","AKR7A3","AKR7A2","FABP3","GADD45A")
my_pseudotime_gene <-A
plot_genes_in_pseudotime(my_cds_subset[my_pseudotime_gene,])

#Calculate the heat map of the top 50 genes
my_pseudotime_de %>% arrange(qval) %>% head(50) %>% select(gene_short_name) -> gene_to_cluster
gene_to_cluster <- gene_to_cluster$gene_short_name
my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds_subset[gene_to_cluster,],num_clusters = 3,cores = 22,show_rownames = TRUE,return_heatmap = TRUE)
