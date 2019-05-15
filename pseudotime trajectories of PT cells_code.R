library(Seurat)

#PT data loading
PT=readRDS(file = "/home/yuzhenyuan/single cell data/kid/PT.rds")

#Seurat data convert to monocle data
library(monocle)
my_cds=importCDS(PT)
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
clustering_DEG_genes <- differentialGeneTest(my_cds_subset,fullModelFormulaStr = '~Cluster',cores = 16)
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
head(pData(my_cds_subset))
my_pseudotime_de <- differentialGeneTest(my_cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 16)
my_pseudotime_de %>% arrange(qval) %>% head()
my_pseudotime_de %>% arrange(qval) %>% head() %>% select(gene_short_name) -> my_pseudotime_gene
plot_cell_trajectory(my_cds_subset, color_by = "Pseudotime")

#"A" stand for top 6 genes of affecting the fate decisions
A=c("RPS8","FAM151A","S100A1","RPS27","ATP1B1","G0S2")
my_pseudotime_gene <-A
plot_genes_in_pseudotime(my_cds_subset[my_pseudotime_gene,])

#Calculate the heat map of the top 50 genes
my_pseudotime_de %>% arrange(qval) %>% head(50) %>% select(gene_short_name) -> gene_to_cluster
gene_to_cluster <- gene_to_cluster$gene_short_name
my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds_subset[gene_to_cluster,],num_clusters = 3,cores = 8,show_rownames = TRUE,return_heatmap = TRUE)
