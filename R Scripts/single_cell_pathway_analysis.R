library(Seurat)
library(dplyr)
library(monocle)
library(edgeR)
library(clusterProfiler)
library(pheatmap)
library(gplots)
setwd("~/Documents/Lund/Vor2018_Haust2014/SingleCell/output_pathway_analysis/")

pbmc <- readRDS("~/Documents/Lund/Vor2018_Haust2014/SingleCell/filtered.idc.rds")
pbmc

DimPlot(object = pbmc,reduction.key = "tSNE_1", group.by="cell_type", cols=c("grey", "darkgoldenrod3", "darkred", "dodgerblue4","darkgreen", "darkorchid4"))
head(x = Idents(object = pbmc), 5)

Idents(object = pbmc) <- pbmc$cell_type ## Run for supervised

pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = T, min.pct = 0.1, logfc.threshold = 0.3)
table(pbmc.markers$cluster)

pbmc.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 2, wt = avg_logFC)
top10 <- pbmc.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_logFC)

pdf("~/Documents/Lund/heatmap_logFC0.3_sup_vol2.pdf")
DoHeatmap(object = pbmc, features = top10$gene, group.by = "cell_type") + NoLegend()
dev.off()

#### Filtering genes:
#Connecting toghter name and cell type:
number.list <- unlist(strsplit(colnames(pbmc@assays$RNA@counts),"-"))[seq(2,2*3693,2)]
HSMM_sample_sheet_name <- data.frame(cell_name = number.list, row.names = colnames(pbmc@assays$RNA@counts))
HSMM_sample_sheet <- data.frame(cell_type = ifelse(number.list == 1, "HEF",
                                                   ifelse(number.list == 2, "idc_3",
                                                          ifelse(number.list==3, "idc_9",
                                                                 ifelse(number.list==4, "cDC1",
                                                                        ifelse(number.list==5, "cDC2","pDC"))))), row.names = colnames(pbmc@assays$RNA@counts))
# Steps in Monocle
pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
HSMM_gene_annotation <- data.frame(gene_short_name = rownames(pbmc@assays$RNA@counts), row.names = rownames(pbmc@assays$RNA@counts))
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
HSMM <- newCellDataSet(as.matrix(pbmc@assays$RNA@counts), phenoData = pd, 
                       featureData = fd, expressionFamily=negbinomial.size())
pData(HSMM)$cell_name = HSMM_sample_sheet_name
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))
pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6] ## no cells above 1e6

# Plot exp of indivitual genes
gene_list_sel <- c("CTSS","CD68", "M6PR", "MCL1", "LAMP2", "LRRK2", "MYD88", "TICAM1", "TFEC", "STAT1")
# gene_list_sel <- bitr(c("54751", "1969", "3339", "55616", "2275", "85440" ,"23032", "91624", "1266", "8394", "9181"),
#                       fromType = "ENTREZID", toType = "SYMBOL",
#                       OrgDb = "org.Hs.eg.db")[,2]

for (i in 1:length(gene_list_sel)) {
  to_be_tested <- row.names(subset(fData(HSMM),
                                   gene_short_name %in% gene_list_sel[i]))
  cds_subset <- HSMM[to_be_tested,]
  cds_subset$cell_type <- factor(cds_subset$cell_typ, levels = c("HEF", "idc_3", "idc_9", "cDC1", "cDC2", "pDC"))
  
  cbPalette=c("grey", "darkgoldenrod3", "darkred", "dodgerblue4","darkgreen", "darkorchid4")
  
  pdf(paste0("~/Downloads/", gene_list_sel[i], "_exp.pdf"))
  print(monocle::plot_genes_jitter(cds_subset=cds_subset,
                             grouping = "cell_type",
                             color_by = "cell_type",
                             nrow= 2,
                             ncol = NULL,
                             plot_trend = TRUE,  panel_order=gene_list_sel)+ scale_color_manual(values = cbPalette))
  dev.off()
  
}
# Violin plot
# monocle::plot_genes_violin(cds_subset,relative_expr=T,
#                            grouping = "cell_type",
#                            color_by = "cell_type",
#                            nrow=2,
#                            ncol= 3,
#                            plot_trend = TRUE, panel_order=gene_list_sel)+ scale_fill_manual(values = cbPalette)

#### Pre-work for edgeR
counts <-  exprs(HSMM[expressed_genes,])
colnames(counts) <- pData(HSMM)[,1]

cds <- DGEList(counts=counts,group=pData(HSMM)[,1])
cds <-calcNormFactors(cds)
cds <- estimateDisp(cds)

###### DE analysis
gene.sign.edgeR = list()
sample.list.edgeR <- levels(pData(HSMM)$cell_type)
for (i in 1:(length(sample.list.edgeR)-1)){
  for (j in (i+1):length(sample.list.edgeR)) {
    compare.sample.edgeR <- paste0(sample.list.edgeR[i], "|", sample.list.edgeR[j])
    print(compare.sample.edgeR)
    compare.edgeR <- exactTest(cds, pair=c(sample.list.edgeR[i], sample.list.edgeR[j]))
    compare.edgeR <- compare.edgeR$table
    compare.edgeR$p.adjust <- p.adjust(compare.edgeR$PValue, "BH")
    
    sig_gene.edgeR<- compare.edgeR[compare.edgeR$p.adjust < 0.05,]
    sig_gene.edgeR <- sig_gene.edgeR[abs(sig_gene.edgeR$logFC) > 0.5,]
    
    write.table(sig_gene.edgeR, file=paste0("edgeR_sig_gene_",
                                            compare.sample.edgeR,".txt"),
                row.names = T, col.names = T, quote = F, sep="\t")
    
    gene.sign.edgeR[[compare.sample.edgeR]] = compare.edgeR
  }
}

all.genes.edgeR <- rownames(pbmc@assays$RNA@counts[expressed_genes,])
names.edgeR <- "all.genes"
for (i in 1:(length(sample.list.edgeR)-1)) {
  for (j in (i+1):length(sample.list.edgeR)) {
    compare.sample.edgeR <- paste0(sample.list.edgeR[i], "|", sample.list.edgeR[j])
    
    all.genes.edgeR <- cbind(all.genes.edgeR, gene.sign.edgeR[[compare.sample.edgeR]][,4], 
                             gene.sign.edgeR[[compare.sample.edgeR]][,1])
    
    names.edgeR <- c(names.edgeR,paste0("qval.", compare.sample.edgeR), paste0("logFC.", compare.sample.edgeR))
  }
}

colnames(all.genes.edgeR) <- names.edgeR
write.table(all.genes.edgeR, file= "edgeR_logFC_qval.txt", sep="\t", quote = F, row.names = F, col.names = T)

####################
# No need to run edgeR again:
gene.sign.edgeR <- read.table("edgeR_output/edgeR_logFC_qval.txt", header = T, row.names = 1)

edgeR_results <- function(sample1, sample2, qval, logFC) {
 results <-  gene.sign.edgeR[, c(sample1, sample2)]
 results <- results[results[,sample1] < qval,]
 results <- results[abs(results[,sample2]) > logFC,]
}

deHEF_idc3 <- edgeR_results("qval.HEF.idc_3", "logFC.HEF.idc_3", 0.05, 0.5)
deHEF_idc9 <- edgeR_results("qval.HEF.idc_9", "logFC.HEF.idc_9", 0.05, 0.5)
deidc3_idc9 <- edgeR_results("qval.idc_3.idc_9", "logFC.idc_3.idc_9", 0.05, 0.5)
ded9_cDC1 <- edgeR_results("qval.cDC1.idc_9", "logFC.cDC1.idc_9", 0.05, 0.5) 


####### 
cpm_counts_sig_hef_idc_upDown <- cpm(cds)
counts_scaled_sig_hef_id_upDown <-t(scale(t(cpm_counts_sig_hef_idc_upDown)))

MinMax <- function(data, min, max){
  data2 <- data
  data2[data2 < min] = min
  data2[data2 > max] = max
  return(data2)
}

counts_scaled_sig_hef_id_upDown2 = MinMax(counts_scaled_sig_hef_id_upDown, -2.5, 2.5)
mycol <- colorpanel(100, "purple", "black", "yellow")

##### Heatmap
make_matrix <- function(x, sample1, sample2) {
  return_matrix <- as.matrix(counts_scaled_sig_hef_id_upDown2[rownames(x),
                                                              c(which(colnames(counts_scaled_sig_hef_id_upDown2) == sample1),
                                                                which(colnames(counts_scaled_sig_hef_id_upDown2) == sample2))])
  colnames(return_matrix) <- rownames(pData(HSMM))[c(which(HSMM$cell_type==sample1), which(HSMM$cell_type==sample2))]
  return(return_matrix)
}

make_datafram <- function(sample1, sample2) {
  d3_hef_colname_dataframe <- data.frame(Cell_type = HSMM$cell_type[c(which(HSMM$cell_type == sample1), which(HSMM$cell_type == sample2))], row.names = rownames(pData(HSMM))[c(which(HSMM$cell_type==sample1), which(HSMM$cell_type== sample2))])
}

plot_heatmap <- function(matrix, df) {
  return_plot <- pheatmap(matrix, col=mycol, show_colnames = F, cutree_rows = 2, annotation_names_row = F,
                          show_rownames = F, cluster_cols = F, annotation_col = df,
                          clustering_distance_rows = "correlation",
                          clustering_method = "ward.D2", annotation_names_col = F )
}

idc3_hef <- make_matrix(deHEF_idc3, "idc_3", "HEF")
idc3_hef_df <- make_datafram("idc_3", "HEF")
png("heatmaps/d3_hef.png")
print(plot_heatmap(idc3_hef, idc3_hef_df))
dev.off()

dc9_hef <- make_matrix(deHEF_idc9, "idc_9", "HEF")
idc9_hef_df <- make_datafram("idc_9", "HEF")
png("heatmaps/d9_hef.png")
print(plot_heatmap(idc9_hef, idc9_hef_df))
dev.off()

idc9_idc3 <- make_matrix(deidc3_idc9, "idc_9", "idc_3")
idc9_idc3_df <- make_datafram("idc_9", "idc_3")
png("heatmaps/d9_d3.png")
print(plot_heatmap(idc9_idc3, idc9_idc3_df))
dev.off()

idc9_cDC1 <- make_matrix(ded9_cDC1, "idc_9", "cDC1")
idc9_cDC1_df <- make_datafram("idc_9", "cDC1")
png("heatmaps/d9_cDC1.png")
print(plot_heatmap(idc9_cDC1, idc9_cDC1_df))
dev.off()

## Pathway analysis:
define_cluster <- function(x, logFC) {
  cluster <- list()
  cluster$Patways_up <- bitr(rownames(x[x[,1] > logFC, ]),
                             fromType = "SYMBOL", toType = "ENTREZID",
                             OrgDb = "org.Hs.eg.db")[,2]
  cluster$Patways_down <-  bitr(rownames(x[x[,1] < -logFC, ]),
                                fromType = "SYMBOL", toType = "ENTREZID",
                                OrgDb = "org.Hs.eg.db")[,2]
  return(cluster)
}

d3_hef_cluster <- define_cluster(deHEF_idc3, 0.5)
d9_hef_cluster <- define_cluster(gdeHEF_idc9, 0.5)
d9_d3_cluster <- define_cluster(deidc3_idc9, 0.5)
d9_cDC1_cluster <- define_cluster(ded9_cDC1, 0.5)

## GO pathway
pathway_analysis <- function(x, ont, name) {
  results <- compareCluster(geneClusters = x, fun="enrichGO", OrgDb="org.Hs.eg.db", ont=ont, qvalueCutoff=0.05)
  saveRDS(results, file=paste0("sign_pathway/Rdata/", name, "_", ont, ".Rdata"))
  write.table(as.data.frame(results), file=paste0("sign_pathway/text_files/", name, "_", ont ,".results.txt"), quote = F, row.names = T, sep="\t")
  return(results)
}

### CC

#d3_hef_CC <- pathway_analysis(d3_hef_cluster, "CC", "d3_hef")
d3_hef_CC <- readRDS("sign_pathway/Rdata/d3_hef_CC.Rdata")
pdf("sign_pathway/plot/d3_hef_CC.pdf", width = 10, height = 12, useDingbats=FALSE)
dotplot(d3_hef_CC, showCategory=5) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#### Day 9 vs Hef
#d9_hef_CC <- pathway_analysis(d9_hef_cluster, "CC", name="d9_hef")
d9_hef_CC <- readRDS("sign_pathway/Rdata/d9_hef_CC.Rdata")
pdf("sign_pathway/d9_hef_CC.pdf", width =9, height = 12, useDingbats=FALSE)
dotplot(d9_hef_CC, showCategory=5) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#### Day 9 vs Day 3
#d9_d3_CC <- pathway_analysis(d9_d3_cluster, "CC", name="d9_d3")
d9_d3_CC <- readRDS("sign_pathway/Rdata/d9_d3_CC.Rdata")
pdf("sign_pathway/d3_d9_CC.pdf", width = 9, height = 12, useDingbats=FALSE)
dotplot(d9_d3_CC, showCategory=5) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#### Day 9 vs cDC1
#d9_cDC1_CC <- pathway_analysis(d9_cDC1_cluster, "CC", name="d9_cDC1")
d9_cDC1_CC <- readRDS("sign_pathway/Rdata/d9_cDC1_CC.Rdata")
pdf("sign_pathway/d9_cDC1_CC.pdf", width = 9, height = 12, useDingbats=FALSE)
dotplot(d9_cDC1_CC, showCategory=5) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

### BP
#d3_hef_BP <- pathway_analysis(d3_hef_cluster, "BP", "d3_hef")
d3_hef_BP <- readRDS("sign_pathway/Rdata/3_hef_BP.Rdata")
pdf("sign_pathway/d3_hef_BP.pdf", width = 16, height = 12, useDingbats=FALSE)
dotplot(d3_hef_BP, showCategory=5) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#### Day 9 vs Hefs 
#d9_hef_BP <- pathway_analysis(d9_hef_cluster, "BP", "d9_hef")
d9_hef_BP <- readRDS("sign_pathway/Rdata/d9_hef_BP.Rdata")
pdf("sign_pathway/d9_hef_BP.pdf", width = 13, height = 12, useDingbats=FALSE)
dotplot(d9_hef_BP, showCategory=5) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#### Day 9 vs Day 3
#d9_d3_BP <- pathway_analysis(d9_d3_cluster, "BP", "d9_d3")
d9_d3_BP <- readRDS("sign_pathway/Rdata/d9_d3_BP.Rdata")
pdf("sign_pathway/d3_d9_BP.pdf", width = 15, height = 12, useDingbats=FALSE)
dotplot(d9_d3_BP, showCategory=5) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#### Day 9 vs cDC1
#d9_cDC1_BP <- pathway_analysis(d9_cDC1_cluster, "BP", "d9_cDC1")
d9_cDC1_BP <- readRDS("sign_pathway/Rdata/d9_cDC1_BP.Rdata")
pdf("sign_pathway/d9_cDC1_BP_vol2.pdf", width = 14, height = 12, useDingbats=FALSE)
dotplot(d9_cDC1_BP, showCategory=5) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

### KEGG pathway analysis
KEGG_pathway_analysis <- function(x, name) {
  results <- compareCluster(geneClusters = x, fun="enrichKEGG", organism="hsa", keyType="kegg", qvalueCutoff=0.05)
  saveRDS(results, file=paste0("sign_pathway/Rdata/", name, "_", "KEGG.Rdata"))
  write.table(as.data.frame(results), file=paste0("sign_pathway/text_files/", namw, "_KEGG_results.txt"), quote = F, row.names = T, sep="\t")
  return(results)
}

### Day 3 vs Hefs
#d3_hef_KEGG <- KEGG_pathway_analysis(d3_hef_cluster, name="d3_hef")
d3_hef_KEGG <- readRDS("sign_pathway/Rdata/d3_hef_KEGG.Rdata")
pdf("sign_pathway/d3_hef_KEGG.pdf", width = 11, height = 14, useDingbats=FALSE)
dotplot(d3_hef_KEGG, showCategory=5) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

### Day 9 vs Hefs
#d9_hef_KEGG <- KEGG_pathway_analysis(d9_hef_cluster, name="d9_hef")
d9_hef_KEGG <- readRDS("sign_pathway/Rdata/d9_hef_KEGG.Rdata")
pdf("sign_pathway/d9_hef_KEGG.pdf", width = 9, height = 14, useDingbats=FALSE)
dotplot(d9_hef_KEGG, showCategory=5) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

### Day 9 vs day 3
#d9_d3_KEGG <- KEGG_pathway_analysis(d9_d3_cluster, name="d9_d3")
d9_d3_KEGG <- readRDS("sign_pathway/Rdata/d9_d3_KEGG.Rdata")
pdf("sign_pathway/d3_d9_KEGG.pdf", width = 9, height = 14, useDingbats=FALSE)
dotplot(d9_d3_KEGG, showCategory=5) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

### Day 9 vs cDC1
#d9_cDC1_KEGG <- KEGG_pathway_analysis(d9_cDC1_cluster, name="d9_cDC1")
d9_cDC1_KEGG <- readRDS("sign_pathway/Rdata/d9_cDC1_KEGG.Rdata")
pdf("sign_pathway/d9_cDC1_KEGG.pdf", width = 8, height = 14, useDingbats=FALSE)
dotplot(d9_cDC1_KEGG, showCategory=5) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


#### KEGG pathways
browsKEGG_mine <- function(compareCluster_results, pathway){
  genes_enriched <- compareCluster_results[,"geneID"][which(compareCluster_results[,"ID"] == pathway)]
  site <- paste0("https://www.kegg.jp/kegg-bin/show_pathway?", pathway,"/", genes_enriched)
  
  return(browseURL(site))
}

browsKEGG_mine(d3_mef_ck_kegg_logFC, "mmu04621")
browsKEGG_mine(d9_mef_ck_kegg_logFC, "hsa03050")
browsKEGG_mine(d3_d9_ck_kegg_logFC, "mmu04142")
browsKEGG_mine(d9_cDC1_ck_kegg_logFC, "hsa")


