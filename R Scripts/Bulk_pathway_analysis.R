library(edgeR)
library(pheatmap)
library(genefilter)
library(clusterProfiler)
library(org.Mm.eg.db)
library(gplots)
library(ggplot2)
setwd("~/Documents/Lund/Vor2018_Haust2014/Data/count_results_3/output_bulk_pathway")
### reading in the data: matrix where rows are genes and colums are sample
counts <- read.delim("~/Documents/Lund/Vor2018_Haust2014/Data/count_results_3/entrenzID.final.matrix.txt", row.names = 1)

### Reorder the matrix, mef/day_5/day_7/day_8/day_9/day_11/day_15/cDC1/cCD2/pDC
counts <- counts[c(38:44,21:37,13:20,1:12,45:48)]

### remove samples that were irregular according to correlation (same samples that were used in paper)
counts <- counts[-c(4, 6,10,14,18,22,34,37,39:40,45)]

sampleCondition <- c(rep("mef", 5), rep("day_5", 4), 
                     rep("day_7", 3), rep("day_8", 3), rep("day_9", 3),
                     rep("day_11", 4), rep("day_15", 4), 
                     rep("cDC1",4), rep("cDC2", 4), rep("pDC", 3))

### filter out lowly expressed gene vol 1:
rv <- rowVars(as.matrix(counts))
counts <- counts[rv >= 1,]

subset.test1 <- c()
for (j in 1:dim(counts)[1]){
  x <- as.numeric(counts[j,])
  list.cell <- unique(as.character(sampleCondition))
  list.all <- c()
  for (i in list.cell){
    list.all <- c(list.all, sum(x[which(sampleCondition %in% i)]!=0))
  }
  if (sum(list.all >= 2) >= 7) {
    subset.test1 <- c(subset.test1, TRUE)
  }
  else{
    subset.test1 <- c(subset.test1, FALSE)}
}

counts = counts[subset.test1,]
### steps in EdgeR:
cds <- DGEList(counts=counts,group=sampleCondition)
cds <-calcNormFactors(cds)
cds <- estimateDisp(cds)

### look at the exp for indivitual gene
genes <- c("Mitf", "Tfe3", "Tfec", "Flcn", "Fnip1", "Fnip2") 
genes <- c("Ctss", "Cd68", "M6pr", "Mcl1", "Lamp2", "Lrrk2", "Myd88", "Ticam1", "Tfec", "Stat1") 

cbPalette=c("grey", "darkgoldenrod3", "green", "aquamarine4", "darkred", "chocolate1", "lightcoral","dodgerblue4","darkgreen", "darkorchid4")

for (i in 1:length(genes)) {
  gene <- genes[i]

  gene2 <- bitr(gene,  fromType = "SYMBOL", toType = "ENTREZID",
                OrgDb = "org.Mm.eg.db")[,2]
  test <- data.frame(counts = log(cds$counts[gene2,]+1), sample = sampleCondition)
  test$sample <- factor(test$sample, levels = c("mef", "day_5", "day_7", "day_8", "day_9", "day_11", "day_15", "cDC1", "cDC2", "pDC"))
  pdf(paste0("knockout_genes//mouse_exp_", gene,".pdf"))
  print(ggboxplot(test, "sample", "counts", ylab = "log counts",
            color = "sample", palette =cbPalette,
            add = "jitter", title=gene)+ labs(title=gene) + theme(plot.title = element_text(size=32, hjust = 0.5,face = "bold"), legend.position = "none"))
  dev.off()
}


###

###################################################
#### Pairwise comparison:

exact_test <- function(x,y, pvalue, logFC) {
  test_control <- exactTest(cds, pair=c(x,y))
  de_test_control <- test_control$table
  de_test_control$p.adjust <- p.adjust(de_test_control$PValue, "BH") 
  de_test_control<- de_test_control[de_test_control$p.adjust < pvalue,]
  de_test_control <- de_test_control[abs(de_test_control$logFC) > logFC,]
}

ded5_mef <- exact_test("mef", "day_5", 0.05, 0.5)
ded9_d5 <- exact_test("day_5", "day_9", 0.05, 0.5)
ded15_d9 <- exact_test("day_9", "day_15", 0.05, 0.5)
ded9_cDC1 <- exact_test("day_9", "cDC1", 0.05, 0.5)



### Save gene lists
save_gene_list <- function(x) {
  x_gene <- x
  x_gene$gene_symbol <- mapIds(org.Mm.eg.db, rownames(x), keytype="ENTREZID",
                                      column="SYMBOL")
  x_gene <- x_gene[c(5,1,2,3,4)]
  return(x_gene[order(abs(x_gene$logFC), decreasing = T),])
}

write.table(save_gene_list(ded5_mef), file="gene_lists_results/d5_mef_results.txt", quote = F, sep="\t", row.names = T, col.names = NA)
write.table(save_gene_list(ded9_d5), file="gene_lists_results/d9_d5_results.txt", quote = F, sep="\t", row.names = T, col.names = NA)
write.table(save_gene_list(ded15_d9), file="gene_lists_results/d15_d9_results.txt", quote = F, sep="\t", row.names = T, col.names = NA)
write.table(save_gene_list(ded9_cDC1), file="gene_lists_results/d9_cDC1_results.txt", quote = F, sep="\t", row.names = T, col.names = NA)


##### Scale the cpm and define max min exp and color
cpm_counts_sig_mef_idc_upDown <- cpm(cds)
counts_scaled_sig_mef_id_upDown <-t(scale(t(cpm_counts_sig_mef_idc_upDown)))

MinMax <- function(data, min, max){
  data2 <- data
  data2[data2 < min] = min
  data2[data2 > max] = max
  return(data2)
}

counts_scaled_sig_mef_id_upDown = MinMax(counts_scaled_sig_mef_id_upDown, -2.5, 2.5)
mycol <- colorpanel(100, "purple", "black", "yellow")

#### Plot heatmaps:
plot_heatmaps <- function(x, x_dataFrame) {
  plot <- pheatmap(as.matrix(counts_scaled_sig_mef_id_upDown[rownames(x),rownames(x_dataFrame)]),
                                      col=mycol, show_colnames = F, annotation_col=x_dataFrame,
                                      labels_col = x_dataframe[,1], cutree_rows = 2, annotation_names_row = F,
                                      show_rownames = F, cluster_cols = F,
                                      clustering_distance_rows = "correlation",
                                      clustering_method = "ward.D2", annotation_names_col=F)
  return(plot)
}

ded5_mef_dataframe <- as.data.frame(as.character(sampleCondition[1:9]), row.names = colnames(cpm_counts_sig_mef_idc_upDown[,1:9]))
colnames(ded5_mef_dataframe) <- "Cell_type"
pdf("d5_mef.pdf")
plot_heatmaps(ded5_mef, ded5_mef_dataframe)
dev.off()

ded9_d5_dataframe <- as.data.frame(as.character(sampleCondition[c(6:9, 16:18)]), row.names = colnames(cpm_counts_sig_mef_idc_upDown[,c(6:9, 16:18)]))
colnames(ded9_d5_dataframe) <- "Cell_type"
pdf("d9_d5.pdf")
plot_heatmaps(ded9_d5, ded9_d5_dataframe)
dev.off()

d15_d9_dataframe <- as.data.frame(as.character(sampleCondition[c(16:18, 23:26)]), row.names = colnames(cpm_counts_sig_mef_idc_upDown[,c(16:18, 23:26)]))
colnames(d15_d9_dataframe) <- "Cell_type"
pdf("d15_d9.pdf")
plot_heatmaps(ded15_d9, d15_d9_dataframe)
dev.off()

d9_cDC1_dataframe <- as.data.frame(as.character(sampleCondition[c(16:18,27:30)]), row.names = colnames(cpm_counts_sig_mef_idc_upDown[,c(16:18,27:30)]))
colnames(d9_cDC1_dataframe) <- "Cell_type"
pdf("d9_cDC1.pdf")
plot_heatmaps(ded9_cDC1, d9_cDC1_dataframe)
dev.off()



### Sett up cluster for pahtway analysis:
define_cluster <- function(x, logFC, name) {
  cluster <- list()
  cluster$Patways_up <- rownames(x[x$logFC > logFC, ])
  cluster$Patways_down<-  rownames(x[x$logFC < -logFC, ])
  write.table(rownames(x[abs(x$logFC) > logFC,]), file=paste0("gene_lists_results/", name ,"_", logFC, ".txt"), sep="\t", quote = F, row.names = F, col.names=F)
  return(cluster)
}

d5_mef_cluster <- define_cluster(ded5_mef, 0.5, name="d5_mef")
d9_d5_cluster <- define_cluster(ded9_d5, 0.5, name="d9_d5")
d15_d9_cluster <- define_cluster(ded15_d9, 0.5, name="d15_d9")
d9_cDC1_cluster <- define_cluster(ded9_cDC1, 0.5, name="d9_cDC1")

### GO Pathway analysis Using clusterProfiler:
pathway_analysis <- function(x, ont, name) {
  results <- compareCluster(geneClusters = x, fun="enrichGO", OrgDb="org.Mm.eg.db", ont=ont, qvalueCutoff=0.05)
  saveRDS(results, file=paste0("compareCluster/", name, "_", ont, ".Rdata"))
  return(results)
}
#### CC
d5_mef_CC <- pathway_analysis(d5_mef_cluster, "CC", name="d5_mef")
d9_d5_CC <- pathway_analysis(d9_d5_cluster, "CC", name="d9_d5")
d15_d9_CC <- pathway_analysis(d15_d9_cluster, "CC", name="d15_d9")
d9_cDC1_CC <- pathway_analysis(d9_cDC1_cluster, "CC ", name="d9_cDC1")


### Plot pathway results (needs to be one by one for optimal size)
pdf("sign_pathways/d5_mef_CC.pdf", width = 11, height = 14, useDingbats=FALSE)
dotplot(d5_mef_CC, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("sign_pathways/d9_d5_CC.pdf", width = 11, height = 14, useDingbats=FALSE)
dotplot(d9_d5_CC, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("sign_pathways/d15_d9_CC.pdf", width = 11, height = 14, useDingbats=FALSE)
dotplot(d15_d9_CC, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("sign_pathways/d9_cDC1_CC.pdf", width = 13, height = 14, useDingbats=FALSE)
dotplot(d9_cDC1_CC, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#### BP
d5_mef_BP <- pathway_analysis(d5_mef_cluster, "BP", name="d5_mef")
d9_d5_BP <- pathway_analysis(d9_d5_cluster, "BP", name="d9_d5")
d15_d9_BP <- pathway_analysis(d15_d9_cluster, "BP", name="d15_d9")
d9_cDC1_BP <- pathway_analysis(d9_cDC1_cluster, "BP ", name="d9_cDC1")

### Plot pathway results (needs to be one by one for optimal size)
pdf("sign_pathways/d5_mef_BP.pdf", width = 12, height = 14, useDingbats=FALSE)
dotplot(d5_mef_BP, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("sign_pathways/d9_d5_BP.pdf", width = 18, height = 14, useDingbats=FALSE)
dotplot(d9_d5_BP, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("sign_pathways/d15_d9_BP.pdf", width = 13, height = 14, useDingbats=FALSE)
dotplot(d15_d9_BP, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("sign_pathways/d9_cDC1_BP.pdf", width = 13, height = 14, useDingbats=FALSE)
dotplot(d9_cDC1_BP, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#### MF
d5_mef_MF <- pathway_analysis(d5_mef_cluster, "MF", name="d5_mef")
d9_d5_MF <- pathway_analysis(d9_d5_cluster, "MF", name="d9_d5")
d15_d9_MF <- pathway_analysis(d15_d9_cluster, "MF", name="d15_d9")
d9_cDC1_MF <- pathway_analysis(d9_cDC1_cluster, "MF ", name="d9_cDC1")

pdf("sign_pathways/d5_mef_MF.pdf", width = 14, height = 14, useDingbats=FALSE)
dotplot(d5_mef_MF, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("sign_pathways/d9_d5__MF.pdf", width = 11, height = 14, useDingbats=FALSE)
dotplot(d9_d5_MF, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("sign_pathways/d15_d9_MF.pdf", width = 11, height = 14, useDingbats=FALSE)
dotplot(d15_d9_MF, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("sign_pathways/d9_cDC1_MF.pdf", width = 11, height = 14, useDingbats=FALSE)
dotplot(d9_cDC1_MF, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


### KEGG Pathway analysis Using clusterProfiler:
KEGG_pathway_analysis <- function(x, name) {
  results <- compareCluster(geneClusters = x, fun="enrichKEGG", organism="mmu", keyType="kegg", qvalueCutoff=0.05)
  saveRDS(results, file=paste0("compareCluster/", name, "_", "KEGG.Rdata"))
  return(results)
}

d5_mef_KEGG <- KEGG_pathway_analysis(d5_mef_cluster, name="d5_mef")
d9_d5_KEGG <- KEGG_pathway_analysis(d9_d5_cluster, name="d9_d5")
d15_d9_KEGG <- KEGG_pathway_analysis(d15_d9_cluster, name="d15_d9")
d9_cDC1_KEGG <- KEGG_pathway_analysis(d9_cDC1_cluster, name="d9_cDC1")

pdf("sign_pathways/d5_mef_KEGG.pdf", width = 13, height = 14, useDingbats=FALSE)
dotplot(d5_mef_KEGG, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("sign_pathways/d9_d5_KEGG.pdf", width = 13, height = 14, useDingbats=FALSE)
dotplot(d9_d5_KEGG, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("sign_pathways/d15_d9_KEGG.pdf", width = 13, height = 14, useDingbats=FALSE)
dotplot(d15_d9_KEGG, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("sign_pathways/d9_cDC1_KEGG.pdf", width = 13, height = 14, useDingbats=FALSE)
dotplot(d9_cDC1_KEGG, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


############################################
## Plot kegg pathways with genes dif exp colored red
# d5_mef_kegg_logFC[which(d5_mef_kegg_logFC[,"Description"] == "Phagosome")]


browsKEGG_mine <- function(compareCluster_results, pathway){
  genes_enriched <- compareCluster_results[,"geneID"][which(compareCluster_results[,"ID"] == pathway)]
  site <- paste0("https://www.kegg.jp/kegg-bin/show_pathway?", pathway,"/", genes_enriched)
  
  return(browseURL(site))
}

browsKEGG_mine(d5_mef_kegg_logFC, "mmu04621")
browsKEGG_mine(d9_d5_kegg_logFC, "mmu04142")
browsKEGG_mine(d9_cDC1_kegg_logFC, "mmu04142")
