library(Seurat)
library(dplyr)
library(monocle)
library(edgeR)
library(clusterProfiler)
setwd("~/Documents/Lund/Vor2018_Haust2014/SingleCell/output_pathway_analysis/")

pbmc <- readRDS("~/Documents/Lund/Vor2018_Haust2014/SingleCell/filtered.idc.rds")
pbmc

DimPlot(object = pbmc,reduction.key = "tSNE_1", group.by="cell_type", cols=c("grey", "darkgoldenrod3", "darkred", "dodgerblue4","darkgreen", "darkorchid4"))
head(x = Idents(object = pbmc), 5)

Idents(object = pbmc) <- pbmc$cell_type ## Run for supervised

pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = T, min.pct = 0.1, logfc.threshold = 0.3)
table(pbmc.markers$cluster)

pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

pdf("~/Documents/Lund/heatmap_logFC0.3_sup_vol2.pdf")
DoHeatmap(object = pbmc, features = top10$gene, group.by = "cell_type") + NoLegend()
dev.off()

#### Filtering genes:
number.list <- unlist(strsplit(colnames(pbmc@assays$RNA@counts),"-"))[seq(2,2*3693,2)]
HSMM_sample_sheet_name <- data.frame(cell_name = number.list)
rownames(HSMM_sample_sheet_name) <- colnames(pbmc@assays$RNA@counts)
colnames(HSMM_sample_sheet_name) <- "cell_name"
HSMM_sample_sheet <- as.data.frame(ifelse(number.list == 1, "HEF",
                                          ifelse(number.list == 2, "idc_3",
                                                 ifelse(number.list==3, "idc_9",
                                                        ifelse(number.list==4, "cDC1",
                                                               ifelse(number.list==5, "cDC2","pDC"))))))
rownames(HSMM_sample_sheet) <- colnames(pbmc@assays$RNA@counts)
colnames(HSMM_sample_sheet) <- "cell_type"
pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
HSMM_gene_annotation <- data.frame(gene_short_name = rownames(pbmc@assays$RNA@counts), row.names = rownames(pbmc@assays$RNA@counts))
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
HSMM <- newCellDataSet(as.matrix(pbmc@assays$RNA@counts), phenoData = pd, 
                       featureData = fd, expressionFamily=negbinomial.size())
pData(HSMM)$cell_name = HSMM_sample_sheet_name
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))
pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6] ## no cells above 1e6

gene_list_sel <- c("CTSS","CD68", "M6PR", "MCL1", "LAMP2", "LRRK2", "MYD88", "TICAM1", "TFEC", "STAT1")
gene_list_sel <- c("TAOK1", "TAOK2", "TAOK3")
gene_list_sel <- bitr(c("54751", "1969", "3339", "55616", "2275", "85440" ,"23032", "91624", "1266", "8394", "9181"),
                      fromType = "ENTREZID", toType = "SYMBOL",
                      OrgDb = "org.Hs.eg.db")[,2]

gene_list_sel <- bitr(c("55920", "912", "911", "6233", "50489", "8673" ,"54918", "309", "3122", "3127"),
                      fromType = "ENTREZID", toType = "SYMBOL",
                      OrgDb = "org.Hs.eg.db")[,2]

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

monocle::plot_genes_violin(cds_subset,relative_expr=T,
                           grouping = "cell_type",
                           color_by = "cell_type",
                           nrow=2,
                           ncol= 3,
                           plot_trend = TRUE, panel_order=gene_list_sel)+ scale_fill_manual(values = cbPalette)

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
gene.sign.edgeR <- read.table("edgeR_logFC_qval.txt", header = T, row.names = 1)

deHEF_idc3 <- gene.sign.edgeR[,c("qval.HEF.idc_3", "logFC.HEF.idc_3")]
deHEF_idc3 <- deHEF_idc3[deHEF_idc3$qval.HEF.idc_3 < 0.05,]
deHEF_idc3 <- deHEF_idc3[abs(deHEF_idc3$logFC.HEF.idc_3) > 0.5,]

deHEF_idc9 <- gene.sign.edgeR[,c("qval.HEF.idc_9", "logFC.HEF.idc_9")]
deHEF_idc9 <- deHEF_idc9[deHEF_idc9$qval.HEF.idc_9 < 0.05,]
deHEF_idc9 <- deHEF_idc9[abs(deHEF_idc9$logFC.HEF.idc_9) > 0.5,]

deidc3_idc9 <- gene.sign.edgeR[,c("qval.idc_3.idc_9", "logFC.idc_3.idc_9")]
deidc3_idc9 <- deidc3_idc9[deidc3_idc9$qval.idc_3.idc_9 < 0.05,]
deidc3_idc9 <- deidc3_idc9[abs(deidc3_idc9$logFC.idc_3.idc_9) > 0.5,]

ded9_cDC1 <- gene.sign.edgeR[,c("qval.cDC1.idc_9", "logFC.cDC1.idc_9")]
ded9_cDC1 <- ded9_cDC1[ded9_cDC1$qval.cDC1.idc_9 < 0.05,]
ded9_cDC1 <- ded9_cDC1[abs(ded9_cDC1$logFC.cDC1.idc_9) > 0.5,]


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
rm(counts_scaled_sig_hef_id_upDown)
mycol <- colorpanel(100, "purple", "black", "yellow")

logFC <- c(0.5)
#####


deHEF_idc3_logFC <- deHEF_idc3[abs(deHEF_idc3$logFC.HEF.idc_3) > logFC,]
d3_hefmatrix <- as.matrix(counts_scaled_sig_hef_id_upDown2[rownames(deHEF_idc3_logFC),
                                                           c(which(colnames(counts_scaled_sig_hef_id_upDown2) == "idc_3"),
                                                             which(colnames(counts_scaled_sig_hef_id_upDown2) == "HEF"))])
colnames(d3_hefmatrix) <- rownames(pData(HSMM))[c(which(HSMM$cell_type=="idc_3"), which(HSMM$cell_type=="HEF"))]
d3_hef_colname_dataframe <- as.data.frame(HSMM$cell_type[c(which(HSMM$cell_type == "idc_3"), which(HSMM$cell_type == "HEF"))], row.names = rownames(pData(HSMM))[c(which(HSMM$cell_type=="idc_3"), which(HSMM$cell_type=="HEF"))])
colnames(d3_hef_colname_dataframe) <- "Cell_type"

out_d3_hef_upDown <- pheatmap(d3_hefmatrix, col=mycol, show_colnames = F, cutree_rows = 2, annotation_names_row = F,
                              show_rownames = F, cluster_cols = F, annotation_col = d3_hef_colname_dataframe,
                              clustering_distance_rows = "correlation",
                              clustering_method = "ward.D2", annotation_names_col = F )
png(paste0("heatmaps/d3_hef_logFC",logFC, ".png"))
print(out_d3_hef_upDown)
dev.off()

d3_hef_cluster <- list()
d3_hef_cluster$Patways_up_day_3_vs_hef <- bitr(rownames(deHEF_idc3[deHEF_idc3$logFC.HEF.idc_3 > logFC, ]),
                                               fromType = "SYMBOL", toType = "ENTREZID",
                                               OrgDb = "org.Hs.eg.db")[,2]
d3_hef_cluster$Patways_down_day_3_vs_hef <- bitr(rownames(deHEF_idc3[deHEF_idc3$logFC.HEF.idc_3 < -logFC, ]),
                                                 fromType = "SYMBOL", toType = "ENTREZID",
                                                 OrgDb = "org.Hs.eg.db")[,2]

d3_hef_ck_CC <- readRDS("~/Documents/Lund/Vor2018_Haust2014/SingleCell/output_pathway_analysis/sign_pathway/d3_hef_ck_CC_logFC0.5.Rdata")
d3_hef_ck_CC <- compareCluster(geneClusters = d3_hef_cluster, OrgDb="org.Hs.eg.db", ont="CC", qvalueCutoff=0.05)
pdf(paste0("sign_pathway/d3_hef_ck_CC","_logFC", logFC,".pdf"), width = 10, height = 12, useDingbats=FALSE)
print(dotplot(d3_hef_ck_CC, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()
write.table(as.data.frame(d3_hef_ck_CC), file=paste0("sign_pathway/d3_hef_ck_CC","_logFC", logFC,".results.txt"), quote = F, row.names = T, sep="\t")
saveRDS(d3_hef_ck_CC, file=paste0("sign_pathway/d3_hef_ck_CC","_logFC", logFC,".Rdata"))
# 
# d3_hef_ck_MF <- compareCluster(geneClusters = d3_hef_cluster, OrgDb="org.Hs.eg.db", ont="MF", qvalueCutoff=0.05)
# pdf(paste0("sign_pathway/d3_hef_ck_MF","_logFC", logFC,".pdf"), width = 15, height = 12, useDingbats=FALSE)
# print(dotplot(d3_hef_ck_MF, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
# dev.off()
# write.table(as.data.frame(d3_hef_ck_MF), file=paste0("sign_pathway/d3_hef_ck_MF","_logFC", logFC,".results.txt"), quote = F, row.names = T, sep="\t")
# saveRDS(d3_hef_ck_MF, file=paste0("sign_pathway/d3_hef_ck_MF","_logFC", logFC,".Rdata"))

d3_hef_ck_BP <- readRDS("~/Documents/Lund/Vor2018_Haust2014/SingleCell/output_pathway_analysis/sign_pathway/Rdata/d3_hef_ck_BP_logFC0.5.Rdata")
d3_hef_ck_BP <- compareCluster(geneClusters = d3_hef_cluster, OrgDb="org.Hs.eg.db", ont="BP", qvalueCutoff=0.05)
pdf(paste0("sign_pathway/d3_hef_ck_BP","_logFC", logFC,".pdf"), width = 16, height = 12, useDingbats=FALSE)
print(dotplot(d3_hef_ck_BP, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()
write.table(as.data.frame(d3_hef_ck_BP), file=paste0("sign_pathway/d3_hef_ck_BP","_logFC", logFC,".results.txt"), quote = F, row.names = T, sep="\t")
saveRDS(d3_hef_ck_BP, file=paste0("sign_pathway/d3_hef_ck_BP","_logFC", logFC,".Rdata"))

#####
deHEF_idc9_logFC <- deHEF_idc9[abs(deHEF_idc9$logFC.HEF.idc_9) > logFC,]
d9_hefmatrix <- as.matrix(counts_scaled_sig_hef_id_upDown2[rownames(deHEF_idc9_logFC),
                                                         c(which(colnames(counts_scaled_sig_hef_id_upDown2) == "idc_9"),
                                                           which(colnames(counts_scaled_sig_hef_id_upDown2) == "HEF"))])
colnames(d9_hefmatrix) <- rownames(pData(HSMM))[c(which(HSMM$cell_type=="idc_9"), which(HSMM$cell_type=="HEF"))]
d9_hef_colname_dataframe <- as.data.frame(HSMM$cell_type[c(which(HSMM$cell_type == "idc_9"), which(HSMM$cell_type == "HEF"))], row.names = rownames(pData(HSMM))[c(which(HSMM$cell_type=="idc_9"), which(HSMM$cell_type=="HEF"))])
colnames(d9_hef_colname_dataframe) <- "Cell_type"

out_d9_hef_upDown <- pheatmap(d9_hefmatrix, col=mycol, show_colnames = F, cutree_rows = 2, annotation_names_row = F,
                            show_rownames = F, cluster_cols = F, annotation_col = d9_hef_colname_dataframe,
                            clustering_distance_rows = "correlation",
                            clustering_method = "ward.D2", annotation_names_col = F )
png(paste0("heatmaps/d9_hef_logFC",logFC, ".png"))
print(out_d9_hef_upDown)
dev.off()

d9_hef_cluster <- list()
d9_hef_cluster$Patways_up_day_9_vs_hef <- bitr(rownames(deHEF_idc9[deHEF_idc9$logFC.HEF.idc_9 > logFC, ]),
                                               fromType = "SYMBOL", toType = "ENTREZID",
                                               OrgDb = "org.Hs.eg.db")[,2]
d9_hef_cluster$Patways_down_day_9_vs_hef <- bitr(rownames(deHEF_idc9[deHEF_idc9$logFC.HEF.idc_9 < -logFC, ]),
                                                 fromType = "SYMBOL", toType = "ENTREZID",
                                                 OrgDb = "org.Hs.eg.db")[,2]

d9_hef_ck_CC <- readRDS("~/Documents/Lund/Vor2018_Haust2014/SingleCell/output_pathway_analysis/sign_pathway/Rdata/d9_hef_ck_CC_logFC0.5.Rdata")
d9_hef_ck_CC <- compareCluster(geneClusters = d9_hef_cluster, OrgDb="org.Hs.eg.db", ont="CC", qvalueCutoff=0.05)
pdf(paste0("sign_pathway/d9_hef_ck_CC","_logFC", logFC,".pdf"), width =9, height = 12, useDingbats=FALSE)
print(dotplot(d9_hef_ck_CC, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()
write.table(as.data.frame(d9_hef_ck_CC), file=paste0("sign_pathway/d9_hef_ck_CC","_logFC", logFC,".results.txt"), quote = F, row.names = T, sep="\t")
saveRDS(d9_hef_ck_CC, file=paste0("sign_pathway/d9_hef_ck_CC","_logFC", logFC,".Rdata"))

# d9_hef_ck_MF <- compareCluster(geneClusters = d9_hef_cluster, OrgDb="org.Hs.eg.db", ont="MF", qvalueCutoff=0.05)
# pdf(paste0("sign_pathway/d9_hef_ck_MF","_logFC", logFC,".pdf"), width = 15, height = 12, useDingbats=FALSE)
# print(dotplot(d9_hef_ck_MF, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
# dev.off()
# write.table(as.data.frame(d9_hef_ck_MF), file=paste0("sign_pathway/d9_hef_ck_MF","_logFC", logFC,".results.txt"), quote = F, row.names = T, sep="\t")
# saveRDS(d9_hef_ck_MF, file=paste0("sign_pathway/d9_hef_ck_MF","_logFC", logFC,".Rdata"))

d9_hef_ck_BP <- readRDS("~/Documents/Lund/Vor2018_Haust2014/SingleCell/output_pathway_analysis/sign_pathway/Rdata/d9_hef_ck_BP_logFC0.5.Rdata")
d9_hef_ck_BP <- compareCluster(geneClusters = d9_hef_cluster, OrgDb="org.Hs.eg.db", ont="BP", qvalueCutoff=0.05)
pdf(paste0("sign_pathway/d9_hef_ck_BP","_logFC", logFC,".pdf"), width = 13, height = 12, useDingbats=FALSE)
print(dotplot(d9_hef_ck_BP, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()
write.table(as.data.frame(d9_hef_ck_BP), file=paste0("sign_pathway/d9_hef_ck_BP","_logFC", logFC,".results.txt"), quote = F, row.names = T, sep="\t")
saveRDS(d9_hef_ck_BP, file=paste0("sign_pathway/d9_hef_ck_BP","_logFC", logFC,".Rdata"))

########
deidc3_idc9_logFC <- deidc3_idc9[abs(deidc3_idc9$logFC.idc_3.idc_9) > logFC,]
d3_day9_matrix <- as.matrix(counts_scaled_sig_hef_id_upDown2[rownames(deidc3_idc9_logFC),
                                                         c(which(colnames(counts_scaled_sig_hef_id_upDown2) == "idc_3"),
                                                           which(colnames(counts_scaled_sig_hef_id_upDown2) == "idc_9"))])
colnames(d3_day9_matrix) <- rownames(pData(HSMM))[c(which(HSMM$cell_type=="idc_3"), which(HSMM$cell_type=="idc_9"))]
d3_d9_colname_dataframe <- as.data.frame(HSMM$cell_type[c(which(HSMM$cell_type == "idc_3"), which(HSMM$cell_type == "idc_9"))], row.names = rownames(pData(HSMM))[c(which(HSMM$cell_type=="idc_3"), which(HSMM$cell_type=="idc_9"))])
colnames(d3_d9_colname_dataframe) <- "Cell_type"

out_d3_d9_upDown <- pheatmap(d3_day9_matrix, col=mycol, show_colnames = F, cutree_rows = 2, annotation_names_row = F,
                           show_rownames = F, cluster_cols = F, annotation_col = d3_d9_colname_dataframe,
                           clustering_distance_rows = "correlation",
                           clustering_method = "ward.D2", annotation_names_col = F )

png(paste0("heatmaps/d9_d3_logFC",logFC, ".png"))
print(out_d3_d9_upDown)
dev.off()

d3_d9_cluster <- list()
d3_d9_cluster$Patways_up_day_9_vs_day_3 <- bitr(rownames(deidc3_idc9[deidc3_idc9$logFC.idc_3.idc_9 > logFC, ]),
                                                fromType = "SYMBOL", toType = "ENTREZID",
                                                OrgDb = "org.Hs.eg.db")[,2]
d3_d9_cluster$Patways_down_day_9_vs_day_3 <- bitr(rownames(deidc3_idc9[deidc3_idc9$logFC.idc_3.idc_9 < -logFC, ]),
                                                  fromType = "SYMBOL", toType = "ENTREZID",
                                                  OrgDb = "org.Hs.eg.db")[,2]

d3_d9_ck_CC <- readRDS("~/Documents/Lund/Vor2018_Haust2014/SingleCell/output_pathway_analysis/sign_pathway/Rdata/d3_d9_ck_CC_logFC0.5.Rdata")
d3_d9_ck_CC <- compareCluster(geneClusters = d3_d9_cluster, OrgDb="org.Hs.eg.db", ont="CC", qvalueCutoff=0.05)
pdf(paste0("sign_pathway/d3_d9_ck_CC","_logFC", logFC,".pdf"), width = 9, height = 12, useDingbats=FALSE)
print(dotplot(d3_d9_ck_CC, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()
write.table(as.data.frame(d3_d9_ck_CC), file=paste0("sign_pathway/d3_d9_ck_CC","_logFC", logFC,".results.txt"), quote = F, row.names = T, sep="\t")
saveRDS(d3_d9_ck_CC, file=paste0("sign_pathway/d3_d9_ck_CC","_logFC", logFC,".Rdata"))

# d3_d9_ck_MF <- compareCluster(geneClusters = d3_d9_cluster, OrgDb="org.Hs.eg.db", ont="MF", qvalueCutoff=0.05)
# pdf(paste0("sign_pathway/d3_d9_ck_MF","_logFC", logFC,".pdf"), width = 15, height = 12, useDingbats=FALSE)
# print(dotplot(d3_d9_ck_MF, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
# dev.off()
# write.table(as.data.frame(d3_d9_ck_MF), file=paste0("sign_pathway/d3_d9_ck_MF","_logFC", logFC,".results.txt"), quote = F, row.names = T, sep="\t")
# saveRDS(d3_d9_ck_MF, file=paste0("sign_pathway/d3_d9_ck_MF","_logFC", logFC,".Rdata"))

d3_d9_ck_BP <- readRDS("~/Documents/Lund/Vor2018_Haust2014/SingleCell/output_pathway_analysis/sign_pathway/Rdata/d3_d9_ck_BP_logFC0.5.Rdata")
d3_d9_ck_BP <- compareCluster(geneClusters = d3_d9_cluster, OrgDb="org.Hs.eg.db", ont="BP", qvalueCutoff=0.05)
pdf(paste0("sign_pathway/d3_d9_ck_BP","_logFC", logFC,".pdf"), width = 15, height = 12, useDingbats=FALSE)
print(dotplot(d3_d9_ck_BP, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()
write.table(as.data.frame(d3_d9_ck_BP), file=paste0("sign_pathway/d3_d9_ck_BP","_logFC", logFC,".results.txt"), quote = F, row.names = T, sep="\t")
saveRDS(d3_d9_ck_BP, file=paste0("sign_pathway/d3_d9_ck_BP","_logFC", logFC,".Rdata"))


#####
ded9_cDC1_logFC <- ded9_cDC1[abs(ded9_cDC1$logFC.cDC1.idc_9) > logFC,]
d9_cDC1_matrix <- as.matrix(counts_scaled_sig_hef_id_upDown2[rownames(ded9_cDC1_logFC),
                                                         c(which(colnames(counts_scaled_sig_hef_id_upDown2) == "idc_9"),
                                                           which(colnames(counts_scaled_sig_hef_id_upDown2) == "cDC1"))])
colnames(d9_cDC1_matrix) <- rownames(pData(HSMM))[c(which(HSMM$cell_type=="idc_9"), which(HSMM$cell_type=="cDC1"))]
d9_cDC1_colname_dataframe <- as.data.frame(HSMM$cell_type[c(which(HSMM$cell_type == "idc_9"), which(HSMM$cell_type == "cDC1"))], row.names = rownames(pData(HSMM))[c(which(HSMM$cell_type=="idc_9"), which(HSMM$cell_type=="cDC1"))])
colnames(d9_cDC1_colname_dataframe) <- "Cell_type"

out_d9_cDC1_upDown <- pheatmap(d9_cDC1_matrix, col=mycol, show_colnames = F, cutree_rows = 2, annotation_names_row = F,
                           show_rownames = F, cluster_cols = F, annotation_col = d9_cDC1_colname_dataframe,
                           clustering_distance_rows = "correlation",
                           clustering_method = "ward.D2", annotation_names_col = F )

png(paste0("heatmaps/d9_cDC1_logFC",logFC, ".png"))
print(out_d9_cDC1_upDown)
dev.off()
d9_cDC1_cluster <- list()
d9_cDC1_cluster$Patways_up_day_9_vs_cDC1 <- bitr(rownames(ded9_cDC1[ded9_cDC1$logFC.cDC1.idc_9 > logFC, ]),
                                                 fromType = "SYMBOL", toType = "ENTREZID",
                                                 OrgDb = "org.Hs.eg.db")[,2]
d9_cDC1_cluster$Patways_down_day_9_vs_cDC1 <- bitr(rownames(ded9_cDC1[ded9_cDC1$logFC.cDC1.idc_9 < -logFC, ]),
                                                   fromType = "SYMBOL", toType = "ENTREZID",
                                                   OrgDb = "org.Hs.eg.db")[,2]

d9_cDC1_ck_CC <- readRDS("~/Documents/Lund/Vor2018_Haust2014/SingleCell/output_pathway_analysis/sign_pathway/Rdata/d9_cDC1_ck_CC_logFC0.5.Rdata")
d9_cDC1_ck_CC <- compareCluster(geneClusters = d9_cDC1_cluster, OrgDb="org.Hs.eg.db", ont="CC", qvalueCutoff=0.05)
pdf(paste0("sign_pathway/d9_cDC1_ck_CC","_logFC", logFC,".pdf"), width = 9, height = 12, useDingbats=FALSE)
print(dotplot(d9_cDC1_ck_CC, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()
write.table(as.data.frame(d9_cDC1_ck_CC), file=paste0("sign_pathway/d9_cDC1_ck_CC","_logFC", logFC,".results.txt"), quote = F, row.names = T, sep="\t")
saveRDS(d9_cDC1_ck_CC, file=paste0("sign_pathway/d9_cDC1_ck_CC","_logFC", logFC,".Rdata"))

# d9_cDC1_ck_MF <- compareCluster(geneClusters = d9_cDC1_cluster, OrgDb="org.Hs.eg.db", ont="MF", qvalueCutoff=0.05)
# pdf(paste0("sign_pathway/d9_cDC1_ck_MF","_logFC", logFC,".pdf"), width = 15, height = 12, useDingbats=FALSE)
# print(dotplot(d9_cDC1_ck_MF, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
# dev.off()
# write.table(as.data.frame(d9_cDC1_ck_MF), file=paste0("sign_pathway/d9_cDC1_ck_MF","_logFC", logFC,".results.txt"), quote = F, row.names = T, sep="\t")
# saveRDS(d9_cDC1_ck_MF, file=paste0("sign_pathway/d9_cDC1_ck_MF","_logFC", logFC,".Rdata"))

d9_cDC1_ck_BP <- readRDS("~/Documents/Lund/Vor2018_Haust2014/SingleCell/output_pathway_analysis/sign_pathway/Rdata/d9_cDC1_ck_BP_logFC0.5.Rdata")
d9_cDC1_ck_BP <- compareCluster(geneClusters = d9_cDC1_cluster, OrgDb="org.Hs.eg.db", ont="ALL", qvalueCutoff=0.05, pool=T)
pdf(paste0("sign_pathway/d9_cDC1_ck_BP","_logFC", logFC,"_vol2.pdf"), width = 14, height = 12, useDingbats=FALSE)
print(dotplot(d9_cDC1_ck_BP, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()
write.table(as.data.frame(d9_cDC1_ck_BP), file=paste0("sign_pathway/d9_cDC1_ck_BP","_logFC", logFC,".results.txt"), quote = F, row.names = T, sep="\t")
saveRDS(d9_cDC1_ck_BP, file=paste0("sign_pathway/d9_cDC1_ck_BP","_logFC", logFC,".Rdata"))

######

d3_mef_ck_kegg_logFC <- readRDS("~/Documents/Lund/Vor2018_Haust2014/SingleCell/output_pathway_analysis/sign_pathway/Rdata/d3_hef_logFC0.5_KEGG.Rdata")
d3_mef_ck_kegg_logFC <- compareCluster(geneClusters = d3_hef_cluster, fun="enrichKEGG", organism="hsa", keyType="kegg", qvalueCutoff=0.05)
pdf(paste0("sign_pathway/d3_hef_logFC", logFC, "_KEGG.pdf"), width = 11, height = 14, useDingbats=FALSE)
print(dotplot(d3_mef_ck_kegg_logFC, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()
write.table(as.data.frame(d3_mef_ck_kegg_logFC), file=paste0("sign_pathway/d3_hef_logFC", logFC, "_KEGG_results.txt"), quote = F, row.names = T, sep="\t")
saveRDS(d3_mef_ck_kegg_logFC, file=paste0("sign_pathway/d3_hef_logFC", logFC, "_","KEGG.Rdata"))

d9_mef_ck_kegg_logFC <- readRDS("~/Documents/Lund/Vor2018_Haust2014/SingleCell/output_pathway_analysis/sign_pathway/Rdata/d9_hef_logFC0.5_KEGG.Rdata")
d9_mef_ck_kegg_logFC <- compareCluster(geneClusters = d9_hef_cluster, fun="enrichKEGG", organism="hsa", keyType="kegg", qvalueCutoff=0.05)
pdf(paste0("sign_pathway/d9_hef_logFC", logFC, "_KEGG.pdf"), width = 9, height = 14, useDingbats=FALSE)
print(dotplot(d9_mef_ck_kegg_logFC, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()
write.table(as.data.frame(d9_mef_ck_kegg_logFC), file=paste0("sign_pathway/d9_hef_logFC", logFC, "_KEGG_results.txt"), quote = F, row.names = T, sep="\t")
saveRDS(d9_mef_ck_kegg_logFC, file=paste0("sign_pathway/d9_hef_logFC", logFC, "_","KEGG.Rdata"))

d3_d9_ck_kegg_logFC <- readRDS("~/Documents/Lund/Vor2018_Haust2014/SingleCell/output_pathway_analysis/sign_pathway/Rdata/d3_d9_logFC0.5_KEGG.Rdata")
d3_d9_ck_kegg_logFC <- compareCluster(geneClusters = d3_d9_cluster, fun="enrichKEGG", organism="hsa", keyType="kegg", qvalueCutoff=0.05)
pdf(paste0("sign_pathway/d3_d9_logFC", logFC, "_KEGG.pdf"), width = 9, height = 14, useDingbats=FALSE)
print(dotplot(d3_d9_ck_kegg_logFC, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()
write.table(as.data.frame(d3_d9_ck_kegg_logFC), file=paste0("sign_pathway/d3_d9_logFC", logFC, "_KEGG_results.txt"), quote = F, row.names = T, sep="\t")
saveRDS(d3_d9_ck_kegg_logFC, file=paste0("sign_pathway/d3_d9_logFC", logFC, "_","KEGG.Rdata"))

d9_cDC1_ck_kegg_logFC <- readRDS("~/Documents/Lund/Vor2018_Haust2014/SingleCell/output_pathway_analysis/sign_pathway/Rdata/d9_cDC1_logFC0.5_KEGG.Rdata")
d9_cDC1_ck_kegg_logFC <- compareCluster(geneClusters = d9_cDC1_cluster, fun="enrichKEGG", organism="hsa", keyType="kegg", qvalueCutoff=0.05)
pdf(paste0("sign_pathway/d9_cDC1_logFC", logFC, "_KEGG.pdf"), width = 8, height = 14, useDingbats=FALSE)
print(dotplot(d9_cDC1_ck_kegg_logFC, showCategory=5, font.size =24) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()
write.table(as.data.frame(d9_cDC1_ck_kegg_logFC), file=paste0("sign_pathway/d9_cDC1_logFC", logFC, "_KEGG_results.txt"), quote = F, row.names = T, sep="\t")
saveRDS(d9_cDC1_ck_kegg_logFC, file=paste0("sign_pathway/d9_cDC1_logFC", logFC, "_","KEGG.Rdata"))


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


