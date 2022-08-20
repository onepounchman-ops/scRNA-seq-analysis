####  Code Description              ####
#---  1. Written by WoLin @ 2019.03.15，last update 19.07.14 ---#
#---  2. Analysis for single sample    ---#
#---  3. Support 10X data & expression matrix ---#
#---  4. Need to change:sample data, sam.name, dims ---#

#### 1. 加载分析使用的工具包 ####
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)

E12.5 <- Read10X_h5(file.path(dataset_loc, ids[1], "GSM4504959_E12.5.h5"), use.names = T)
E15.5 <- Read10X_h5(file.path(dataset_loc, ids[1], "GSM4504960_E15.5.h5"), use.names = T)
E17.5 <- Read10X_h5(file.path(dataset_loc, ids[1], "GSM4504961_E17.5.h5"), use.names = T)
colnames(E12.5) <- paste(colnames(E12.5),"E12.5",sep = "_")
colnames(E15.5) <- paste(colnames(E15.5),"E15.5",sep = "_")
colnames(E17.5) <- paste(colnames(E17.5),"E17.5",sep = "_")

experiment.data <- cbind(E12.5,E15.5)
experiment.data <- cbind(experiment.data,E17.5)

sam.name <- "multi"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

experiment.aggregate <- CreateSeuratObject(
  experiment.data,
  project = "multi", 
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "_")
save(experiment.aggregate,file=paste0("./",sam.name,"/",sam.name,"_raw_SeuratObject.RData"))

slotNames(experiment.aggregate)
experiment.aggregate@assays
dim(experiment.aggregate@meta.data)
View(experiment.aggregate@meta.data)

experiment.aggregate.matrix <- as.matrix(experiment.aggregate@assays$RNA@counts)

experiment.aggregate[["percent.mt"]] <- PercentageFeatureSet(experiment.aggregate, 
                                                             pattern = "^mt-")
pdf(paste0("./",sam.name,"/QC-VlnPlot.pdf"),width = 8,height = 4.5)
VlnPlot(experiment.aggregate, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)
dev.off()

gene.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$nFeature_RNA,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
rna.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$nCount_RNA,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
mt.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$percent.mt,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
freq.combine <- as.data.frame(cbind(gene.freq,rna.freq,mt.freq))
colnames(freq.combine) <- c(paste(colnames(gene.freq),"Gene",sep = "_"),
                            paste(colnames(rna.freq),"RNA",sep = "_"),
                            paste(colnames(mt.freq),"mt",sep = "_"))
write.table(freq.combine,file = paste0(sam.name,"/QC-gene_frequency.txt"),quote = F,sep = "\t")
rm(gene.freq,rna.freq,mt.freq)
View(freq.combine)

plot1 <- FeatureScatter(experiment.aggregate, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(experiment.aggregate, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf(paste0("./",sam.name,"/QC-FeatureScatter.pdf"),width = 8,height = 4.5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()
rm(plot1,plot2)

cat("Before filter :",nrow(experiment.aggregate@meta.data),"cells\n")
experiment.aggregate <- subset(experiment.aggregate, 
                               subset = 
                                 nFeature_RNA > 500 & 
                                 nCount_RNA > 1000 & 
                                 nCount_RNA < 20000 &
                                 percent.mt < 5)
cat("After filter :",nrow(experiment.aggregate@meta.data),"cells\n")

experiment.aggregate <- NormalizeData(experiment.aggregate, 
                                      normalization.method = "LogNormalize",
                                      scale.factor = 10000)

experiment.aggregate <- FindVariableFeatures(experiment.aggregate, 
                                             selection.method = "vst",
                                             nfeatures = 1000)

top10 <- head(x = VariableFeatures(experiment.aggregate), 10)
plot1 <- VariableFeaturePlot(experiment.aggregate)
plot2 <- LabelPoints(plot = plot1, points = top10)
pdf(file = paste0(sam.name,"/Norm-feature_variable_plot.pdf"),width = 8,height = 5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()

experiment.aggregate <- ScaleData(
  object = experiment.aggregate,
  do.scale = FALSE,
  do.center = FALSE,
  vars.to.regress = c("orig.ident","percent.mt"))

experiment.aggregate <- RunPCA(object = experiment.aggregate, 
                               features = VariableFeatures(experiment.aggregate),
                               verbose = F,npcs = 50)

pdf(paste0("./",sam.name,"/PCA-VizDimLoadings.pdf"),width = 7,height = 5)
VizDimLoadings(experiment.aggregate, dims = 1:2, reduction = "pca")
dev.off()

pdf(paste0("./",sam.name,"/PCA-DimPlot.pdf"),width = 5,height = 4)
DimPlot(experiment.aggregate, reduction = "pca")
dev.off()

pdf(paste0("./",sam.name,"/PCA-DimHeatmap.pdf"))
DimHeatmap(experiment.aggregate, dims = 1:6, cells = 500, balanced = TRUE)
dev.off()

experiment.aggregate <- JackStraw(experiment.aggregate, num.replicate = 100,dims = 40)
experiment.aggregate <- ScoreJackStraw(experiment.aggregate, dims = 1:40)
pdf(paste0("./",sam.name,"/PCA-JackStrawPlot_40.pdf"),width = 6,height = 5)
JackStrawPlot(object = experiment.aggregate, dims = 1:40)
dev.off()

pdf(paste0("./",sam.name,"/PCA-ElbowPlot.pdf"),width = 8,height = 5)
ElbowPlot(experiment.aggregate,ndims = 40)
dev.off()

dim.use <- 1:14

experiment.aggregate <- FindNeighbors(experiment.aggregate, dims = dim.use)
experiment.aggregate <- FindClusters(experiment.aggregate, resolution = 0.5)

experiment.aggregate <- RunTSNE(experiment.aggregate, dims = dim.use, 
                                do.fast = TRUE)
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_res0.5_",max(dim.use),"PC.pdf"),width = 8.5,height = 7)
DimPlot(object = experiment.aggregate, pt.size=0.5,label = T) 
dev.off()

pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_SamGroup_",max(dim.use),"PC.pdf"),width = 8,height = 7)
DimPlot(object = experiment.aggregate, 
        group.by="orig.ident", 
        pt.size=0.5,reduction = "tsne")
dev.off()

pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_SamGroup_slipt_",max(dim.use),"PC.pdf"),width = 10,height = 4)
DimPlot(object = experiment.aggregate, 
        split.by ="orig.ident", 
        pt.size=0.5,reduction = "tsne")
dev.off()

table(experiment.aggregate@meta.data$orig.ident)

all.markers <- FindAllMarkers(experiment.aggregate, only.pos = TRUE, 
                              min.pct = 0.3, logfc.threshold = 0.25)
write.table(all.markers,
            file=paste0("./",sam.name,"/",sam.name,"_total_marker_genes_tsne_",max(dim.use),"PC.txt"),
            sep="\t",quote = F,row.names = F)

marker.sig <- all.markers %>% 
  mutate(Ratio = round(pct.1/pct.2,3)) %>%
  filter(p_val_adj <= 0.05)  # 本条件为过滤统计学不显著的基因

for(cluster_id in unique(marker.sig$cluster)){
  # cluster.markers <- FindMarkers(experiment.aggregate, ident.1 = cluster, min.pct = 0.3)
  # cluster.markers <- as.data.frame(cluster.markers) %>% 
  #   mutate(Gene = rownames(cluster.markers))
  cl4.genes <- marker.sig %>% 
    filter(cluster == cluster_id) %>%
    arrange(desc(avg_logFC))
  cl4.genes <- cl4.genes[1:min(nrow(cl4.genes),4),"gene"]
  Selective.marker.genes<-c("Epcam","Ager","Hopx","Sftpb","Sftpc","S100a8",
                            "S100a9","Hba-a2","Hbb-bt","Tgfbi","Actg2",
                            "Acta2","Myh11","Tagln","Col1a1","Col13a1",
                            "Col14a1","Scgb3a2","Scgb1a1","Cyp2f2",
                            "Upk3b","Wt1","Ube2c","Hmmr","Cldn5",
                            "Emcn","Pecam1","Cd68","C1qa","Pdgfrb",
                            "Higd1b","Gucy1a3","Cox4i2","Postn",
                            "Pdzd2","Cspg4","Ptn","Lum","Prrx1")
  Epithelium.marker.genes<-c("Epcam","Ager","Hopx","Sftpb","Sftpc")  
}

experiment.merged <- RenameIdents(
  object = experiment.aggregate,
  "0" = "Mesenchymal Progenitors",
  "1" = "Epithelial Progenitors",
  "2" = "Matrix Fibroblasts",
  "3" = "Intermediate Fibroblasts",
  "4" = "Endothelial Cells",
  "5" = "Mesenchymal Progenitors",
  "6" = "AT1 Cells",
  "7" = "Neturophils",
  "8" = "AT2 Cells",
  "9" = "Matrix Fibroblasts",
  "10" = "Myofibroblasts/Smooth Muscel Cells",
  "11" = "Endothelial Cells",
  "12" = "Erythrocyte",
  "13" = "Epithelial Progenitors",
  "14" = "Myofibroblasts/Smooth Muscel Cells",
  "15" = "Pericyte1",
  "16" = "Endothelial Cells",
  "17" = "Mesothelial Cells",
  "18" = "Pericyte2",
  "19" = "Macrophage/monocyte",
  "20" = "Myofibroblasts/Smooth Muscel Cells",
  "21" = "Mesothelial Cells",
  "22" = "Erythrocyte",
  "23" = "Club Cells",
  "24" = "Matrix Fibroblasts"
)
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_Rename_",max(dim.use),"PC.pdf"),width = 10,height = 7)
DimPlot(object = experiment.merged, pt.size=0.5,
        reduction = "tsne",label = F) +
  ggsci::scale_color_igv()
dev.off()

pdf(paste0("./",sam.name,"/MarkerGene-DotPlot_all_cluster_tsne_",max(dim.use),"PC.pdf"),width = 50,height = 5)
DotPlot(experiment.aggregate, features = Selective.marker.genes)+
  RotatedAxis()
dev.off()

meta_data <- experiment.aggregate@meta.data 
plot_data <- data.frame(table(meta_data$orig.ident,meta_data$seurat_clusters))
plot_data$Total <- apply(plot_data,1,function(x)sum(plot_data[plot_data$Var1 == x[1],3]))
plot_data <- plot_data %>% mutate(Percentage = round(Freq/Total,3) * 100)
ggplot(plot_data,aes(x = Var1,y = Percentage,fill = Var2)) +
  geom_bar(stat = "identity",position = "stack") +
  theme_classic() + 
  theme(axis.title.x = element_blank()) + labs(fill = "Cluster")
