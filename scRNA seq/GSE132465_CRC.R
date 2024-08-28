library(Seurat)
library(SeuratData)
library(SeuratObject)
library(ggplot2)
library(patchwork)
library(dplyr)
####Step 1: Preprocessing + Epithelial cells annotation####
path <- paste0("data/CRC_GSE132465/","matrix",".txt")
ref <- read.table(file=path,head=TRUE,sep="\t")
a<-ref$Index
rownames(ref)<-a
ref$Index<-NULL
ref <- CreateSeuratObject(counts = ref, project = "CRC_GSE144735" , 
                          min.cells = 3, min.features = 200)
ref[["percent.mt"]] <- PercentageFeatureSet(ref, pattern = "^MT-")
VlnPlot(ref, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
library(readr)
path <- paste0("data/CRC_GSE132465/","annotation",".txt")
metadata <- read.table(file=path,head=TRUE,sep="\t")
metadata_cell_type_major<-metadata$Cell_type
metadata_cell_type_minor<-metadata$Cell_subtype
metadata_cell_patient<-metadata$Patient
metadata_cell_classr<-metadata$Class
metadata_cell_sample<-metadata$Sample
ref<-AddMetaData(ref,metadata_cell_type_major,col.name = "celltype_major")
ref<-AddMetaData(ref,metadata_cell_type_minor,col.name = "celltype_minor")
ref<-AddMetaData(ref,metadata_cell_patient,col.name = "patient")
ref<-AddMetaData(ref,metadata_cell_classr,col.name = "class")
ref<-AddMetaData(ref,metadata_cell_sample,col.name = "sample")
Idents(ref)<-"celltype_major"
ref<-subset(ref,idents = c("Epithelial cells"))
Idents(ref)<-"class"
ref<-subset(ref,idents=c("Tumor"))
saveRDS(ref,"data/CRC_GSE132465/ref.rds")
ref <- NormalizeData(ref,normalization.method = "LogNormalize", scale.factor = 1e4)
ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(ref)
ref <- ScaleData(ref, features = all.genes)
ref <- RunPCA(ref, features = VariableFeatures(object = ref))
ref <- FindNeighbors(ref, dims = 1:20)
ref <- FindClusters(ref, resolution = 1.5)
ref <- RunUMAP(ref, dims = 1:20)
Idents(ref)<-"seurat_clusters"
DimPlot(ref, reduction = "umap",label=T)
Idents(ref)<-"patient"
DimPlot(ref, reduction = "umap")
Idents(ref)<-"celltype_minor"
DimPlot(ref, reduction = "umap")
Idents(ref)<-"class"
DimPlot(ref, reduction = "umap")
Idents(ref)<-"sample"
DimPlot(ref, reduction = "umap")
Idents(ref)<-"seurat_clusters"
#cluster in previous research
gene<-c("IFI6","TM4SF4","CA2","CA7","PYY","CHGA","ITLN1","MUC2","MKI67","PTTG1","PLCG2","SH2D6")
gene<-c("KRT8","KRT18","AGR2","TFF3")
DotPlot(object = ref,
        features = gene)+RotatedAxis()
ggsave(file=paste0("plot/CRC_GSE144735/","dot_plot_epi_MAIN",".pdf"),width=15,height=10)
Idents(ref)<-"seurat_clusters"
DoHeatmap(ref , features = gene, size = 3)
ggsave(file=paste0("plot/CRC_GSE144735/","heatmap_epi_MAIN.pdf"),height=10,width=25,dpi=600)
FeaturePlot(object = ref,label.size = 2,
            features = gene)
ggsave(file=paste0("plot/CRC_GSE144735/","feature_plot_epi_MAIN_2",".pdf"),width=15,height=15)
#epi subtype
gene<-c("IFI6","TM4SF4","SLC26A3","CA2","PIGR","PHGR1","BEST4","CA7",
        "PYY","CHGA","ITLN1","MUC2","RND3","DNAJB1","SCD","FDPS","SLC2A1","GPRC5A","S100A11",
        "IL32","UBE2C","PTTG1","PCNA","MKI67","ASCL2","PTPRO","PLCG2","SH2D6")
gene<-c("ASCL2","PTPRO")
DotPlot(object = ref,
        features = gene)+RotatedAxis()
ggsave(file=paste0("plot/CRC_GSE132465/","dot_plot_epi",".pdf"),width=15,height=10)
Idents(ref)<-"seurat_clusters"
DoHeatmap(ref , features = gene, size = 3)
ggsave(file=paste0("plot/CRC_GSE132465/","heatmap_epi.pdf"),height=10,width=25,dpi=600)
FeaturePlot(object = ref,label.size = 2,
            features = gene)
ggsave(file=paste0("plot/CRC_GSE132465/","feature_plot_epi",".pdf"),width=15,height=20)
#cancer stem cell annotation by marker expression in all samples
Idents(ref)<-"seurat_clusters"
new.cluster.ids <- c("other","other","cancer stem cell","other","other","other",
                     "other","other","other","other","other",
                     "cancer stem cell","other","other","other","other",
                     "cancer stem cell","other","other","other","other",
                     "other","other","other","other","other",
                     "other","other","other","other","cancer stem cell",
                     "other")
names(new.cluster.ids) <- levels(ref)
ref <- RenameIdents(ref, new.cluster.ids)
DimPlot(ref, reduction = "umap",label=F)
ref$csc_origin <- ref@active.ident
table(ref$celltype_minor,ref$csc_origin)
table(ref$patient,ref$csc_origin)
####Step2: CSC annotation####
sample<-c("SMC02","SMC03","SMC07","SMC08","SMC09","SMC18","SMC22")
id<-"SMC22"
Idents(ref)<-"patient"
ref_sub<-subset(ref,idents=id)
ref_sub <- NormalizeData(ref_sub,normalization.method = "LogNormalize", scale.factor = 1e4)
ref_sub <- FindVariableFeatures(ref_sub, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(ref_sub)
ref_sub <- ScaleData(ref_sub, features = all.genes)
ref_sub <- RunPCA(ref_sub, features = VariableFeatures(object = ref_sub))
ref_sub <- FindNeighbors(ref_sub, dims = 1:20)
ref_sub <- FindClusters(ref_sub, resolution = 1)
ref_sub <- RunUMAP(ref_sub, dims = 1:20)
DimPlot(ref_sub, reduction = "umap",label=T)
ref_sub<-ref_list[["SMC22"]]
Idents(ref_sub)<-"seurat_clusters"
library(readr)
colorectal_cancer_stem_cell_markers <- read_csv("data/colorectal cancer stem cell markers.csv")
marker<-list(colorectal_cancer_stem_cell_markers$combined)
name <-"csc_combined"
ref_sub<-AddModuleScore(ref_sub,features =marker,name =name)
ggplot(ref_sub@meta.data,aes(x=seurat_clusters,y=csc_combined1,fill=seurat_clusters))+
  geom_boxplot()+labs(title=name)
ref.markers <- FindAllMarkers(ref_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
df<-ref.markers[ref.markers$p_val_adj<0.05,]
df<-df[df$gene%in%marker[[1]],]
print(df)
table(df$cluster)
for(i in sort(unique(ref_sub$seurat_clusters))){
  print(paste0("cluster",i,":"))
  print(df[df$cluster==i,"gene"])
}
#P2
new.cluster.ids <- c("other","cancer stem cell","other","other","other")
#P3
new.cluster.ids <- c("cancer stem cell","other","cancer stem cell","other","other","other")
#P7
new.cluster.ids <- c("cancer stem cell","other","cancer stem cell","other","other","other","other")
#P8
new.cluster.ids <- c("cancer stem cell","other","cancer stem cell","other")
#P9
new.cluster.ids <- c("cancer stem cell","other","cancer stem cell","other","other","other","other","other","other","other")
#P18
new.cluster.ids <- c("other","cancer stem cell","other","other","other","cancer stem cell","other","other","other")
#P22
new.cluster.ids <- c("cancer stem cell","other","other","other","other","other","other","other","other","other")
Idents(ref_sub)<-"seurat_clusters"
names(new.cluster.ids) <- levels(ref_sub)
ref_sub<- RenameIdents(ref_sub, new.cluster.ids)
ref_sub$csc<- ref_sub@active.ident
table(ref_sub$csc)
Idents(ref_sub)<-"csc"
DimPlot(ref_sub, reduction = "umap",label=F)
ref_list[[id]]<-ref_sub

Idents(ref_sub)<-"seurat_clusters"
gene<-c("ASCL2","PTPRO")
DotPlot(object = ref_sub,
        features = gene)+RotatedAxis()
#P22
new.cluster.ids <- c("other","cancer stem cell","cancer stem cell","other","cancer stem cell","cancer stem cell","other","cancer stem cell","cancer stem cell","other")
Idents(ref)<-"seurat_clusters"
names(new.cluster.ids) <- levels(ref)
ref<- RenameIdents(ref, new.cluster.ids)
ref$csc_marker_individual<- ref@active.ident
table(ref$csc_marker_individual)
Idents(ref)<-"csc_marker_individual"
DimPlot(ref, reduction = "umap",label=F)
ref_list[[id]]<-ref
saveRDS(ref_list,"data/CRC_GSE132465/ref_list.rds")

####Step3: genertate tsv file####
sample<-c("SMC02","SMC03","SMC07","SMC08","SMC09","SMC18")
for(id in sample){
  ref<-ref_list[[id]]
  # Specify the filename for the TSV file
  filename <- paste0("data/Compass/CRC_GSE132465/",id,".tsv")
  # Write the data frame to a TSV file
  data<-ref@assays$RNA@layers$counts
  rownames(data)<-rownames(ref)[ref@assays$RNA@features@.Data[,"counts"]]
  colnames(data)<-colnames(ref)[ref@assays$RNA@cells@.Data[,"counts"]]
  write.table(data, file = filename, sep = "\t", quote = FALSE, row.names = T)
  filename <- paste0("data/Compass/CRC_GSE132465/",id,"_cluster_meta.csv")
  data<-as.data.frame(ref@meta.data[,c("csc_origin","csc")])
  rownames(data)<-rownames(ref@meta.data)
  colnames(data)<-c("csc_marker_all","csc_addmodulescore")
  write.csv(data, file = filename)
}
id<-"SMC22"
ref<-ref_list[[id]]
# Specify the filename for the TSV file
filename <- paste0("data/Compass/CRC_GSE132465/",id,".tsv")
# Write the data frame to a TSV file
data<-ref@assays$RNA@layers$counts
rownames(data)<-rownames(ref)[ref@assays$RNA@features@.Data[,"counts"]]
colnames(data)<-colnames(ref)[ref@assays$RNA@cells@.Data[,"counts"]]
write.table(data, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
filename <- paste0("data/Compass/CRC_GSE132465/",id,"_cell_meta.csv")
data<-ref@meta.data[,c("csc_origin","csc","csc_marker_individual")]
rownames(data)<-rownames(ref@meta.data)
colnames(data)<-c("csc_marker_all","csc_addmodulescore","csc_marker_individual")
write.csv(data, file = filename)

####Step4: cell cycle score####
library(Seurat)
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
id<-names(ref_list)[1]
ref_sub<-ref_list[[id]]
library(readr)
colorectal_cancer_stem_cell_markers <- read_csv("data/colorectal cancer stem cell markers.csv")
marker<-list(colorectal_cancer_stem_cell_markers$combined)
name <-"csc_combined"
ref_sub<-AddModuleScore(ref_sub,features =marker,name =name)
DimPlot(ref_sub, reduction = "umap",label=T)
ref_sub <- CellCycleScoring(ref_sub, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
cellcycle <- data.frame(G2M=ref_sub$G2M.Score,S = ref_sub$S.Score,seurat_clusters = ref_sub$seurat_clusters)
rownames(cellcycle)<-rownames(ref_sub@meta.data)
library(ggplot2)
P1<-ggplot(cellcycle,aes(x=seurat_clusters,y=G2M,fill=seurat_clusters))+
  geom_boxplot()+labs(title=paste0(id,"_G2M"))
P2<-ggplot(cellcycle,aes(x=seurat_clusters,y=S,fill=seurat_clusters))+
  geom_boxplot()+labs(title=paste0(id,"_S"))
P1+P2
#correlation with csc score
for(id in names(ref_list)){
  ref_sub<-ref_list[[id]]
  colorectal_cancer_stem_cell_markers <- read_csv("data/colorectal cancer stem cell markers.csv")
  marker<-list(colorectal_cancer_stem_cell_markers$combined)
  name <-"csc_combined"
  ref_sub<-AddModuleScore(ref_sub,features =marker,name =name)
  ref_sub <- CellCycleScoring(ref_sub, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  pearson_corr <- cor(ref_sub$csc_combined1, ref_sub$S.Score, method = "pearson")
  spearman_corr <- cor(ref_sub$csc_combined1, ref_sub$S.Score, method = "spearman")
  kendall_corr <- cor(ref_sub$csc_combined1, ref_sub$S.Score, method = "kendall")
  s_score<-c(pearson_corr,spearman_corr,kendall_corr)
  
  pearson_corr <- cor(ref_sub$csc_combined1, ref_sub$G2M.Score, method = "pearson")
  spearman_corr <- cor(ref_sub$csc_combined1, ref_sub$G2M.Score, method = "spearman")
  kendall_corr <- cor(ref_sub$csc_combined1, ref_sub$G2M.Score, method = "kendall")
  g2m_score<-c(pearson_corr,spearman_corr,kendall_corr)
  df<-data.frame("S.Score"=s_score,"G2M.Score"=g2m_score)
  rownames(df)<-c("pearson_corr","spearman_corr","kendall_corr")
  print(df)
  print(id)
}
