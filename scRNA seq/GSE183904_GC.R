library(Seurat)
library(SeuratData)
library(SeuratObject)
library(ggplot2)
library(patchwork)
library(dplyr)
library(readr)
####Step 1: preprocessing +Epithelial cells annotation####
ref_index_list<-c(1:40)
ref_list <- sapply(ref_index_list,function(id){
  path <- paste0("data/GC_GSE183904/GSE183904_RAW/","sample",id,".csv")
  ref <- read_csv(path)
  a<-ref$...1
  ref$...1<-NULL
  rownames(ref)<-a
  ref <- CreateSeuratObject(counts = ref, project = paste0("sample",id) , min.cells = 3, min.features = 500)
  ref[["percent.mt"]] <- PercentageFeatureSet(ref, pattern = "^MT-")
  ref <- subset(ref, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20)
  return(ref)
})
saveRDS(ref_list,"data/GC_GSE183904/ref_list.rds")#save

#since there are 30 samples, we classifed 10-20 of them as a group(1-15,15-28,29-40)
#normalization
library(harmony)
ref_index_list<-c(1:40)
normal <- c(6,21,25,4,23,1,9,11,35,37)
# Remove specified numbers
ref_index_list <- ref_index_list[!ref_index_list%in%normal]
ref <- merge(ref_list[[ref_index_list[21]]], y = ref_list[ref_index_list[c(22:30)]],add.cell.ids = ref_index_list[c(1:10)], merge.data = TRUE)
unique(sapply(X = strsplit(colnames(ref), split = "_"), FUN = "[", 1))
#saveRDS(ref,file="data/GC_GSE183904/ref.rds")
#remove(ref_list)
ref <- NormalizeData(ref,normalization.method = "LogNormalize", scale.factor = 1e4)
#HVG
ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(ref)
#saveRDS(ref,file="data/BC_GSE180286/ref.rds")
ref <- ScaleData(ref, features = all.genes)
#saveRDS(ref,file="data/BC_GSE180286/ref.rds")
ref <- RunPCA(ref, features = VariableFeatures(object = ref))
ref <- RunHarmony(ref, group.by.vars = "orig.ident")
ref <- IntegrateLayers(
  object = ref, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
#saveRDS(ref,file="data/BC_GSE180286/ref.rds")
ref <- FindNeighbors(ref, dims = 1:20,reduction = "harmony")
ref <- FindClusters(ref, resolution = 0.8,reduction = "harmony")
ref <- RunUMAP(ref, dims = 1:20,reduction = "harmony")
Idents(ref)<-"seurat_clusters"
DimPlot(ref, reduction = "umap",label=T)
Idents(ref)<-"orig.ident"
DimPlot(ref, reduction = "umap",label=F)
Idents(ref)<-"seurat_clusters"
gene1<-c("CDH1","MUC5AC","TFF1","LIPF","PGA3","REG4","FN1",
         "RGS5","NOTCH3","LUM","DCN","PLVAP","ACKR1")
gene2<-c("CD8A","IL2RA","KLRD1","MS4A1","TNFRSF17","KIT","PLD4","CD163")
FeaturePlot(object = ref,label.size = 2,
            features = gene1)
ggsave(file=paste0("plot/GC_GSE183904/","feature_plot_gene1_3",".pdf"),width=15,height=15)
FeaturePlot(object = ref,label.size = 2,
            features = gene2)
ggsave(file=paste0("plot/GC_GSE183904/","feature_plot_gene2_3",".pdf"),width=15,height=15)
DotPlot(object = ref,features = gene1)
ggsave(file=paste0("plot/GC_GSE183904/","dot_plot_gene1_3",".pdf"),width=15,height=10)
DotPlot(object = ref,features = gene2)
ggsave(file=paste0("plot/GC_GSE183904/","dot_plot_gene2_3",".pdf"),width=15,height=10)
saveRDS(ref,file="data/GC_GSE183904/ref_3.rds")
Idents(ref)<-"seurat_clusters"
#1
new.cluster.ids <- c("T cells","T cells","Plasma cells","T cells","Epithelial Cells","Epithelial Cells","B cells",
                     "Fibroblasts","Endothelial Cells","NK cells","Macrophages","T regs","Epithelial Cells",
                     "Macrophages","Macrophages","Mast cells","T cells","Pericytes","Epithelial Cells",
                     "Fibroblasts","Fibroblasts","Plasma cells","T cells","Plasma cells","Pericytes","Epithelial Cells")
#2
new.cluster.ids <- c("T cells","T cells","Plasma cells","Fibroblasts","Macrophages","Endothelial Cells",
                     "Epithelial Cells","Epithelial Cells","Plasma cells","Macrophages","B cells",
                     "Mast cells","Pericytes","Macrophages","Fibroblasts","T cells",
                     "T reg","Fibroblasts","Plasma cells","Epithelial Cells","Pericytes",
                     "Fibroblasts","Epithelial Cells","T cells","Fibroblasts","Fibroblasts","T cells")
#3
new.cluster.ids <- c("T cells","T cells","Epithelial Cells","Epithelial Cells","NK cells","Plasma cells",
                     "Plasma cells","B cells","Endothelial Cells","T reg","Fibroblasts",
                     "Macrophages","Epithelial Cells","Macrophages","Pericytes","T cells",
                     "Mast cells","DCs","Epithelial Cells","Epithelial Cells","Epithelial Cells",
                     "Fibroblasts","Plasma cells","Epithelial Cells","Epithelial Cells","Pericytes",
                     "Plasma cells","Mast cells")
Idents(ref)<-"seurat_clusters"
names(new.cluster.ids) <- levels(ref)
ref <- RenameIdents(ref, new.cluster.ids)
ref$primary_cell_types <- ref@active.ident
DimPlot(ref, reduction = "umap")
Idents(ref)<-"primary_cell_types"
ref<-subset(ref, idents = c("Epithelial Cells"))
saveRDS(ref,file="data/GC_GSE183904/ref_3.rds")

####Step 2: CSC annotaion####
Idents(ref_1)<-"orig.ident"
sample<-unique(ref_1$orig.ident)[10]
ref_sub<-subset(ref_1,idents = c(sample))
ref_sub <- NormalizeData(ref_sub,normalization.method = "LogNormalize", scale.factor = 1e4)
ref_sub <- FindVariableFeatures(ref_sub, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(ref_sub)
ref_sub <- ScaleData(ref_sub, features = all.genes)
ref_sub <- RunPCA(ref_sub, features = VariableFeatures(object = ref_sub))
ref_sub <- FindNeighbors(ref_sub, dims = 1:20)
ref_sub <- FindClusters(ref_sub, resolution = 0.8)
ref_sub <- RunUMAP(ref_sub, dims = 1:20)
DimPlot(ref_sub, reduction = "umap",label=T)
gene<-c("CEACAM6","EGFR","MET","CCND1","KRAS")
DotPlot(object = ref_sub,
        features = gene)+RotatedAxis()
ref_sub<-AddModuleScore(ref_sub,features =gene,name = "malignancy")
ggplot(ref_sub@meta.data,aes(x=seurat_clusters,y=malignancy1,fill=seurat_clusters))+
  geom_boxplot()+labs(title="malignancy")
colorectal_cancer_stem_cell_markers <- read_csv("data/Gastric cancer stem cell markers.csv")
marker<-list(colorectal_cancer_stem_cell_markers$combined)
name <-"csc_combined"
ref_sub<-AddModuleScore(ref_sub,features =marker,name =name,assay="RNA",slot="scale.data")
ggplot(ref_sub@meta.data,aes(x=seurat_clusters,y=csc_combined1,fill=seurat_clusters))+
  geom_boxplot()+labs(title=name)
#probably better for individual,serve batch effect!
#dominant gene expressionz
FeaturePlot(object = ref_sub,label.size = 2,
            features = colorectal_cancer_stem_cell_markers$combined)
ggsave(file=paste0("plot/GC_GSE183904/",sample,".pdf"),width=15,height=15)
ref.markers <- FindAllMarkers(ref_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
df<-ref.markers[ref.markers$gene%in%marker[[1]],]
print(df)
table(df$cluster)
for(i in sort(unique(ref_sub$seurat_clusters))){
  print(paste0("cluster",i,":"))
  print(df[df$cluster==i,"gene"])
}
#29
new.cluster.ids <- c("other","other","other","other","other","other","other","cancer stem cell","other","other","cancer stem cell","other","other")
#30
new.cluster.ids <- c("other","other","cancer stem cell","other","other","cancer stem cell","other","other","other","other","other")
#31
new.cluster.ids <- c("other","other","other","other","other","other","cancer stem cell","other","other","other")
#36
new.cluster.ids <- c("other","other","other","other","other","other","cancer stem cell")
#39
new.cluster.ids <- c("other","cancer stem cell","other","cancer stem cell","other","cancer stem cell","other","other","other","cancer stem cell","other","other","other")
#17
new.cluster.ids <- c("other","cancer stem cell","other","other","other")
#18
new.cluster.ids <- c("other","other","other","other","other","cancer stem cell","other","other","other","other","other")
#24
new.cluster.ids <- c("other","other","other","other","cancer stem cell","cancer stem cell","other")
#26
new.cluster.ids <- c("cancer stem cell","other","other","cancer stem cell","other","other","other")
#27
new.cluster.ids <- c("cancer stem cell","other","other","cancer stem cell")
#28
new.cluster.ids <- c("cancer stem cell","other","cancer stem cell","other","other","other","other","other","other","other")
#2
new.cluster.ids <- c("cancer stem cell","other","cancer stem cell","other","other","other","other")
#5
new.cluster.ids <- c("cancer stem cell","other","other","other","other","cancer stem cell","other","other","other")
#7
new.cluster.ids <- c("other","other","cancer stem cell","other","cancer stem cell","cancer stem cell","other","other")
#8
new.cluster.ids <- c("other","cancer stem cell","other","other","cancer stem cell")
#10
new.cluster.ids <- c("other","other","cancer stem cell","other","other","other","other")
#14
new.cluster.ids <- c("other","other","other","cancer stem cell","other","cancer stem cell","other","other","cancer stem cell","cancer stem cell","other")
Idents(ref_sub)<-"seurat_clusters"
names(new.cluster.ids) <- levels(ref_sub)
ref_sub<- RenameIdents(ref_sub, new.cluster.ids)
ref_sub$csc<- ref_sub@active.ident
table(ref_sub$csc)
Idents(ref_sub)<-"csc"
DimPlot(ref_sub, reduction = "umap",label=F)

####Step3: generate tsv files####
sample_list<-names(ref_list)
for(id in sample_list){
  ref<-ref_list[[id]]
  # Specify the filename for the TSV file
  filename <- paste0("data/Compass/GC_GSE183904/",id,".tsv")
  # Write the data frame to a TSV file
  data<-ref@assays$RNA@layers$counts
  rownames(data)<-rownames(ref)[ref@assays$RNA@features@.Data[,paste0("counts.",id)]]
  colnames(data)<-colnames(ref)[ref@assays$RNA@cells@.Data[,paste0("counts.",id)]]
  write.table(data, file = filename, sep = "\t", quote = FALSE, row.names = T)
  filename <- paste0("data/Compass/GC_GSE183904/",id,"_cluster_meta.csv")
  data<-as.data.frame(ref@meta.data[,"seurat_clusters"])
  rownames(data)<-rownames(ref@meta.data)
  colnames(data)<-"seurat_cluster"
  write.csv(data, file = filename)
}
