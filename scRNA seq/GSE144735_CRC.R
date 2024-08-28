library(Seurat)
library(SeuratData)
library(SeuratObject)
library(ggplot2)
library(patchwork)
library(dplyr)
####Step 1:pre-processing + Separate Epithelial Cells####
path <- paste0("data/CRC_GSE144735/","count_matrix",".txt")
ref <- read.table(file=path,head=TRUE,sep="\t")
a<-ref$Index
rownames(ref)<-a
ref$Index<-NULL
ref <- CreateSeuratObject(counts = ref, project = "CRC_GSE144735" , 
                          min.cells = 3, min.features = 200)
ref[["percent.mt"]] <- PercentageFeatureSet(ref, pattern = "^MT-")
VlnPlot(ref, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
library(readr)
path <- paste0("data/CRC_GSE144735/","annotation",".txt")
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
ref<-subset(ref,idents=c("Border","Tumor"))
saveRDS(ref,"data/CRC_GSE144735/ref.rds")
data<-as.data.frame(ref@assays$RNA@layers$counts)
rownames(data)<-rownames(ref)
colnames(data)<-colnames(ref)
cluster<-as.data.frame(ref@meta.data[,"seurat_clusters"])
rownames(cluster)<-colnames(ref)
colnames(cluster)<-NULL
write.csv(data,file="data/CRC_GSE144735/count.csv")
write.csv(cluster,file="data/CRC_GSE144735/cluster.csv")
ref <- NormalizeData(ref,normalization.method = "LogNormalize", scale.factor = 1e4)
ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(ref)
ref <- ScaleData(ref, features = all.genes)
ref <- RunPCA(ref, features = VariableFeatures(object = ref))
ref <- FindNeighbors(ref, dims = 1:20)
ref <- FindClusters(ref, resolution = 0.8)
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

#####Step 2: CSC annotation####
sample<-unique(ref$patient)
id<-sample[6]
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

#KUL01
new.cluster.ids <- c("cancer stem cell","other","cancer stem cell","cancer stem cell","other","other","other","cancer stem cell","cancer stem cell","cancer stem cell")
#KUL19
new.cluster.ids <- c("other","other","other","other","cancer stem cell","other","cancer stem cell")
#KUL31
new.cluster.ids <- c("cancer stem cell","other","cancer stem cell","other","other","other")
Idents(ref_sub)<-"seurat_clusters"
names(new.cluster.ids) <- levels(ref_sub)
ref_sub<- RenameIdents(ref_sub, new.cluster.ids)
ref_sub$csc<- ref_sub@active.ident
table(ref_sub$csc)
Idents(ref_sub)<-"csc"
DimPlot(ref_sub, reduction = "umap",label=F)
ref_list[[id]]<-ref_sub
saveRDS(ref_list,"data/CRC_GSE144735/ref_list.rds")

#####Step 3: generate tsv file####
sample<-c("KUL01","KUL19","KUL31")
for(id in sample){
  ref<-ref_list[[id]]
  #Specify the filename for the TSV file
  filename <- paste0("data/Compass/CRC_GSE144735/",id,".tsv")
  # Write the data frame to a TSV file
  data<-ref@assays$RNA@layers$counts
  rownames(data)<-rownames(ref)[ref@assays$RNA@features@.Data[,"counts"]]
  colnames(data)<-colnames(ref)[ref@assays$RNA@cells@.Data[,"counts"]]
  write.table(data, file = filename, sep = "\t", quote = FALSE, row.names = T)
  filename <- paste0("data/Compass/CRC_GSE144735/",id,"_cluster_meta.csv")
  data<-as.data.frame(ref@meta.data[,"csc"])
  rownames(data)<-rownames(ref@meta.data)
  colnames(data)<-c("celltype")
  write.csv(data, file = filename)
}
