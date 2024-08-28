library(Seurat)
library(SeuratData)
library(SeuratObject)
library(ggplot2)
library(patchwork)
library(dplyr)
sample_list<-c("CID4066","CID45171","CID4067","CID4290A","CID4530N","CID4535",
               "CID4495","CID4513","CID4515","CID4523","CID44971","CID44991")

####Step 1:pre-processing + Separate Epithelial Cells####
ref_list <- sapply(sample_list,function(name){
  expression_matrix <- ReadMtx(mtx = paste0("data\\BC_GSE176078\\",name,"\\count_matrix_sparse.mtx"),
                               cells = paste0("data\\BC_GSE176078\\",name,"\\count_matrix_barcodes.tsv"),
                               features = paste0("data\\BC_GSE176078\\",name,"\\count_matrix_genes.tsv"),
                               feature.column = 1)
  pbmc <- CreateSeuratObject(counts = expression_matrix, project = name, 
                             min.cells = 3, min.features = 200)
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2
  pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA > 250 & percent.mt < 20)
  #metadata
  library(readr)
  metadata <- read_csv(paste0("data\\BC_GSE176078\\",name,"\\metadata.csv"))
  metadata_cell_type_major<-metadata$celltype_major
  names(metadata_cell_type_major)<-metadata$...1
  metadata_cell_type_minor<-metadata$celltype_minor
  names(metadata_cell_type_minor)<-metadata$...1
  metadata_cell_type_subtype<-metadata$celltype_subset
  names(metadata_cell_type_subtype)<-metadata$...1
  pbmc<-AddMetaData(pbmc,metadata_cell_type_major,col.name = "celltype_major")
  pbmc<-AddMetaData(pbmc,metadata_cell_type_minor,col.name = "celltype_minor")
  pbmc<-AddMetaData(pbmc,metadata_cell_type_subtype,col.name = "celltype_subset")
  return(pbmc)
})
id<-sample_list[1]
ref<-ref_list[[id]]
Idents(ref)<-"celltype_major"
ref_sub<-subset(ref, idents = c("Cancer Epithelial"))

#####Step 2: CSC annotation####
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
colorectal_cancer_stem_cell_markers <- read_csv("data/breast cancer stem cell markers.csv")
marker<-list(colorectal_cancer_stem_cell_markers$`Walcher, L et.al`)
name <-"csc_combined"
ref_sub<-AddModuleScore(ref_sub,features =marker,name =name)
ggplot(ref_sub@meta.data,aes(x=seurat_clusters,y=csc_combined1,fill=seurat_clusters))+
  geom_boxplot()+labs(title=name)
ref_sub<-AddModuleScore(ref_sub,features =marker,name =paste0(name,"1"))
ggplot(ref_sub@meta.data,aes(x=seurat_clusters,y=csc_combined11,fill=seurat_clusters))+
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
#CID4066
new.cluster.ids <- c("cancer stem cell","other","other","other","other","other")
#CID4513
new.cluster.ids <- c("cancer stem cell","cancer stem cell","other","other","other","other","other","other")
#CID4523
new.cluster.ids <- c("other","other","other","other","cancer stem cell","cancer stem cell","other","other","other")
#CID4067
new.cluster.ids <- c("other","cancer stem cell","other","other","other","other","other","other","other","other","other")
#CID4495
new.cluster.ids <- c("other","other","other","other","other","cancer stem cell","other","other")
#CID4515
new.cluster.ids <- c("other","cancer stem cell","other","other","other","other","other","other","cancer stem cell","other","other")
#CID44971
new.cluster.ids <- c("cancer stem cell","other","other","other","other","cancer stem cell","cancer stem cell","other","other","other","cancer stem cell","other")
Idents(ref_sub)<-"seurat_clusters"
names(new.cluster.ids) <- levels(ref_sub)
ref_sub<- RenameIdents(ref_sub, new.cluster.ids)
ref_sub$csc<- ref_sub@active.ident
table(ref_sub$csc)
Idents(ref_sub)<-"csc"
DimPlot(ref_sub, reduction = "umap",label=F)
ref_list[[id]]<-ref_sub

#####Step 3: generate tsv file####
for(id in sample_list){
  ref<-ref_list[[id]]
  Idents(ref)<-"celltype_major"
  ref<-subset(ref,idents = c("Cancer Epithelial"))
  # Specify the filename for the TSV file
  filename <- paste0("data/Compass/ST_BC_GSE176078/",id,".tsv")
  # Write the data frame to a TSV file
  data<-ref@assays$RNA@layers$counts
  rownames(data)<-rownames(ref)[ref@assays$RNA@features@.Data[,"counts"]]
  colnames(data)<-colnames(ref)[ref@assays$RNA@cells@.Data[,"counts"]]
  write.table(data, file = filename, sep = "\t", quote = FALSE, row.names = T,col.names = TRUE)
  filename <- paste0("data/Compass/ST_BC_GSE176078/",id,"_cell_meta.csv")
  data<-as.data.frame(ref@meta.data[,"csc"])
  rownames(data)<-rownames(ref@meta.data)
  colnames(data)<-"celltype"
  write.csv(data, file = filename)
}

####Step 4: cell cycle score#####
library(Seurat)
sample_list<-c("CID4515","CID4067","CID44971","CID4495","CID4513","CID4066","CID4523")
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# Read in the expression matrix The first row is a header row, the first column is rownames
id<-sample_list[1]
ref<-ref_list[[id]]
Idents(ref)<-"celltype_major"
ref_sub<-subset(ref, idents = c("Cancer Epithelial"))
ref_sub <- NormalizeData(ref_sub,normalization.method = "LogNormalize", scale.factor = 1e4)
ref_sub <- FindVariableFeatures(ref_sub, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(ref_sub)
ref_sub <- ScaleData(ref_sub, features = all.genes)
ref_sub <- RunPCA(ref_sub, features = VariableFeatures(object = ref_sub))
ref_sub <- FindNeighbors(ref_sub, dims = 1:20)
ref_sub <- FindClusters(ref_sub, resolution = 1)
ref_sub <- RunUMAP(ref_sub, dims = 1:20)
Idents(ref_sub)<-"seurat_clusters"
DimPlot(ref_sub, reduction = "umap",label=T)
ref_sub <- CellCycleScoring(ref_sub, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
cellcycle <- data.frame(,G2M=ref_sub$G2M.Score,S = ref_sub$S.Score,seurat_clusters = ref_sub$seurat_clusters)
rownames(cellcycle)<-rownames(ref_sub@meta.data)
library(ggplot2)
P1<-ggplot(cellcycle,aes(x=seurat_clusters,y=G2M,fill=seurat_clusters))+
  geom_boxplot()+labs(title=paste0(id,"_G2M"))
P2<-ggplot(cellcycle,aes(x=seurat_clusters,y=S,fill=seurat_clusters))+
  geom_boxplot()+labs(title=paste0(id,"_S"))
P1+P2
