library(Seurat)
library(SeuratData)
library(SeuratObject)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Seurat)
library(SeuratData)
library(SeuratObject)
library(ggplot2)
library(patchwork)
library(dplyr)
####Step 1:pre-processing + Separate Epithelial Cells####
ref_index_list <- c("A2019_1","A2019_2","A2019_3",
                    "B2019_1","B2019_2","B2019_3",
                    "C2020_1","C2020_2","C2020_3",
                    "D2020_1","D2020_2","D2020_3",
                    "E2020_1","E2020_2","E2020_3")
ref_list <- sapply(ref_index_list,function(id){
  path <- paste0("data/BC_GSE180286/",id,".txt")
  ref <- read.table(file=path,head=TRUE,sep="\t")
  ref <- CreateSeuratObject(counts = ref, project = id , min.cells = 3, min.features = 500)
  ref[["percent.mt"]] <- PercentageFeatureSet(ref, pattern = "^MT-")
  VlnPlot(ref, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
  plot1 <- FeatureScatter(ref, feature1 = "nCount_RNA", feature2 = "percent.mt") 
  plot2 <- FeatureScatter(ref, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
  plot1 + plot2
  ref <- subset(ref, subset = nFeature_RNA > 100 & nCount_RNA > 200 & percent.mt < 50)
  return(ref)
})
saveRDS(ref_list,"data/BC_GSE180286/ref_list.rds")#save

#normalization
library(harmony)
ref <- merge(ref_list[['A2019_1']], y = ref_list[ref_index_list[c(2:15)]],add.cell.ids = ref_index_list, merge.data = TRUE)
unique(sapply(X = strsplit(colnames(ref), split = "_"), FUN = "[", 1))
saveRDS(ref,file="data/BC_GSE180286/ref.rds")
remove(ref_list)
ref <- NormalizeData(ref,normalization.method = "LogNormalize", scale.factor = 1e4)
#HVG
ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(ref)
#saveRDS(ref,file="data/BC_GSE180286/ref.rds")
ref <- ScaleData(ref, features = all.genes)
#saveRDS(ref,file="data/BC_GSE180286/ref.rds")
ref <- RunPCA(ref, features = VariableFeatures(object = ref))
ref <- RunHarmony(ref, group.by.vars = "patient")
ref <- IntegrateLayers(
  object = ref, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
#saveRDS(ref,file="data/BC_GSE180286/ref.rds")
ref <- FindNeighbors(ref, dims = 1:20,reduction = "harmony")
ref <- FindClusters(ref, resolution = 0.8,reduction = "harmony")
ref <- RunUMAP(ref, dims = 1:20,reduction = "harmony")
DimPlot(ref, reduction = "umap")
saveRDS(ref,file="data/BC_GSE180286/ref.rds")
for(id in unique(ref@meta.data$orig.ident)){
  if(id == "A2019_1"||id == "A2019_2"||id == "A2019_3"){
    ref@meta.data[ref$orig.ident==id,"patient"] = "P1"
  }else if(id == "B2019_1"||id == "B2019_2"||id == "B2019_3"){
    ref@meta.data[ref@meta.data$orig.ident==id,"patient"] = "P2"
  }else if(id == "C2020_1"||id == "C2020_2"||id == "C2020_3"){
    ref@meta.data[ref@meta.data$orig.ident==id,"patient"] = "P3"
  }else if(id == "D2020_1"||id == "D2020_2"||id == "D2020_3"){
    ref@meta.data[ref@meta.data$orig.ident==id,"patient"] = "P4"
  }else if(id == "E2020_1"||id == "E2020_2"||id == "E2020_3"){
    ref@meta.data[ref@meta.data$orig.ident==id,"patient"] = "P5"
  }
}
Idents(ref)<-"patient"
DimPlot(ref, reduction = "umap")
for(id in c("P1","P2","P3","P4","P5")){
  if(id == "P5"|| id=="P4"){
    ref@meta.data[ref@meta.data$patient==id,"TNBC"] = "TNBC"
  }else{
    ref@meta.data[ref@meta.data$patient==id,"TNBC"] = "non_TNBC"
  }
}
Idents(ref)<-"TNBC"
DimPlot(ref, reduction = "umap")
marker_gene <- c("EPCAM","ERBB2", "ESR1", "PGR", "MKI67", "CD44", "CD24", "ALDH1A1",
                 "DCN", "COL1A2", "ACTA2", "VWF", "PECAM1", "KRT14", "MCAM", "MYH11", "TAGLN",
                 "CD3G","CD3D","CD8A","NKG7","FOXP3","IL2RA","KLRD1", "KLRB1", "CD79A", "MS4A1", 
                 "JCHAIN","IGHG4","CD163","CD86","MSR1", "TPSD1", "KIT", "LAMP3", "CCR7", "IL3RA", "CLEC4C")
#dot plot
DotPlot(ref, features = marker_gene) + RotatedAxis()
ggsave(file=paste0("plot/TNBC/BC_GSE180286/","dot_anno_plot.pdf"),height=10,width=20,dpi=600)
DotPlot(object = ref,
        features = c("EPCAM","MKI67","CD3D","CD86","MS4A1","JCHAIN","PECAM1","PDGFRB"))+RotatedAxis()
ggsave(file=paste0("plot/TNBC/BC_GSE180286/","dot_plot_new_marker",".pdf"),width=15,height=15)
DotPlot(object = ref,
        features = c("DCN", "COL1A2", "ACTA2", "VWF", "PECAM1", "KRT14", "MCAM", "MYH11", "TAGLN"))+RotatedAxis()
ggsave(file=paste0("plot/TNBC/BC_GSE180286/","dot_plot_new_marker_2",".pdf"),width=15,height=15)
#heatmap
Idents(ref)<-"secondary_cell_types"
DoHeatmap(ref , features = marker_gene, size = 3)
ggsave(file=paste0("plot/TNBC/BC_GSE180286/","heatmap.pdf"),height=10,width=25,dpi=600)

#violin plot
VlnPlot(object = ref,ncol=3,
        features = marker_gene)
ggsave(file=paste0("plot/TNBC/BC_GSE180286/","vln_anno_plot",".pdf"),width=30,height=20)

#feature plot
#tumor marker
FeaturePlot(object = ref,label.size = 2,
            features = c("MKI67"))
#tumor marker
FeaturePlot(object = ref,label.size = 2,
            features = c("EPCAM","ERBB2", "ESR1", "PGR", "MKI67", "CD44", "CD24", "ALDH1A1",
                         "DCN", "COL1A2", "ACTA2", "VWF", "PECAM1", "KRT14", "MCAM", "MYH11", "TAGLN"))
ggsave(file=paste0("plot/TNBC/BC_GSE180286/","feature_plot_1",".pdf"),width=32,height=32)
FeaturePlot(object = ref,label.size = 2,
            features = c("CD3G","CD3D","CD8A","NKG7","FOXP3","IL2RA","KLRD1", "KLRB1", "CD79A", "MS4A1", 
                         "JCHAIN","IGHG4","CD163","CD86","MSR1", "TPSD1", "KIT", "LAMP3", "CCR7", "IL3RA", "CLEC4C"))
ggsave(file=paste0("plot/TNBC/BC_GSE180286/","feature_plot_2",".pdf"),width=32,height=32)
FeaturePlot(object = ref,label.size = 2,
            features = c("EPCAM","MKI67","CD3D","CD86","MS4A1","JCHAIN","PECAM1","PDGFRB"))
ggsave(file=paste0("plot/TNBC/BC_GSE180286/","feature_plot_new_marker",".pdf"),width=15,height=15)
Idents(ref)<-"seurat_clusters"
new.cluster.ids<-c("Epithelial cells","Fibroblasts","CD4 naive T cells","B cells","Epithelial cells","CD4 naive T cells","Epithelial cells","Cycling cells","Myeloid cells","NK","CD4 naive T cells",
                   "Cycling cells","Epithelial cells","Plasma cells","Myeloid cells","Endothelial cells","B cells","mature DC","CD4 naive T cells","Fibroblasts","pDC",
                   "CD4 naive T cells","Epithelial cells","Fibroblasts","Epithelial cells","Myeloid cells","Epithelial cells")
names(new.cluster.ids) <- levels(ref)
ref <- RenameIdents(ref, new.cluster.ids)
ref$secondary_cell_types <- ref@active.ident
Idents(ref)<-"secondary_cell_types"
DimPlot(ref, reduction = "umap",label=T)

#####Step 2: CSC annotation####
id<-"P5"
ref<-ref_sub_list[[id]]
DimPlot(ref, reduction = "umap",label=T)
library(readr)
colorectal_cancer_stem_cell_markers <- read_csv("data/breast cancer stem cell markers.csv")
marker<-list(colorectal_cancer_stem_cell_markers$`Walcher, L et.al`)
name <-"csc_combined"
ref<-JoinLayers(ref)
ref<-AddModuleScore(ref,features =marker,name =name)
ggplot(ref@meta.data,aes(x=seurat_clusters,y=csc_combined1,fill=seurat_clusters))+
  geom_boxplot()+labs(title=name)
ref.markers <- FindAllMarkers(ref, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
df<-ref.markers[ref.markers$p_val_adj<0.05,]
df<-df[df$gene%in%marker[[1]],]
print(df)
table(df$cluster)
for(i in sort(unique(ref$seurat_clusters))){
  print(paste0("cluster",i,":"))
  print(df[df$cluster==i,"gene"])
}
#P1
new.cluster.ids <- c("cancer stem cell","other","other","other","other")
Idents(ref)<-"seurat_clusters"
names(new.cluster.ids) <- levels(ref)
ref<- RenameIdents(ref, new.cluster.ids)
ref$csc_add_module_score<- ref@active.ident
Idents(ref)<-"csc_add_module_score"
DimPlot(ref, reduction = "umap",label=F)

#####Step 3: generate tsv file####
data<-ref@assays$RNA@layers$counts
rownames(data)<-rownames(ref)[ref@assays$RNA@features@.Data[,"counts"]]
colnames(data)<-colnames(ref)[ref@assays$RNA@cells@.Data[,"counts"]]
filename <- paste0("data/Compass/BC_GSE180286/",id,".tsv")
write.table(data, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
filename <- paste0("data/Compass/BC_GSE180286/",id,"_cell_meta.csv")
data<-ref@meta.data[,c("csc")]
write.csv(data, file = filename)
