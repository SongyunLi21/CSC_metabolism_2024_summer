library(Seurat)
#library(SeuratData)
library(SeuratObject)
library(ggplot2)
library(patchwork)
library(dplyr)
library(RColorBrewer)
####Step 1: preprocessing####
# Load the expression data
st_index_list<-c("C1","C2","C3","C4","L1","L2")
st_list <- sapply(st_index_list,function(id){
  # Load the expression data
  expr.url <- paste0("data/CRC_ST_GSE225857/", id)
  expr.data <- Read10X(data.dir =  expr.url)
  stref <- Seurat::CreateSeuratObject(counts = expr.data, project = id, assay = 'Spatial')
  # Load the image data
  path <- paste0("data/CRC_ST_GSE225857/",id,'//spatial')
  img <- Seurat::Read10X_Image(image.dir = path)
  Seurat::DefaultAssay(object = img) <- 'Spatial'
  img <- img[colnames(x = stref)]
  stref[['image']] <- img
  return(stref)
})
#saveRDS(st_list,"data/st_list_h.rds")#save
marker<-c("ALDH1A1","PROM1","CD44","EPCAM","LGR5","SOX2")
i=1
for(id in st_index_list){
  st_list[[id]]<- NormalizeData(st_list[[id]],normalization.method = "LogNormalize", scale.factor = 1e4)
  #SpatialFeaturePlot(st_list[[id]], features = marker,
                     #image.alpha=0,ncol=4,pt.size.factor =c(1.8,1.5,1.5,1.5,1.5,1.5)[i],slot = "data",stroke=0)
  #ggsave(file=paste0("plot/CRC_ST_GSE225857/",id,"csc_marker_expression_normalized",".pdf"),width=18,height=15)
  i =i+1
}


####Step 2: marker expression####
for(top_percent in c(0.03,0.05,0.08,0.10)){
  i=1
  for(id in st_index_list){
    ident_list<-c()
    st_list[[id]]@meta.data[,paste0("top_normalized_",top_percent*100,"%","_positive")]<-rep(0,length(rownames(st_list[[id]]@meta.data)))
    for(gene in marker){
      st_list[[id]]@meta.data[,paste0(gene,"_normalized")]<-st_list[[id]]@assays$Spatial@layers$data[rownames(st_list[[id]])==gene,]
      num_top <- ceiling(length(st_list[[id]]@meta.data[,paste0(gene,"_normalized")]) * top_percent)
      top_indices <- order(st_list[[id]]@meta.data[,paste0(gene,"_normalized")], decreasing = TRUE)[1:num_top]
      st_list[[id]]@meta.data[,paste0(gene,"top_normalized_",top_percent*100,"%")]<-rep(0,length(st_list[[id]]@meta.data[,paste0(gene,"_normalized")]))
      top_indices<-top_indices[st_list[[id]]@meta.data[top_indices,paste0(gene,"_normalized")]>0]
      st_list[[id]]@meta.data[top_indices,paste0(gene,"top_normalized_",top_percent*100,"%")]<-1
      Idents(st_list[[id]])<-paste0(gene,"top_normalized_",top_percent*100,"%")
      st_list[[id]]@meta.data[,paste0("top_normalized_",top_percent*100,"%","_positive")]<-st_list[[id]]@meta.data[,paste0("top_normalized_",top_percent*100,"%","_positive")]+st_list[[id]]@meta.data[,paste0(gene,"top_normalized_",top_percent*100,"%")]
      ident_list <- append(ident_list,paste0(gene,"top_normalized_",top_percent*100,"%"))
    }
    SpatialDimPlot(st_list[[id]],group.by = ident_list,stroke=0.01,
                   facet.highlight = TRUE,label = FALSE, label.size = 1,ncol=3,
                   image.alpha=0,pt.size.factor=c(1.8,1.5,1.5,1.5,1.5,1.5)[i])
    ggsave(file=paste0("plot/CRC_ST_GSE225857/",id,"normalized_",top_percent,"_expression",".pdf"),width=15,height=13)
    SpatialDimPlot(st_list[[id]],group.by = paste0("top_normalized_",top_percent*100,"%","_positive"),
                   facet.highlight = TRUE,label = FALSE, label.size = 1,ncol=1,
                   image.alpha=0,stroke=0.01,
                   pt.size.factor =c(1.8,1.5,1.5,1.5,1.5,1.5)[i])+scale_fill_manual(values=brewer.pal(name="YlGnBu", n = length(marker)+1))
    ggsave(file=paste0("plot/CRC_ST_GSE225857/",id,"normalized_",top_percent,"_csc_niche",".pdf"),width=5,height=5)
    i=i+1
  }
}

####Step3: Scoring method####
i=1
for(id in st_index_list){
  vectors <- list(st_list[[id]]$CD44_normalized, st_list[[id]]$ALDH1A1_normalized, st_list[[id]]$PROM1_normalized,st_list[[id]]$EPCAM_normalized, st_list[[id]]$LGR5_normalized,st_list[[id]]$SOX2_normalized)
  # Normalize function
  normalize <- function(vec) {
    (vec - min(vec)) / (max(vec) - min(vec))
  }
  # Normalize each vector to [0, 1]
  normalized_vectors <- lapply(vectors, normalize)
  # Sum normalized vectors element-wise
  result_vector <- Reduce(`+`, normalized_vectors)
  st_list[[id]]<-AddMetaData(st_list[[id]],result_vector,col="sum_of_csc_marker")
  SpatialFeaturePlot(st_list[[id]], features = "sum_of_csc_marker",
                     image.alpha=0,pt.size.factor =c(1.8,1.5,1.5,1.5,1.5,1.5)[i],slot = "data")
  ggsave(file=paste0("plot/CRC_ST_GSE225857/",id,"sum_of_csc_marker_normalized",".pdf"),width=5,height=5)
  i =i+1
}

library(UCell)
i=1
for(id in st_index_list){
  gene.sets <- list(cance_stem_cell = marker)
  sample.matrix<-st_list[[id]]@assays$Spatial@layers$data
  colnames(sample.matrix)<-colnames(st_list[[id]])[st_list[[id]]@assays$Spatial@cells@.Data[,"data"]]
  rownames(sample.matrix)<-rownames(st_list[[id]])[st_list[[id]]@assays$Spatial@features@.Data[,"data"]]
  scores <- ScoreSignatures_UCell(sample.matrix, features=gene.sets)
  st_list[[id]]<-AddMetaData(st_list[[id]],scores,col="cance_stem_cell_UCell")
  SpatialFeaturePlot(st_list[[id]], features = "cance_stem_cell_UCell",
                     image.alpha=0,pt.size.factor =c(1.8,1.5,1.5,1.5,1.5,1.5)[i],slot = "data")
  ggsave(file=paste0("plot/CRC_ST_GSE225857/",id,"cance_stem_cell_UCell",".pdf"),width=5,height=5)
  i =i+1
}
i=1
for(id in st_index_list){
  ref<-st_list[[id]]
  ref<-AddModuleScore(ref,features=marker,name="csc_add_module_score",slot="data",nbin=24)
  SpatialFeaturePlot(ref, features = "csc_add_module_score1",
                     image.alpha=0,pt.size.factor =c(1.8,1.5,1.5,1.5,1.5,1.5)[i],slot = "count")
  ggsave(file=paste0("plot/CRC_ST_GSE225857/",id,"csc_add_modules_core",".pdf"),width=5,height=5)
  i =i+1
}
write.csv(ref$celltype,"data/CRC_ST_GSE225857/celltype.csv")

####Step4: CARD####
#scRNA seq data
#This dataset contain two files: immune cells and unimmune cells, my current reuslt only use non-immmune cells
#non immune cells
path<- "data/CRC_ST_GSE225857/non_immune/non_immune_counts.txt"
ref <- read.table(file=path,head=TRUE,sep="\t")
a<-ref$X
rownames(ref)<-a
ref$X<-NULL
ref <- CreateSeuratObject(counts = ref, project = "non_immune" , 
                          min.cells = 3, min.features = 200)
ref[["percent.mt"]] <- PercentageFeatureSet(ref, pattern = "^MT-")
VlnPlot(ref, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
metadata <- read.table(paste0("data/CRC_ST_GSE225857/non_immune/non_immune_meta.txt"),header = TRUE, sep = "\t")
sample_id_info<-as.data.frame(metadata$orig.ident)
library(do)
a<-Replace(metadata$X,from="-",to=".")
rownames(sample_id_info)<-a
ref<-AddMetaData(ref,sample_id_info,col.name = "sample_id")
table(ref$sample_id)
sample_info<-as.data.frame(metadata$Sample)
cell_type_info<-as.data.frame(metadata$cluster)
rownames(cell_type_info)<-a
ref<-AddMetaData(ref,cell_type_info,col.name = "cell_type")
table(ref$cell_type)
cluster_info<-as.data.frame(metadata$seurat_clusters)
rownames(cluster_info)<-a
ref<-AddMetaData(ref,cluster_info,col.name = "seurat_clusters")
table(ref$seurat_clusters)
table(ref$seurat_clusters,ref$cell_type)
Idents(ref)<-"cell_type"
new.cluster.ids <- c("Tumor cells","Tumor cells","Tumor cells","fibrblast","Tumor cells",
                     "Tumor cells","Tumor cells","Tumor cells","fibrblast","Tumor cells",
                     "Tumor cells","fibrblast","fibrblast","endothelial","fibrblast",
                     "fibrblast","endothelial","endothelial","endothelial","endothelial",
                     "endothelial","Tumor cells","Tumor cells")
#immune cells
path<- "data/CRC_ST_GSE225857/immune/immune_counts.txt"
ref_1 <- read.table(file=path,head=TRUE,sep="\t")
a<-ref_1$X
rownames(ref_1)<-a
ref_1$X<-NULL
ref_1 <- CreateSeuratObject(counts = ref_1, project = "immune" , 
                          min.cells = 3, min.features = 200)
ref_1[["percent.mt"]] <- PercentageFeatureSet(ref_1, pattern = "^MT-")
VlnPlot(ref_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
metadata <- read.table(paste0("data/CRC_ST_GSE225857/immune/immune_meta.txt"),header = TRUE, sep = "\t")
sample_id_info<-as.data.frame(metadata$orig.ident)
library(do)
a<-Replace(metadata$X,from="-",to=".")
rownames(sample_id_info)<-a
ref_1<-AddMetaData(ref_1,sample_id_info,col.name = "sample_id")
table(ref_1$sample_id)
sample_info<-as.data.frame(metadata$Sample)
cell_type_info<-as.data.frame(metadata$cluster)
rownames(cell_type_info)<-a
ref_1<-AddMetaData(ref_1,cell_type_info,col.name = "cell_type")
table(ref_1$cell_type)
cluster_info<-as.data.frame(metadata$seurat_clusters)
rownames(cluster_info)<-a
ref_1<-AddMetaData(ref_1,cluster_info,col.name = "seurat_clusters")
table(ref_1$seurat_clusters)
table(ref_1$seurat_clusters,ref_1$cell_type)
Idents(ref_1)<-"cell_type"
names(new.cluster.ids) <- levels(ref_1)
ref_1 <- RenameIdents(ref_1, new.cluster.ids)
ref_1$cell_type_major<-ref@active.ident
Idents(ref_1)<-"cell_type_major"
ref_1$celltype<-ref_1$cell_type_major
celltype <- c()
table(ref@active.ident)

#csc annotation
ref_index_list<-unique(ref$sample_id)
id <- ref_index_list[5]
Idents(ref)<-"sample_id"
ref_sub<-subset(ref,idents = c(id))
Idents(ref_sub)<-"cell_type_major"
ref_sub<-subset(ref_sub,idents = c("Tumor cells"))
ref_sub <- NormalizeData(ref_sub,normalization.method = "LogNormalize", scale.factor = 1e4)
ref_sub <- FindVariableFeatures(ref_sub, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(ref_sub)
ref_sub <- ScaleData(ref_sub, features = all.genes)
ref_sub <- RunPCA(ref_sub, features = VariableFeatures(object = ref_sub))
ref_sub <- FindNeighbors(ref_sub, dims = 1:20)
ref_sub <- FindClusters(ref_sub, resolution = 0.8)
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
df<-ref.markers[ref.markers$gene%in%marker[[1]],]
print(df)
table(df$cluster)
for(i in sort(unique(ref_sub$seurat_clusters))){
  print(paste0("cluster",i,":"))
  print(df[df$cluster==i,"gene"])
}
new.cluster.ids <- c("other","other","cancer stem cell","cancer stem cell","cancer stem cell","other","cancer stem cell","other","cancer stem cell")
new.cluster.ids <- c("cancer stem cell","other","other","cancer stem cell","other","cancer stem cell","other","cancer stem cell","other","other","other","other","cancer stem cell","other","cancer stem cell","cancer stem cell")
new.cluster.ids <- c("other","cancer stem cell","cancer stem cell","other","other","other","other","other","other","other","cancer stem cell","other")
new.cluster.ids <- c("other","cancer stem cell","other")
new.cluster.ids <- c("cancer stem cell","cancer stem cell","other","cancer stem cell","cancer stem cell","cancer stem cell","other","other","other","other","other","other")
Idents(ref_sub)<-"seurat_clusters"
names(new.cluster.ids) <- levels(ref_sub)
ref_sub<- RenameIdents(ref_sub, new.cluster.ids)
ref_sub$csc<- ref_sub@active.ident
table(ref_sub$csc)
Idents(ref_sub)<-"csc"
DimPlot(ref_sub, reduction = "umap",label=F)
celltype<-append(celltype,ref_sub$csc)
celltype<-append(celltype,ref$cell_type_major[ref$cell_type_major!="Tumor cells"])
ref<-AddMetaData(ref,celltype,col="celltype")
ref$celltype<-as.character(ref$celltype)
ref<-NormalizeData(ref,normalization.method = "LogNormalize", scale.factor = 1e4)
#saveRDS(ref,file="data/ST_BC_GSE176078/ref.rds")

#CARD
for(id in names(st_list)){
  spatial_count<-st_list[[id]]@assays$Spatial@layers$data
  rownames(spatial_count)<-rownames(st_list[[id]])[st_list[[id]]@assays$Spatial@features@.Data[,"data"]]
  colnames(spatial_count)<-colnames(st_list[[id]])[st_list[[id]]@assays$Spatial@cells@.Data[,"data"]]
  spatial_count[1:4,1:4]
  spatial_location<-st_list[[id]]@images$image@coordinates[st_list[[id]]@assays$Spatial@cells@.Data[,"data"],]
  spatial_location$tissue<-NULL
  library(CARD)
  spatial_location$row<-NULL
  spatial_location$col<-NULL
  library(do)
  name<-rownames(spatial_location)
  name<-Replace(name,from=paste0(id,"_"),to="")
  rownames(spatial_location)<-name
  colnames(spatial_location)<-c("x","y")
  spatial_location[1:4,]
  spatial_count<-spatial_count[,colnames(spatial_count)%in%rownames(spatial_location)]
  spatial_count[1:4,1:4]
  spatial_location<-spatial_location[rownames(spatial_location)%in%colnames(spatial_count),]
  spatial_location[1:4,]
  sc_count<-ref@assays$RNA@layers$data
  rownames(sc_count)<-rownames(ref)[ref@assays$RNA@features@.Data[,"data"]]
  colnames(sc_count)<-colnames(ref)[ref@assays$RNA@cells@.Data[,"data"]]
  #write.csv(sc_count,"data/CRC_ST_GSE225857/sc_count.csv")
  #sc_count<-ref_list[["CID44971"]]@assays$RNA@layers$counts
  #rownames(sc_count)<-rownames(ref_list[["CID44971"]])[ref_list[["CID44971"]]@assays$RNA@features@.Data[,"counts"]]
  #colnames(sc_count)<-colnames(ref_list[["CID44971"]])[ref_list[["CID44971"]]@assays$RNA@cells@.Data[,"counts"]]
  sc_count[1:4,1:4]
  sc_meta<-as.data.frame(ref@meta.data[,c("celltype","orig.ident")])
  #sc_meta<-as.data.frame(ref_list[["CID44971"]]@meta.data[,c("celltype","orig.ident")])
  colnames(sc_meta)<-c("cellType","sampleInfo")
  sc_meta[1:4,]
  #write.csv(sc_meta,"data/CRC_ST_GSE225857/sc_meta.csv")
  CARD_obj = createCARDObject(
    sc_count = sc_count,
    sc_meta = sc_meta,
    spatial_count = spatial_count,
    spatial_location = spatial_location,
    ct.varname = "cellType",
    ct.select = unique(sc_meta$cellType),
    sample.varname = "sampleInfo",
    minCountGene = 100,
    minCountSpot = 5)
  CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
  print(CARD_obj@Proportion_CARD[1:2,])
  write.csv(CARD_obj@Proportion_CARD,paste0("data/ST_BC_GSE176078/",id,"_","CARD_score_all_normalized",".csv"))
  ct.visualize = unique(sc_meta$cellType[!is.na(sc_meta$cellType)])
  ## visualize the spatial distribution of the cell type proportion
  p2 <- CARD.visualize.prop(
    proportion = CARD_obj@Proportion_CARD,
    spatial_location = CARD_obj@spatial_location,
    ct.visualize = ct.visualize,                 ### selected cell types to visualize
    colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
    NumCols = 4,                                 ### number of columns in the figure panel
    pointSize = 1.2)                             ### point size in ggplot2 scatterplot
  print(p2)
  ggsave(file=paste0("plot/CRC_ST_GSE225857/CARD/",id,"_","all_cells","_single_harmony_normalized_normalized",".pdf"),width=15,height=13)
  ## visualize the spatial distribution of two cell types on the same plot
  p3 = CARD.visualize.prop.2CT(
    proportion = CARD_obj@Proportion_CARD,                             ### Cell type proportion estimated by CARD
    spatial_location = CARD_obj@spatial_location,                      ### spatial location information
    ct2.visualize = c("cancer stem cell","other"),              ### two cell types you want to visualize
    colors = list(c("lightblue","lightyellow","red"),c("lightblue","lightyellow","black")))       ### two color scales
  print(p3)
  ggsave(file=paste0("plot/CRC_ST_GSE225857/CARD/",id,"_","all_cells","_double_all_normalized",".pdf"),width=12,height=12)
  p4 <- CARD.visualize.Cor(CARD_obj@Proportion_CARD,colors = NULL) # if not provide, we will use the default colors
  print(p4)
  ggsave(file=paste0("plot/CRC_ST_GSE225857/CARD/",id,"_","portion_corr_all_normalized",".pdf"),width=10,height=10)
}
