library(Seurat)
#library(SeuratData)
library(SeuratObject)
library(ggplot2)
library(patchwork)
library(dplyr)
library(RColorBrewer)
####Step1:preprocessing####
# Load the expression data
st_index_list<-c("1142243F","1160920F","CID4290","CID4465","CID4535","CID44971")
st_list <- sapply(st_index_list,function(id){
  # Load the expression data
  expr.url <- paste0("data/ST_BC_GSE176078/", id)
  expr.data <- Read10X(data.dir =  expr.url,gene.column = 1)
  stref <- Seurat::CreateSeuratObject(counts = expr.data, project = id, assay = 'Spatial')
  # Load the image data
  path <- paste0("data/ST_BC_GSE176078/",id,'//spatial')
  img <- Seurat::Read10X_Image(image.dir = path)
  Seurat::DefaultAssay(object = img) <- 'Spatial'
  img <- img[colnames(x = stref)]
  stref[['image']] <- img
  #calculate the percentage of MT genes。
  stref[["percent.mt"]] <- PercentageFeatureSet(stref, pattern = "^MT-")
  #stref <- SCTransform(stref, assay = "Spatial", verbose = T)
  #plot the portion of MT and feature scores
  #VlnPlot(stref, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)
  #SpatialFeaturePlot(stref, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"),image.alpha = 0,pt.size.factor = 2.5)
  return(stref)
})
saveRDS(st_list,"data/st_list_h.rds")#save

####Step2:marker expression####
marker<-c("CD44", "CD24","ALDH1A1","PROM1","ITGA6")
i=1
for(id in st_index_list){
  st_list[[id]]<- NormalizeData(st_list[[id]],normalization.method = "LogNormalize", scale.factor = 1e4)
  SpatialFeaturePlot(st_list[[id]], features = marker,
                     image.alpha=0,ncol=3,pt.size.factor =c(2.5,2.5,1.8,1.8,1.8,1.8)[i],slot = "data",stroke=0)
  ggsave(file=paste0("plot/ST_BC_GSE176078/",id,"csc_marker_expression_normalized",".pdf"),width=12,height=15)
  SpatialFeaturePlot(st_list[[id]], features = marker,
                     image.alpha=0,ncol=3,pt.size.factor =c(2.5,2.5,1.8,1.8,1.8,1.8)[i],slot = "counts",stroke=0)
  ggsave(file=paste0("plot/ST_BC_GSE176078/",id,"csc_marker_expression",".pdf"),width=12,height=15)
  i =i+1
}

for(top_percent in c(0.03,0.05,0.08,0.10)){
  i=1
  for(id in st_index_list){
    ident_list<-c()
    st_list[[id]]@meta.data[,paste0("top_",top_percent*100,"%","_positive")]<-rep(0,length(rownames(st_list[[id]]@meta.data)))
    for(gene in marker){
      st_list[[id]]@meta.data[,gene]<-st_list[[id]]@assays$Spatial@layers$counts[rownames(st_list[[id]])==gene,]
      num_top <- ceiling(length(st_list[[id]]@meta.data[,gene]) * top_percent)
      top_indices <- order(st_list[[id]]@meta.data[,gene], decreasing = TRUE)[1:num_top]
      st_list[[id]]@meta.data[,paste0(gene,"top_",top_percent*100,"%")]<-rep(0,length(st_list[[id]]@meta.data[,gene]))
      top_indices<-top_indices[st_list[[id]]@meta.data[top_indices,gene]>0]
      st_list[[id]]@meta.data[top_indices,paste0(gene,"top_",top_percent*100,"%")]<-1
      Idents(st_list[[id]])<-paste0(gene,"top_",top_percent*100,"%")
      st_list[[id]]@meta.data[,paste0("top_",top_percent*100,"%","_positive")]<-st_list[[id]]@meta.data[,paste0("top_",top_percent*100,"%","_positive")]+st_list[[id]]@meta.data[,paste0(gene,"top_",top_percent*100,"%")]
      ident_list <- append(ident_list,paste0(gene,"top_",top_percent*100,"%"))
    }
    SpatialDimPlot(st_list[[id]],group.by = ident_list,stroke=0,
                   facet.highlight = TRUE,label = FALSE, label.size = 3,ncol=3,
                   image.alpha=0,pt.size.factor =c(2.9,2.9,2.2,2.2,2.2,2.2)[i])
    ggsave(file=paste0("plot/ST_BC_GSE176078/",id,"_",top_percent,"_expression",".pdf"),width=20,height=16)
    SpatialDimPlot(st_list[[id]],group.by = paste0("top_",top_percent*100,"%","_positive"),
                   facet.highlight = TRUE,label = FALSE, label.size = 3,ncol=1,
                   image.alpha=0,stroke=0,
                   pt.size.factor =c(2.9,2.9,2.2,2.2,2.2,2.2)[i])+scale_fill_manual(values=brewer.pal(name="YlGnBu", n = length(marker)+1))
    i=i+1
    ggsave(file=paste0("plot/ST_BC_GSE176078/",id,"_",top_percent,"_csc_niche",".pdf"),width=7,height=7)
  }
}

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
                   image.alpha=0,pt.size.factor =c(2.7,2.7,2,2,2,2)[i])
    ggsave(file=paste0("plot/ST_BC_GSE176078/",id,"normalized_",top_percent,"_expression",".pdf"),width=15,height=13)
    SpatialDimPlot(st_list[[id]],group.by = paste0("top_normalized_",top_percent*100,"%","_positive"),
                   facet.highlight = TRUE,label = FALSE, label.size = 1,ncol=1,
                   image.alpha=0,stroke=0.01,
                   pt.size.factor =c(2.7,2.7,2,2,2,2)[i])+scale_fill_manual(values=brewer.pal(name="YlGnBu", n = length(marker)+1))
    i=i+1
    ggsave(file=paste0("plot/ST_BC_GSE176078/",id,"normalized_",top_percent,"_csc_niche",".pdf"),width=5,height=5)
  }
}
ident_list<-c("CD44top_normalized_3%","CD44top_3%",
              "CD24top_normalized_3%","CD24top_3%",
              "ALDH1A1top_normalized_3%","ALDH1A1top_3%",
              "PROM1top_normalized_3%","PROM1top_3%",
              "ITGA6top_normalized_3%","ITGA6top_3%")
i=1
for(id in st_index_list){
  SpatialDimPlot(st_list[[id]],group.by = ident_list,stroke=0.01,
                 facet.highlight = TRUE,label = FALSE, label.size = 0.5,ncol=4,
                 image.alpha=0,pt.size.factor =c(2.7,2.7,2,2,2,2)[i])
  i=i+1
  ggsave(file=paste0("plot/ST_BC_GSE176078/",id,"_3_compare",".pdf"),width=20,height=13)
}
saveRDS(st_list,"data/ST_BC_GSE176078/st_list.rds")
####Step3: scoring based method#####
#histogram
for(id in st_index_list){
  for(gene in marker){
    obj<-st_list[[id]]
    cts_values <- obj@meta.data[,paste0(gene,"_normalized")] # Replace with your actual vector of values
    # Gene name (replace with your actual gene name)
    gene_name <- gene
    # Calculate percentile values
    percentiles <- c(97, 95, 92, 90) 
    percentile_values <- quantile(cts_values, probs = percentiles / 100)
    pdf(paste0("plot/ST_BC_GSE176078/",id,"normalize_",gene,"_histogram",".pdf"), width = 7, height = 5)
    # Plot histogram with annotations
    hist(cts_values, breaks = 30, col = "lightblue", main = gene_name,
         xlab = "Continuous Values", ylab = "Frequency")
    # Add percentile annotation
    abline(v = percentile_values, col = c("red", "blue", "green", "orange"), lty = 2)
    text(percentile_values, par("usr")[4], labels = paste(percentiles, "%"), col = c("red", "blue", "green", "orange"), pos = 3)
    # Add legend
    legend("topright", legend = paste(percentiles, "% Percentile"), col = c("red", "blue", "green", "orange"), lty = 2, bty = "n")
    dev.off()
  }
}
#addmodulescore
i=1
for(id in st_index_list){
  ref<-st_list[[id]]
  ref<-AddModuleScore(ref,features=marker,name="csc_score",slot="data",nbin=24)
  SpatialFeaturePlot(ref, features = "csc_score1",
                     image.alpha=0,pt.size.factor =c(2.5,2.5,1.8,1.8,1.8,1.8)[i],slot = "count")
  ggsave(file=paste0("plot/ST_BC_GSE176078/",id,"csc_add_modules_core",".pdf"),width=5,height=5)
  i =i+1
}
#sum of normalized score
i=1
for(id in st_index_list){
  vectors <- list(st_list[[id]]$CD44_normalized, st_list[[id]]$CD24_normalized, st_list[[id]]$ALDH1A1_normalized, st_list[[id]]$PROM1_normalized, st_list[[id]]$ITGA6_normalized)
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
                     image.alpha=0,pt.size.factor =c(2.5,2.5,1.8,1.8,1.8,1.8)[i],slot = "count")
  ggsave(file=paste0("plot/ST_BC_GSE176078/",id,"sum_of_csc_marker_normalized",".pdf"),width=5,height=5)
  i =i+1
}
#ucell
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
                     image.alpha=0,pt.size.factor =c(2.5,2.5,1.8,1.8,1.8,1.8)[i],slot = "count")
  ggsave(file=paste0("plot/ST_BC_GSE176078/",id,"cance_stem_cell_UCell",".pdf"),width=5,height=5)
  i =i+1
}
####Step4: deconvolution medthod--CARD####
#scRNA seq data of GSE176078, preprocessing code is in the scRNA seq folder
ref_index_list<-names(ref_list)
ref <- merge(ref_list[[ref_index_list[1]]], y = ref_list[ref_index_list[c(2:length(ref_index_list))]],add.cell.ids = ref_index_list, merge.data = TRUE)
unique(sapply(X = strsplit(colnames(ref), split = "_"), FUN = "[", 1))
remove(ref_list)
ref <- NormalizeData(ref,normalization.method = "LogNormalize", scale.factor = 1e4)
ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(ref)
#saveRDS(ref,file="data/BC_GSE180286/ref.rds")
ref <- ScaleData(ref, features = all.genes)
#saveRDS(ref,file="data/BC_GSE180286/ref.rds")
ref <- RunPCA(ref, features = VariableFeatures(object = ref))
ref <- RunHarmony(ref, group.by.vars = "orig.ident")
ref <- IntegrateLayers(object = ref, method = HarmonyIntegration,
                       orig.reduction = "pca", new.reduction = "harmony",
                       verbose = FALSE)
ref <- JoinLayers(ref)
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
  #write.csv(sc_count,"data/ST_BC_GSE176078/sc_count.csv")
  #sc_count<-ref_list[["CID44971"]]@assays$RNA@layers$counts
  #rownames(sc_count)<-rownames(ref_list[["CID44971"]])[ref_list[["CID44971"]]@assays$RNA@features@.Data[,"counts"]]
  #colnames(sc_count)<-colnames(ref_list[["CID44971"]])[ref_list[["CID44971"]]@assays$RNA@cells@.Data[,"counts"]]
  sc_count[1:4,1:4]
  sc_meta<-as.data.frame(ref@meta.data[,c("celltype","orig.ident")])
  #sc_meta<-as.data.frame(ref_list[["CID44971"]]@meta.data[,c("celltype","orig.ident")])
  colnames(sc_meta)<-c("cellType","sampleInfo")
  sc_meta[1:4,]
  #write.csv(sc_meta,"data/ST_BC_GSE176078/sc_meta.csv")
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
  write.csv(CARD_obj@Proportion_CARD,paste0("data/ST_BC_GSE176078/",id,"_","CARD_score_harmony_normalized_normalized",".csv"))
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
  ggsave(file=paste0("plot/ST_BC_GSE176078/CARD/",id,"_","all_cells","_single_harmony_normalized_normalized",".pdf"),width=15,height=13)
  ## visualize the spatial distribution of two cell types on the same plot
  p3 = CARD.visualize.prop.2CT(
    proportion = CARD_obj@Proportion_CARD,                             ### Cell type proportion estimated by CARD
    spatial_location = CARD_obj@spatial_location,                      ### spatial location information
    ct2.visualize = c("cancer stem cell","other"),              ### two cell types you want to visualize
    colors = list(c("lightblue","lightyellow","red"),c("lightblue","lightyellow","black")))       ### two color scales
  print(p3)
  ggsave(file=paste0("plot/ST_BC_GSE176078/CARD/",id,"_","all_cells","_double_harmony_normalized_normalized",".pdf"),width=12,height=12)
  p4 <- CARD.visualize.Cor(CARD_obj@Proportion_CARD,colors = NULL) # if not provide, we will use the default colors
  print(p4)
  ggsave(file=paste0("plot/ST_BC_GSE176078/CARD/",id,"_","portion_corr_harmony_normalized_normalized",".pdf"),width=10,height=10)
}
