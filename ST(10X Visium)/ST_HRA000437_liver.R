library(Seurat)
library(SeuratData)
library(SeuratObject)
library(ggplot2)
library(patchwork)
library(dplyr)
library(RColorBrewer)
####Step 1:preprocessing####
# Load the expression data
st_index_list<-c("c1T","1T","2T","3T","4T")
st_list <- sapply(st_index_list,function(id){
  # Load the expression data
  expr.url <- paste0("data/ST_liver_HRA000437/", id)
  expr.data <- Read10X(data.dir =  expr.url )
  stref <- Seurat::CreateSeuratObject(counts = expr.data, project = 'anterior1', assay = 'Spatial')
  # Load the image data
  path <- paste0("data/ST_liver_HRA000437/",id,'//spatial')
  img <- Seurat::Read10X_Image(image.dir = path)
  Seurat::DefaultAssay(object = img) <- 'Spatial'
  img <- img[colnames(x = stref)]
  stref[['image']] <- img
  #calculate the percentage of MT genes。
  stref[["percent.mt"]] <- PercentageFeatureSet(stref, pattern = "^MT-")
  #plot the portion of MT and feature scores
  #VlnPlot(stref, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)
  #SpatialFeaturePlot(stref, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"),image.alpha = 0,pt.size.factor = 2.5)
  return(stref)
})
#saveRDS(st_list,"data/st_list_h.rds")#save
marker<-c("CD47", "EPCAM","KRT19","PROM1","SOX9")
#i=1
#for(id in st_index_list){
  #SpatialFeaturePlot(st_list[[id]], features = marker,
  #                   image.alpha=0,ncol=3,pt.size.factor =c(3.1,3.2,2.5,2.7,2.5)[i],slot = "count")
  #ggsave(file=paste0("plot/ST_liver_HRA000437/",id,"csc_marker_expression",".pdf"),width=12,height=15)
  #i =i+1
#}
####Step2:marker expression####
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
      st_list[[id]]@meta.data[top_indices,paste0(gene,"top_",top_percent*100,"%")]<-1
      Idents(st_list[[id]])<-paste0(gene,"top_",top_percent*100,"%")
      st_list[[id]]@meta.data[,paste0("top_",top_percent*100,"%","_positive")]<-st_list[[id]]@meta.data[,paste0("top_",top_percent*100,"%","_positive")]+st_list[[id]]@meta.data[,paste0(gene,"top_",top_percent*100,"%")]
      ident_list <- append(ident_list,paste0(gene,"top_",top_percent*100,"%"))
    }
    #SpatialDimPlot(st_list[[id]],group.by = ident_list,
    #               facet.highlight = TRUE,label = FALSE, label.size = 3,ncol=3,
    #               image.alpha=0,alpha=c(0.8,0.8,0.8,0.8,0.8),
    #               pt.size.factor = c(3.1,3.2,2.5,2.7,2.5)[i])
    #ggsave(file=paste0("plot/ST_liver_HRA000437/",id,"_",top_percent,"_expression",".pdf"),width=15,height=13)
    #SpatialDimPlot(st_list[[id]],group.by = paste0("top_",top_percent*100,"%","_positive"),
    #              facet.highlight = TRUE,label = FALSE, label.size = 3,ncol=1,
    #              image.alpha=0,
    #              pt.size.factor = c(3.1,3.2,2.5,2.7,2.5)[i])+scale_fill_manual(values=brewer.pal(name="YlGnBu", n = length(marker)+1))
    #i=i+1
    #ggsave(file=paste0("plot/ST_liver_HRA000437/",id,"_",top_percent,"_csc_niche",".pdf"),width=5,height=5)
  }
}
