library(Seurat)
#library(SeuratData)
library(SeuratObject)
library(ggplot2)
library(patchwork)
library(dplyr)
library(RColorBrewer)
####Step1: add donor information###
Idents(total_harmony)<-"celltype"
ref<-subset(total_harmony,idents = c("CD4_T"))
remove(total_harmony)
table(ref$disease)
table(total_harmony$donor,total_harmony$disease)
saveRDS(ref,"data/CD4/ref.rds")
ref<-subset(total_harmony,idents = c("CD4_T"))
ref$donor_info<-ref$donor
sample_list<-c("1037_V4_BLD", "1037_v4_CSF", "1044_V5_BLD", "1044_V5_CSF", 
               "1070_base_BLD", "1070_base_CSF", "1071_base_BLD", "1071_base_CSF", 
               "1072_base_BLD", "1072_base_CSF", "1073_base_BLD", "1073_base_CSF", 
               "1075_base_BLD_GEX", "1075_base_CSF_GEX", "3009_V3_BLD", "3009_V3_CSF", 
               "3013_v4_BLD", "3013_v4_CSF", "3029_base_BLD", "3029_base_CSF", 
               "3031_base_BLD", "3031_base_CSF")
donor<-c(1037,1044,1070,1071,1072,1073,1075,3009,3013,3029,3031)
k<-1
i<-1
while(i < length(sample_list)){
  j <-i+1
  ref@meta.data[ref$sample==sample_list[i]|ref$sample==sample_list[j],"donor_info"]<-donor[k]
  i<-i+2
  k<-k+1
}
ref@meta.data[ref$dataset=="gse117937","donor_info"]<-ref@meta.data[ref$dataset=="ggse117937","sample"]
sample_list<-c("hiv1_bld","hiv1_csf","hiv2_bld","hiv2_csf")
donor<-c("hiv1","hiv2")
k<-1
i<-1
while(i < length(sample_list)){
  j <-i+1
  ref@meta.data[ref$sample==sample_list[i]|ref$sample==sample_list[j],"donor_info"]<-donor[k]
  i<-i+2
  k<-k+1
}

ref@meta.data[ref$dataset=="gse157829","donor_info"]<-ref@meta.data[ref$dataset=="gse157829","sample"]
ref@meta.data[ref$dataset=="gse239916","donor_info"]<-ref@meta.data[ref$dataset=="gse239916","sample"]
ref@meta.data[ref$dataset=="gse128879","donor_info"]<-ref@meta.data[ref$dataset=="gse128879","sample"]
sample<-unique(ref$donor)
sample<-sample[!is.na(sample)]
####Step2:generate tsv file####
sample<-c("1073","1075","1037","pid471","pid630","pid876","YW1_suppressed","YW5_suppressed","YW8_viremic","3013","3009","hd1","hd1_unf")
Idents(ref)<-"donor_info"
for(id in sample){
  ref_sub<-subset(ref,idents=c(id))
  filename <- paste0("data/Compass/CD4/",id,".tsv")
  # Write the data frame to a TSV file
  data<-ref_sub@assays$RNA@layers$data
  rownames(data)<-rownames(ref_sub)[ref_sub@assays$RNA@features@.Data[,"data"]]
  colnames(data)<-colnames(ref_sub)[ref_sub@assays$RNA@cells@.Data[,"data"]]
  write.table(data, file = filename, sep = "\t", quote = FALSE, row.names = T)
  filename <- paste0("data/Compass/CD4/",id,"_cell_meta.csv")
  data<-as.data.frame(ref_sub@meta.data[,"disease"])
  rownames(data)<-rownames(ref_sub@meta.data)
  colnames(data)<-c("celltype")
  write.csv(data, file = filename)
}
Idents(ref)<-"donor_info"
ref<-subset(ref,idents = c("1073","1075","1037","pid630","pid876","YW1_suppressed","YW8_viremic","3013","3009","hd1","hd1_unf","C3","C4","P1"))
####Step3: use a different way to calculate metabolism socre###
library(scMetabolism)
library(ggplot2)
library(rsvd)
ref<-sc.metabolism.Seurat(obj = ref, method = "AUCell", imputation = F, ncores = 2, metabolism.type = "KEGG")
