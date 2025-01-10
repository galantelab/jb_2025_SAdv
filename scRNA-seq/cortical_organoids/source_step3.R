#Integrate all samples ###########################################################################################

#Load library
library(Seurat)

#Create directory to store integrated samples
dir.create("3_integrated_samples", showWarnings=T)

#List input files
input_files = list.files(path="2_results_norm", pattern="rds", full.names=T)

#Read each input file into Seurat object
seurat_objects=list()
for(i in 1:length(input_files)){
  seurat_objects[[i]] = LoadSeuratRds(input_files[i])
}

#Find anchors between individual samples
integration_anchors = FindIntegrationAnchors(object.list=seurat_objects, anchor.features=2000, dims = 1:20)

#Integrate samples to construct a batch-corrected expression matrix encompassing all cells
integrated_data = IntegrateData(anchorset=integration_anchors, dims = 1:20)

#Save integrated Seurat object
SaveSeuratRds(integrated_data, "3_integrated_samples/integrated.rds")
