#Integrate reference datasets (1.5 + 2 months) ######################################################################################

#Load library
library(Seurat)

#Create Seurat object with input references
seurat_objects=list()
seurat_objects[[1]] = LoadSeuratRds("1.5M_harmonizedOBJ.rds")
seurat_objects[[2]] = LoadSeuratRds("2M_harmonizedOBJ.rds")

#Find anchors between individual samples
integration_anchors = FindIntegrationAnchors(object.list=seurat_objects, anchor.features=2000, dims = 1:20)

#Integrate samples to construct a batch-corrected expression matrix encompassing all cells
integrated_data = IntegrateData(anchorset=integration_anchors, dims = 1:20)

#Save integrated Seurat object
SaveSeuratRds(integrated_data, "integrated_1.5M_2M_harmonizedOBJ.rds")

