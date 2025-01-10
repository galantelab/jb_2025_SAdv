#Map and annotate integrated dataset with reference cell types; run PCA + UMAP analyses ###############

#Load library
library(Seurat)

#Read integrated Seurat object
df = LoadSeuratRds("3_integrated_samples/integrated.rds")

#Load Seurat object from reference paper
ref = LoadSeuratRds("integrated_1.5M_2M_harmonizedOBJ.rds")

#Run ScaleData to perform dimensionality reduction by PCA in reference data
ref = ScaleData(ref, features=rownames(ref))
ref = RunPCA(ref, npcs=30, verbose=F)

#Find transfer anchors, transfer data and add predictions to metadata of integrated data
transfer_anchors = FindTransferAnchors(reference=ref, query=df, dims=1:20, reference.reduction="pca")
predictions = TransferData(anchorset=transfer_anchors, refdata=ref$FinalName, dims=1:20)
df = AddMetaData(df, metadata=predictions)

#Scale integrated data to enable running PCA and UMAP analyses
df = ScaleData(df, features=rownames(df))

#Perform dimensionality reduction by PCA
df = RunPCA(df, npcs=30, verbose=F)

#Perform dimensionality reduction by UMAP
df = RunUMAP(df, dims=1:20, verbose=F)

#Save Seurat object
SaveSeuratRds(df, "3_integrated_samples/integrated_mixed_1.5_2M_with_reductions.rds")
