#Map and annotate integrated dataset with reference cell types from paper; run PCA and UMAP analyses ##########################################

#Load library
library(Seurat)

#Read integrated Seurat object
df = LoadSeuratRds("3_integrated_samples/integrated.rds")

#Load Seurat object from reference paper and run FindVariableFeatures
ref = LoadSeuratRds("1stTRI_6-10Weeks_Thalamus-norm.rds")
ref = FindVariableFeatures(ref, selection.method="vst", nfeatures=2000)

#Find transfer anchors, transfer data and add predictions to metadata of integrated data
transfer_anchors = FindTransferAnchors(reference=ref, query=df, dims=1:20, reference.reduction="pca")
predictions = TransferData(anchorset=transfer_anchors, refdata=ref$clusters, dims=1:20)
df = AddMetaData(df, metadata=predictions)

#Scale integrated data to enable running PCA and UMAP analyses
df = ScaleData(df, features=rownames(df))

#Perform dimensionality reduction by PCA
df = RunPCA(df, npcs=30, verbose=F)

#Perform dimensionality reduction by UMAP
df = RunUMAP(df, dims=1:20, verbose=F)

#Save Seurat object
SaveSeuratRds(df, "3_integrated_samples/integrated_with_reductions.rds")
