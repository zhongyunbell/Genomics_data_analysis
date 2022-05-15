#####################
### May 15, 2022
### Julie Huang

'''
The analysis is using the Nov 2020 publiation <<VEGF-B promotes endocardium-derived coronary vessel development and cardiac regenration>>
https://www.ahajournals.org/doi/10.1161/CIRCULATIONAHA.120.050635

Data:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128509

- Single cell changes in postnatal P7 and in adult cardiac endothelium upon VEGF-B overexpression

'''

setwd("/Users/huangz36/Documents/Julie_learn_genomics/GSE128509_RAW/")

###### STEP 0: Restructure input files ######
# fs=list.files('/Users/huangz36/Documents/Julie_learn_genomics/GSE128509_RAW/','genes.tsv.gz')
fs=list.files('./','genes.tsv.gz')
fs
samples1=gsub('_genes.tsv.gz', '',fs)
samples1
library(stringr)
samples2=str_split(fs,'_',simplify = T)[,2]
samples2
samples2 = samples1
samples2

lapply(1:length(samples2), function(i){
  x=samples2[i]
  y=samples1[i]
  dir.create(x,recursive = T)
  file.copy(from=paste0(y,'_genes.tsv.gz'),
            to=file.path(x, 'features.tsv.gz'))
  file.copy(from=paste0(y,'_matrix.mtx.gz'),
            to=file.path(x,'matrix.mtx.gz'))
  file.copy(from=paste0(y,'_barcodes.tsv.gz'),
            to=file.path(x,'barcodes.tsv.gz'))
})

###### STEP1: read in data ######
install.packages('Seurat')
library(Seurat)
library(data.table)
dir='/Users/huangz36/Documents/Julie_learn_genomics/GSE128509_RAW/outputs'
samples=list.files( dir)
samples
sceList = lapply(samples, function(pro){
  # pro=samples[1]
  folder=file.path( dir, pro)
  print(pro)
  print(folder)
  print(list.files(folder))
  sce=CreateSeuratObject(counts = Read10X(folder),
                         project = pro, 
                         min.cels = 5, 
                         min.features = 300)
  return(sce)
})

names(sceList)

samples = gsub('^GSM[1-9]*_','',samples)
samples
names(sceList) = samples
names(sceList)

sce.all <- merge(sceList[[1]], y=sceList[ -1 ], 
                 add.cell.ids = samples)

as.data.frame(sce.all@assays$RNA@counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all@meta.data$orig.ident)


###### STEP2: reduce dimensions, cluster ######
## Since data coming from 4 samples, need to use harmony to treat samples
sce.all <- NormalizeData(sce.all, 
                     normalization.method = "LogNormalize",
                     scale.factor = 1e4)
sce.all <- FindVariableFeatures(sce.all)
sce.all <- ScaleData(sce.all)
## Run PCA
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all))
sce.all2 <- RunTSNE(sce.all, features = VariableFeatures(object = sce.all)) # add-on
DimPlot(sce.all2, reduction="tsne", label=T)

install.packages('harmony')
library(harmony)

seuratObj <- RunHarmony(sce.all, "orig.ident")
names(seuratObj@reductions)
seuratObj <- RunUMAP(seuratObj, dims = 1:15, 
                     reduction = "harmony")
DimPlot(seuratObj, reduction="umap", label=T )

sce=seuratObj
# compute nearest neighbor graph
sce <- FindNeighbors(sce, reduction = "harmony", 
                     dims = 1:15)
sce.all = sce

# set different resolutions, observe the clustering result, pick the proper one
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1)) {
  sce.all=FindClusters(sce.all, resolution = res, algorithm = 1)
}
colnames(sce.all@meta.data)
apply(sce.all@meta.data[,grep("RNA_snn", colnames(sce.all@meta.data))],2,table)

# Use clustree to visualize
install.packages("clustree")
library(clustree)

p2_tree=clustree(sce.all@meta.data,prefix="RNA_snn_res.")
ggsave(plot=p2_tree, filename="Tree_diff_resolution.pdf")
# For following analysis, use res as 0.8
sel.clust = "RNA_snn_res.0.8"
sce.all <- SetIdent(sce.all, value=sel.clust)
table(sce.all@active.ident)
#saveRDS(sce.all, "sce.all_int.rds")
