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
sceList
sceList[[1]]
sce.all <- merge(sceList[[1]], y=sceList[ -1 ], 
                 add.cell.ids = samples)
sce.all
as.data.frame(sce.all@assays$RNA@counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all@meta.data$orig.ident)

###### STEP2: reduce dimensions, cluster ######
