############### step 1 cell-type annotation ###############
library(Seurat)
library(ggplot2)
library(scales)
library(harmony)


data = read.table("Combined_P18-SMK-Index_L4_811D01_RSEC_MolsPerCell.csv",sep = ",",header = T,check.names = F,stringsAsFactors = F)


rownames(data) = paste("Cell_",data$Cell_Index)
data = data[,-1]


sample_to_cell = read.table("P18-SMK-Index_L4_811D01_Sample_Tag_Calls.csv",sep = ",",header = T,check.names = F,stringsAsFactors = F)
sample_to_cell$Cell_Index = paste("Cell_",sample_to_cell$Cell_Index)
sample_to_cell = sample_to_cell[,-3]
rownames(sample_to_cell) = sample_to_cell$Cell_Index


unique(sample_to_cell$Sample_Tag)


sample_to_cell$Patient = NA
sample_to_cell$Patient[ sample_to_cell$Sample_Tag ==  "SampleTag04_hs"] = "Patient18" 
sample_to_cell$Patient[ sample_to_cell$Sample_Tag ==  "SampleTag05_hs"] = "Patient18" 
sample_to_cell$Patient[ sample_to_cell$Sample_Tag ==  "SampleTag06_hs"] = "Patient18" 


sample_to_cell$Sample_Tag[ sample_to_cell$Sample_Tag ==  "SampleTag04_hs"] = "ccRcc_AT" 
sample_to_cell$Sample_Tag[ sample_to_cell$Sample_Tag ==  "SampleTag05_hs"] = "ccRcc_CT" 
sample_to_cell$Sample_Tag[ sample_to_cell$Sample_Tag ==  "SampleTag06_hs"] = "ccRcc_TT" 


rt = CreateSeuratObject(counts =  t(data),project = "ccrcc",meta.data = sample_to_cell,min.cells = 3,min.features = 200) 


rt$percent.mt = PercentageFeatureSet(object = rt,pattern = "^MT-")
rt$percent.rp = PercentageFeatureSet(object = rt,pattern = "^RP[SL]")
rt$percent.hsp = PercentageFeatureSet(object = rt,pattern = "HSP")

VlnPlot(object = rt,features = "nCount_RNA",group.by = "Sample_Tag",pt.size = 0)
VlnPlot(object = rt,features = "nFeature_RNA",group.by = "Sample_Tag",pt.size = 0)
VlnPlot(object = rt,features = "percent.mt",group.by = "Sample_Tag",pt.size = 0)
VlnPlot(object = rt,features = "percent.rp",group.by = "Sample_Tag",pt.size = 0)
VlnPlot(object = rt,features = "percent.hsp",group.by = "Sample_Tag",pt.size = 0)


FeatureScatter(object = rt,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",group.by = "Sample_Tag")
FeatureScatter(object = rt,feature1 = "nCount_RNA",feature2 = "percent.mt",group.by = "Sample_Tag")
FeatureScatter(object = rt,feature1 = "nCount_RNA",feature2 = "percent.rp",group.by = "Sample_Tag")
FeatureScatter(object = rt,feature1 = "nCount_RNA",feature2 = "percent.hsp",group.by = "Sample_Tag")


rt = rt[,! rt$Sample_Tag %in% c("Multiplet","Undetermined")]


summary(rt$nCount_RNA)
summary(rt$nFeature_RNA)
summary(rt$percent.mt)
summary(rt$percent.rp)
summary(rt$percent.hsp)


hist(rt$nCount_RNA,breaks = 50)
hist(rt$nFeature_RNA,breaks = 50)
hist(rt$percent.mt,breaks = 50)
hist(rt$percent.rp,breaks = 50)
hist(rt$percent.hsp,breaks = 50)

ncol(rt)


median(rt$percent.mt)- 3*mad(rt$percent.mt,constant = 1)
median(rt$percent.mt)+ 3*mad(rt$percent.mt,constant = 1)

quantile(rt$nCount_RNA,probs = 0.99)
quantile(rt$nCount_RNA,probs = 0.01)
quantile(rt$nFeature_RNA,probs = 0.99)
quantile(rt$nFeature_RNA,probs = 0.01)


rt = subset(rt, subset = percent.mt < 20 & nFeature_RNA<5000 & nCount_RNA < 20000 & nCount_RNA > 1000)


hist(rt$nCount_RNA,breaks = 50)
hist(rt$nFeature_RNA,breaks = 50)
hist(rt$percent.mt,breaks = 50)
hist(rt$percent.rp,breaks = 50)
hist(rt$percent.hsp,breaks = 50)

ncol(rt)


rt = NormalizeData(object = rt,normalization.method = "LogNormalize",scale.factor = 10000)

rt = FindVariableFeatures(object = rt,selection.method = "vst",nfeatures = 3000)

rt = ScaleData(object = rt,features = rownames(rt))

rt = RunPCA(object = rt,features = VariableFeatures(rt),npcs = 100)

ElbowPlot(rt,ndims = 100)

rt = RunHarmony(object = rt,group.by.vars = "Sample_Tag",reduction = "pca",dims.use = 1:50,theta = 2,lambda = 1,max.iter.harmony = 10,verbose = T)

rt = FindNeighbors(object = rt,reduction = "harmony",dims = c(1:50),k.param = 20)
rt = FindClusters(object = rt,resolution = 1.5)

rt = RunUMAP(object = rt,reduction = "harmony",dims = c(1:50))

UMAPPlot(object = rt,label = T)
UMAPPlot(object = rt,group.by = "Sample_Tag")
UMAPPlot(object = rt,group.by = "Patient")

FeaturePlot(rt,"nCount_RNA")
FeaturePlot(rt,"nFeature_RNA")
FeaturePlot(rt,"percent.mt")
FeaturePlot(rt,"percent.rp")
FeaturePlot(rt,"percent.hsp")

options(repr.plot.height = 12, repr.plot.width = 12)
FeaturePlot(object = rt,features = c("EPCAM","PECAM1","PTPRC","ACTA2"))


FeaturePlot(object = rt,features = c("NDUFA4L2","CA9","SLC17A3"))
FeaturePlot(object = rt,features = c("CDH1","EPCAM","VIM"))


FeaturePlot(object = rt,features = c("SLC13A3","SLC34A1","SLC7A13","SLC16A9"))
FeaturePlot(object = rt,features = c("SLC22A7","SLC17A3","SLC22A8"))



FeaturePlot(object = rt,features = c("CLCNKB","ATP6V0D2","SLC4A1","SLC26A4"))

FeaturePlot(object = rt,features = c("KCJN1","SLC8A1","AVPR2","CLDN8"))

FeaturePlot(object = rt,features = c("AQP2","SLC12A1","CLDN16"))


FeaturePlot(object = rt,features = c("WT1","PODXL","PTPRO"))




FeaturePlot(object = rt,features = c("PECAM1","PTPRB","KDR","PLVAP"))
FeaturePlot(object = rt,features = c("SLC14A1","AQP1","SEMA3G","CLDN5"))


FeaturePlot(object = rt,features = c("CD14","LYZ","CD68","CD163"))
FeaturePlot(object = rt,features = c("CD1C","FCER1A","RNASE1","APOE"))
FeaturePlot(object = rt,features = c("VCAN","FCGR3A","FCN1","FN1"))



FeaturePlot(object = rt,features = c("CD3D","CD3E","CD8A","CD8B"))
FeaturePlot(object = rt,features = c("GZMA","GZMB","GZMH","GZMK"))
FeaturePlot(object = rt,features = c("GZMM","PDCD1","CTLA4","LAG3"))
FeaturePlot(object = rt,features = c("TMIGD2","TIGIT","HAVCR2"))

FeaturePlot(object = rt,features = c("CD4","CD40LG","IL7R","FOXP3"))


FeaturePlot(object = rt,features = c("KLRD1","GNLY","NCAM1","FCGR3A"))
FeaturePlot(object = rt,features = "NKG7")


FeaturePlot(object = rt,features = c("IGHM","MS4A1","IGHD","CCR7"))


FeaturePlot(object = rt,features = c("CD79A","CD79B","MZB1","JCHAIN"))
FeaturePlot(object = rt,features = c("IGHG3","IGHA2","IGKC"))


FeaturePlot(object = rt,features = c("ACTA2","DCN","PDGFRB"))

FeaturePlot(object = rt,features = c("RERGL","COX4L2","GATA3","POSTN"))

FeaturePlot(object = rt,features = c("RGS5","NOTCH3","AGTR1","HIGD1B"))

FeaturePlot(object = rt,features = c("ADH1B","FBLN1","C1S","MFAP4"))

FeaturePlot(object = rt,features = c("POSTN","COL1A1","THBS2","CDH11"))


FeaturePlot(object = rt,features = c("KIT","SLC18A2","FCER1A"))


FeaturePlot(object = rt,features = c("SELL","CXCL8", "FCGR3B", "MNDA"))

UMAPPlot(rt,label = T)

rt$All_first_cluster = NA
rt$All_first_cluster[ rt$seurat_clusters %in% c(26,11,22,23,21) ] = "Epithelial cell"
rt$All_first_cluster[ rt$seurat_clusters %in% c(12,17,8,15,1) ] = "Endothelial cell"
rt$All_first_cluster[ rt$seurat_clusters %in% c(6,10) ] = "Myeloid cell"
rt$All_first_cluster[ rt$seurat_clusters %in% c(9,5,16,3,2,14,19) ] = "T/NK cell"
rt$All_first_cluster[ rt$seurat_clusters %in% c(18,20) ] = "B/Plasma cell"
rt$All_first_cluster[ rt$seurat_clusters %in% c(13,0,27,24) ] = "Mesenchymal cell"
rt$All_first_cluster[ rt$seurat_clusters %in% c(25) ] = "Mast cell"
rt$All_first_cluster[ rt$seurat_clusters %in% c(4,7) ] = "Neutrophil"


rt$new = rt$All_first_cluster
legend_name = paste(seq(from = 1,to = length(levels(as.factor(rt$new)))),levels(as.factor(rt$new)))
rt$new = as.numeric(as.factor(rt$new))


nn = length(legend_name)
tmp = rt@reductions$umap@cell.embeddings
tmp = cbind(tmp,rt$new,rt$Sample_Tag)
colnames(tmp) = c("UMAP_1","UMAP_2","Cell_labels","Sample_Tag")
tmp = as.data.frame(tmp,stringsAsFactors = F)
tmp$UMAP_1 = as.numeric(tmp$UMAP_1)
tmp$UMAP_2 = as.numeric(tmp$UMAP_2)
tmp$Cell_labels = as.factor(tmp$Cell_labels)
p = ggplot(data = tmp)+geom_point(mapping = aes(x = UMAP_1,y=UMAP_2,color = Cell_labels),size=0.5)+scale_color_manual(values = hue_pal()(nn),labels = legend_name)+
    guides(color = guide_legend(override.aes = list(size=5)))
for(i in 1:nn){
    
    x = median(tmp[tmp$Cell_labels == i,1])
    y = median(tmp[tmp$Cell_labels == i,2])
    p = p+annotate("text",x,y,label = i) 
}

p = p+
theme_bw()+
theme(axis.text = element_blank(),
      axis.line  = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank())

options(repr.plot.height = 12, repr.plot.width = 12)
p

saveRDS(object = rt,file = "Patient18_1st.rds")

############### step 2 removal of low-quality cells and doublets ###############

library(Seurat)
library(ggplot2)
library(scales)
library(harmony)


data = readRDS("Patient18_1st.rds")

options(repr.plot.height = 12, repr.plot.width = 12)
DimPlot(object = data,reduction = "umap",group.by = "All_first_cluster",label = T)

levels(as.factor(data$All_first_cluster))


raw_matrix = data@assays$RNA@counts


name = colnames(data)[ data$All_first_cluster == "B/Plasma cell"]
one_cell_matrix = raw_matrix[,name]
one_cell_meta = data.frame(data$Sample_Tag[name])
colnames(one_cell_meta) = "Sample_Tag"

data2 = CreateSeuratObject(one_cell_matrix,min.cells = 0,min.features = 0,meta.data = one_cell_meta)


data2$percent.mt = PercentageFeatureSet(object = data2,pattern = "^MT-")
data2$percent.rp = PercentageFeatureSet(object = data2,pattern = "^RP[SL]")
data2$percent.hsp = PercentageFeatureSet(object = data2,pattern = "HSP")


data2 = NormalizeData(object = data2,normalization.method = "LogNormalize",scale.factor = 10000)

data2 = FindVariableFeatures(object = data2,selection.method = "vst",nfeatures = 2000)
data2 = ScaleData(object = data2,features = rownames(data2))
npcs = 20
data2 = RunPCA(object = data2,features = VariableFeatures(data2),npcs = npcs)
data2 = RunHarmony(object = data2,group.by.vars = "Sample_Tag",reduction = "pca",dims.use = 1:npcs,theta = 1,lambda = 1,max.iter.harmony = 2,verbose = T)

data2 = FindNeighbors(object = data2,reduction = "harmony",dims = c(1:npcs),k.param = 10)
data2 = FindClusters(object = data2,resolution = 1)

data2 = RunUMAP(object = data2,reduction = "harmony",dims = c(1:npcs))

options(repr.plot.height = 8, repr.plot.width = 8)
UMAPPlot(object = data2,label = T)
UMAPPlot(object = data2,group.by = "Sample_Tag")

options(repr.plot.height = 12, repr.plot.width = 12)

FeaturePlot(object = data2,features = c("MS4A1","IGHD","IGHM","CCR7"))


FeaturePlot(object = data2,features = c("MZB1","JCHAIN","IGHG3","IGHA2"))

options(repr.plot.height = 8, repr.plot.width = 8)
FeaturePlot(data2,c("PTPRC","PTPRB","EPCAM","ACTA2"))
FeaturePlot(data2,c("CA9","SLC13A3","PECAM1","KDR"))
FeaturePlot(data2,c("PDGFRB","DCN","RGS5","NOTCH3"))
FeaturePlot(data2,c("GATA3","RERGL","ADH1B","POSTN"))
FeaturePlot(data2,c("CD3D","CD3E","CD8A","CD8B"))
FeaturePlot(data2,c("CD14","LYZ","FCGR3B","C1QC"))
FeaturePlot(data2,c("CD4","CD40LG","KLRD1","GNLY"))
FeaturePlot(data2,c("KIT","SLC18A2","FCGR3B","SELL"))
FeaturePlot(data2,c("nCount_RNA","nFeature_RNA","percent.mt","percent.rp"))
FeaturePlot(data2,c("nCount_RNA","nFeature_RNA","percent.mt","percent.hsp"))

DEG = FindMarkers(object = data2,ident.1 = c(11),slot = "data",logfc.threshold = 0.5,min.pct = 0.1,test.use = "wilcox",only.pos = T)
FeaturePlot(data2,rownames(DEG)[1:4])
DEG

data2 = FindNeighbors(object = data2,reduction = "harmony",dims = c(1:npcs),k.param = 5)
data2 = FindClusters(object = data2,resolution = 2)

options(repr.plot.height = 8, repr.plot.width = 8)
UMAPPlot(object = data2,label = T)
UMAPPlot(object = data2,group.by = "Sample_Tag")

cell_name = colnames(data2)[ data2$seurat_clusters %in% c(0)]
data$All_first_cluster[cell_name] = "B cell"
cell_name = colnames(data2)[ data2$seurat_clusters %in% c(1,2)]
data$All_first_cluster[cell_name] = "Plasma cell"



options(repr.plot.height = 8, repr.plot.width = 8)
DimPlot(data,group.by = "All_first_cluster")


name = colnames(data)[ data$All_first_cluster %in% c("Epithelial cell","Podocyte")]
one_cell_matrix = raw_matrix[,name]
one_cell_meta = data.frame(data$Sample_Tag[name])
colnames(one_cell_meta) = "Sample_Tag"

data2 = CreateSeuratObject(one_cell_matrix,min.cells = 0,min.features = 0,meta.data = one_cell_meta)


data2$percent.mt = PercentageFeatureSet(object = data2,pattern = "^MT-")
data2$percent.rp = PercentageFeatureSet(object = data2,pattern = "^RP[SL]")
data2$percent.hsp = PercentageFeatureSet(object = data2,pattern = "HSP")


data2 = NormalizeData(object = data2,normalization.method = "LogNormalize",scale.factor = 10000)

data2 = FindVariableFeatures(object = data2,selection.method = "vst",nfeatures = 2000)
data2 = ScaleData(object = data2,features = rownames(data2))
npcs = 20
data2 = RunPCA(object = data2,features = VariableFeatures(data2),npcs = npcs)
data2 = RunHarmony(object = data2,group.by.vars = "Sample_Tag",reduction = "pca",dims.use = 1:npcs,theta = 2,lambda = 1,max.iter.harmony = 2,verbose = T)

data2 = FindNeighbors(object = data2,reduction = "harmony",dims = c(1:npcs),k.param = 10)
data2 = FindClusters(object = data2,resolution = 1)

data2 = RunUMAP(object = data2,reduction = "harmony",dims = c(1:npcs))

UMAPPlot(object = data2,label = T)
UMAPPlot(object = data2,group.by = "Sample_Tag")

options(repr.plot.height = 8, repr.plot.width = 8)
FeaturePlot(data2,c("PTPRC","PTPRB","EPCAM","ACTA2"))
FeaturePlot(data2,c("CA9","SLC13A3","PECAM1","KDR"))
FeaturePlot(data2,c("PDGFRB","DCN","RGS5","NOTCH3"))
FeaturePlot(data2,c("GATA3","RERGL","ADH1B","POSTN"))
FeaturePlot(data2,c("CD3D","CD3E","CD8A","CD8B"))
FeaturePlot(data2,c("CD14","LYZ","FCGR3B","C1QC"))
FeaturePlot(data2,c("CD4","CD40LG","KLRD1","GNLY"))
FeaturePlot(data2,c("KIT","SLC18A2","FCGR3B","SELL"))
FeaturePlot(data2,c("nCount_RNA","nFeature_RNA","percent.mt","percent.rp"))
FeaturePlot(data2,c("nCount_RNA","nFeature_RNA","percent.mt","percent.hsp"))

options(repr.plot.height = 12, repr.plot.width = 12)

FeaturePlot(object = data2,features = c("NDUFA4L2","CA9","SLC17A3"))
FeaturePlot(object = data2,features = c("CDH1","EPCAM","VIM"))


FeaturePlot(object = data2,features = c("SLC13A3","SLC34A1","SLC7A13","SLC16A9"))
FeaturePlot(object = data2,features = c("SLC22A7","SLC17A3","SLC22A8"))



FeaturePlot(object = data2,features = c("CLCNKB","ATP6V0D2","SLC4A1","SLC26A4"))

FeaturePlot(object = data2,features = c("KCJN1","SLC8A1","AVPR2","CLDN8"))

FeaturePlot(object = data2,features = c("AQP2","SLC12A1","CLDN16"))


FeaturePlot(object = data2,features = c("WT1","PODXL","PTPRO"))

DEG = FindMarkers(object = data2,ident.1 = c(11),slot = "data",logfc.threshold = 0.5,min.pct = 0.1,test.use = "wilcox",only.pos = T)
FeaturePlot(data2,rownames(DEG)[1:4])
DEG

data2 = FindNeighbors(object = data2,reduction = "harmony",dims = c(1:npcs),k.param = 5)
data2 = FindClusters(object = data2,resolution = 1.5)

DimPlot(data2,label = T)

doublets_name = colnames(data2)[ data2$seurat_clusters %in% c(15,13)]
data$All_first_cluster[doublets_name] = "Noise"
cell_name = colnames(data2)[ data2$seurat_clusters %in% c(12,10,3,2,8,16,4,6,9,11)]
data$All_first_cluster[cell_name] = "Cancer cell"
cell_name = colnames(data2)[ data2$seurat_clusters %in% c(0,1,14,7,5)]
data$All_first_cluster[cell_name] = "Epithelial cell"
cell_name = colnames(data2)[ data2$seurat_clusters %in% c(5,8)]
data$All_first_cluster[cell_name] = "Podocyte"

options(repr.plot.height = 12, repr.plot.width = 12)
DimPlot(data,group.by = "All_first_cluster",label = T)


name = colnames(data)[ data$All_first_cluster == "T/NK cell"]
one_cell_matrix = raw_matrix[,name]
one_cell_meta = data.frame(data$Sample_Tag[name])
colnames(one_cell_meta) = "Sample_Tag"

data2 = CreateSeuratObject(one_cell_matrix,min.cells = 0,min.features = 0,meta.data = one_cell_meta)


data2$percent.mt = PercentageFeatureSet(object = data2,pattern = "^MT-")
data2$percent.rp = PercentageFeatureSet(object = data2,pattern = "^RP[SL]")
data2$percent.hsp = PercentageFeatureSet(object = data2,pattern = "HSP")


data2 = NormalizeData(object = data2,normalization.method = "LogNormalize",scale.factor = 10000)

data2 = FindVariableFeatures(object = data2,selection.method = "vst",nfeatures = 2000)
data2 = ScaleData(object = data2,features = rownames(data2))
npcs = 30
data2 = RunPCA(object = data2,features = VariableFeatures(data2),npcs = npcs)
data2 = RunHarmony(object = data2,group.by.vars = "Sample_Tag",reduction = "pca",dims.use = 1:npcs,theta = 2,lambda = 1,max.iter.harmony = 3,verbose = T)

data2 = FindNeighbors(object = data2,reduction = "harmony",dims = c(1:npcs),k.param = 20)
data2 = FindClusters(object = data2,resolution = 1)

data2 = RunUMAP(object = data2,reduction = "harmony",dims = c(1:npcs))

options(repr.plot.height = 12, repr.plot.width = 12)
UMAPPlot(object = data2,label = T)
UMAPPlot(object = data2,group.by = "Sample_Tag")

options(repr.plot.height = 12, repr.plot.width = 12)

FeaturePlot(object = data2,features = c("CD3D","CD3E","CD8A","CD8B"))
FeaturePlot(object = data2,features = c("GZMA","GZMB","GZMH","GZMK"))
FeaturePlot(object = data2,features = c("GZMM","PDCD1","CTLA4","LAG3"))
FeaturePlot(object = data2,features = c("TMIGD2","TIGIT","HAVCR2"))

FeaturePlot(object = data2,features = c("CD4","CD40LG","IL7R","FOXP3"))


FeaturePlot(object = data2,features = c("KLRD1","GNLY","NCAM1","FCGR3A"))
FeaturePlot(object = data2,features = "NKG7")
FeaturePlot(object = data2,features = "ZNF683")

options(repr.plot.height = 8, repr.plot.width = 8)
FeaturePlot(data2,c("PTPRC","PTPRB","EPCAM","ACTA2"))
FeaturePlot(data2,c("CA9","SLC13A3","PECAM1","KDR"))
FeaturePlot(data2,c("PDGFRB","DCN","RGS5","NOTCH3"))
FeaturePlot(data2,c("GATA3","RERGL","ADH1B","POSTN"))
FeaturePlot(data2,c("CD3D","CD3E","CD8A","CD8B"))
FeaturePlot(data2,c("CD14","LYZ","FCGR3B","C1QC"))
FeaturePlot(data2,c("CD4","CD40LG","KLRD1","GNLY"))
FeaturePlot(data2,c("KIT","SLC18A2","FCGR3B","SELL"))
FeaturePlot(data2,c("nCount_RNA","nFeature_RNA","percent.mt","percent.rp"))
FeaturePlot(data2,c("nCount_RNA","nFeature_RNA","percent.mt","percent.hsp"))

data2 = FindNeighbors(object = data2,reduction = "harmony",dims = c(1:npcs),k.param = 10)
data2 = FindClusters(object = data2,resolution = 2)

DimPlot(data2,label = T)

DEG = FindMarkers(object = data2,ident.1 = c(16),slot = "data",logfc.threshold = 0.5,min.pct = 0.1,test.use = "wilcox",only.pos = T)
FeaturePlot(data2,rownames(DEG)[1:4])
DEG

library(scHCL)
hcl_result <- scHCL(scdata = data2[,data2$seurat_clusters %in% c(16)]@assays$RNA@counts, numbers_plot = 3)

sort(table(hcl_result$scHCL),decreasing = T)

doublets_name = colnames(data2)[ data2$seurat_clusters %in% c(20,19)]
data$All_first_cluster[doublets_name] = "Noise"
cell_name = colnames(data2)[ data2$seurat_clusters %in% c(14)]
data$All_first_cluster[cell_name] = "CD4+ Treg cell"
cell_name = colnames(data2)[ data2$seurat_clusters %in% c(11,9,5,8,10,18)]
data$All_first_cluster[cell_name] = "CD4+ T cell"
cell_name = colnames(data2)[ data2$seurat_clusters %in% c(13,15)]
data$All_first_cluster[cell_name] = "Cycling CD8+ T cell"


cell_name = colnames(data2)[ data2$seurat_clusters %in% c(3,4,7,17)]
data$All_first_cluster[cell_name] = "Exhausted CD8+ T cell"
cell_name = colnames(data2)[ data2$seurat_clusters %in% c(0,16)]
data$All_first_cluster[cell_name] = "NKT cell"
cell_name = colnames(data2)[ data2$seurat_clusters %in% c(12,6,1,2)]
data$All_first_cluster[cell_name] = "NK cell"

DimPlot(data,group.by = "All_first_cluster",label = T)


name = colnames(data)[ data$All_first_cluster == "Endothelial cell"]
one_cell_matrix = raw_matrix[,name]
one_cell_meta = data@meta.data[name,]

data2 = CreateSeuratObject(one_cell_matrix,min.cells = 0,min.features = 0,meta.data = one_cell_meta)


data2$percent.mt = PercentageFeatureSet(object = data2,pattern = "^MT-")
data2$percent.rp = PercentageFeatureSet(object = data2,pattern = "^RP[SL]")
data2$percent.hsp = PercentageFeatureSet(object = data2,pattern = "HSP")


data2 = NormalizeData(object = data2,normalization.method = "LogNormalize",scale.factor = 10000)

data2 = FindVariableFeatures(object = data2,selection.method = "vst",nfeatures = 2000)
data2 = ScaleData(object = data2,features = rownames(data2))
npcs = 30
data2 = RunPCA(object = data2,features = VariableFeatures(data2),npcs = npcs)
data2 = RunHarmony(object = data2,group.by.vars = "Sample_Tag",reduction = "pca",dims.use = 1:npcs,theta = 2,lambda = 1,max.iter.harmony = 2,verbose = T)

data2 = FindNeighbors(object = data2,reduction = "harmony",dims = c(1:npcs),k.param = 20)
data2 = FindClusters(object = data2,resolution = 1)

data2 = RunUMAP(object = data2,reduction = "harmony",dims = c(1:npcs))

UMAPPlot(object = data2,label = T)
UMAPPlot(object = data2,group.by = "Sample_Tag")





FeaturePlot(object = data2,features = c("PECAM1","PTPRB","KDR","PLVAP"))
FeaturePlot(object = data2,features = c("SLC14A1","AQP1","SEMA3G","CLDN5"))

options(repr.plot.height = 8, repr.plot.width = 8)
FeaturePlot(data2,c("PTPRC","PTPRB","EPCAM","ACTA2"))
FeaturePlot(data2,c("CA9","SLC13A3","PECAM1","KDR"))
FeaturePlot(data2,c("PDGFRB","DCN","RGS5","NOTCH3"))
FeaturePlot(data2,c("GATA3","RERGL","ADH1B","POSTN"))
FeaturePlot(data2,c("CD3D","CD3E","CD8A","CD8B"))
FeaturePlot(data2,c("CD14","LYZ","FCGR3B","C1QC"))
FeaturePlot(data2,c("CD4","CD40LG","KLRD1","GNLY"))
FeaturePlot(data2,c("KIT","SLC18A2","FCGR3B","SELL"))
FeaturePlot(data2,c("nCount_RNA","nFeature_RNA","percent.mt","percent.rp"))
FeaturePlot(data2,c("nCount_RNA","nFeature_RNA","percent.mt","percent.hsp"))

DEG = FindMarkers(object = data2,ident.1 = c(7),slot = "data",logfc.threshold = 0.5,min.pct = 0.1,test.use = "wilcox",only.pos = T)
FeaturePlot(data2,rownames(DEG)[1:4])
DEG

data2 = FindNeighbors(object = data2,reduction = "harmony",dims = c(1:npcs),k.param = 10)
data2 = FindClusters(object = data2,resolution = 0.1)

UMAPPlot(object = data2,label = T)
UMAPPlot(object = data2,group.by = "Sample_Tag")

doublets_name = colnames(data2)[ data2$seurat_clusters %in% c(8)]
data$All_first_cluster[doublets_name] = "Noise"

DimPlot(data,group.by = "All_first_cluster")


name = colnames(data)[ data$All_first_cluster %in% c("Mast cell","Neutrophil","Myeloid cell")]
one_cell_matrix = raw_matrix[,name]
one_cell_meta = data@meta.data[name,]

data2 = CreateSeuratObject(one_cell_matrix,min.cells = 0,min.features = 0,meta.data = one_cell_meta)


data2$percent.mt = PercentageFeatureSet(object = data2,pattern = "^MT-")
data2$percent.rp = PercentageFeatureSet(object = data2,pattern = "^RP[SL]")
data2$percent.hsp = PercentageFeatureSet(object = data2,pattern = "HSP")


data2 = NormalizeData(object = data2,normalization.method = "LogNormalize",scale.factor = 10000)

data2 = FindVariableFeatures(object = data2,selection.method = "vst",nfeatures = 2000)
data2 = ScaleData(object = data2,features = rownames(data2))
npcs = 30
data2 = RunPCA(object = data2,features = VariableFeatures(data2),npcs = npcs)
data2 = RunHarmony(object = data2,group.by.vars = "Sample_Tag",reduction = "pca",dims.use = 1:npcs,theta = 2,lambda = 1,max.iter.harmony = 2,verbose = T)

data2 = FindNeighbors(object = data2,reduction = "harmony",dims = c(1:npcs),k.param = 20)
data2 = FindClusters(object = data2,resolution = 1)

data2 = RunUMAP(object = data2,reduction = "harmony",dims = c(1:npcs))

UMAPPlot(object = data2,label = T)
UMAPPlot(object = data2,group.by = "Sample_Tag")

options(repr.plot.height = 12, repr.plot.width = 12)

FeaturePlot(object = data2,features = c("CD14","LYZ","CD68","CD163"))
FeaturePlot(object = data2,features = c("CD1C","FCER1A","RNASE1","APOE"))
FeaturePlot(object = data2,features = c("VCAN","FCGR3A","FCN1","FN1"))


FeaturePlot(object = data2,features = c("KIT","SLC18A2","FCER1A"))


FeaturePlot(object = data2,features = c("SELL","CXCL8", "FCGR3B", "MNDA"))

options(repr.plot.height = 8, repr.plot.width = 8)
FeaturePlot(data2,c("PTPRC","PTPRB","EPCAM","ACTA2"))
FeaturePlot(data2,c("CA9","SLC13A3","PECAM1","KDR"))
FeaturePlot(data2,c("PDGFRB","DCN","RGS5","NOTCH3"))
FeaturePlot(data2,c("GATA3","RERGL","ADH1B","POSTN"))
FeaturePlot(data2,c("CD3D","CD3E","CD8A","CD8B"))
FeaturePlot(data2,c("CD14","LYZ","FCGR3B","C1QC"))
FeaturePlot(data2,c("CD4","CD40LG","KLRD1","GNLY"))
FeaturePlot(data2,c("KIT","SLC18A2","FCGR3B","SELL"))
FeaturePlot(data2,c("nCount_RNA","nFeature_RNA","percent.mt","percent.rp"))
FeaturePlot(data2,c("nCount_RNA","nFeature_RNA","percent.mt","percent.hsp"))

DEG = FindMarkers(object = data2,ident.1 = c(8),slot = "data",logfc.threshold = 0.5,min.pct = 0.1,test.use = "wilcox",only.pos = T)
FeaturePlot(data2,rownames(DEG)[1:4])
DEG

data2 = FindNeighbors(object = data2,reduction = "harmony",dims = c(1:npcs),k.param = 10)
data2 = FindClusters(object = data2,resolution = 1)

DimPlot(data2,label = T)

re_name = colnames(data2)[ data2$seurat_clusters %in% c(15,2,9,8,4,13)]
data$All_first_cluster[re_name] = "Macrophage"
re_name = colnames(data2)[ data2$seurat_clusters %in% c(7)]
data$All_first_cluster[re_name] = "Dendritic cell"
re_name = colnames(data2)[ data2$seurat_clusters %in% c(3,1,0,6,10)]
data$All_first_cluster[re_name] = "Neutrophil"
re_name = colnames(data2)[ data2$seurat_clusters %in% c(14)]
data$All_first_cluster[re_name] = "Mast cell"
re_name = colnames(data2)[ data2$seurat_clusters %in% c(5,11)]
data$All_first_cluster[re_name] = "Monocyte"
doublets_name = colnames(data2)[ data2$seurat_clusters %in% c(12)]
data$All_first_cluster[doublets_name] = "Noise"

DimPlot(data,group.by = "All_first_cluster",label = T)


name = colnames(data)[ data$All_first_cluster %in% c("Dendritic cell","Macrophage","Monocyte")]
one_cell_matrix = raw_matrix[,name]
one_cell_meta = data@meta.data[name,]

data2 = CreateSeuratObject(one_cell_matrix,min.cells = 0,min.features = 0,meta.data = one_cell_meta)


data2$percent.mt = PercentageFeatureSet(object = data2,pattern = "^MT-")
data2$percent.rp = PercentageFeatureSet(object = data2,pattern = "^RP[SL]")
data2$percent.hsp = PercentageFeatureSet(object = data2,pattern = "HSP")


data2 = NormalizeData(object = data2,normalization.method = "LogNormalize",scale.factor = 10000)

data2 = FindVariableFeatures(object = data2,selection.method = "vst",nfeatures = 2000)
data2 = ScaleData(object = data2,features = rownames(data2))
npcs = 30
data2 = RunPCA(object = data2,features = VariableFeatures(data2),npcs = npcs)
data2 = RunHarmony(object = data2,group.by.vars = "Sample_Tag",reduction = "pca",dims.use = 1:npcs,theta = 2,lambda = 1,max.iter.harmony = 2,verbose = T)

data2 = FindNeighbors(object = data2,reduction = "harmony",dims = c(1:npcs),k.param = 20)
data2 = FindClusters(object = data2,resolution = 1)

data2 = RunUMAP(object = data2,reduction = "harmony",dims = c(1:npcs))

UMAPPlot(object = data2,label = T)
UMAPPlot(object = data2,group.by = "Sample_Tag")

options(repr.plot.height = 12, repr.plot.width = 12)

FeaturePlot(object = data2,features = c("CD14","LYZ","CD68","CD163"))
FeaturePlot(object = data2,features = c("CD1C","FCER1A","RNASE1","APOE"))
FeaturePlot(object = data2,features = c("VCAN","FCGR3A","FCN1","FN1"))


FeaturePlot(object = data2,features = c("KIT","SLC18A2","FCER1A"))


FeaturePlot(object = data2,features = c("SELL","CXCL8", "FCGR3B", "MNDA"))

options(repr.plot.height = 8, repr.plot.width = 8)
FeaturePlot(data2,c("PTPRC","PTPRB","EPCAM","ACTA2"))
FeaturePlot(data2,c("CA9","SLC13A3","PECAM1","KDR"))
FeaturePlot(data2,c("PDGFRB","DCN","RGS5","NOTCH3"))
FeaturePlot(data2,c("GATA3","RERGL","ADH1B","POSTN"))
FeaturePlot(data2,c("CD3D","CD3E","CD8A","CD8B"))
FeaturePlot(data2,c("CD14","LYZ","FCGR3B","C1QC"))
FeaturePlot(data2,c("CD4","CD40LG","KLRD1","GNLY"))
FeaturePlot(data2,c("KIT","SLC18A2","FCGR3B","SELL"))
FeaturePlot(data2,c("nCount_RNA","nFeature_RNA","percent.mt","percent.rp"))
FeaturePlot(data2,c("nCount_RNA","nFeature_RNA","percent.mt","percent.hsp"))

DEG = FindMarkers(object = data2,ident.1 = c(10),slot = "data",logfc.threshold = 0.5,min.pct = 0.1,test.use = "wilcox",only.pos = T)
FeaturePlot(data2,rownames(DEG)[1:4])
DEG

data2 = FindNeighbors(object = data2,reduction = "harmony",dims = c(1:npcs),k.param = 3)
data2 = FindClusters(object = data2,resolution = 2)

DimPlot(data2,label = T)
DimPlot(data2,label = T,group.by = "All_first_cluster")







re_name = colnames(data2)[ data2$seurat_clusters %in% c(22)]
data$All_first_cluster[re_name] = "Noise"

DimPlot(data,group.by = "All_first_cluster",label = T)


name = colnames(data)[ data$All_first_cluster %in% c("Mesenchymal cell")]
one_cell_matrix = raw_matrix[,name]
one_cell_meta = data@meta.data[name,]

data2 = CreateSeuratObject(one_cell_matrix,min.cells = 0,min.features = 0,meta.data = one_cell_meta)


data2$percent.mt = PercentageFeatureSet(object = data2,pattern = "^MT-")
data2$percent.rp = PercentageFeatureSet(object = data2,pattern = "^RP[SL]")
data2$percent.hsp = PercentageFeatureSet(object = data2,pattern = "HSP")


data2 = NormalizeData(object = data2,normalization.method = "LogNormalize",scale.factor = 10000)

data2 = FindVariableFeatures(object = data2,selection.method = "vst",nfeatures = 2000)
data2 = ScaleData(object = data2,features = rownames(data2))
npcs = 20
data2 = RunPCA(object = data2,features = VariableFeatures(data2),npcs = npcs)
data2 = RunHarmony(object = data2,group.by.vars = "Sample_Tag",reduction = "pca",dims.use = 1:npcs,theta = 2,lambda = 1,max.iter.harmony = 2,verbose = T)

data2 = FindNeighbors(object = data2,reduction = "harmony",dims = c(1:npcs),k.param = 20)
data2 = FindClusters(object = data2,resolution = 1)

data2 = RunUMAP(object = data2,reduction = "harmony",dims = c(1:npcs))

UMAPPlot(object = data2,label = T)
UMAPPlot(object = data2,group.by = "Sample_Tag")

options(repr.plot.height = 12, repr.plot.width = 12)

FeaturePlot(object = data2,features = c("ACTA2","DCN","PDGFRB"))

FeaturePlot(object = data2,features = c("RERGL","COX4L2","GATA3","POSTN"))

FeaturePlot(object = data2,features = c("RGS5","NOTCH3","AGTR1","HIGD1B"))

FeaturePlot(object = data2,features = c("ADH1B","FBLN1","C1S","MFAP4"))

FeaturePlot(object = data2,features = c("POSTN","COL1A1","THBS2","CDH11"))

options(repr.plot.height = 8, repr.plot.width = 8)
FeaturePlot(data2,c("PTPRC","PTPRB","EPCAM","ACTA2"))
FeaturePlot(data2,c("CA9","SLC13A3","PECAM1","KDR"))
FeaturePlot(data2,c("PDGFRB","DCN","RGS5","NOTCH3"))
FeaturePlot(data2,c("GATA3","RERGL","ADH1B","POSTN"))
FeaturePlot(data2,c("CD3D","CD3E","CD8A","CD8B"))
FeaturePlot(data2,c("CD14","LYZ","FCGR3B","C1QC"))
FeaturePlot(data2,c("CD4","CD40LG","KLRD1","GNLY"))
FeaturePlot(data2,c("KIT","SLC18A2","FCGR3B","SELL"))
FeaturePlot(data2,c("nCount_RNA","nFeature_RNA","percent.mt","percent.rp"))
FeaturePlot(data2,c("nCount_RNA","nFeature_RNA","percent.mt","percent.hsp"))

DEG = FindMarkers(object = data2,ident.1 = c(10),slot = "data",logfc.threshold = 0.5,min.pct = 0.1,test.use = "wilcox",only.pos = T)
FeaturePlot(data2,rownames(DEG)[1:4])
DEG

data2 = FindNeighbors(object = data2,reduction = "harmony",dims = c(1:npcs),k.param = 5)
data2 = FindClusters(object = data2,resolution = 1)

DimPlot(data2,label = T)

re_name = colnames(data2)[ data2$seurat_clusters %in% c(16,14)]
data$All_first_cluster[re_name] = "Noise"

DimPlot(data,group.by = "All_first_cluster",label = T)

table(data$All_first_cluster)


name = colnames(data)[  !data$All_first_cluster %in% c("Noise")]
all_matrix = raw_matrix[,name]
all_meta = data@meta.data[name,]


data2 = CreateSeuratObject(all_matrix,min.cells = 0,min.features = 0,meta.data = all_meta)


data2$percent.mt = PercentageFeatureSet(object = data2,pattern = "^MT-")
data2$percent.rp = PercentageFeatureSet(object = data2,pattern = "^RP[SL]")
data2$percent.hsp = PercentageFeatureSet(object = data2,pattern = "HSP")


data2 = NormalizeData(object = data2,normalization.method = "LogNormalize",scale.factor = 10000)

data2 = FindVariableFeatures(object = data2,selection.method = "vst",nfeatures = 3000)
data2 = ScaleData(object = data2,features = rownames(data2))
npcs = 50
data2 = RunPCA(object = data2,features = VariableFeatures(data2),npcs = npcs)
data2 = RunHarmony(object = data2,group.by.vars = "Sample_Tag",reduction = "pca",dims.use = 1:npcs,theta = 2,lambda = 1,max.iter.harmony = 2,verbose = T)

data2 = RunUMAP(object = data2,reduction = "harmony",dims = c(1:npcs))

options(repr.plot.height = 12, repr.plot.width = 12)
UMAPPlot(object = data2,label = T,group.by = "All_first_cluster",repel = T)
UMAPPlot(object = data2,group.by = "Sample_Tag")

saveRDS(data,"Patient18_2ed.rds")

############### step 3 discrimination of normal and malignant epithelial cells ###############
library(Seurat)
library(harmony)

data = readRDS("Patient18_2ed.rds")

options(repr.plot.height = 12, repr.plot.width = 12)
UMAPPlot(object = data,label = T,group.by = "All_first_cluster",repel = T)

cell_type = unique(data$All_first_cluster)

cell_type = cell_type[!cell_type %in% c("Cancer cell","Epithelial cell","Noise")]

all_cell_names = colnames(data)
subset_cells1 = vector()
for(i in cell_type){
    tmp_cell_names = all_cell_names[ data$All_first_cluster == i ]
    if(length(tmp_cell_names) <= 50){
        subset_cells1 = c(subset_cells1,tmp_cell_names)
    }else if(length(tmp_cell_names) > 50){
        subset_cells1 = c(subset_cells1,sample(x = tmp_cell_names,size = 50,replace = FALSE))
    }
}

subset_cells2 = all_cell_names[data$All_first_cluster%in% c("Cancer cell","Epithelial cell")]

subset_cells = c(subset_cells1,subset_cells2)

options(repr.plot.height = 12, repr.plot.width = 12)
UMAPPlot(object = data[,subset_cells],label = T,group.by = "All_first_cluster",repel = T)

data = data[,subset_cells]

data$infercnv = data$All_first_cluster

data$infercnv[ ! data$infercnv %in% c("Cancer cell","Epithelial cell")] = "Reference cell"

UMAPPlot(object = data[,subset_cells],label = T,group.by = "infercnv",repel = T)


name = colnames(data)[  !data$infercnv %in% c("Reference cell")]
raw_matrix = data@assays$RNA@counts
all_matrix = raw_matrix[,name]
all_meta = data@meta.data[name,]


data2 = CreateSeuratObject(all_matrix,min.cells = 0,min.features = 0,meta.data = all_meta)


data2$percent.mt = PercentageFeatureSet(object = data2,pattern = "^MT-")
data2$percent.rp = PercentageFeatureSet(object = data2,pattern = "^RP[SL]")
data2$percent.hsp = PercentageFeatureSet(object = data2,pattern = "HSP")


data2 = NormalizeData(object = data2,normalization.method = "LogNormalize",scale.factor = 10000)

data2 = FindVariableFeatures(object = data2,selection.method = "vst",nfeatures = 2000)
data2 = ScaleData(object = data2,features = rownames(data2))
npcs = 30
data2 = RunPCA(object = data2,features = VariableFeatures(data2),npcs = npcs)
data2 = RunHarmony(object = data2,group.by.vars = "Sample_Tag",reduction = "pca",dims.use = 1:npcs,theta = 2,lambda = 1,max.iter.harmony = 5,verbose = T)

data2 = FindNeighbors(object = data2,reduction = "harmony",dims = c(1:npcs),k.param = 20)
data2 = FindClusters(object = data2,resolution = 1)

data2 = RunUMAP(object = data2,reduction = "harmony",dims = c(1:npcs))

data2 = FindNeighbors(object = data2,reduction = "harmony",dims = c(1:npcs),k.param = 5)
data2 = FindClusters(object = data2,resolution = 3)

UMAPPlot(object = data2,label = T)
UMAPPlot(object = data2,group.by = "Sample_Tag")
UMAPPlot(object = data2,group.by = "All_first_cluster")

options(repr.plot.height = 12, repr.plot.width = 12)

FeaturePlot(object = data2,features = c("NDUFA4L2","CA9","SLC17A3"))
FeaturePlot(object = data2,features = c("CDH1","EPCAM","VIM"))


FeaturePlot(object = data2,features = c("SLC13A3","SLC34A1","SLC7A13","SLC16A9"))
FeaturePlot(object = data2,features = c("SLC22A7","SLC17A3","SLC22A8"))



FeaturePlot(object = data2,features = c("CLCNKB","ATP6V0D2","SLC4A1","SLC26A4"))

FeaturePlot(object = data2,features = c("KCJN1","SLC8A1","AVPR2","CLDN8"))

FeaturePlot(object = data2,features = c("AQP2","SLC12A1","CLDN16"))


FeaturePlot(object = data2,features = c("WT1","PODXL","PTPRO"))

options(repr.plot.height = 8, repr.plot.width = 8)
FeaturePlot(data2,c("PTPRC","PTPRB","EPCAM","ACTA2"))
FeaturePlot(data2,c("CA9","SLC13A3","PECAM1","KDR"))
FeaturePlot(data2,c("PDGFRB","DCN","RGS5","NOTCH3"))
FeaturePlot(data2,c("GATA3","RERGL","ADH1B","POSTN"))
FeaturePlot(data2,c("CD3D","CD3E","CD8A","CD8B"))
FeaturePlot(data2,c("CD14","LYZ","FCGR3B","C1QC"))
FeaturePlot(data2,c("CD4","CD40LG","KLRD1","GNLY"))
FeaturePlot(data2,c("KIT","SLC18A2","FCGR3B","SELL"))
FeaturePlot(data2,c("nCount_RNA","nFeature_RNA","percent.mt","percent.rp"))
FeaturePlot(data2,c("nCount_RNA","nFeature_RNA","percent.mt","percent.hsp"))

UMAPPlot(object = data2,label = T)
















UMAPPlot(object = data,label = T,group.by = "infercnv",repel = T)

saveRDS(data,"P18_infercnv_data.rds")

system("Rscript infercnv.R ./P18_infercnv_data.rds P18")

CNV_observation = read.table("./infercnv_result/P18/infercnv.observations.txt",header = T,row.names = 1,check.names = F,stringsAsFactors = F,sep = " ")

meta = data@meta.data[,c("infercnv","Sample_Tag")]

pheatmap::pheatmap(CNV_observation,cluster_rows = F,cluster_cols = T,show_rownames = F,show_colnames = F,annotation_col = meta,cutree_cols = 4,clustering_method = "ward.D")

h = hclust(dist(t(CNV_observation)),method = "ward.D")

ch = cutree(h,k = 4)

ch = sort(ch)

pheatmap::pheatmap(CNV_observation[names(ch)],cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,annotation_col = meta,gaps_col = cumsum(table(ch))[1:3])

rt = readRDS("Patient18_2ed.rds")

rt$All_first_cluster[names(ch[ch %in% c(1,2,3)])] = "Cancer cell"
rt$All_first_cluster[names(ch[ch %in% c(4)])] = "Epithelial cell"

data2$All_first_cluster[names(ch[ch %in% c(1,2,3)])] = "Cancer cell"
data2$All_first_cluster[names(ch[ch %in% c(4)])] = "Epithelial cell"

DimPlot(data2,group.by = "All_first_cluster")

saveRDS(rt,"Patient18_3rd.rds")
