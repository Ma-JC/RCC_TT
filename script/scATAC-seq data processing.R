library(Signac)
library(Seurat)
library(SeuratWrappers)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86) #---GRCh38 (hg38)
library(patchwork)
library(GenomicRanges)
library(ggpubr)
library(harmony)
library(clustree)

library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)
library(TFBSTools)

library(chromVARmotifs)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

library(slingshot)
library(SummarizedExperiment)
library(destiny)
library(Biobase)

library(RColorBrewer)

library(dplyr)
suppressPackageStartupMessages(library(PseudotimeDE))
suppressPackageStartupMessages(library(scales))
library(parallel)

######## scRNA-seq ############
plot_genebypseudotime = function(Pseudotime,gene_expression,gene_name,Clusters,legend.position=c(0.8,0.8),col2){
    pseu_matrix = cbind(Pseudotime,gene_expression)
    pseu_matrix = as.data.frame(pseu_matrix)
    colnames(pseu_matrix) = c("Pseudotime","Gene_expression")
    pseu_matrix$Clusters = Clusters
    pseu_matrix = pseu_matrix[ order(pseu_matrix$Pseudotime),]
    
    ggplot(pseu_matrix,aes(x = Pseudotime,y=Gene_expression))+
    scale_color_manual(values = col2)+
    geom_point(aes(color=Clusters))+
    geom_smooth(method = "gam")+
    ylab(gene_name)+
    theme_bw()+
    theme(panel.grid.major = element_line(linewidth=0.5,linetype = "dashed",colour = "grey60"),
      legend.margin = margin(16, 16, 16, 16),
      legend.title = element_text(size=16,face="bold"),
      legend.text = element_text(size=14,face="bold"),
      legend.position = legend.position,
      legend.background = element_rect(fill = "lightgreen",linewidth = 0.5,color="black"),
      axis.title = element_text(size=18,face = "bold"),
      axis.text = element_text(size=14),
      panel.border = element_rect(linewidth = 1.5),
     panel.grid.minor = element_line(linewidth=0.5,linetype = "dashed",colour = "grey80"),
     
     
     )
}

scRNA <- readRDS("/home/ncpsb/new_mnt/ccRCC_data/Validation/scATAC_seq_ccRCC/All_scRNA.rds")

raw_matrix = scRNA@assays$RNA@counts

name = colnames(scRNA)[ scRNA$All_first_cluster  %in% c("Fibroblast","Pericyte")]
one_cell_matrix = raw_matrix[,name]
one_cell_meta = scRNA@meta.data[name,]

data2 = CreateSeuratObject(one_cell_matrix,min.cells = 0,min.features = 0,meta.data = one_cell_meta)

data2$percent.mt = PercentageFeatureSet(object = data2,pattern = "^MT-")
data2$percent.rp = PercentageFeatureSet(object = data2,pattern = "^RP[SL]")
data2$percent.hsp = PercentageFeatureSet(object = data2,pattern = "HSP")


data2 = NormalizeData(object = data2,normalization.method = "LogNormalize",scale.factor = 10000)
data2 = FindVariableFeatures(object = data2,selection.method = "vst",nfeatures = 3000)
data2 = ScaleData(object = data2,features = rownames(data2))
npcs = 30
data2 = RunPCA(object = data2,features = VariableFeatures(data2),npcs = npcs)

data2 = FindNeighbors(object = data2,reduction = "pca",dims = c(1:npcs),k.param = 10)
data2 = FindClusters(object = data2,resolution = 1)

data2 = RunUMAP(object = data2,reduction = "pca",dims = c(1:npcs))
data2$All_second_cluster = data2$All_first_cluster

options(repr.plot.height = 8, repr.plot.width = 8)
UMAPPlot(object = data2,label = T)

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
FeaturePlot(object = data2,features = c("THBS2","FAP","CD36","S100A4"))
FeaturePlot(object = data2,features = c("FBLN5","MFAP5","ADH1B","CYSLTR2"))
FeaturePlot(object = data2,features = c("CREB3L1","TBX2","MAF","FOXO1"))

re_name = colnames(data2)[ data2$seurat_clusters %in% c(8,9)]
data2$All_second_cluster[re_name] = "FAP+ Fibroblast"
re_name = colnames(data2)[ data2$seurat_clusters %in% c(5)]
data2$All_second_cluster[re_name] = "CYSLTR2+ Fibroblast"
re_name = colnames(data2)[ data2$seurat_clusters %in% c(0,1,4,7)]
data2$All_second_cluster[re_name] = "CD36+ Pericyte"
re_name = colnames(data2)[ data2$seurat_clusters %in% c(2,3,6,10)]
data2$All_second_cluster[re_name] = "S100A4+ Pericyte"

scRNA_Mesen = data2

options(repr.plot.height = 8, repr.plot.width = 8)
UMAPPlot(object = scRNA_Mesen,label = T,group.by="All_second_cluster")

name = colnames(scRNA_Mesen)[ scRNA_Mesen$All_second_cluster  %in% c("FAP+ Fibroblast","CYSLTR2+ Fibroblast")]
one_cell_matrix = raw_matrix[,name]
one_cell_meta = scRNA_Mesen@meta.data[name,]

data2 = CreateSeuratObject(one_cell_matrix,min.cells = 0,min.features = 0,meta.data = one_cell_meta)
data2$percent.mt = PercentageFeatureSet(object = data2,pattern = "^MT-")
data2$percent.rp = PercentageFeatureSet(object = data2,pattern = "^RP[SL]")
data2$percent.hsp = PercentageFeatureSet(object = data2,pattern = "HSP")


data2 = NormalizeData(object = data2,normalization.method = "LogNormalize",scale.factor = 10000)
data2 = FindVariableFeatures(object = data2,selection.method = "vst",nfeatures = 2000)
data2 = ScaleData(object = data2,features = rownames(data2))
npcs = 20
data2 = RunPCA(object = data2,features = VariableFeatures(data2),npcs = npcs)

data2 = FindNeighbors(object = data2,reduction = "pca",dims = c(1:npcs),k.param = 10)
data2 = FindClusters(object = data2,resolution = 1)

data2 = RunUMAP(object = data2,reduction = "pca",dims = c(1:npcs))
data2$All_second_cluster = data2$All_first_cluster

data2 = FindNeighbors(object = data2,reduction = "pca",dims = c(1:npcs),k.param = 10)
data2 = FindClusters(object = data2,resolution = 2)

options(repr.plot.height = 8, repr.plot.width = 8)
UMAPPlot(object = data2,label = T)
UMAPPlot(object = data2,label = T,group.by="All_second_cluster")

options(repr.plot.height = 12, repr.plot.width = 12)
FeaturePlot(object = data2,features = c("THBS2","FAP","CD36","S100A4"))
FeaturePlot(object = data2,features = c("FBLN5","MFAP5","ADH1B","CYSLTR2"))
FeaturePlot(object = data2,features = c("CREB3L1","TBX2","MAF","FOXO1"))

UMAPPlot(object = data2,label = T)

re_name = colnames(data2)[ data2$seurat_clusters %in% c(7)]
scRNA_Mesen$All_second_cluster[re_name] = "MFAP5+ Fibroblast"
re_name = colnames(data2)[ data2$seurat_clusters %in% c(3,6,5,2)]
scRNA_Mesen$All_second_cluster[re_name] = "CYSLTR2+ Fibroblast"
re_name = colnames(data2)[ data2$seurat_clusters %in% c(1,4,0,8)]
scRNA_Mesen$All_second_cluster[re_name] = "FAP+ Fibroblast"

options(repr.plot.height = 8, repr.plot.width = 8)
UMAPPlot(object = scRNA_Mesen,label = T,group.by="All_second_cluster")

options(repr.plot.height = 5, repr.plot.width = 6)
n = length(unique(scRNA_Mesen$All_second_cluster))
col2 = colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(n)
names(col2) = unique(scRNA_Mesen$All_second_cluster)
UMAPPlot(object = scRNA_Mesen,label = F,group.by="All_second_cluster",label.box=FALSE,cols=col2)+    
    ggtitle("scRNA-seq")

clean_scRNA_Mesen = scRNA_Mesen


######## scATAC-seq ############

scATAC = readRDS("/home/ncpsb/new_mnt/ccRCC_data/Validation/scATAC_seq_ccRCC/All_scATAC.rds")

scATAC@assays$peaks@fragments[[1]]@path = "/home/ncpsb/new_mnt/ccRCC_data/Validation/scATAC_seq_ccRCC/data/scATAC_seq/fragments.tsv.gz"

macs2_peaks = readRDS("/home/ncpsb/new_mnt/ccRCC_data/Validation/scATAC_seq_ccRCC/MACS2/cellType.peak.rds")

mesenchymal_peaks = macs2_peaks[ macs2_peaks$peak_called_in %in% c("Fibroblast","Pericyte"),] 

mesenchymal_peaks <- keepStandardChromosomes(mesenchymal_peaks, pruning.mode = "coarse")
mesenchymal_peaks <- subsetByOverlaps(x = mesenchymal_peaks, ranges = blacklist_hg38_unified, invert = TRUE)


macs2_counts <- FeatureMatrix(
  fragments = Fragments(scATAC), # from cellranger fragment result
  features = mesenchymal_peaks,
  cells = colnames(scATAC)[ scATAC$All_first_cluster %in% c("Fibroblast","Pericyte") ]
)

annotations = readRDS("/home/ncpsb/new_mnt/ccRCC_data/Validation/scATAC_seq_ccRCC/data/scATAC_seq/annotations.rds")
chrom_assay <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(scATAC),
  annotation = annotations, 
)

meta = scATAC@meta.data[colnames(scATAC)[ scATAC$All_first_cluster %in% c("Fibroblast","Pericyte") ] ,]

scATAC_Mesen <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "MACS2",
  meta.data = meta
)

DefaultAssay(scATAC_Mesen) <- "MACS2"
gene.activities <- GeneActivity(scATAC_Mesen)

scATAC_Mesen[['Macs2ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
scATAC_Mesen <- NormalizeData(
  object = scATAC_Mesen,
  assay = 'Macs2ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(scATAC_Mesen$nCount_MACS2)
)

clean_scATAC_Mesen = scATAC_Mesen

scATAC_Mesen = clean_scATAC_Mesen
scATAC_Mesen = RunTFIDF(object = scATAC_Mesen,assay = "MACS2",method = 1)
scATAC_Mesen = FindTopFeatures(object = scATAC_Mesen,assay = "MACS2",min.cutoff = "q10")
scATAC_Mesen = RunSVD(object = scATAC_Mesen,assay = "MACS2",n = 20)
scATAC_Mesen <- RunUMAP(scATAC_Mesen, dims = 2:20, reduction = 'lsi')
scATAC_Mesen <- FindNeighbors(object = scATAC_Mesen, reduction = 'lsi', dims = 2:20,k.param = 10)
scATAC_Mesen <- FindClusters(object = scATAC_Mesen, verbose = FALSE, algorithm = 3,resolution = 1)

scATAC_Mesen <- FindNeighbors(object = scATAC_Mesen, reduction = 'lsi', dims = 2:20,k.param = 5)
scATAC_Mesen <- FindClusters(object = scATAC_Mesen, verbose = FALSE, algorithm = 3,resolution = 3)

options(repr.plot.height = 8, repr.plot.width = 8)
DimPlot(scATAC_Mesen,label = T)

options(repr.plot.height = 8, repr.plot.width = 8)
DimPlot(scATAC_Mesen, group.by = 'All_first_cluster', pt.size = 1, reduction = 'umap',label = T)
DimPlot(scATAC_Mesen, group.by = 'predict_cellType', pt.size = 1, reduction = 'umap',label = F)

options(repr.plot.height = 12, repr.plot.width = 12)
DefaultAssay(scATAC_Mesen) <- "Macs2ACTIVITY"
scATAC_Mesen <- ScaleData(scATAC_Mesen, features = rownames(scATAC_Mesen))

FeaturePlot(object = scATAC_Mesen,features = c("THBS2","FAP","CD36","S100A4"))
FeaturePlot(object = scATAC_Mesen,features = c("FBLN5","MFAP5","ADH1B","CYSLTR2"))
FeaturePlot(object = scATAC_Mesen,features = c("CREB3L1","TBX2","MAF","FOXO1"))

DefaultAssay(scRNA_Mesen) <- "RNA"
var.genes <- VariableFeatures(scRNA_Mesen)

transfer.anchors <- FindTransferAnchors(    reference = scRNA_Mesen, 
                                            query = scATAC_Mesen, 
                                            features = var.genes, #SCT model：scRNA.merge@assays$SCT@var.features
                                            reference.assay = "RNA",
                                            query.assay = "Macs2ACTIVITY", 
                                            reduction = "cca",
                                            npcs = 20, # Number of PCs to compute on reference
                                            dims = 1:20) # Which dimensions to use from the reduction to specify the neighbor search space

celltype.predictions <- TransferData(  anchorset = transfer.anchors, 
                                        refdata = scRNA_Mesen@meta.data[,"All_second_cluster"], 
                                        weight.reduction = scATAC_Mesen[["lsi"]], 
                                        dims = 2:20)  

scATAC_Mesen <- AddMetaData(scATAC_Mesen, metadata = celltype.predictions)

options(repr.plot.height = 12, repr.plot.width = 12)
n = length(unique(scATAC_Mesen$predict_cellType))
col2 = colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(n)
col2 = c(col2[seq(1,n,by = 2)],col2[seq(2,n,by = 2)])
DimPlot(scATAC_Mesen, group.by = "predicted.id", label = TRUE, repel = TRUE, reduction = "umap",cols = col2)

scATAC_Mesen$predicted.id2 = scATAC_Mesen$predicted.id
scATAC_Mesen$predicted.id2[ scATAC_Mesen$prediction.score.max < 0.5 ] = NA
DimPlot(scATAC_Mesen, group.by = "predicted.id2", label = TRUE, repel = TRUE, reduction = "umap",cols = col2)

for(i in unique(scATAC_Mesen$predicted.id)){
    scATAC_Mesen$tmp = NA
    scATAC_Mesen$tmp[ scATAC_Mesen$predicted.id == i ] = i
    print(DimPlot(scATAC_Mesen,group.by = "tmp",pt.size = 2))
}
scATAC_Mesen$tmp = NULL

for(i in unique(scATAC_Mesen$seurat_clusters)){
    scATAC_Mesen$tmp = NA
    scATAC_Mesen$tmp[ scATAC_Mesen$seurat_clusters == i ] = i
    print(DimPlot(scATAC_Mesen,group.by = "tmp",pt.size = 2))
}
scATAC_Mesen$tmp = NULL

scATAC_Mesen$All_second_cluster = scATAC_Mesen$All_first_cluster

cell_name = colnames(scATAC_Mesen)[ scATAC_Mesen$seurat_clusters %in% c(5)]
scATAC_Mesen$All_second_cluster[cell_name] = "MFAP5+ Fibroblast"

cell_name = colnames(scATAC_Mesen)[ scATAC_Mesen$seurat_clusters %in% c(1)]
scATAC_Mesen$All_second_cluster[cell_name] = "FAP+ Fibroblast"

cell_name = colnames(scATAC_Mesen)[ scATAC_Mesen$seurat_clusters %in% c(8,0,3,6,4,2,16,19,12,15,18,17)]
scATAC_Mesen$All_second_cluster[cell_name] = "CD36+ Pericyte"

cell_name = colnames(scATAC_Mesen)[ scATAC_Mesen$seurat_clusters %in% c(7,11,9,10,14,13)]
scATAC_Mesen$All_second_cluster[cell_name] = "S100A4+ Pericyte"

options(repr.plot.height = 12, repr.plot.width = 12)
n = length(unique(scATAC_Mesen$All_second_cluster))
col2 = colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(n)
col2 = c(col2[seq(1,n,by = 2)],col2[seq(2,n,by = 2)])
DimPlot(scATAC_Mesen, group.by = "All_second_cluster", label = TRUE, repel = TRUE, reduction = "umap",cols = col2)

DefaultAssay(scATAC_Mesen) <- "MACS2"

if( !dir.exists("./MACS2") ){
    dir.create("./MACS2")

}
peaks <- CallPeaks(
                          object = scATAC_Mesen,
                          group.by = "All_second_cluster",
                          macs2.path = "/home/ncpsb/anaconda3/envs/scATAC/bin/macs2",
                          outdir = "./MACS2"
                        )

saveRDS(peaks, "./MACS2/Mesenchymal_cellType.peak.rds") 


peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)


Mesen_macs2_counts <- FeatureMatrix(
  fragments = Fragments(scATAC_Mesen), # from cellranger fragment result
  features = peaks,
  cells = colnames(scATAC_Mesen)
)

annotations = readRDS("/home/ncpsb/new_mnt/ccRCC_data/Validation/scATAC_seq_ccRCC/data/scATAC_seq/annotations.rds")

chrom_assay <- CreateChromatinAssay(
  counts = Mesen_macs2_counts,
  fragments = Fragments(scATAC_Mesen),
  annotation = annotations, 
)

meta = scATAC_Mesen@meta.data

scATAC_Mesen_2nd <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "MACS2",
  meta.data = meta
)

DefaultAssay(scATAC_Mesen_2nd) <- "MACS2"
gene.activities <- GeneActivity(scATAC_Mesen_2nd)

scATAC_Mesen_2nd[['Macs2ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
scATAC_Mesen_2nd <- NormalizeData(
  object = scATAC_Mesen_2nd,
  assay = 'Macs2ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(scATAC_Mesen_2nd$nCount_MACS2)
)

clean_scATAC_Mesen_2nd = scATAC_Mesen_2nd

scATAC_Mesen_2nd = clean_scATAC_Mesen_2nd
scATAC_Mesen_2nd = RunTFIDF(object = scATAC_Mesen_2nd,assay = "MACS2",method = 1)
scATAC_Mesen_2nd = FindTopFeatures(object = scATAC_Mesen_2nd,assay = "MACS2",min.cutoff = "q10")
scATAC_Mesen_2nd = RunSVD(object = scATAC_Mesen_2nd,assay = "MACS2",n = 20)
scATAC_Mesen_2nd <- RunUMAP(scATAC_Mesen_2nd, dims = 2:20, reduction = 'lsi')
scATAC_Mesen_2nd <- FindNeighbors(object = scATAC_Mesen_2nd, reduction = 'lsi', dims = 2:20,k.param = 10)
scATAC_Mesen_2nd <- FindClusters(object = scATAC_Mesen_2nd, verbose = FALSE, algorithm = 3,resolution = 1)

options(repr.plot.height = 8, repr.plot.width = 8)
DimPlot(scATAC_Mesen_2nd, group.by = 'All_second_cluster', pt.size = 1, reduction = 'umap',label = T)
DimPlot(scATAC_Mesen_2nd, group.by = 'predicted.id', pt.size = 1, reduction = 'umap',label = F)

scATAC_Mesen_2nd <- FindNeighbors(object = scATAC_Mesen_2nd, reduction = 'lsi', dims = 2:20,k.param = 10)
scATAC_Mesen_2nd <- FindClusters(object = scATAC_Mesen_2nd, verbose = FALSE, algorithm = 3,resolution = 2)

options(repr.plot.height = 8, repr.plot.width = 8)
DimPlot(scATAC_Mesen_2nd,label = T)

options(repr.plot.height = 12, repr.plot.width = 12)

DefaultAssay(scATAC_Mesen_2nd) <- "Macs2ACTIVITY"
scATAC_Mesen_2nd <- ScaleData(scATAC_Mesen_2nd, features = rownames(scATAC_Mesen_2nd))

FeaturePlot(object = scATAC_Mesen_2nd,features = c("THBS2","FAP","CD36","S100A4"),max.cutoff = 0.1,pt.size = 2)
FeaturePlot(object = scATAC_Mesen_2nd,features = c("FBLN5","MFAP5","ADH1B","CYSLTR2"))
FeaturePlot(object = scATAC_Mesen_2nd,features = c("CREB3L1","TBX2","MAF","FOXO1"),max.cutoff = 0.15)

cell_name = colnames(scATAC_Mesen_2nd)[ scATAC_Mesen_2nd$seurat_clusters %in% c(10)]
scATAC_Mesen_2nd$All_second_cluster[cell_name] = "MFAP5+ Fibroblast"

cell_name = colnames(scATAC_Mesen_2nd)[ scATAC_Mesen_2nd$seurat_clusters %in% c(5,3)]
scATAC_Mesen_2nd$All_second_cluster[cell_name] = "FAP+ Fibroblast"

cell_name = colnames(scATAC_Mesen_2nd)[ scATAC_Mesen_2nd$seurat_clusters %in% c(0,1,4,7,9)]
scATAC_Mesen_2nd$All_second_cluster[cell_name] = "CD36+ Pericyte"

cell_name = colnames(scATAC_Mesen_2nd)[ scATAC_Mesen_2nd$seurat_clusters %in% c(6,2,8)]
scATAC_Mesen_2nd$All_second_cluster[cell_name] = "S100A4+ Pericyte"

cell_name = colnames(scATAC_Mesen_2nd)[ scATAC_Mesen_2nd$seurat_clusters %in% c(10)]
clean_scATAC_Mesen_2nd$All_second_cluster[cell_name] = "MFAP5+ Fibroblast"

cell_name = colnames(scATAC_Mesen_2nd)[ scATAC_Mesen_2nd$seurat_clusters %in% c(5,3)]
clean_scATAC_Mesen_2nd$All_second_cluster[cell_name] = "FAP+ Fibroblast"

cell_name = colnames(scATAC_Mesen_2nd)[ scATAC_Mesen_2nd$seurat_clusters %in% c(0,1,4,7,9)]
clean_scATAC_Mesen_2nd$All_second_cluster[cell_name] = "CD36+ Pericyte"

cell_name = colnames(scATAC_Mesen_2nd)[ scATAC_Mesen_2nd$seurat_clusters %in% c(6,2,8)]
clean_scATAC_Mesen_2nd$All_second_cluster[cell_name] = "S100A4+ Pericyte"

options(repr.plot.height = 12, repr.plot.width = 12)
n = length(unique(scATAC_Mesen_2nd$All_second_cluster))
col2 = colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(n)
col2 = c(col2[seq(1,n,by = 2)],col2[seq(2,n,by = 2)])
DimPlot(scATAC_Mesen_2nd, group.by = "All_second_cluster", label = TRUE, repel = TRUE, reduction = "umap",cols = col2,pt.size = 2)

options(repr.plot.height = 5, repr.plot.width = 6)
# n = length(unique(scATAC_Mesen_2nd$All_second_cluster))
# col2 = colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(n)
# col2 = c(col2[seq(1,n,by = 2)],col2[seq(2,n,by = 2)])
DimPlot(scATAC_Mesen_2nd, group.by = "All_second_cluster",label.box = F, label = F, repel = FALSE, reduction = "umap",cols = col2,pt.size = 2)+
#     NoLegend()+
    ggtitle("scATAC-seq")

saveRDS(list("clean"=clean_scATAC_Mesen_2nd,"final"=scATAC_Mesen_2nd),"./data/scATAC_seq_Mesenchymal_list.rds")

############### Trajectory inference based on scATAC-seq data ###########################

peaks = readRDS("./MACS2/Mesenchymal_cellType.peak.rds") 
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

feature = paste(peaks@seqnames[!grepl("S100A4\\+_Pericyte|MFAP5\\+_Fibroblast",peaks$peak_called_in),],
                peaks@ranges[!grepl("S100A4\\+_Pericyte|MFAP5\\+_Fibroblast",peaks$peak_called_in),],sep = "-")


clean_scATAC_Mesen_2nd = readRDS("./data/scATAC_seq_Mesenchymal_list.rds")$clean

Mesen_pseudo_macs2_counts <- FeatureMatrix(
  fragments = Fragments(clean_scATAC_Mesen_2nd), # from cellranger fragment result
  features = peaks,
  cells = colnames(clean_scATAC_Mesen_2nd)[ clean_scATAC_Mesen_2nd$All_second_cluster %in% c("FAP+ Fibroblast","CD36+ Pericyte")]
)

annotations = readRDS("/home/ncpsb/new_mnt/ccRCC_data/Validation/scATAC_seq_ccRCC/data/scATAC_seq/annotations.rds")

chrom_assay <- CreateChromatinAssay(
  counts = Mesen_pseudo_macs2_counts,
  fragments = Fragments(clean_scATAC_Mesen_2nd),
  annotation = annotations, 
)

meta = clean_scATAC_Mesen_2nd@meta.data[ clean_scATAC_Mesen_2nd$All_second_cluster %in% c("FAP+ Fibroblast","CD36+ Pericyte"),]

scATAC_Mesen_3rd <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "MACS2",
  meta.data = meta
)

DefaultAssay(scATAC_Mesen_3rd) <- "MACS2"
gene.activities <- GeneActivity(scATAC_Mesen_3rd)

scATAC_Mesen_3rd[['Macs2ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
scATAC_Mesen_3rd <- NormalizeData(
  object = scATAC_Mesen_3rd,
  assay = 'Macs2ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(scATAC_Mesen_3rd$nCount_MACS2)
)

clean_scATAC_Mesen_3rd = scATAC_Mesen_3rd

scATAC_Mesen_3rd = clean_scATAC_Mesen_3rd
scATAC_Mesen_3rd = RunTFIDF(object = scATAC_Mesen_3rd,assay = "MACS2",method = 1)
scATAC_Mesen_3rd = RunSVD(object = scATAC_Mesen_3rd,assay = "MACS2",n = 10,features = feature)
scATAC_Mesen_3rd <- RunUMAP(scATAC_Mesen_3rd, dims = 2:10, reduction = 'lsi')
scATAC_Mesen_3rd <- FindNeighbors(object = scATAC_Mesen_3rd, reduction = 'lsi', dims = 2:10,k.param = 10)
scATAC_Mesen_3rd <- FindClusters(object = scATAC_Mesen_3rd, verbose = FALSE, algorithm = 3,resolution = 1)

options(repr.plot.height = 8, repr.plot.width = 8)
DimPlot(scATAC_Mesen_3rd, group.by = 'All_second_cluster', pt.size = 1, reduction = 'umap',label = T)
DimPlot(scATAC_Mesen_3rd, group.by = 'predicted.id', pt.size = 1, reduction = 'umap',label = F)

saveRDS(list("clean"=clean_scATAC_Mesen_3rd,"final"=scATAC_Mesen_3rd),"./data/scATAC_seq_Mesenchymal_pseudo_list.rds")

tmp_data = readRDS("./data/scATAC_seq_Mesenchymal_pseudo_list.rds")
clean_scATAC_Mesen_3rd = tmp_data$clean
scATAC_Mesen_3rd = tmp_data$final
remove(tmp_data)

 CoveragePlot(
      object = scATAC_Mesen_3rd,
      group.by = "All_second_cluster",
      region = "CREB3L1",
      ranges.title = "MACS2",
      links = F,
      peaks = T,
      extend.upstream = 1000,
      extend.downstream = 1000,
 )

ct = as.ExpressionSet(as.data.frame(scATAC_Mesen_3rd@reductions$lsi@cell.embeddings)[,2:10])
ct$celltype <- as.character(scATAC_Mesen_3rd@meta.data[,c("All_second_cluster")])

dm <- DiffusionMap(ct,sigma = "local",n_pcs = NA)

plot(dm,c(1,2), pch = 20,col_by = "celltype")

scATAC_Mesen_3rd[["DC"]] = CreateDimReducObject(embeddings = dm@eigenvectors[,1:2],key = "DC_",assay = "MACS2")
st = as.SingleCellExperiment(scATAC_Mesen_3rd,assay = "MACS2")

ss = slingshot(data = st,clusterLabels = "All_second_cluster",reducedDim = "DC",start.clus = c("CD36+ Pericyte"))

n = length(unique(scATAC_Mesen_3rd$All_second_cluster))
col3 = colorRampPalette(brewer.pal(9,"Set1"))(n)
col3 = c(col3[seq(1,n,by = 2)],col3[seq(2,n,by = 2)])
names(col3) = unique(scATAC_Mesen_3rd$All_second_cluster)

plot(reducedDims(ss)$DC,col = col3[ ss$All_second_cluster ], pch=16, asp = 1)
lines(SlingshotDataSet(ss), lwd=2, col='black')

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(ss$slingPseudotime_1, breaks=100)]

plot(reducedDims(ss)$DC, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(ss), lwd=2, col='black')

DC_Matrix = as.data.frame(cbind(as.data.frame(reducedDims(ss)$DC),pseudotime = ss$slingPseudotime_1,Cell_type=as.character(ss$All_second_cluster)))

library(viridis)

options(repr.plot.height = 8, repr.plot.width = 8)
ggplot(DC_Matrix)+geom_point(aes(x=DC_1,y=DC_2,color=Cell_type))+theme_bw()+
scale_color_manual(values = c("red","blue"))+
theme(panel.grid.major = element_line(linewidth=0.5,linetype = "dashed",colour = "grey60"),
      legend.margin = margin(16, 16, 16, 16),
      legend.title = element_text(size=16,face="bold"),
      legend.text = element_text(size=14,face="bold"),
      legend.position = c(0.5,0.8),
      legend.background = element_rect(fill = "lightgreen",linewidth = 0.5,color="black"),
      axis.title = element_text(size=18,face = "bold"),
      axis.text = element_text(size=14),
      panel.border = element_rect(linewidth = 1.5),
     panel.grid.minor = element_line(linewidth=0.5,linetype = "dashed",colour = "grey80"),
     
     
     )

options(repr.plot.height = 8, repr.plot.width = 8)
ggplot(DC_Matrix)+geom_point(aes(x=DC_1,y=DC_2,color=pseudotime))+theme_bw()+
# scale_color_manual(values = col)+
scale_color_viridis()+
theme(panel.grid.major = element_line(linewidth=0.5,linetype = "dashed",colour = "grey60"),
      legend.margin = margin(16, 16, 16, 16),
      legend.title = element_text(size=16,face="bold"),
      legend.text = element_text(size=14,face="bold"),
      legend.position = c(0.5,0.75),
      legend.background = element_rect(fill = "lightgreen",linewidth = 0.5,color="black"),
      axis.title = element_text(size=18,face = "bold"),
      axis.text = element_text(size=14),
      panel.border = element_rect(linewidth = 1.5),
     panel.grid.minor = element_line(linewidth=0.5,linetype = "dashed",colour = "grey80"),
     
     
     )

plot_genebypseudotime = function(Pseudotime,gene_expression,gene_name,Clusters,legend.position=c(0.8,0.8),col2){
    pseu_matrix = cbind(Pseudotime,gene_expression)
    pseu_matrix = as.data.frame(pseu_matrix)
    colnames(pseu_matrix) = c("Pseudotime","Gene_expression")
    pseu_matrix$Clusters = Clusters
    pseu_matrix = pseu_matrix[ order(pseu_matrix$Pseudotime),]
    
    ggplot(pseu_matrix,aes(x = Pseudotime,y=Gene_expression))+
    scale_color_manual(values = col2)+
    geom_point(aes(color=Clusters))+
    geom_smooth(method = "gam")+
    ylab(gene_name)+
    theme_bw()+
    theme(panel.grid.major = element_line(linewidth=0.5,linetype = "dashed",colour = "grey60"),
      legend.margin = margin(16, 16, 16, 16),
      legend.title = element_text(size=16,face="bold"),
      legend.text = element_text(size=14,face="bold"),
      legend.position = legend.position,
      legend.background = element_rect(fill = "lightgreen",linewidth = 0.5,color="black"),
      axis.title = element_text(size=18,face = "bold"),
      axis.text = element_text(size=14),
      panel.border = element_rect(linewidth = 1.5),
     panel.grid.minor = element_line(linewidth=0.5,linetype = "dashed",colour = "grey80"),
     
     
     )
}


options(repr.plot.height = 8, repr.plot.width = 8)
plot_genebypseudotime(Pseudotime = ss$slingPseudotime_1[!is.na(ss$slingPseudotime_1)],legend.position = c(0.2,0.8),col2=c("red","blue"),
                      gene_expression = as.numeric(scATAC_Mesen_3rd@assays$Macs2ACTIVITY@data["MCAM",])[!is.na(ss$slingPseudotime_1)],
                      gene_name = "MCAM(+)",
                      Clusters = ss$All_second_cluster[!is.na(ss$slingPseudotime_1)]
                     )
plot_genebypseudotime(Pseudotime = ss$slingPseudotime_1[!is.na(ss$slingPseudotime_1)],legend.position = c(0.2,0.8),col2=c("red","blue"),
                      gene_expression = as.numeric(scATAC_Mesen_3rd@assays$Macs2ACTIVITY@data["DCN",])[!is.na(ss$slingPseudotime_1)],
                      gene_name = "DCN(+)",
                      Clusters = ss$All_second_cluster[!is.na(ss$slingPseudotime_1)]
                     )


options(repr.plot.height = 10, repr.plot.width = 10)
tmp_ss = ss
tmp_st = st[,colnames(tmp_ss)]

plot(reducedDims(tmp_ss)$DC, col =col3[ tmp_ss$All_second_cluster], pch=16, asp = 1)
lines(SlingshotDataSet(tmp_ss), lwd=3, type = 'lineages', col = 'black')

ss_ori_tbl <- tibble(cell = colnames(tmp_st), pseudotime = rescale(colData(tmp_ss)$slingPseudotime_1))
tmp_name = colnames(tmp_st)
set.seed( 1234, kind = "L'Ecuyer-CMRG" );
options(mc.cores =12)
n = 100
ss_index <- mclapply(seq_len(n), function(x) {
  sample(x = c(1:dim(tmp_st)[2]), size = 0.8*dim(tmp_st)[2], replace = FALSE)
})

clean_Fibro <- clean_scATAC_Mesen_3rd
DefaultAssay(clean_Fibro) <- "MACS2"

ss_sub_tbl <- mclapply(ss_index, function(x) {
    clean_Fibro2 <- clean_Fibro[,tmp_name[x]]
    clean_Fibro2 <- RunTFIDF(object = clean_Fibro2,assay = "MACS2",method = 1)
    clean_Fibro2 <- FindTopFeatures(object = clean_Fibro2,assay = "MACS2",min.cutoff = "q10")
    clean_Fibro2 <- RunSVD(object = clean_Fibro2,assay = "MACS2",n = 10)

    ct2 = as.ExpressionSet(as.data.frame(clean_Fibro2@reductions$lsi@cell.embeddings)[,2:10])
    ct2$celltype <- as.character(clean_Fibro2@meta.data[,c("All_second_cluster")])

    dm2 <- DiffusionMap(ct2)

    clean_Fibro2[["DC"]] = CreateDimReducObject(embeddings = dm2@eigenvectors[,1:2],key = "DC_",assay = "MACS2")

    sce = as.SingleCellExperiment(clean_Fibro2,assay = "MACS2")
    fit = slingshot(data = sce,clusterLabels = "All_second_cluster",reducedDim = "DC",start.clus = c("CD36+ Pericyte"))

    if( sum(grepl("slingPseudotime",colnames(sce@colData))) > 1){
            pseu_M = cbind(fit$slingPseudotime_1,fit$slingPseudotime_2,fit$slingPseudotime_3,fit$slingPseudotime_4,fit$slingPseudotime_5,fit$slingPseudotime_6,fit$slingPseudotime_7)
            fit$all_pseudotime = rowMeans(pseu_M,na.rm = T)
            tbl <- tibble(cell = colnames(sce), pseudotime = rescale(colData(fit)$all_pseudotime))
    }else{
         tbl <- tibble(cell = colnames(sce), pseudotime = rescale(colData(fit)$slingPseudotime_1))
    }


      ## Make sure the direction of pseudotime is the same as the original pseudotime
      merge.tbl <- left_join(tbl, ss_ori_tbl, by = "cell")

      if(cor(merge.tbl$pseudotime.x, merge.tbl$pseudotime.y) < 0) {
        tbl <- dplyr::mutate(tbl, pseudotime = 1-pseudotime)
      }
      tbl
})

PseudotimeDE::plotUncertainty(ss_ori_tbl, ss_sub_tbl)

TF = read.table("/home/ncpsb/new_mnt/Other/Reference_data/SCENIC/allTFs_hg38.txt")$V1

system.time(res <- PseudotimeDE::runPseudotimeDE(gene.vec = rownames(scATAC_Mesen_3rd@assays$Macs2ACTIVITY@data),
                                                 assay.use = "logcounts",
                                                 ori.tbl = ss_ori_tbl,
                                                 sub.tbl = ss_sub_tbl, ## To save time, use 100 subsamples
                                                 mat = scATAC_Mesen_3rd@assays$Macs2ACTIVITY@data, ## You can also use a matrix or SeuratObj as the input
                                                 model = "gaussian",
                                                 mc.cores = 16))

saveRDS(list(ss,st,ss_ori_tbl,ss_sub_tbl,res),"./data/scATAC_seq_Mesenchymal_pseudotime_list.rds")

######### chromVAR ###############

data("human_pwms_v2") #filtered collection of human motifs from cisBP database


clean_Fibro <- clean_scATAC_Mesen_3rd
DefaultAssay(clean_Fibro) <- "MACS2"

clean_Fibro <- AddMotifs(object = clean_Fibro,genome = BSgenome.Hsapiens.UCSC.hg38,pfm = human_pwms_v2)

clean_Fibro <- RunChromVAR(object = clean_Fibro,genome = BSgenome.Hsapiens.UCSC.hg38)

DefaultAssay(clean_Fibro) <- "chromvar"

data("human_pwms_v2")
TF_name = vector()
for(i in human_pwms_v2@listData){
    TF_name = c(TF_name,i@name)
}

system.time(res <- PseudotimeDE::runPseudotimeDE(gene.vec = rownames(clean_Fibro),
                                                 assay.use = "logcounts",
                                                 ori.tbl = ss_ori_tbl,
                                                 sub.tbl = ss_sub_tbl, ## To save time, use 100 subsamples
                                                 mat = clean_Fibro@assays$chromvar@data, ## You can also use a matrix or SeuratObj as the input
                                                 model = "gaussian",
                                                 mc.cores = 12))

res$gene = TF_name

saveRDS(list(clean_Fibro,ss,res),"./data/scATAC_seq_Mesenchymal_chromvar_list.rds")

######## optical match ##############
scATAC = readRDS("./data/scATAC_seq_Mesenchymal_list.rds")$final
scRNA <- readRDS("/home/ncpsb/new_mnt/ccRCC_data/Validation/scATAC_seq_ccRCC/All_scRNA.rds")

ene.use = VariableFeatures(scRNA)

transfer.anchors <- FindTransferAnchors(reference = scRNA,
                                        query = scATAC,
                                        features = gene.use,
                                        reference.assay = "RNA",
                                        query.assay = "Macs2ACTIVITY",
                                        reduction = "cca")

refdata <- GetAssayData(scRNA, assay = "RNA", slot = "data")[gene.use, ]

imputation <- TransferData(anchorset = transfer.anchors, 
                           refdata = refdata,  
                           weight.reduction = scATAC[["lsi"]], # 因为只有一个患者，没有harmony
                           dims = 1:20)
scATAC[["RNA"]] <- imputation

DefaultAssay(scATAC) <- "RNA"
scRNA$tech <- "RNA"
scATAC$tech <- "ATAC"
coembed <- merge(x = scATAC, y = scRNA)

coembed <- coembed %>%
    ScaleData(features = gene.use, do.scale = FALSE) %>%
    RunPCA(features = gene.use) %>%
    RunHarmony(group.by.vars = "tech",reduction = "pca",dims.use = 1:30,max.iter.harmony = 2)%>%
    RunUMAP(dims = 1:30, verbose = FALSE,reduction = "harmony")

DimPlot(coembed,group.by = "tech")

DimPlot(coembed,group.by = "All_second_cluster")

options(repr.plot.height = 5, repr.plot.width = 10)
DimPlot(coembed, reduction = "umap", group.by = "All_second_cluster", split.by = "tech")

obj.atac <- subset(coembed, tech == "ATAC")
obj.rna <- subset(coembed, tech == "RNA")
obj.atac
obj.rna

cca_umap_df <- as.data.frame(coembed@reductions$umap@cell.embeddings)
colnames(cca_umap_df) <- c("UMAP1", "UMAP2")
head(cca_umap_df)

library("optmatch")

options(repr.plot.height = 6, repr.plot.width = 6)

df_cell_pairing <- cell_pairing(ATACpcs = obj.atac@reductions$umap@cell.embeddings,
                                RNApcs = obj.rna@reductions$umap@cell.embeddings,
                                max_multimatch=1,
                                cca_umap_df = cca_umap_df,
                                nCores = 24
                               )
     

rownames(df_cell_pairing) = df_cell_pairing$ATAC

df_cell_pairing

sel_cells <- c(df_cell_pairing$ATAC,df_cell_pairing$RNA)
coembed.sub <- coembed[, sel_cells]

options(repr.plot.height = 5, repr.plot.width = 10)
DimPlot(coembed.sub, reduction = "umap", group.by = "All_second_cluster", split.by = "tech")

coembed$pairs = NA
coembed$pairs[ "AAACGAATCCGGCTGA-1" ] = "Pair_1"
coembed$pairs[ "T1_GTGCTGGCATCCGTGG-1" ] = "Pair_2"
# coembed$pairs[ "T1_AGGCATTAGCGACTAG-1" ] = "Pair_3"
# coembed$pairs[ "T1_ATATCCTCAATGTCTG-1" ] = "Pair_4"
# coembed$pairs[ "T1_TCGTCCAAGACGACGT-1" ] = "Pair_5"
# coembed$pairs[ "T1_CCTTCAGCAAATTGCC-1" ] = "Pair_6"
DimPlot(coembed,group.by = "pairs",pt.size = 0.1)

#############  correlation among four TF-associated metrics ###############
library(slingshot)
library(SummarizedExperiment)
library(destiny)
library(Biobase)

Regulons1 = read.csv("SCENIC/OUTPUT/RESULT/Regulons_scRNA_1.csv",row.names = 1)
Regulons1 = as.data.frame(t(Regulons1))

rownames(Regulons1) = sapply(rownames(Regulons1),function(x){gsub("\\.\\.\\.","",x)})

Regulons2 = read.csv("SCENIC/OUTPUT/RESULT/Regulons_scRNA_2.csv",row.names = 1)
Regulons2 = as.data.frame(t(Regulons2))

rownames(Regulons2) = sapply(rownames(Regulons2),function(x){gsub("\\.\\.\\.","",x)})

intersect(colnames(Regulons1)[Regulons1["CREB3L1",,drop=T] != 0],colnames(Regulons2)[Regulons2["CREB3L1",,drop=T] != 0])

intersect(colnames(Regulons1)[Regulons1["MEF2C",,drop=T] != 0],colnames(Regulons2)[Regulons2["MEF2C",,drop=T] != 0])

RegulonsAUC = read.csv("SCENIC/OUTPUT/RESULT/RegulonsAUC_scRNA_1.csv",row.names = 1)
RegulonsAUC = as.data.frame(t(RegulonsAUC))

rownames(RegulonsAUC) = sapply(rownames(RegulonsAUC),function(x){gsub("\\.\\.\\.","",x)})

tmp_data = readRDS("./data/scATAC_seq_Mesenchymal_chromvar_list.rds")
clean_Fibro = tmp_data[[1]]
ss = tmp_data[[2]]
# res = tmp_data[[3]]
remove(tmp_data)

RNA_data = as.matrix(scRNA@assays$RNA@data[,df_cell_pairing[colnames(ss),"RNA"]])
Regulons_data = as.matrix(RegulonsAUC)[,df_cell_pairing[colnames(ss),"RNA"]]
ATAC_data = as.matrix(scATAC@assays$Macs2ACTIVITY@data)[,colnames(ss)]


RNA_data = RNA_data[rowSums(RNA_data > 0) >= length(colnames(ss))*0.05,]
Regulons_data = Regulons_data[rowSums(Regulons_data > 0) >= length(colnames(ss))*0.05,]
ATAC_data = ATAC_data[rowSums(ATAC_data > 0) >= length(colnames(ss))*0.05,]

share_genes = intersect(rownames(RNA_data),rownames(Regulons_data))
share_genes = intersect(share_genes,rownames(ATAC_data))

RNA_data = RNA_data[share_genes,]
Regulons_data = Regulons_data[share_genes,]
ATAC_data = ATAC_data[share_genes,]


result = vector()
for(i in share_genes){
    
    result = rbind(result,c(cor(RNA_data[i,],Regulons_data[i,],method = "spearman"),
                            cor(RNA_data[i,],ATAC_data[i,],method = "spearman"),
                            cor(Regulons_data[i,],ATAC_data[i,],method = "spearman")
                       
                       
                            )
                  )
}

rownames(result) = share_genes
colnames(result) = c("RNA_Regulons","RNA_ATAC","Regulons_ATAC")

result[ order(rowSums(result),decreasing = T),]

result = as.data.frame(result)

result$gene_name = rownames(result)
result$gene_name[ rowSums(result < 0.1) >= 1 ] = NA

library(plot3D)

options(repr.plot.height = 7, repr.plot.width = 7)
scatter3D(result$RNA_Regulons,result$RNA_ATAC,result$Regulons_ATAC,ticktype = "detailed",
          xlab = "Cor: RNA and Regulon",
          ylab = "Cor: RNA and ATAC",
          zlab = "Cor: Regulon and ATAC",
          pch = 19, cex = 1,bty = "u",phi = 10,col="purple",theta = 60,col.panel ="white", expand =1, col.grid = "grey80",lwd.grid=3)
text3D(result$RNA_Regulons,result$RNA_ATAC,result$Regulons_ATAC,labels = result$gene_name,pch = 19, cex = 1,add = TRUE,bty = "b2",phi = 2)

ThreeD1 = result

DC_Matrix = as.data.frame(cbind(as.data.frame(reducedDims(ss)$DC),pseudotime = ss$slingPseudotime_1,Cell_type=as.character(ss$All_second_cluster)))

DC_Matrix = DC_Matrix[ order(DC_Matrix$pseudotime),]
# clean_Fibro2 = clean_Fibro[,rownames(DC_Matrix)]
DC_Matrix$Cell_type2 = scRNA$All_second_cluster[df_cell_pairing[rownames(DC_Matrix),"RNA"]]

gene_name = "CREB3L1"
DC_Matrix[["Regulon"]] = as.numeric(RegulonsAUC[gene_name,df_cell_pairing[rownames(DC_Matrix),"RNA"]])[!is.na(ss$slingPseudotime_1)]
DC_Matrix[["RNA"]] = scRNA@assays$RNA@data[gene_name,df_cell_pairing[rownames(DC_Matrix),"RNA"]]
DC_Matrix[["ATAC"]] = scATAC@assays$Macs2ACTIVITY@data[gene_name,rownames(DC_Matrix)]


options(repr.plot.height = 8, repr.plot.width = 8)
ggplot(data = DC_Matrix,aes(x=pseudotime,y=Regulon))+geom_point(color="red")+geom_smooth(method = "gam")
ggplot(data = DC_Matrix,aes(x=pseudotime,y=RNA))+geom_point(color="blue")+geom_smooth(method = "gam")
ggplot(data = DC_Matrix,aes(x=pseudotime,y=ATAC))+geom_point(color="green")+geom_smooth(method = "gam")

options(repr.plot.height = 2.1, repr.plot.width = 10)
p1 = ggplot(data = DC_Matrix)+
    geom_smooth(aes(x=pseudotime,y=RNA),method = "gam",color="red")+theme_bw()+
    ylab("RNA\nExpression")+
    theme(panel.grid.major = element_line(linewidth=0.5,linetype = "dashed",colour = "grey60"),
          axis.title = element_text(size=18,face = "bold"),
          axis.text = element_text(size=14),
          panel.border = element_rect(linewidth = 1.5),
         panel.grid.minor = element_line(linewidth=0.5,linetype = "dashed",colour = "grey80"),
         )
p2 = ggplot(data = DC_Matrix)+
    geom_smooth(aes(x=pseudotime,y=Regulon),method = "gam",color="blue")+theme_bw()+
    ylab("Regulon\nActivity")+
    theme(panel.grid.major = element_line(linewidth=0.5,linetype = "dashed",colour = "grey60"),
          axis.title = element_text(size=18,face = "bold"),
          axis.text = element_text(size=14),
          panel.border = element_rect(linewidth = 1.5),
         panel.grid.minor = element_line(linewidth=0.5,linetype = "dashed",colour = "grey80"),
         )
p3 = ggplot(data = DC_Matrix)+
    geom_smooth(aes(x=pseudotime,y=ATAC),method = "gam",color="green")+theme_bw()+
    ylab("Chromatin\nAccessibility")+
    theme(panel.grid.major = element_line(linewidth=0.5,linetype = "dashed",colour = "grey60"),
          axis.title = element_text(size=18,face = "bold"),
          axis.text = element_text(size=14),
          panel.border = element_rect(linewidth = 1.5),
         panel.grid.minor = element_line(linewidth=0.5,linetype = "dashed",colour = "grey80"),
         )
p1 + p2 + p3


library(slingshot)
library(SummarizedExperiment)
library(destiny)
library(Biobase)

RegulonsAUC = read.csv("SCENIC/OUTPUT/RESULT/RegulonsAUC_scRNA_1.csv",row.names = 1)
RegulonsAUC = as.data.frame(t(RegulonsAUC))

rownames(RegulonsAUC) = sapply(rownames(RegulonsAUC),function(x){gsub("\\.\\.\\.","",x)})

tmp_data = readRDS("./data/scATAC_seq_Mesenchymal_chromvar_list.rds")
clean_Fibro = tmp_data[[1]]
ss = tmp_data[[2]]
res = tmp_data[[3]]
remove(tmp_data)

RNA_data = as.matrix(scRNA@assays$RNA@data[,df_cell_pairing[colnames(ss),"RNA"]])
Regulons_data = as.matrix(RegulonsAUC)[,df_cell_pairing[colnames(ss),"RNA"]]
Chrom_data = as.matrix(clean_Fibro@assays$chromvar@data)[,colnames(ss)]

rownames(Chrom_data) = res$gene


RNA_data = RNA_data[rowSums(RNA_data > 0) >= length(colnames(ss))*0.05,]
Regulons_data = Regulons_data[rowSums(Regulons_data > 0) >= length(colnames(ss))*0.05,]
Chrom_data = Chrom_data[rowSums(ATAC_data > 0) >= length(colnames(ss))*0.05,]

share_genes = intersect(rownames(RNA_data),rownames(Regulons_data))
share_genes = intersect(share_genes,rownames(Chrom_data))

RNA_data = RNA_data[share_genes,]
Regulons_data = Regulons_data[share_genes,]
Chrom_data = Chrom_data[share_genes,]

result = vector()
for(i in share_genes){
    
    result = rbind(result,c(cor(RNA_data[i,],Regulons_data[i,],method = "spearman"),
                            cor(RNA_data[i,],Chrom_data[i,],method = "spearman"),
                            cor(Regulons_data[i,],Chrom_data[i,],method = "spearman")
                       
                       
                            )
                  )
}

rownames(result) = share_genes
colnames(result) = c("RNA_Regulons","RNA_Chrom","Regulons_Chrom")

result[ order(rowSums(result),decreasing = T),]

result = as.data.frame(result)

result$gene_name = rownames(result)
result$gene_name[ rowSums(result < 0.1) >= 1 ] = NA

library(plot3D)

options(repr.plot.height = 7, repr.plot.width = 7)
scatter3D(result$RNA_Regulons,result$RNA_Chrom,result$Regulons_Chrom,ticktype = "detailed",
          xlab = "Cor: RNA and Regulon",
          ylab = "Cor: RNA and ChromVAR",
          zlab = "Cor: Regulon and ChromVAR",
          pch = 19, cex = 1,bty = "u",phi = 10,col="purple",theta = 60,col.panel ="white", expand =1, col.grid = "grey80",lwd.grid=3)
text3D(result$RNA_Regulons,result$RNA_Chrom,result$Regulons_Chrom,labels = result$gene_name,pch = 19, cex = 1,add = TRUE,bty = "b2",phi = 2)

ThreeD2 = result

saveRDS(list(f=ThreeD1,s=ThreeD2),"./data/ThreeD_1.rds")

gene_name = "MEF2C"
DC_Matrix[["Regulon"]] = as.numeric(RegulonsAUC[gene_name,df_cell_pairing[rownames(DC_Matrix),"RNA"]])[!is.na(ss$slingPseudotime_1)]
DC_Matrix[["RNA"]] = scRNA@assays$RNA@data[gene_name,df_cell_pairing[rownames(DC_Matrix),"RNA"]]
DC_Matrix[["Chrom"]] = Chrom_matrix[gene_name,rownames(DC_Matrix)]

options(repr.plot.height = 2.1, repr.plot.width = 10)
p1 = ggplot(data = DC_Matrix)+
    geom_smooth(aes(x=pseudotime,y=RNA),method = "gam",color="red")+theme_bw()+
    ylab("RNA\nExpression")+
    theme(panel.grid.major = element_line(linewidth=0.5,linetype = "dashed",colour = "grey60"),
          axis.title = element_text(size=18,face = "bold"),
          axis.text = element_text(size=14),
          panel.border = element_rect(linewidth = 1.5),
         panel.grid.minor = element_line(linewidth=0.5,linetype = "dashed",colour = "grey80"),
         )
p2 = ggplot(data = DC_Matrix)+
    geom_smooth(aes(x=pseudotime,y=Regulon),method = "gam",color="blue")+theme_bw()+
    ylab("Regulon\nActivity")+
    theme(panel.grid.major = element_line(linewidth=0.5,linetype = "dashed",colour = "grey60"),
          axis.title = element_text(size=18,face = "bold"),
          axis.text = element_text(size=14),
          panel.border = element_rect(linewidth = 1.5),
         panel.grid.minor = element_line(linewidth=0.5,linetype = "dashed",colour = "grey80"),
         )
p3 = ggplot(data = DC_Matrix)+
    geom_smooth(aes(x=pseudotime,y=Chrom),method = "gam",color="yellow")+theme_bw()+
    ylab("TF Binding\nActivity")+
    theme(panel.grid.major = element_line(linewidth=0.5,linetype = "dashed",colour = "grey60"),
          axis.title = element_text(size=18,face = "bold"),
          axis.text = element_text(size=14),
          panel.border = element_rect(linewidth = 1.5),
         panel.grid.minor = element_line(linewidth=0.5,linetype = "dashed",colour = "grey80"),
         )
p1 + p2 + p3


############# pseudotimDE pseudotime-dependent TFs ####################

library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)
library(RColorBrewer)

library(slingshot)
library(SummarizedExperiment)
library(destiny)
library(Biobase)
library(optmatch)
library(Matrix)
library(FNN)
library(dplyr)


tmp_data = readRDS("../04. Discovery in Mesenchymal Cell//data/scATAC_seq_Mesenchymal_chromvar_list.rds")
# clean_Fibro = tmp_data[[1]]
# ss = tmp_data[[2]]
chromVar_res = tmp_data[[3]]
remove(tmp_data)

chromVar_res

RegulonsAUC = read.csv("../04. Discovery in Mesenchymal Cell/SCENIC/OUTPUT/RESULT/RegulonsAUC_scRNA_1.csv",row.names = 1)
RegulonsAUC = as.data.frame(t(RegulonsAUC))

rownames(RegulonsAUC) = sapply(rownames(RegulonsAUC),function(x){gsub("\\.\\.\\.","",x)})

RegulonsAUC

tmp_data = readRDS("../04. Discovery in Mesenchymal Cell//data/scATAC_seq_Mesenchymal_pseudotime_list.rds")
ss = tmp_data[[1]]
st = tmp_data[[2]]
ss_ori_tbl = tmp_data[[3]]
ss_sub_tbl = tmp_data[[4]]
ATAC_res = tmp_data[[5]]
remove(tmp_data)

RNA_data = as.matrix(scRNA@assays$RNA@data[,df_cell_pairing[colnames(ss),"RNA"]])
colnames(RNA_data) = colnames(ss)

TF_candidates = intersect(intersect(intersect(rownames(RNA_data),
                                              rownames(RegulonsAUC)),
                                    ATAC_res$gene),
                          chromVar_res$gene)

length(TF_candidates)

system.time(res <- PseudotimeDE::runPseudotimeDE(gene.vec = TF_candidates,
                                                 assay.use = "logcounts",
                                                 ori.tbl = ss_ori_tbl,
                                                 sub.tbl = ss_sub_tbl, ## To save time, use 100 subsamples
                                                 mat = RNA_data, ## You can also use a matrix or SeuratObj as the input
                                                 model = "gaussian",
                                                 mc.cores = 12))

fp = res$para.pv
names(fp) = res$gene

GE = as.data.frame(-log10(sort(fp) + 10^(-314)))
GE$gene_name = rownames(GE)
colnames(GE)[1] = "-log10(P_value)"
GE$gene_name = factor(GE$gene_name,levels = rev(GE$gene_name))

GE

RegulonsAUC = RegulonsAUC[,df_cell_pairing[colnames(ss),"RNA"]]
colnames(RegulonsAUC) = colnames(ss)

system.time(res <- PseudotimeDE::runPseudotimeDE(gene.vec = TF_candidates,
                                                 assay.use = "logcounts",
                                                 ori.tbl = ss_ori_tbl,
                                                 sub.tbl = ss_sub_tbl, ## To save time, use 100 subsamples
                                                 mat = as.matrix(RegulonsAUC), ## You can also use a matrix or SeuratObj as the input
                                                 model = "gaussian",
                                                 mc.cores = 12))

fp = res$para.pv
names(fp) = res$gene

regulon = as.data.frame(-log10(sort(fp) + 10^(-314)))
regulon$gene_name = rownames(regulon)
colnames(regulon)[1] = "-log10(P_value)"
regulon$gene_name = factor(regulon$gene_name,levels = rev(regulon$gene_name))

regulon

fp = chromVar_res$para.pv
names(fp) = chromVar_res$gene

ChromVar = as.data.frame(-log10(sort(fp) + 10^(-314)))
ChromVar$gene_name = rownames(ChromVar)
colnames(ChromVar)[1] = "-log10(P_value)"
ChromVar$gene_name = factor(ChromVar$gene_name,levels = rev(ChromVar$gene_name))

ChromVar

fp = ATAC_res$para.pv
names(fp) = ATAC_res$gene

ATAC = as.data.frame(-log10(sort(fp) + 10^(-314)))
ATAC$gene_name = rownames(ATAC)
colnames(ATAC)[1] = "-log10(P_value)"
ATAC$gene_name = factor(ATAC$gene_name,levels = rev(ATAC$gene_name))

ATAC

ER = merge(merge(merge(GE,regulon,by="gene_name"),ATAC,by="gene_name"),ChromVar,by="gene_name")

colnames(ER) = c("gene_name","RNA_expression","Regulon_activity","Chromatin accessibility","TF_binding_activity")

ER2 = ER
ER2[,2:4] = ER[,2:4] > -log10(0.05 + 10^(-314))

ER2$gene_name2 = ER2$gene_name
ER2$gene_name2[ rowSums(ER2[,2:4]) < 3 ] = NA

library(plot3D)

options(repr.plot.height = 10, repr.plot.width = 10)
pdf("1.PseudotimeDE_TF_1.pdf",width = 10,height = 10)
scatter3D(ER$RNA_expression,ER$Regulon_activity,ER$`Chromatin accessibility`,pch = 19, cex = 2,bty = "u",phi = 10,col="purple",theta = 30,col.panel ="white", expand =1, col.grid = "grey80",lwd.grid=3)
text3D(ER$RNA_expression,ER$Regulon_activity,ER$`Chromatin accessibility`,labels = ER2$gene_name2,pch = 19, cex = 2,add = TRUE,bty = "b2",phi = 2)
dev.off()

ER2 = ER
ER2[,c(2,3,5)] = ER[,c(2,3,5)] > -log10(0.05 + 10^(-314))

ER2$gene_name2 = ER2$gene_name
ER2$gene_name2[ rowSums(ER2[,c(2,3,5)]) < 3 ] = NA

library(plot3D)

options(repr.plot.height = 15, repr.plot.width = 15)
pdf("1.PseudotimeDE_TF_2.pdf",width = 10,height = 10)
scatter3D(ER$RNA_expression,ER$Regulon_activity,ER$TF_binding_activity,pch = 19, cex = 2,bty = "u",phi = 10,col="purple",theta = 30,col.panel ="white", expand =1, col.grid = "grey80",lwd.grid=3)
text3D(ER$RNA_expression,ER$Regulon_activity,ER$TF_binding_activity,labels = ER2$gene_name2,pch = 19, cex = 2,add = TRUE,bty = "b2",phi = 2)
dev.off()