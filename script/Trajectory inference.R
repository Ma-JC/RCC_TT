library(Seurat)
library(harmony)
library(RColorBrewer)
library(ggplot2)
library(destiny)
library(Biobase)
library(slingshot)
library(SummarizedExperiment)

data = readRDS("..//04. Discovery in Mesenchymal Cell/data/Mesenchymal_clean.rds")

data$All_third_cluster = factor(data$All_third_cluster,
                                labels = c("CD36+ Pericyte 1","CD36+ Pericyte 2","IFI+ Pericyte","Cycling Pericyte","Pericyte",
                                                                    "S100A4+ vSMC 1","S100A4+ vSMC 2","RRAD+ vSMC","ITGA8+ vSMC","FOSB+ vSMC",
                                                                    "ADH1B+ Fibroblast","FBLN5+ Fibroblast","MFAP5+ Fibroblast","CYSLTR2+ Fibroblast","FAP+ Fibroblast",
                                                                    "Mesangial cell"),
                                levels = c("CD36+ Pericyte 1","CD36+ Pericyte 2","IFI+ Pericyte","Cycling Mesenchymal cell","vSMC/Pericyte",
                                                                    "S100A4+ vSMC 1","S100A4+ vSMC 2","RRAD+ vSMC","ITGA8+ vSMC","FOSB+ vSMC",
                                                                    "ADH1B+ Fibroblast","FBLN5+ Fibroblast","MFAP5+ Fibroblast","CYSLTR2+ Fibroblast","FAP+ Fibroblast",
                                                                    "Mesangial cell")
                                )

# data$All_fourth_cluster = as.character(data$All_third_cluster)
data$All_fifth_cluster = ifelse(data$All_third_cluster %in% c("S100A4+ vSMC 1","S100A4+ vSMC 2"),"vSMC 2",
                                ifelse(data$All_third_cluster %in% c("CD36+ Pericyte 1","CD36+ Pericyte 2"),"Pericyte 2",
                                             ifelse(data$All_third_cluster %in% c("RRAD+ vSMC","ITGA8+ vSMC","FOSB+ vSMC"),"vSMC 1",
                                                   ifelse(data$All_third_cluster %in% c("Pericyte"),"Pericyte 1",NA
                                                         )
                                                   )
                                             
                                      )
                                )

data$All_fifth_cluster[ is.na(data$All_fifth_cluster) ] = as.character(data$All_third_cluster)[ is.na(data$All_fifth_cluster) ]

target_cell = c("Mesangial cell","vSMC 1","vSMC 2","Cycling Pericyte","IFI+ Pericyte","Pericyte 1","Pericyte 2",
                "CYSLTR2+ Fibroblast","FAP+ Fibroblast","FBLN5+ Fibroblast","MFAP5+ Fibroblast","ADH1B+ Fibroblast")

data$All_fifth_cluster = factor(data$All_fifth_cluster,levels = target_cell)

data = RunUMAP(object = data,reduction = "harmony",dims = c(1:50),n.components = 2)

n = length(unique(data$All_fifth_cluster))
col = colorRampPalette(brewer.pal(9,"Set1"))(n)
col = c(col[seq(1,n,by = 2)],col[seq(2,n,by = 2)])

names(col) = levels(data$All_fifth_cluster)

options(repr.plot.height = 12, repr.plot.width = 14)
DimPlot(object = data,group.by = "All_fifth_cluster",pt.size = 1,raster = FALSE,label = F,cols = col,label.box = F)+
NoAxes()+
theme(plot.title = element_blank(),legend.position = "top",legend.text = element_text(size = 18,face = "bold"))+
guides(color=guide_legend(override.aes = list(size = 10,text = 1)))

name = colnames(data)[ data$All_fifth_cluster %in% c("Pericyte 2","FAP+ Fibroblast","CYSLTR2+ Fibroblast") ]
raw_matrix = data@assays$RNA@counts
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
data3 = data2 
data2 = RunHarmony(object = data2,group.by.vars = "Patient", plot_convergence = TRUE,max.iter.harmony = 2,dims.use = 1:npcs,theta = 2,lambda =1,reduction = "pca")
data2 = RunHarmony(object = data2,group.by.vars = "Sample_Tag", plot_convergence = TRUE,reduction = "harmony",dims.use = 1:npcs,theta = 2,lambda = 1,max.iter.harmony = 2,verbose = T)

data2 = FindNeighbors(object = data2,reduction = "harmony",dims = c(1:npcs),k.param = 20)
data2 = FindClusters(object = data2,resolution = 2)

data2 = RunUMAP(object = data2,reduction = "harmony",dims = c(1:npcs))

options(repr.plot.height = 8, repr.plot.width = 8)
DimPlot(data2,label=T)
DimPlot(data2,group.by="All_third_cluster")

tmp_data = data2[ ,data2$seurat_clusters %in% c(9,13,3,10,8,5,7,18,16,12,6,15)]

ct = as.ExpressionSet(as.data.frame(tmp_data@reductions$harmony@cell.embeddings))
ct$celltype <- as.character(tmp_data@meta.data[,c("All_fifth_cluster")])
ct$TT_status <- as.character(tmp_data@meta.data[,c("newmeta")])

sigmas <- find_sigmas(ct, verbose = FALSE)
optimal_sigma(sigmas)

dm <- DiffusionMap(ct,sigma = "local",n_eigs = round(optimal_sigma(sigmas), 2))

options(repr.plot.height = 8, repr.plot.width = 8)
palette(col[unique(tmp_data$All_fifth_cluster)])
plot(dm,c(1,2), pch = 20,col_by = "celltype")

tmp_data[["DC"]] = CreateDimReducObject(embeddings = dm@eigenvectors[,1:2],key = "DC_",assay = "RNA")

st = as.SingleCellExperiment(tmp_data)
ss = slingshot(data = st,clusterLabels = "All_fifth_cluster",reducedDim = "DC")

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(ss$slingPseudotime_1, breaks=100)]

plot(reducedDims(ss)$DC, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(ss), lwd=2, col='black')

scatterplot3d::scatterplot3d(dm$DC1,dm$DC2,dm$DC3,plotcol,pch = 20,cex.symbols = 0.1,angle = 30)

pseu_M = cbind(ss$slingPseudotime_1,ss$slingPseudotime_2,ss$slingPseudotime_3,ss$slingPseudotime_4,ss$slingPseudotime_5,ss$slingPseudotime_6,ss$slingPseudotime_7)
ss$all_pseudotime = rowMeans(pseu_M,na.rm = T)

library(dplyr)
group_by(data.frame(clusters = as.character(ss$All_third_cluster),pseudotime = ss$all_pseudotime),clusters) %>% summarise(mean(pseudotime))

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(ss$all_pseudotime, breaks=100)]

plot(reducedDims(ss)$DC, col = plotcol, pch=16, asp = 1)
# lines(SlingshotDataSet(ss), lwd=2, col='black')

scatterplot3d::scatterplot3d(dm$DC1,dm$DC2,dm$DC3,plotcol,pch = 20,cex.symbols = 1,angle = 30)

DC_Matrix = as.data.frame(cbind(as.data.frame(reducedDims(ss)$DC),pseudotime = ss$all_pseudotime,Cell_type=as.character(ss$All_fifth_cluster)))

library(viridis)

options(repr.plot.height = 8, repr.plot.width = 8)
ggplot(DC_Matrix)+geom_point(aes(x=DC_1,y=DC_2,color=Cell_type))+theme_bw()+
scale_color_manual(values = col)+
theme(panel.grid.major = element_line(linewidth=0.5,linetype = "dashed",colour = "grey60"),
      legend.margin = margin(16, 16, 16, 16),
      legend.title = element_text(size=16,face="bold"),
      legend.text = element_text(size=14,face="bold"),
      legend.position = c(0.75,0.85),
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
      legend.position = c(0.8,0.8),
      legend.background = element_rect(fill = "lightgreen",linewidth = 0.5,color="black"),
      axis.title = element_text(size=18,face = "bold"),
      axis.text = element_text(size=14),
      panel.border = element_rect(linewidth = 1.5),
     panel.grid.minor = element_line(linewidth=0.5,linetype = "dashed",colour = "grey80"),
     
     
     )

plot_genebypseudotime = function(Pseudotime,gene_expression,gene_name,Clusters,legend.position=c(0.8,0.8)){
    pseu_matrix = cbind(Pseudotime,gene_expression)
    pseu_matrix = as.data.frame(pseu_matrix)
    colnames(pseu_matrix) = c("Pseudotime","Gene_expression")
    pseu_matrix$Clusters = Clusters
    pseu_matrix = pseu_matrix[ order(pseu_matrix$Pseudotime),]
    
    ggplot(pseu_matrix,aes(x = Pseudotime,y=Gene_expression))+
    scale_color_manual(values = col)+
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

plot_genebypseudotime(Pseudotime = ss$slingPseudotime_1[!is.na(ss$slingPseudotime_1)],
                      gene_expression = as.numeric(tmp_data["CYGB",]@assays$RNA@data)[!is.na(ss$slingPseudotime_1)],
                      gene_name = "CYGB",
                      Clusters = ss$All_fifth_cluster[!is.na(ss$slingPseudotime_1)]
                     )
plot_genebypseudotime(Pseudotime = ss$slingPseudotime_1[!is.na(ss$slingPseudotime_1)],
                      gene_expression = as.numeric(tmp_data["CD248",]@assays$RNA@data)[!is.na(ss$slingPseudotime_1)],
                      gene_name = "CD248",
                      Clusters = ss$All_fifth_cluster[!is.na(ss$slingPseudotime_1)]
                     )
plot_genebypseudotime(Pseudotime = ss$slingPseudotime_1[!is.na(ss$slingPseudotime_1)],
                      gene_expression = as.numeric(tmp_data["MCAM",]@assays$RNA@data)[!is.na(ss$slingPseudotime_1)],
                      gene_name = "MCAM",
                      Clusters = ss$All_fifth_cluster[!is.na(ss$slingPseudotime_1)],
                      legend.position=c(0.8,0.8)
                     )

plot_genebypseudotime(Pseudotime = ss$slingPseudotime_1[!is.na(ss$slingPseudotime_1)],
                      gene_expression = as.numeric(tmp_data["DCN",]@assays$RNA@data)[!is.na(ss$slingPseudotime_1)],
                      gene_name = "DCN",
                      Clusters = ss$All_fifth_cluster[!is.na(ss$slingPseudotime_1)],
                      legend.position=c(0.8,0.2)
                     )

plot_genebypseudotime(Pseudotime = ss$slingPseudotime_1[!is.na(ss$slingPseudotime_1)],
                      gene_expression = as.numeric(tmp_data["LUM",]@assays$RNA@data)[!is.na(ss$slingPseudotime_1)],
                      gene_name = "LUM",
                      Clusters = ss$All_fifth_cluster[!is.na(ss$slingPseudotime_1)]
                     )

options(repr.plot.height = 10, repr.plot.width = 10)
tmp_ss = ss
tmp_st = st[,colnames(tmp_ss)]

plot(reducedDims(tmp_ss)$DC, col =col[ tmp_ss$All_fifth_cluster], pch=16, asp = 1)
lines(SlingshotDataSet(tmp_ss), lwd=3, type = 'lineages', col = 'black')

library(dplyr)
suppressPackageStartupMessages(library(PseudotimeDE))
suppressPackageStartupMessages(library(scales))

ss_ori_tbl <- tibble(cell = colnames(tmp_st), pseudotime = rescale(colData(tmp_ss)$slingPseudotime_1))
tmp_name = colnames(tmp_st)
set.seed( 1234, kind = "L'Ecuyer-CMRG" );
## Set the cores for parallelization. Note that mclapply doesnot work on Windows.
options(mc.cores =6)

## Number of subsmaples
n = 100

## Ganerate random subsamples
ss_index <- mclapply(seq_len(n), function(x) {
  sample(x = c(1:dim(tmp_st)[2]), size = 0.8*dim(tmp_st)[2], replace = FALSE)
})

print(npcs)
ss_sub_tbl <- mclapply(ss_index, function(x) {
    
    tmp_data2 = data3[,tmp_name[x]]
    tmp_data2 = RunHarmony(object = tmp_data2,group.by.vars = "Patient", plot_convergence = TRUE,max.iter.harmony = 2,dims.use = 1:npcs,theta = 2,lambda =1,reduction = "pca")
    tmp_data2 = RunHarmony(object = tmp_data2,group.by.vars = "Sample_Tag", plot_convergence = TRUE,reduction = "harmony",dims.use = 1:npcs,theta = 2,lambda = 1,max.iter.harmony = 2,verbose = T)

    ct2 = as.ExpressionSet(as.data.frame(tmp_data2@reductions$harmony@cell.embeddings))
    ct2$celltype <- tmp_data2@meta.data[,c("All_fifth_cluster")]

    sigmas <- find_sigmas(ct2, verbose = FALSE)


    dm2 <- DiffusionMap(ct2,sigma = "local",n_eigs = round(optimal_sigma(sigmas), 2))

    tmp_data2[["DC"]] = CreateDimReducObject(embeddings = dm2@eigenvectors[,1:2],key = "DC_",assay = "RNA")

    sce = as.SingleCellExperiment(tmp_data2)
    fit = slingshot(data = sce,clusterLabels = "All_fifth_cluster",reducedDim = "DC")
    
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

system.time(res <- PseudotimeDE::runPseudotimeDE(gene.vec = rownames(tmp_st),
                                                 assay.use = "logcounts",
                                                 ori.tbl = ss_ori_tbl,
                                                 sub.tbl = ss_sub_tbl, ## To save time, use 100 subsamples
                                                 mat = tmp_st, ## You can also use a matrix or SeuratObj as the input
                                                 model = "gaussian",
                                                 mc.cores = 6))

print(res)

saveRDS(res,"./Disc_scRNA_TimeGene_new.rds")

res = readRDS("./Disc_scRNA_TimeGene_new.rds")

fp = res$para.pv
names(fp) = res$gene
sort(fp)

library(ComplexHeatmap)
library(circlize)

tep = as.matrix(tmp_data@assays$RNA@data)[names(sort(fp))[1:200],colnames(ss)[order(ss$slingPseudotime_1)]]
tep = t(apply(tep,1,function(x){(x - mean(x))/sd(x)}))
tep[ tep > 2] =2
tep[ tep < -2] = -2

highlight_genes <- c("FAP", "RGS5","CYSLTR2","COL14A1")
highlight_idx <- match(highlight_genes, rownames(tep)) 

row_anno <- rowAnnotation(
  mark = anno_mark(
    at = highlight_idx,
    labels = highlight_genes,
    labels_gp = gpar(fontsize = 10),
    link_width = unit(1, "cm"),
    padding = unit(1, "mm")
  )
)

options(repr.plot.height = 10, repr.plot.width = 10)
hc <- HeatmapAnnotation(show_legend = F,
     annotation_name_side = "right",
     which = "column",
     Pseudotime = meta$slingPseudotime,
     Cell = meta$All_fifth_cluster,
     col = list(Cell = c("Pericyte 2"="#FFAD12","CYSLTR2+ Fibroblast"="#419486","FAP+ Fibroblast"="#B6742A"),
                Pseudotime=colorRamp2(c(min(meta$slingPseudotime),(max(meta$slingPseudotime)-min(meta$slingPseudotime))/2,max(meta$slingPseudotime)),c("#440154", "#21908C", "#FDE725"))
               )
)

ht = Heatmap(tep,
             top_annotation = hc,
             right_annotation = row_anno,
             clustering_method_rows = "ward.D",
             row_names_gp = gpar(fontface = "bold",cex = 1),
             show_column_names = F,
             show_row_names = F,
             cluster_columns = FALSE,
             cluster_rows = TRUE,
             heatmap_legend_param = list(title = "Z-score", legend_height = unit(4, "cm"),
                                    title_gp = gpar(fontsize = 18,fontface = "bold"),
                                    labels_gp = gpar(fontsize = 16,fontface = "bold"),
                                    title_position = "leftcenter-rot",
                                    at = c( -2,0,2))
       
       )

draw(ht,
         heatmap_legend_side ="left"      )
		 
png("PFT_heatmap_dis.png",width = 1000,height = 1000,res = 100)
draw(ht,
         heatmap_legend_side ="left")
dev.off()


plot_genebypseudotime2 = function(Pseudotime,gene_expression,gene_name,Clusters,legend.position=c(0.8,0.8)){
    pseu_matrix = cbind(Pseudotime,gene_expression)
    pseu_matrix = as.data.frame(pseu_matrix)
    colnames(pseu_matrix) = c("Pseudotime","Gene_expression")
    pseu_matrix$Cell_subset = Clusters
    pseu_matrix = pseu_matrix[ order(pseu_matrix$Pseudotime),]
    
    y_range <- max(pseu_matrix$Gene_expression, na.rm = TRUE) - min(pseu_matrix$Gene_expression, na.rm = TRUE)

    ggplot(pseu_matrix,aes(x = Pseudotime,y=Gene_expression))+
    scale_color_manual(values = col)+
    geom_point(aes(color=Cell_subset))+
    geom_smooth(method = "gam")+
    geom_tile(aes(x = Pseudotime,y=-0.5,fill=Cell_subset), height = 0.05*y_range,show.legend = TRUE) +
    ylab(gene_name)+
    theme_bw()+
    scale_fill_manual(values = c("Pericyte 2"="#FFAD12","CYSLTR2+ Fibroblast"="#419486","FAP+ Fibroblast"="#B6742A"))+
    theme(panel.grid.major = element_line(linewidth=0.5,linetype = "dashed",colour = "grey60"),
#       legend.margin = margin(16, 16, 16, 16),
#       legend.title = element_text(size=16,face="bold"),
#       legend.text = element_text(size=14,face="bold"),
#       legend.position = legend.position,
#       legend.background = element_rect(fill = "lightgreen",linewidth = 0.5,color="black"),
          
      axis.title = element_text(size=18,face = "bold"),
      axis.text.y = element_text(size=14),
      axis.text.x = element_text(size=11),
      panel.border = element_rect(linewidth = 1.5),
     panel.grid.minor = element_line(linewidth=0.5,linetype = "dashed",colour = "grey80"),
     
     
     )
}

options(repr.plot.height = 5, repr.plot.width = 5)
plot_list <- list()
for(i in highlight_genes){
    
    p = plot_genebypseudotime2(Pseudotime = ss$slingPseudotime_1[!is.na(ss$slingPseudotime_1)],
                      gene_expression = as.numeric(tmp_data[i,]@assays$RNA@data)[!is.na(ss$slingPseudotime_1)],
                      gene_name = i,
                      Clusters = ss$All_fifth_cluster[!is.na(ss$slingPseudotime_1)]
                     )
    plot_list[[i]] = p
}


library(patchwork)

options(repr.plot.height = 15/2, repr.plot.width = 25/2)
combined_plot <- wrap_plots(plot_list, ncol = 5)+ 
                 plot_layout(guides = "collect")& 
                 theme(legend.position = "bottom",
                       legend.text = element_text(face = "bold",size = 12),
                       legend.title = element_text(face = "bold",size = 18))  # 每行 5 张图
print(combined_plot)

png("PFT_detialed_genes_pseudotime.png",width = 2500/2,height = 1500/2,res = 100)
print(combined_plot)
dev.off()

plot_genebypseudotime2 = function(Pseudotime,gene_expression,gene_name,Clusters,legend.position=c(0.8,0.8)){
    pseu_matrix = cbind(Pseudotime,gene_expression)
    pseu_matrix = as.data.frame(pseu_matrix)
    colnames(pseu_matrix) = c("Pseudotime","Gene_expression")
    pseu_matrix$Cell_subset = Clusters
    pseu_matrix = pseu_matrix[ order(pseu_matrix$Pseudotime),]
    
    y_range <- max(pseu_matrix$Gene_expression, na.rm = TRUE) - min(pseu_matrix$Gene_expression, na.rm = TRUE)

    ggplot(pseu_matrix,aes(x = Pseudotime,y=Gene_expression))+
    scale_color_manual(values = col)+
#     geom_point(aes(color=Cell_subset))+
    geom_smooth(method = "gam")+
    geom_tile(aes(x = Pseudotime,y=-0.5,fill=Cell_subset), height = 0.05*y_range,show.legend = TRUE) +
    ylab(gene_name)+
    theme_bw()+
    scale_fill_manual(values = c("Pericyte 2"="#FFAD12","CYSLTR2+ Fibroblast"="#419486","FAP+ Fibroblast"="#B6742A"))+
    theme(panel.grid.major = element_line(linewidth=0.5,linetype = "dashed",colour = "grey60"),
#       legend.margin = margin(16, 16, 16, 16),
#       legend.title = element_text(size=16,face="bold"),
#       legend.text = element_text(size=14,face="bold"),
#       legend.position = legend.position,
#       legend.background = element_rect(fill = "lightgreen",linewidth = 0.5,color="black"),
          
      axis.title = element_text(size=18,face = "bold"),
      axis.text.y = element_text(size=14),
      axis.text.x = element_text(size=11),
      panel.border = element_rect(linewidth = 1.5),
     panel.grid.minor = element_line(linewidth=0.5,linetype = "dashed",colour = "grey80"),
     
     
     )
}

options(repr.plot.height = 5, repr.plot.width = 5)
i = "CREB3L1"
P1 = plot_genebypseudotime2(Pseudotime = ss$slingPseudotime_1[!is.na(ss$slingPseudotime_1)],
                      gene_expression = as.numeric(tmp_data[i,]@assays$RNA@data)[!is.na(ss$slingPseudotime_1)],
                      gene_name = i,
                      Clusters = ss$All_fifth_cluster[!is.na(ss$slingPseudotime_1)],
                       legend.position= NULL
                     )
i = "MEF2C"
P2 = plot_genebypseudotime2(Pseudotime = ss$slingPseudotime_1[!is.na(ss$slingPseudotime_1)],
                      gene_expression = as.numeric(tmp_data[i,]@assays$RNA@data)[!is.na(ss$slingPseudotime_1)],
                      gene_name = i,
                      Clusters = ss$All_fifth_cluster[!is.na(ss$slingPseudotime_1)],
                       legend.position=NULL
                     )

P1+P2+plot_layout(guides = "collect")& theme(legend.position = "bottom",
                       legend.text = element_text(face = "bold",size = 12),
                       legend.title = element_text(face = "bold",size = 18))