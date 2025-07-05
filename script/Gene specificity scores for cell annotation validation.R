##### our data #########
data = readRDS("/home/ncpsb/mnt/Analysis_for_each_Patient/02.Cell annotation and remove doublets in all patients/discovery_datasets/data/AllData/merge_postprocess.2.rds")
tmp_data = data@assays$RNA@data
tmp_data = tmp_data[rowSums(tmp_data > 0) >= ncol(tmp_data)*0.001,]
gs = sortGenes(tmp_data, data$All_second_cluster,binarizeMethod = "adaptiveMedian")
generank_ourdata = as.data.frame(gs$specScore)
colnames(generank_ourdata) = paste(colnames(generank_ourdata),"(our data)",sep = "")
generank_ourdata = as.data.frame(gs$specScore)

##### Kuppe et al., 2021 Renal fibrosis ######
generank_kidney_fibrosis = read.csv("/home/ncpsb/mnt/Analysis_for_each_Patient/02.Cell annotation and remove doublets in all patients/discovery_datasets/kidney_fibrosis.csv",check.names = F)
generank_kidney_fibrosis = generank_kidney_fibrosis[! duplicated(generank_kidney_fibrosis$gene_name),]
rownames(generank_kidney_fibrosis) = generank_kidney_fibrosis$gene_name
generank_kidney_fibrosis$gene_name = NULL

share_genes = intersect(rownames(generank_kidney_fibrosis),rownames(generank_ourdata))
combind_data = cbind(generank_kidney_fibrosis[share_genes,],generank_ourdata[share_genes,])
cor_matrix1 = cor(combind_data,method = "pearson")

##### Bi et al., 2021RCC ######

generank_ccRCC_immune = read.csv("/home/ncpsb/new_mnt/ccRCC_data/Validation/ccRCC_immune_cancer_cell2021/Bi.generank.csv",check.names = F)
share_genes = intersect(rownames(generank_ccRCC_immune),rownames(generank_ourdata))
combind_data = cbind(generank_ccRCC_immune[share_genes,],generank_ourdata[share_genes,])
cor_matrix2 = cor(combind_data,method = "pearson")

# heatmap 
my_color_range <- colorRamp2(c(1, 0),colors = c("red", "white"))
options(repr.plot.height = 8, repr.plot.width = 15)
ht = Heatmap(cor_matrix1[colnames(generank_ourdata),colnames(generank_kidney_fibrosis)],column_names_rot = 45,
             cluster_rows = F,
             col=my_color_range,
            heatmap_legend_param = list(title = "spearman", 
                                        legend_height = unit(3, "cm"),
#                                     title_gp = gpar(fontsize = 28,fontface = "bold"),
#                                     labels_gp = gpar(fontsize = 24,fontface = "bold"),
                                    title_position = "leftcenter-rot",
                                    at = c( 0,0.5, 1))
            ) + 
     Heatmap(cor_matrix2[colnames(generank_ourdata),colnames(generank_ccRCC_immune)],column_names_rot = 45,
             cluster_rows = F,
#               col=my_color_range,
            heatmap_legend_param = list(title = "spearman", 
                                        legend_height = unit(3, "cm"),
#                                     title_gp = gpar(fontsize = 28,fontface = "bold"),
#                                     labels_gp = gpar(fontsize = 24,fontface = "bold"),
                                    title_position = "leftcenter-rot",
                                    at = c( 0,0.5, 1))
            
            
            )
			
p = draw(ht,heatmap_legend_side ="left")