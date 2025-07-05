library(Seurat)
library(harmony)
library(RColorBrewer)
library(reshape2)
library(RColorBrewer)
library(ggplot2)

data = readRDS("~/mnt/Analysis_for_each_Patient/02.Cell annotation and remove doublets in all patients/discovery_datasets/data/AllData/merge_postprocess.2.rds")


# All data
cpdb_name = vector()
for(i in unique(data$All_second_cluster)){
  tmp_name = colnames(data)[data$All_second_cluster == i]
  if(length(tmp_name) > 500){
    cpdb_name = c(cpdb_name,sample(tmp_name,size = 500,replace = F))
  }else{
    cpdb_name = c(cpdb_name,tmp_name)
  }
}

cpdb_data = data[,cpdb_name]

output = as.data.frame(cpdb_data@assays$RNA@counts)

output = cbind("Gene" = rownames(output),output)

meta = data.frame("Cell"=colnames(cpdb_data),"cell_type"=cpdb_data$All_second_cluster)

if(!dir.exists("CellPhoneDB")){
    dir.create("CellPhoneDB")
}

if(!dir.exists("CellPhoneDB/All")){
    dir.create("CellPhoneDB/All")
}

if(!dir.exists("CellPhoneDB/All/INPUT")){
    dir.create("CellPhoneDB/All/INPUT")
}

if(!dir.exists("CellPhoneDB/All/OUTPUT")){
    dir.create("CellPhoneDB/All/OUTPUT")
}

write.table(output,"CellPhoneDB/All/INPUT/Cellphonedb_All_data.txt",row.names = F,col.names = T,quote = F,sep = "\t")
write.table(meta,"CellPhoneDB/All/INPUT/Cellphonedb_All_meta.txt",row.names = F,col.names = T,quote = F,sep = "\t")

system("/home/ncpsb/anaconda3/envs/cpdb/bin/cellphonedb method statistical_analysis ./CellPhoneDB/All/INPUT/Cellphonedb_All_meta.txt ./CellPhoneDB/All/INPUT/Cellphonedb_All_data.txt --counts-data hgnc_symbol --threads 24 --output-path ./CellPhoneDB/All/OUTPUT > ./CellPhoneDB/All/OUTPUT/cpdb.log")

system("/home/ncpsb/anaconda3/envs/cpdb/bin/cellphonedb plot heatmap_plot ./CellPhoneDB/All/INPUT/Cellphonedb_All_meta.txt --pvalues-path ./CellPhoneDB/All/OUTPUT/pvalues.txt --output-path ./CellPhoneDB/All/OUTPUT/")

# AT from patients with TT

cpdb_name = vector()
for(i in unique(data$All_second_cluster)){
  tmp_name = colnames(data)[data$All_second_cluster == i & data$Sample_Tag == "ccRcc_AT" & data$newmeta == "with_thrombus"]
  if(length(tmp_name) > 500){
    cpdb_name = c(cpdb_name,sample(tmp_name,size = 500,replace = F))
  }else{
    print(i)
    print(length(tmp_name))
    cpdb_name = c(cpdb_name,tmp_name)
  }
}

cpdb_data = data[,cpdb_name]

output = as.data.frame(cpdb_data@assays$RNA@counts)

output = cbind("Gene" = rownames(output),output)

meta = data.frame("Cell"=colnames(cpdb_data),"cell_type"=cpdb_data$All_second_cluster)

if(!dir.exists("CellPhoneDB")){
    dir.create("CellPhoneDB")
}

if(!dir.exists("CellPhoneDB/AT")){
    dir.create("CellPhoneDB/AT")
}

if(!dir.exists("CellPhoneDB/AT/INPUT")){
    dir.create("CellPhoneDB/AT/INPUT")
}

if(!dir.exists("CellPhoneDB/AT/OUTPUT")){
    dir.create("CellPhoneDB/AT/OUTPUT")
}

write.table(output,"CellPhoneDB/AT/INPUT/Cellphonedb_AT_data.txt",row.names = F,col.names = T,quote = F,sep = "\t")
write.table(meta,"CellPhoneDB/AT/INPUT/Cellphonedb_AT_meta.txt",row.names = F,col.names = T,quote = F,sep = "\t")

system("/home/ncpsb/anaconda3/envs/cpdb/bin/cellphonedb method statistical_analysis ./CellPhoneDB/AT/INPUT/Cellphonedb_AT_meta.txt ./CellPhoneDB/AT/INPUT/Cellphonedb_AT_data.txt --counts-data hgnc_symbol --threads 24 --output-path ./CellPhoneDB/AT/OUTPUT > ./CellPhoneDB/AT/OUTPUT/cpdb.log")

system("/home/ncpsb/anaconda3/envs/cpdb/bin/cellphonedb plot heatmap_plot ./CellPhoneDB/AT/INPUT/Cellphonedb_AT_meta.txt --pvalues-path ./CellPhoneDB/AT/OUTPUT/pvalues.txt --output-path ./CellPhoneDB/AT/OUTPUT/")

# PT from patients with TT

cpdb_name = vector()
for(i in unique(data$All_second_cluster)){
  tmp_name = colnames(data)[data$All_second_cluster == i & data$Sample_Tag == "ccRcc_CT" & data$newmeta == "with_thrombus"]
  if(length(tmp_name) > 500){
    cpdb_name = c(cpdb_name,sample(tmp_name,size = 500,replace = F))
  }else{
    cpdb_name = c(cpdb_name,tmp_name)
  }
}

cpdb_data = data[,cpdb_name]

output = as.data.frame(cpdb_data@assays$RNA@counts)

output = cbind("Gene" = rownames(output),output)

meta = data.frame("Cell"=colnames(cpdb_data),"cell_type"=cpdb_data$All_second_cluster)

if(!dir.exists("CellPhoneDB")){
    dir.create("CellPhoneDB")
}

if(!dir.exists("CellPhoneDB/CT")){
    dir.create("CellPhoneDB/CT")
}

if(!dir.exists("CellPhoneDB/CT/INPUT")){
    dir.create("CellPhoneDB/CT/INPUT")
}

if(!dir.exists("CellPhoneDB/CT/OUTPUT")){
    dir.create("CellPhoneDB/CT/OUTPUT")
}

write.table(output,"CellPhoneDB/CT/INPUT/Cellphonedb_CT_data.txt",row.names = F,col.names = T,quote = F,sep = "\t")
write.table(meta,"CellPhoneDB/CT/INPUT/Cellphonedb_CT_meta.txt",row.names = F,col.names = T,quote = F,sep = "\t")

system("/home/ncpsb/anaconda3/envs/cpdb/bin/cellphonedb method statistical_analysis ./CellPhoneDB/CT/INPUT/Cellphonedb_CT_meta.txt ./CellPhoneDB/CT/INPUT/Cellphonedb_CT_data.txt --counts-data hgnc_symbol --threads 24 --output-path ./CellPhoneDB/CT/OUTPUT > ./CellPhoneDB/CT/OUTPUT/cpdb.log")

system("/home/ncpsb/anaconda3/envs/cpdb/bin/cellphonedb plot heatmap_plot ./CellPhoneDB/CT/INPUT/Cellphonedb_CT_meta.txt --pvalues-path ./CellPhoneDB/CT/OUTPUT/pvalues.txt --output-path ./CellPhoneDB/CT/OUTPUT/")

# TT
cpdb_name = vector()
for(i in unique(data$All_second_cluster)){
  tmp_name = colnames(data)[data$All_second_cluster == i & data$Sample_Tag == "ccRcc_TT" & data$newmeta == "with_thrombus"]
  if(length(tmp_name) > 500){
    cpdb_name = c(cpdb_name,sample(tmp_name,size = 500,replace = F))
  }else{
    cpdb_name = c(cpdb_name,tmp_name)
  }
}

cpdb_data = data[,cpdb_name]

output = as.data.frame(cpdb_data@assays$RNA@counts)

output = cbind("Gene" = rownames(output),output)

meta = data.frame("Cell"=colnames(cpdb_data),"cell_type"=cpdb_data$All_second_cluster)

if(!dir.exists("CellPhoneDB")){
    dir.create("CellPhoneDB")
}

if(!dir.exists("CellPhoneDB/TT")){
    dir.create("CellPhoneDB/TT")
}

if(!dir.exists("CellPhoneDB/TT/INPUT")){
    dir.create("CellPhoneDB/TT/INPUT")
}

if(!dir.exists("CellPhoneDB/TT/OUTPUT")){
    dir.create("CellPhoneDB/TT/OUTPUT")
}

write.table(output,"CellPhoneDB/TT/INPUT/Cellphonedb_TT_data.txt",row.names = F,col.names = T,quote = F,sep = "\t")
write.table(meta,"CellPhoneDB/TT/INPUT/Cellphonedb_TT_meta.txt",row.names = F,col.names = T,quote = F,sep = "\t")

system("/home/ncpsb/anaconda3/envs/cpdb/bin/cellphonedb method statistical_analysis ./CellPhoneDB/TT/INPUT/Cellphonedb_TT_meta.txt ./CellPhoneDB/TT/INPUT/Cellphonedb_TT_data.txt --counts-data hgnc_symbol --threads 24 --output-path ./CellPhoneDB/TT/OUTPUT > ./CellPhoneDB/TT/OUTPUT/cpdb.log")

system("/home/ncpsb/anaconda3/envs/cpdb/bin/cellphonedb plot heatmap_plot ./CellPhoneDB/TT/INPUT/Cellphonedb_TT_meta.txt --pvalues-path ./CellPhoneDB/TT/OUTPUT/pvalues.txt --output-path ./CellPhoneDB/TT/OUTPUT/")

# PT from patients without TT

cpdb_name = vector()
for(i in unique(data$All_second_cluster)){
  tmp_name = colnames(data)[data$All_second_cluster == i & data$Sample_Tag == "ccRcc_CT" & data$newmeta == "without_thrombus"]
  if(length(tmp_name) > 500){
    cpdb_name = c(cpdb_name,sample(tmp_name,size = 500,replace = F))
  }else{
    cpdb_name = c(cpdb_name,tmp_name)
  }
}

cpdb_data = data[,cpdb_name]

output = as.data.frame(cpdb_data@assays$RNA@counts)

output = cbind("Gene" = rownames(output),output)

meta = data.frame("Cell"=colnames(cpdb_data),"cell_type"=cpdb_data$All_second_cluster)

if(!dir.exists("CellPhoneDB")){
    dir.create("CellPhoneDB")
}

if(!dir.exists("CellPhoneDB/CT_withouTT")){
    dir.create("CellPhoneDB/CT_withouTT")
}

if(!dir.exists("CellPhoneDB/CT_withouTT/INPUT")){
    dir.create("CellPhoneDB/CT_withouTT/INPUT")
}

if(!dir.exists("CellPhoneDB/CT_withouTT/OUTPUT")){
    dir.create("CellPhoneDB/CT_withouTT/OUTPUT")
}

write.table(output,"CellPhoneDB/CT_withouTT/INPUT/Cellphonedb_CT_withouTT_data.txt",row.names = F,col.names = T,quote = F,sep = "\t")
write.table(meta,"CellPhoneDB/CT_withouTT/INPUT/Cellphonedb_CT_withouTT_meta.txt",row.names = F,col.names = T,quote = F,sep = "\t")

system("/home/ncpsb/anaconda3/envs/cpdb/bin/cellphonedb method statistical_analysis ./CellPhoneDB/CT_withouTT/INPUT/Cellphonedb_CT_withouTT_meta.txt ./CellPhoneDB/CT_withouTT/INPUT/Cellphonedb_CT_withouTT_data.txt --counts-data hgnc_symbol --threads 24 --output-path ./CellPhoneDB/CT_withouTT/OUTPUT > ./CellPhoneDB/CT_withouTT/OUTPUT/cpdb.log")

system("/home/ncpsb/anaconda3/envs/cpdb/bin/cellphonedb plot heatmap_plot ./CellPhoneDB/CT_withouTT/INPUT/Cellphonedb_CT_withouTT_meta.txt --pvalues-path ./CellPhoneDB/CT_withouTT/OUTPUT/pvalues.txt --output-path ./CellPhoneDB/CT_withouTT/OUTPUT/")

# Interaction ranking
library(dplyr)
library(ComplexHeatmap)

tmp = read.table("./CellPhoneDB/All/OUTPUT/count_network.txt",header = T,sep = "\t")
cm = reshape2::dcast(data = tmp,formula = SOURCE~TARGET,value.var = "count") %>% tibble::column_to_rownames("SOURCE")

m = names(sort(rowSums(cm),decreasing = T))
cm = cm[m,m]

hc <- HeatmapAnnotation(
     which = "column",
  total = anno_barplot(
    rowSums(cm),
    gp = gpar(
      fill = col[colnames(cm)]
    ))
)

hr <- HeatmapAnnotation(
    which = "row",
  total = anno_barplot(
    rowSums(cm),
    gp = gpar(
      fill = col[rownames(cm)]
    ))
)

options(repr.plot.height = 8, repr.plot.width = 10)
Heatmap(cm,top_annotation = hc,right_annotation = hr,cluster_rows = F,cluster_columns = F)

########## differential cellular crosstalk analysis ############

cpdb_AT_w = read.table("./CellPhoneDB/AT/OUTPUT/significant_means.txt",header = T,check.names = F,stringsAsFactors = F,quote = "",sep = "\t")
cpdb_CT_w = read.table("./CellPhoneDB/CT//OUTPUT/significant_means.txt",header = T,check.names = F,stringsAsFactors = F,quote = "",sep = "\t")
cpdb_TT_w = read.table("./CellPhoneDB/TT/OUTPUT/significant_means.txt",header = T,check.names = F,stringsAsFactors = F,quote = "",sep = "\t")
cpdb_CT_w0 = read.table("./CellPhoneDB/CT_withouTT//OUTPUT/significant_means.txt",header = T,check.names = F,stringsAsFactors = F,quote = "",sep = "\t")

cpdb_AT_w = cpdb_AT_w[ !duplicated(cpdb_AT_w$interacting_pair),]
cpdb_CT_w = cpdb_CT_w[ !duplicated(cpdb_CT_w$interacting_pair),]
cpdb_TT_w = cpdb_TT_w[ !duplicated(cpdb_TT_w$interacting_pair),]
cpdb_CT_w0 = cpdb_CT_w0[ !duplicated(cpdb_CT_w0$interacting_pair),]

cpdb_AT_w = cpdb_AT_w[,!grepl("Epithelial|Mesan|Vasa|Glomer|Vein",colnames(cpdb_AT_w))]
cpdb_CT_w = cpdb_CT_w[,!grepl("Epithelial|Mesan|Vasa|Glomer|Vein",colnames(cpdb_CT_w))]
cpdb_TT_w = cpdb_TT_w[,!grepl("Epithelial|Mesan|Vasa|Glomer|Vein",colnames(cpdb_TT_w))]
cpdb_CT_w0 = cpdb_CT_w0[,!grepl("Epithelial|Mesan|Vasa|Glomer|Vein",colnames(cpdb_CT_w0))]

colnames(cpdb_AT_w) = paste(colnames(cpdb_AT_w),"_AT_w",sep = "")
colnames(cpdb_CT_w) = paste(colnames(cpdb_CT_w),"_CT_w",sep = "")
colnames(cpdb_TT_w) = paste(colnames(cpdb_TT_w),"_TT_w",sep = "")
colnames(cpdb_CT_w0) = paste(colnames(cpdb_CT_w0),"_CT_wo",sep = "")

cpdb = merge(cpdb_AT_w,cpdb_CT_w,by.x = "interacting_pair_AT_w",by.y = "interacting_pair_CT_w",all = TRUE)
cpdb = merge(cpdb,cpdb_TT_w,by.x = "interacting_pair_AT_w",by.y = "interacting_pair_TT_w",all = TRUE)
cpdb = merge(cpdb,cpdb_CT_w0,by.x = "interacting_pair_AT_w",by.y = "interacting_pair_CT_wo",all = TRUE)

cpdb_AT_w_2 = cpdb[,colnames(cpdb_AT_w)[13:ncol(cpdb_AT_w)]]
cpdb_CT_w_2 = cpdb[,colnames(cpdb_CT_w)[13:ncol(cpdb_CT_w)]]
cpdb_TT_w_2 = cpdb[,colnames(cpdb_TT_w)[13:ncol(cpdb_TT_w)]]
cpdb_CT_w0_2 = cpdb[,colnames(cpdb_CT_w0)[13:ncol(cpdb_CT_w0)]]

res = vector()
for(i in 1:ncol(cpdb_CT_w0_2)){
    res = c(res,sum((is.na(cpdb_CT_w0_2[,i])+is.na(cpdb_CT_w_2[,i])) == 1 & (is.na(cpdb_CT_w0_2[,i])-is.na(cpdb_CT_w_2[,i])) == 1))
}
names(res) = colnames(cpdb_AT_w_2)
sort(res,decreasing = T)

cell_type_re = c('Cycling T cell','CD8\\+ T cell','NK cell','CD4\\+ T cell','CD4\\+ Treg cell',
                 'Macrophage','Monocyte','Dendritic cell','Cycling Myeloid cell','Neutrophil','Mast cell','B cell','Plasma cell',
                 'Endothelial cell 1\\(Glomerular Capillaries\\)','Endothelial cell 2\\(Arterioles\\)',
                 'Endothelial cell 3\\(Vein\\)','Endothelial cell 4\\(Vasa Recta\\)','Endothelial cell 5\\(Cancer\\)','Endothelial cell 6\\(Cancer\\)',
                 'Pericyte','vSMC','Fibroblast','Mesangial cell','Epithelial cell','Cancer cell')

res2 = vector()
for(i in cell_type_re){
    res2 = c(res2,sum(res[grepl(paste(i,"",sep = ""),names(res))]))
}

cell_type = c('Cycling T cell','CD8+ T cell','NK cell','CD4+ T cell','CD4+ Treg cell',
                 'Macrophage','Monocyte','Dendritic cell','Cycling Myeloid cell','Neutrophil','Mast cell','B cell','Plasma cell',
                 'Endothelial cell 1(Glomerular Capillaries)','Endothelial cell 2(Arterioles)',
                 'Endothelial cell 3(Vein)','Endothelial cell 4(Vasa Recta)','Endothelial cell 5(Cancer)','Endothelial cell 6(Cancer)',
                 'Pericyte','vSMC','Fibroblast','Mesangial cell','Epithelial cell','Cancer cell')
names(res2) = cell_type

sort(res2)

library(igraph)
cross_data = cbind(as.data.frame(do.call(rbind,sapply(names(res),function(x){strsplit(x = x,"\\||_AT_w")}))),res)
colnames(cross_data) = c("from","to","weight")
cross_data = cross_data[ cross_data$weight !=0,]

net_pc<-graph_from_data_frame(d = cross_data,directed=TRUE)

weighted_degrees <- numeric(vcount(net_pc))

# 遍历每个节点，计算加权度
for (i in 1:vcount(net_pc)) {

  edges <- incident(net_pc, i, mode = "all")

  weighted_degrees[i] <- sum(E(net_pc)$weight[edges[-i]])
}

names(weighted_degrees) = V(net_pc)$name

cross_data$from = factor(cross_data$from,names(sort(weighted_degrees,decreasing = T)))
cross_data = cross_data[ order(cross_data$from),]
cross_data$from = as.character(cross_data$from)

net_pc<-graph_from_data_frame(d = cross_data,directed=TRUE)

weighted_degrees <- numeric(vcount(net_pc))

for (i in 1:vcount(net_pc)) {
  
  edges <- incident(net_pc, i, mode = "all")
  
  weighted_degrees[i] <- sum(E(net_pc)$weight[edges[-i]])
}

names(weighted_degrees) = V(net_pc)$name

weighted_degrees

options(repr.plot.height = 10, repr.plot.width = 13)
plot(net_pc,vertex.size=weighted_degrees/25,vertex.label.cex=1.2,vertex.label.dist=0,vertex.label.color="black",
     layout = layout_in_circle(net_pc),
     edge.width = E(net_pc)$weight*0.2,edge.curved=0.2,
     edge.arrow.size = 0.2,
     edge.color=col[cross_data$from],
     vertex.size = 5,
     vertex.color=col[V(net_pc)$name],
    )

title("Differential cellular crosstalk in primary tumor\n(with TT vs without TT)", line = 0.5,adj = 0.5, cex.main = 2, font.main = 2)

options(repr.plot.height = 5, repr.plot.width = 10)
par(mar = c(15, 6, 1, 2.1))
barplot(weighted_degrees,las = 2,col = col[names(weighted_degrees)],ylab = "Weighted degree")