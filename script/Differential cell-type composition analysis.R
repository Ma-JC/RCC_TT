############## Figure 1D & Figure S2E###################
library(Seurat)
library(harmony)
library(reshape2)
library(RColorBrewer)
library(ggplot2)

# PT(with TT) vs PT(without TT)
data = readRDS("~/mnt/Analysis_for_each_Patient/02.Cell annotation and remove doublets in all patients/discovery_datasets/data/AllData/merge_postprocess.2.rds")

sub_data = data[,data$Sample_Tag == "ccRcc_CT"]

tmp_res = vector()

group = list("CAF" = levels(data$All_second_cluster)[20:23],
    "EC" = levels(data$All_second_cluster)[14:19],
    "Epi_cancer"=levels(data$All_second_cluster)[24:25],
    "immune" =c(levels(data$All_second_cluster)[c(1:13)]))

for(i in names(group)){
   

    res = reshape2::melt(table(sub_data$Patient,sub_data$All_second_cluster)[,group[[i]]]/rowSums(table(sub_data$Patient,sub_data$All_second_cluster)[,group[[i]]]))
    res$Var1 = as.character(res$Var1)
    res$Var2 = as.character(res$Var2)

    tmp = reshape2::melt(table(sub_data$Patient,sub_data$newmeta))
    tmp = tmp[tmp$value != 0,]
    tmp$value = NULL

    tmp$Var1 = as.character(tmp$Var1)
    tmp$Var2 = as.character(tmp$Var2)

    rownames(tmp) = tmp$Var1

    res$group = tmp[ res$Var1,"Var2"]
    
    tmp_res = rbind(tmp_res,res)
}

res = tmp_res
res[is.na(res$value),"value"] = 0

pr = vector()
for(i in unique(res$Var2)){
    tmp = res[ res$Var2 == i,]
    pr = c(pr,t.test(value~group,data = tmp)$p.value)
}

pr[is.na(pr)] = 1

names(pr) = unique(res$Var2)
res$Var2 = factor(res$Var2,levels = names(pr))



#################### a refined analysis within a strictly defined clear cell RCC cohort #######################
new_res = res
new_res$mayo = ifelse(new_res$Var1 %in% c("with_thrombus_Patient14","without_thrombus_Patient06","without_thrombus_Patient07"),"0_level","II_level")
new_res$mayo[ new_res$group == "without_thrombus" ] = "CANCER"
new_res$mayo = factor(new_res$mayo,levels=c("CANCER","0_level","II_level"))

new_res$cancer_type = ifelse(new_res$Var1 %in% c("with_thrombus_Patient15","without_thrombus_Patient07"),"pRCC",ifelse(new_res$Var1 %in% c("with_thrombus_Patient17"),"sRCC","ccRCC"))
new_res$cancer_type2 = ifelse(new_res$Var1 %in% c("with_thrombus_Patient15","without_thrombus_Patient07","with_thrombus_Patient17"),"no-ccRCC","ccRCC")

new_res$metastasis = ifelse(new_res$Var1 %in% c("with_thrombus_Patient14"),"N1M1",ifelse(new_res$Var1 %in% c("with_thrombus_Patient17"),"N1M0","N0M0"))
new_res$metastasis2 = ifelse(new_res$Var1 %in% c("with_thrombus_Patient14","with_thrombus_Patient17"),"no-N0M0","N0M0")

new_res$sex = ifelse(new_res$Var1 %in% c("with_thrombus_Patient08","with_thrombus_Patient18"),"Female","Male")

summary(lm(value~group+mayo+cancer_type2+metastasis2+sex,data = new_res[new_res$Var2 == "NK cell" ,] ))
summary(lm(value~group+mayo+cancer_type2+metastasis2+sex,data = new_res[new_res$Var2 == "Fibroblast" ,]))

new_res[new_res$Var2 == "NK cell" & 
        new_res$mayo != "0_level" & 
        new_res$cancer_type == "ccRCC" & 
        new_res$metastasis == "N0M0" & 
        new_res$sex == "Male",]

library(ggplot2)
MyData_NK_clean = new_res[new_res$Var2 == "NK cell" & 
                      new_res$mayo != "0_level" & 
                      new_res$cancer_type == "ccRCC" & 
                      new_res$metastasis == "N0M0" & 
                      new_res$sex == "Male",]
MyData_NK_clean$group = factor(MyData_NK_clean$group,levels = c("without_thrombus","with_thrombus"),labels = c("without_TT","with_TT"))
options(repr.plot.height = 5, repr.plot.width = 6.5)
p1 = ggplot(data = MyData_NK_clean)+
    geom_boxplot(aes(x = group,y = value,color = group),outlier.colour = NA)+
    ggtitle("scRNA-seq")+
    ylab("Cell proportionof NK cell")+
    labs(color="TT status")+
    geom_point(aes(x = group,y = value,color=group),pch = 20,size = 3)+
    theme_bw()+
    scale_color_manual(values = c("blue","red"))+
    theme(panel.border = element_rect(linewidth = 2),
                      legend.text = element_text(size = 15),
                      legend.title = element_text(size = 18,face="bold"),
                      plot.title = element_text(size=28,face = "bold",hjust = 0.5),
                      axis.text.y = element_text(size = 16),
                      axis.text.x = element_text(size = 18,face = "bold"),
                      axis.title = element_text(face = "bold",size = 18),
                      axis.title.x = element_blank(),
                      panel.grid = element_line(colour = "grey66",linetype = "dashed"),

                     )+
    annotate(geom = "text",x = 1.5,y = 0.3,size=6,fontface="bold",
             label = paste("P.value:",round(t.test(value~group,data = MyData_NK_clean)$p.value,3))
            )
p1

library(ggplot2)
MyData_Fibro_clean = new_res[new_res$Var2 == "Fibroblast" & 
                      new_res$mayo != "0_level" & 
                      new_res$cancer_type == "ccRCC" & 
                      new_res$metastasis == "N0M0" & 
                      new_res$sex == "Male",]
MyData_Fibro_clean$group = factor(MyData_Fibro_clean$group,levels = c("without_thrombus","with_thrombus"),labels = c("without_TT","with_TT"))
options(repr.plot.height = 5, repr.plot.width = 6.5)
p2 = ggplot(data = MyData_Fibro_clean)+
    geom_boxplot(aes(x = group,y = value,color = group),outlier.colour = NA)+
    ggtitle("scRNA-seq")+
    ylab("Cell proportionof Fibroblast")+
    labs(color="TT status")+
    geom_point(aes(x = group,y = value,color=group),pch = 20,size = 3)+
    theme_bw()+
    scale_color_manual(values = c("blue","red"))+
#     scale_fill_manual(values = c("blue","red"))+
    theme(panel.border = element_rect(linewidth = 2),
                      legend.text = element_text(size = 15),
                      legend.title = element_text(size = 18,face="bold"),
#                       legend.position = "none",
                      plot.title = element_text(size=28,face = "bold",hjust = 0.5),
                      axis.text.y = element_text(size = 16),
                      axis.text.x = element_text(size = 18,face = "bold"),
                      axis.title = element_text(face = "bold",size = 18),
                      axis.title.x = element_blank(),
                      panel.grid = element_line(colour = "grey66",linetype = "dashed"),

                     )+
#     theme(axis.title = element_blank(),axis.text = element_text(size = 16,face="bold"),legend.text = element_text(size=16),legend.title = element_text(face="bold",size=18))+
    annotate(geom = "text",x = 1.5,y = 0.6,size=6,fontface="bold",
             label = paste("P.value:",round(t.test(value~group,data = MyData_Fibro_clean)$p.value,3))
            )
p2

options(repr.plot.height = 5, repr.plot.width = 11)
library(patchwork)
p1+p2+plot_layout(guides = "collect")
pdf(file = "clean_comparison.pdf",width = 11,height = 5)
p1+p2+plot_layout(guides = "collect")
dev.off()
####################################################################################################

col = scales::brewer_pal(palette = "Set1")(9)

library(ggplot2)
p = ggplot(data = res)+geom_boxplot(aes(x = Var2,y = value,color = group))+
geom_point(aes(x = Var2,y = value,color = group),position = position_jitterdodge(0.1))+
    theme_bw()+
    scale_color_manual(values = c(col[1],col[2]))+
    labs(x ="Cell type",y="The proportion of each cell type in each patient",color = " ")+
    theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
          legend.position = "top",
          axis.text.x = element_text(face = "bold",size = 16,angle = 45,hjust = 1),
          axis.text.y = element_text(face = "bold",size = 15),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18),
          axis.title = element_text(size = 18))

 for(i in  1:length(pr)){
     tmp = res[res$Var2 == names(pr)[i],]
     gene_max = max(tmp$value)
     p = p + annotate("text",i-0.5,gene_max,label = paste("p.value:",round(pr[i],4)),size = 5,hjust = 0) 
 }

options(repr.plot.height = 10, repr.plot.width = 30)
p

res2 = vector()
for(i in names(pr)){
    a = mean(res[ res$Var2 == i & res$group =="with_thrombus","value"])
    b = mean(res[ res$Var2 == i & res$group =="without_thrombus","value"])
    res2 = c(res2,(a-b)/(a+b) )
}
pr[is.na(pr)] = 1
res3 = as.data.frame(cbind(pr,res2))
colnames(res3) = c("p.value","FC")
res3$FC[is.na(res3$FC)] = 0
res3$p.value = -log10(res3$p.value)
res3$label = rownames(res3)
res3$color = c(rep("Mesenchymal cell",length(group[[1]])),rep("Endothelial cell",length(group[[2]])),rep("Epithelial cell",length(group[[3]])),rep("Immune cell",length(group[[4]])))

res3$color[ res3$p.value < -log10(0.05)] = "No.significance"

options(repr.plot.height = 6, repr.plot.width = 12)
ggplot(data = res3)+
ylab("-log10(p.value)")+
theme_classic()+
ylim(0,2.5)+
xlim(-1,1)+
geom_hline(aes(yintercept=-log10(0.05)),linetype = "dashed")+
geom_vline(aes(xintercept=0),linetype = "dashed")+
geom_label(aes(x = FC,y = p.value,label = label,fill = color,size = p.value),color = "white")+
scale_fill_manual(values = c(scales::brewer_pal(palette = "Set1")(4)[1],scales::brewer_pal(palette = "Set1")(4)[4],"grey"))+
guides(size=guide_legend(label = 0))+
theme(axis.title = element_text(size = 18,face = "bold"),
      legend.position = "top",
      axis.text = element_text(size = 12,face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 18,face = "bold"))

tmp = res[ res$Var2 %in% c("NK cell","Fibroblast"),]
tmp$Var2 = as.character(tmp$Var2)

p = ggplot(data = tmp)+geom_boxplot(aes(x = Var2,y = value,color = group))+
geom_point(aes(x = Var2,y = value,color = group),position = position_jitterdodge(0.1))+
    theme_bw()+
    scale_color_manual(values = c(col[1],col[2]))+
    labs(x ="Cell type",y="The proportion of each cell type in each patient",color = " ")+
    theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
          legend.position = "top",
          axis.text.x = element_text(face = "bold",size = 16),
          axis.text.y = element_text(face = "bold",size = 15),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18),
          axis.title = element_text(size = 18))
n = 1
 for(i in  c("Fibroblast","NK cell")){
    
     tmp = res[res$Var2 == i,]
     gene_max = max(tmp$value)
     p = p + annotate("text",n-0.5,gene_max+0.05,label = paste("p.value:",round(pr[i],4)),size = 5,hjust = 0) 
     n = n + 1
 }

options(repr.plot.height = 6, repr.plot.width = 8)
p

res3$VS = "with_TT VS without TT in CT"
WvsWOinCT = res3


# PT(with TT) vs PT(without TT)
sub_data = data[,data$Sample_Tag == "ccRcc_AT"]

tmp_res = vector()

group = list("CAF" = levels(data$All_second_cluster)[20:23],
    "EC" = levels(data$All_second_cluster)[14:19],
    "Epi_cancer"=levels(data$All_second_cluster)[24:25],
    "immune" =c(levels(data$All_second_cluster)[1:13]))

for(i in names(group)){
   

    res = reshape2::melt(table(sub_data$Patient,sub_data$All_second_cluster)[,group[[i]]]/rowSums(table(sub_data$Patient,sub_data$All_second_cluster)[,group[[i]]]))
    res$Var1 = as.character(res$Var1)
    res$Var2 = as.character(res$Var2)

    tmp = reshape2::melt(table(sub_data$Patient,sub_data$newmeta))
    tmp = tmp[tmp$value != 0,]
    tmp$value = NULL

    tmp$Var1 = as.character(tmp$Var1)
    tmp$Var2 = as.character(tmp$Var2)

    rownames(tmp) = tmp$Var1

    res$group = tmp[ res$Var1,"Var2"]
    
    tmp_res = rbind(tmp_res,res)
}

res = tmp_res
res[is.na(res$value),"value"] = 0

pr = vector()
for(i in unique(res$Var2)){
    tmp = res[ res$Var2 == i,]
    pr = c(pr,t.test(value~group,data = tmp)$p.value)
}

pr[is.na(pr)] = 1

names(pr) = unique(res$Var2)
res$Var2 = factor(res$Var2,levels = names(pr))


col = scales::brewer_pal(palette = "Set1")(9)

library(ggplot2)
p = ggplot(data = res)+geom_boxplot(aes(x = Var2,y = value,color = group))+
geom_point(aes(x = Var2,y = value,color = group),position = position_jitterdodge(0.1))+
    theme_bw()+
    scale_color_manual(values = c(col[1],col[2]))+
    labs(x ="Cell type",y="The proportion of each cell type in each patient",color = " ")+
    theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
          legend.position = "top",
          axis.text.x = element_text(face = "bold",size = 16,angle = 45,hjust = 1),
          axis.text.y = element_text(face = "bold",size = 15),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18),
          axis.title = element_text(size = 18))

 for(i in  1:length(pr)){
     tmp = res[res$Var2 == names(pr)[i],]
     gene_max = max(tmp$value)
     p = p + annotate("text",i-0.5,gene_max,label = paste("p.value:",round(pr[i],4)),size = 5,hjust = 0) 
 }

options(repr.plot.height = 10, repr.plot.width = 30)
p

res2 = vector()
for(i in names(pr)){
    a = mean(res[ res$Var2 == i & res$group =="with_thrombus","value"])
    b = mean(res[ res$Var2 == i & res$group =="without_thrombus","value"])
    res2 = c(res2,(a-b)/(a+b) )
}
pr[is.na(pr)] = 1
res3 = as.data.frame(cbind(pr,res2))
colnames(res3) = c("p.value","FC")
res3$FC[is.na(res3$FC)] = 0
res3$p.value = -log10(res3$p.value)
res3$label = rownames(res3)
res3$color = c(rep("Mesenchymal cell",length(group[[1]])),rep("Endothelial cell",length(group[[2]])),rep("Epithelial cell",length(group[[3]])),rep("Immune cell",length(group[[4]])))

res3$color[ res3$p.value < -log10(0.05)] = "No.significance"

options(repr.plot.height = 6, repr.plot.width = 12)
ggplot(data = res3)+
ylab("-log10(p.value)")+
theme_classic()+
ylim(0,2.5)+
xlim(-1,1)+
geom_hline(aes(yintercept=-log10(0.05)),linetype = "dashed")+
geom_vline(aes(xintercept=0),linetype = "dashed")+
geom_label(aes(x = FC,y = p.value,label = label,fill = color,size = p.value),color = "white")+
scale_fill_manual(values = c("grey"))+
guides(size=guide_legend(label = 0))+
theme(axis.title = element_text(size = 18,face = "bold"),
      legend.position = "top",
      axis.text = element_text(size = 12,face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 18,face = "bold"))

res3$VS = "with_TT VS without TT in AT"
WvsWOinAT = res3

# PT vs TT from patients with TT
sub_data = data[,data$newmeta == "with_thrombus" & data$Sample_Tag != "ccRcc_AT" & !data$Patient  %in% c("without_thrombus_Patient06","without_thrombus_Patient07")]
subtmp_CT = data[,data$newmeta == "with_thrombus" & data$Sample_Tag == "ccRcc_CT" & !data$Patient  %in% c("without_thrombus_Patient06","without_thrombus_Patient07")]
subtmp_TT = data[,data$newmeta == "with_thrombus" & data$Sample_Tag == "ccRcc_TT" & !data$Patient  %in% c("without_thrombus_Patient06","without_thrombus_Patient07")]

tmp_res = vector()

group = list("CAF" = levels(data$All_second_cluster)[20:23],
    "EC" = levels(data$All_second_cluster)[14:19],
    "Epi_cancer"=levels(data$All_second_cluster)[24:25],
    "immune" =c(levels(data$All_second_cluster)[1:13]))

for(i in names(group)){
   
    for(j in c("CT","TT")){
            if(j == "CT"){
                res = reshape2::melt(table(subtmp_CT$Patient,subtmp_CT$All_second_cluster)[,group[[i]]]/rowSums(table(subtmp_CT$Patient,subtmp_CT$All_second_cluster)[,group[[i]]]))
                res$group = "CT"
            }else{
                res = reshape2::melt(table(subtmp_TT$Patient,subtmp_TT$All_second_cluster)[,group[[i]]]/rowSums(table(subtmp_TT$Patient,subtmp_TT$All_second_cluster)[,group[[i]]]))
                res$group = "TT"
            }
            res$Var1 = as.character(res$Var1)
            res$Var2 = as.character(res$Var2)

            tmp_res = rbind(tmp_res,res)
    }

}


res = tmp_res
res[is.na(res$value),"value"] = 0

pr = vector()
for(i in unique(res$Var2)){
    tmp = res[ res$Var2 == i,]
    tmp_CT = tmp[tmp$group == "CT","value"]
    tmp_TT = tmp[tmp$group == "TT","value"]
    pr = c(pr,t.test(tmp_CT,tmp_TT,paired = TRUE,)$p.value)
}

pr[is.na(pr)] = 1

names(pr) = unique(res$Var2)
res$Var2 = factor(res$Var2,levels = names(pr))

col = scales::brewer_pal(palette = "Set1")(9)

library(ggplot2)
p = ggplot(data = res)+geom_boxplot(aes(x = Var2,y = value,color = group))+
geom_point(aes(x = Var2,y = value,color = group),position = position_jitterdodge(0.1))+
    theme_bw()+
    scale_color_manual(values = c(col[1],col[2]))+
    labs(x ="Cell type",y="The proportion of each cell type in each patient",color = " ")+
    theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
          legend.position = "top",
          axis.text.x = element_text(face = "bold",size = 16,angle = 45,hjust = 1),
          axis.text.y = element_text(face = "bold",size = 15),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18),
          axis.title = element_text(size = 18))

 for(i in  1:length(pr)){
     tmp = res[res$Var2 == names(pr)[i],]
     gene_max = max(tmp$value)
     p = p + annotate("text",i-0.5,gene_max,label = paste("p.value:",round(pr[i],4)),size = 5,hjust = 0) 
 }

options(repr.plot.height = 10, repr.plot.width = 30)
p

res2 = vector()
for(i in names(pr)){
    a = mean(res[ res$Var2 == i & res$group =="TT","value"])
    b = mean(res[ res$Var2 == i & res$group =="CT","value"])
    res2 = c(res2 ,(a-b)/(a+b))
}
pr[is.na(pr)] = 1
res3 = as.data.frame(cbind(pr,res2))
colnames(res3) = c("p.value","FC")
res3$FC[is.na(res3$FC)] = 0
res3$p.value = -log10(res3$p.value)
res3$label = rownames(res3)
res3$color = c(rep("Mesenchymal cell",length(group[[1]])),rep("Endothelial cell",length(group[[2]])),rep("Epithelial cell",length(group[[3]])),rep("Immune cell",length(group[[4]])))

res3$color[ res3$p.value < -log10(0.05)] = "No.significance"

options(repr.plot.height = 6, repr.plot.width = 12)
ggplot(data = res3)+
# geom_point(aes(x = FC,y = p.value))+
ylab("-log10(p.value)")+
theme_classic()+
ylim(0,2.5)+
xlim(-1,1)+
geom_hline(aes(yintercept=-log10(0.05)),linetype = "dashed")+
geom_vline(aes(xintercept=0),linetype = "dashed")+
geom_label(aes(x = FC,y = p.value,label = label,fill = color,size = p.value),color = "white")+
scale_fill_manual(values = c(scales::brewer_pal(palette = "Set1")(4)[1],scales::brewer_pal(palette = "Set1")(4)[4],"grey"))+
guides(size=guide_legend(label = 0))+
theme(axis.title = element_text(size = 18,face = "bold"),
      legend.position = "top",
      axis.text = element_text(size = 12,face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 18,face = "bold"))

res3$VS = "CT VS  TT in Patient with TT"
CTvsTTinWTT = res3

# AT vs TT from patients with TT
sub_data = data[,data$newmeta == "with_thrombus" & data$Sample_Tag != "ccRcc_CT" & !data$Patient  %in% c("without_thrombus_Patient06","without_thrombus_Patient07")]
subtmp_AT = data[,data$newmeta == "with_thrombus" & data$Sample_Tag == "ccRcc_AT" & !data$Patient  %in% c("without_thrombus_Patient06","without_thrombus_Patient07")]
subtmp_TT = data[,data$newmeta == "with_thrombus" & data$Sample_Tag == "ccRcc_TT" & !data$Patient  %in% c("without_thrombus_Patient06","without_thrombus_Patient07")]

tmp_res = vector()

group = list("CAF" = levels(data$All_second_cluster)[20:23],
    "EC" = levels(data$All_second_cluster)[14:19],
    "Epi_cancer"=levels(data$All_second_cluster)[24:25],
    "immune" =c(levels(data$All_second_cluster)[1:13]))

for(i in names(group)){
   
    for(j in c("AT","TT")){
            if(j == "AT"){
                res = reshape2::melt(table(subtmp_AT$Patient,subtmp_AT$All_second_cluster)[,group[[i]]]/rowSums(table(subtmp_AT$Patient,subtmp_AT$All_second_cluster)[,group[[i]]]))
                res$group = "AT"
            }else{
                res = reshape2::melt(table(subtmp_TT$Patient,subtmp_TT$All_second_cluster)[,group[[i]]]/rowSums(table(subtmp_TT$Patient,subtmp_TT$All_second_cluster)[,group[[i]]]))
                res$group = "TT"
            }
            res$Var1 = as.character(res$Var1)
            res$Var2 = as.character(res$Var2)

            tmp_res = rbind(tmp_res,res)
    }

}


res = tmp_res
res[is.na(res$value),"value"] = 0

pr = vector()
for(i in unique(res$Var2)){
    tmp = res[ res$Var2 == i,]
    tmp_AT = tmp[tmp$group == "AT","value"]
    tmp_TT = tmp[tmp$group == "TT","value"]
    pr = c(pr,t.test(tmp_AT,tmp_TT,paired = TRUE,)$p.value)
}

pr[is.na(pr)] = 1

names(pr) = unique(res$Var2)
res$Var2 = factor(res$Var2,levels = names(pr))


col = scales::brewer_pal(palette = "Set1")(9)

library(ggplot2)
p = ggplot(data = res)+geom_boxplot(aes(x = Var2,y = value,color = group))+
geom_point(aes(x = Var2,y = value,color = group),position = position_jitterdodge(0.1))+
    theme_bw()+
    scale_color_manual(values = c(col[1],col[2]))+
    labs(x ="Cell type",y="The proportion of each cell type in each patient",color = " ")+
    theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
          legend.position = "top",
          axis.text.x = element_text(face = "bold",size = 16,angle = 45,hjust = 1),
          axis.text.y = element_text(face = "bold",size = 15),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18),
          axis.title = element_text(size = 18))

 for(i in  1:length(pr)){
     tmp = res[res$Var2 == names(pr)[i],]
     gene_max = max(tmp$value)
     p = p + annotate("text",i-0.5,gene_max,label = paste("p.value:",round(pr[i],4)),size = 5,hjust = 0) 
 }

options(repr.plot.height = 10, repr.plot.width = 30)
p

res2 = vector()
for(i in names(pr)){
    a = mean(res[ res$Var2 == i & res$group =="TT","value"])
    b = mean(res[ res$Var2 == i & res$group =="AT","value"])
    res2 = c(res2 ,(a-b)/(a+b))
}
pr[is.na(pr)] = 1
res3 = as.data.frame(cbind(pr,res2))
colnames(res3) = c("p.value","FC")
res3$FC[is.na(res3$FC)] = 0
res3$p.value = -log10(res3$p.value)
res3$label = rownames(res3)
res3$color = c(rep("Mesenchymal cell",length(group[[1]])),rep("Endothelial cell",length(group[[2]])),rep("Epithelial cell",length(group[[3]])),rep("Immune cell",length(group[[4]])))

res3$color[ res3$p.value < -log10(0.05)] = "No.significance"

options(repr.plot.height = 6, repr.plot.width = 12)
ggplot(data = res3)+
# geom_point(aes(x = FC,y = p.value))+
ylab("-log10(p.value)")+
theme_classic()+
# ylim(0,2.5)+
xlim(-1,1)+
geom_hline(aes(yintercept=-log10(0.05)),linetype = "dashed")+
geom_vline(aes(xintercept=0),linetype = "dashed")+
geom_label(aes(x = FC,y = p.value,label = label,fill = color,size = p.value),color = "white")+
scale_fill_manual(values = c(scales::brewer_pal(palette = "Set1")(4),"grey"))+
guides(size=guide_legend(label = 0))+
theme(axis.title = element_text(size = 18,face = "bold"),
      legend.position = "top",
      axis.text = element_text(size = 12,face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 18,face = "bold"))

res3$VS = "AT VS  TT in Patient with TT"
ATvsTTinWTT = res3

# AT vs PT from patients with TT
sub_data = data[,data$newmeta == "with_thrombus" & data$Sample_Tag != "ccRcc_TT" & !data$Patient  %in% c("without_thrombus_Patient06","without_thrombus_Patient07","with_thrombus_Patient07")]
subtmp_AT = data[,data$newmeta == "with_thrombus" & data$Sample_Tag == "ccRcc_AT" & !data$Patient  %in% c("without_thrombus_Patient06","without_thrombus_Patient07","with_thrombus_Patient07")]
subtmp_CT = data[,data$newmeta == "with_thrombus" & data$Sample_Tag == "ccRcc_CT" & !data$Patient  %in% c("without_thrombus_Patient06","without_thrombus_Patient07","with_thrombus_Patient07")]

tmp_res = vector()

group = list("CAF" = levels(data$All_second_cluster)[20:23],
    "EC" = levels(data$All_second_cluster)[14:19],
    "Epi_cancer"=levels(data$All_second_cluster)[24:25],
    "immune" =c(levels(data$All_second_cluster)[1:13]))

for(i in names(group)){
   
    for(j in c("AT","CT")){
            if(j == "AT"){
                res = reshape2::melt(table(subtmp_AT$Patient,subtmp_AT$All_second_cluster)[,group[[i]]]/rowSums(table(subtmp_AT$Patient,subtmp_AT$All_second_cluster)[,group[[i]]]))
                res$group = "AT"
            }else{
                res = reshape2::melt(table(subtmp_CT$Patient,subtmp_CT$All_second_cluster)[,group[[i]]]/rowSums(table(subtmp_CT$Patient,subtmp_CT$All_second_cluster)[,group[[i]]]))
                res$group = "CT"
            }
            res$Var1 = as.character(res$Var1)
            res$Var2 = as.character(res$Var2)

            tmp_res = rbind(tmp_res,res)
    }

}


res = tmp_res
res[is.na(res$value),"value"] = 0

pr = vector()
for(i in unique(res$Var2)){
    tmp = res[ res$Var2 == i,]
    tmp_AT = tmp[tmp$group == "AT","value"]
    tmp_CT = tmp[tmp$group == "CT","value"]
    pr = c(pr,t.test(tmp_AT,tmp_CT,paired = TRUE,)$p.value)
}

pr[is.na(pr)] = 1

names(pr) = unique(res$Var2)
res$Var2 = factor(res$Var2,levels = names(pr))


col = scales::brewer_pal(palette = "Set1")(9)

library(ggplot2)
p = ggplot(data = res)+geom_boxplot(aes(x = Var2,y = value,color = group))+
geom_point(aes(x = Var2,y = value,color = group),position = position_jitterdodge(0.1))+
    theme_bw()+
    scale_color_manual(values = c(col[1],col[2]))+
    labs(x ="Cell type",y="The proportion of each cell type in each patient",color = " ")+
    theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
          legend.position = "top",
          axis.text.x = element_text(face = "bold",size = 16,angle = 45,hjust = 1),
          axis.text.y = element_text(face = "bold",size = 15),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18),
          axis.title = element_text(size = 18))

 for(i in  1:length(pr)){
     tmp = res[res$Var2 == names(pr)[i],]
     gene_max = max(tmp$value)
     p = p + annotate("text",i-0.5,gene_max,label = paste("p.value:",round(pr[i],4)),size = 5,hjust = 0) 
 }

options(repr.plot.height = 10, repr.plot.width = 30)
p

res2 = vector()
for(i in names(pr)){
    a = mean(res[ res$Var2 == i & res$group =="CT","value"])
    b = mean(res[ res$Var2 == i & res$group =="AT","value"])
    res2 = c(res2 ,(a-b)/(a+b))
}
pr[is.na(pr)] = 1
res3 = as.data.frame(cbind(pr,res2))
colnames(res3) = c("p.value","FC")
res3$FC[is.na(res3$FC)] = 0
res3$p.value = -log10(res3$p.value)
res3$label = rownames(res3)
res3$color = c(rep("Mesenchymal cell",length(group[[1]])),rep("Endothelial cell",length(group[[2]])),rep("Epithelial cell",length(group[[3]])),rep("Immune cell",length(group[[4]])))

res3$color[ res3$p.value < -log10(0.05)] = "No.significance"

options(repr.plot.height = 6, repr.plot.width = 12)
ggplot(data = res3)+
# geom_point(aes(x = FC,y = p.value))+
ylab("-log10(p.value)")+
theme_classic()+
# ylim(0,2.5)+
xlim(-1,1)+
geom_hline(aes(yintercept=-log10(0.05)),linetype = "dashed")+
geom_vline(aes(xintercept=0),linetype = "dashed")+
geom_label(aes(x = FC,y = p.value,label = label,fill = color,size = p.value),color = "white")+
scale_fill_manual(values = c(scales::brewer_pal(palette = "Set1")(4),"grey"))+
guides(size=guide_legend(label = 0))+
theme(axis.title = element_text(size = 18,face = "bold"),
      legend.position = "top",
      axis.text = element_text(size = 12,face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 18,face = "bold"))

res3$VS = "AT VS  CT in Patient with TT"
ATvsCTinWTT = res3



VS = rbind(ATvsCTinWTT,ATvsTTinWTT,CTvsTTinWTT,WvsWOinAT,WvsWOinCT)


VS$p.value[ VS$VS %in% c("AT VS  CT in Patient with TT","AT VS  TT in Patient with TT","CT VS  TT in Patient with TT")] = -1 * VS$p.value[ VS$VS %in% c("AT VS  CT in Patient with TT","AT VS  TT in Patient with TT","CT VS  TT in Patient with TT")]

VS1 = VS[VS$VS %in% c("with_TT VS without TT in CT","with_TT VS without TT in AT"),]

VS1$FC_label = ifelse(VS1$FC > 0,"Enrich in with_TT","Enrich in without_TT")
VS1$FC2 = abs(VS1$FC)
VS1$label = factor(VS1$label,levels = rev(c(
  "Cycling T cell", "CD8+ T cell", "NK cell", "CD4+ T cell", "CD4+ Treg cell",
  "Macrophage", "Monocyte", "Dendritic cell", "Cycling myeloid cell", "Neutrophil",
  "Mast cell", "B cell", "Plasma cell",
  "Glomerular capillaries EC", "Arterioles EC", "Vein EC", "Vasa recta EC",
  "Tumor EC 1", "Tumor EC 2",
  "Pericyte", "vSMC", "Fibroblast", "Mesangial cell", "Epithelial cell", "Cancer cell"
)))

VS1$VS[VS1$VS  == "with_TT VS without TT in CT"] = "with_TT VS without TT in PT"


VS2 = VS[!VS$VS %in% c("with_TT VS without TT in CT","with_TT VS without TT in AT"),]

VS2$FC_label = ifelse(VS2$FC > 0 & VS2$VS == "AT VS  CT in Patient with TT","Enrich in CT",
                      ifelse(VS2$FC < 0 & VS2$VS == "AT VS  CT in Patient with TT","Enrich in AT",
                            ifelse(VS2$FC > 0 & VS2$VS == "AT VS  TT in Patient with TT","Enrich in TT",
                                   ifelse(VS2$FC < 0 & VS2$VS == "AT VS  TT in Patient with TT","Enrich in AT",
                                          ifelse(VS2$FC > 0 & VS2$VS == "CT VS  TT in Patient with TT","Enrich in TT","Enrich in CT")
                                                ))))
VS2$FC2 = abs(VS2$FC)
VS2$label = factor(VS2$label,levels = rev(c(
  "Cycling T cell", "CD8+ T cell", "NK cell", "CD4+ T cell", "CD4+ Treg cell",
  "Macrophage", "Monocyte", "Dendritic cell", "Cycling myeloid cell", "Neutrophil",
  "Mast cell", "B cell", "Plasma cell",
  "Glomerular capillaries EC", "Arterioles EC", "Vein EC", "Vasa recta EC",
  "Tumor EC 1", "Tumor EC 2",
  "Pericyte", "vSMC", "Fibroblast", "Mesangial cell", "Epithelial cell", "Cancer cell"
)))

VS2$VS[ VS2$VS == "AT VS  CT in Patient with TT" ] = "AT VS  PT in Patient with TT"
VS2$VS[ VS2$VS == "CT VS  TT in Patient with TT" ] = "PT VS  TT in Patient with TT"
VS2$FC_label[ VS2$FC_label == "Enrich in CT" ] = "Enrich in PT"


VS1$p.value = 10^(-VS1$p.value)
VS2$p.value = 10^(VS2$p.value)
VS = rbind(VS1,VS2)
VS$fdr_P = p.adjust(VS$p.value,method = "BH")

VS1 = VS[1:50,]
VS2 = VS[51:125,]

options(repr.plot.height = 9, repr.plot.width = 10)
p1 = ggplot(VS1)+geom_point(aes(x = label,y = -log10(fdr_P),color = VS,shape = FC_label,size = FC2),position = position_jitter(width = 0.01))+
ylab("-log10(FDR)")+
ggtitle(label = "with TT VS without TT")+
geom_hline(aes(yintercept=-log10(0.2)),linetype = "dashed")+
scale_color_manual(values = scales::brewer_pal(palette = "Set1")(9)[1:2])+
coord_flip()+
theme_classic()+
theme(axis.text.y = element_text(hjust = 0.5,size = 18,face = "bold"),
      legend.text = element_text(size = 14),
      legend.title = element_text(face = "bold",size = 16),
      plot.title = element_text(hjust = 0.5,size = 22,face = "bold"),
      axis.title.x = element_text(face = "bold",size = 18),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 14,face = "bold"))+
guides(color = guide_legend(order = 1,title = "Comparison",ncol = 1,override.aes = list(size=5)),
         shape = guide_legend(order = 2,title = "Enrich",ncol = 1,override.aes = list(size=5)),
         size = guide_legend(order = 3,title = "Fold Change",ncol = 1))


options(repr.plot.height = 9, repr.plot.width = 10)
p2 = ggplot(VS2)+geom_point(aes(x = label,y = log10(fdr_P),color = VS,shape = FC_label,size = FC2),position = position_jitter(width = 0.01))+
ylab("-log10(FDR)")+
ggtitle(label = "AT VS CT VS TT")+
geom_hline(aes(yintercept=log10(0.2)),linetype = "dashed")+
scale_color_manual(values = scales::brewer_pal(palette = "Set1")(9)[3:5])+
coord_flip()+
theme_classic()+
scale_x_discrete(position = "top")+
scale_y_continuous(labels = c(6,5,4,3,2,1,0),breaks = c(-6,-5,-4,-3,-2,-1,0))+
theme(axis.text.y = element_blank(),
      axis.text.x = element_text(size = 14,face = "bold"),
      legend.box.just = "right",
      legend.title = element_text(hjust = 1,face = "bold",size = 16),
      legend.text = element_text(size = 14),
      legend.position = "left",
      plot.title = element_text(hjust = 0.5,size = 22,face = "bold"),
      axis.title.x = element_text(face = "bold",size = 18),
      axis.title.y = element_blank())+
guides(color = guide_legend(order = 1,title = "Comparison",ncol = 1,override.aes = list(size=5)),
         shape = guide_legend(order = 2,title = "Enrich",ncol = 1,override.aes = list(size=5)),
         size = guide_legend(order = 3,title = "Fold Change",ncol = 1))


library(patchwork)
options(repr.plot.height = 9, repr.plot.width = 20)
p2 + p1

pdf("Paired_comparison_Figure_1D_FDR.pdf",height = 9,width = 20)
p2 + p1
dev.off()

write.table(VS,file = "Paired_comparison_adjusted_P_Figure_1D.txt",sep = "\t",col.names = T,row.names = F)

############## Figure 3C ###################

library(Seurat)
library(harmony)
library(reshape2)
library(RColorBrewer)
library(ggplot2)

data = readRDS("../04. Discovery in Mesenchymal Cell/data/Mesenchymal_clean.rds")

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

data@active.ident = data$All_third_cluster

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



# PT (with TT) vs PT (without TT)
sub_data = data[,data$Sample_Tag == "ccRcc_CT"]
tmp_res = vector()
subtmp = sub_data[,!sub_data$Patient %in% names(table(sub_data$Patient)[table(sub_data$Patient) < 50])]
group = list("CAF" = unique(data$All_fifth_cluster))

for(i in names(group)){
   

    res = reshape2::melt(table(subtmp$Patient,subtmp$All_fifth_cluster)[,group[[i]]]/rowSums(table(subtmp$Patient,subtmp$All_fifth_cluster)[,group[[i]]]))
    res$Var1 = as.character(res$Var1)
    res$Var2 = as.character(res$Var2)

    tmp = reshape2::melt(table(subtmp$Patient,subtmp$newmeta))
    tmp = tmp[tmp$value != 0,]
    tmp$value = NULL

    tmp$Var1 = as.character(tmp$Var1)
    tmp$Var2 = as.character(tmp$Var2)

    rownames(tmp) = tmp$Var1

    res$group = tmp[ res$Var1,"Var2"]
    
    tmp_res = rbind(tmp_res,res)
}

res = tmp_res
res[is.na(res$value),"value"] = 0

pr = vector()
for(i in unique(res$Var2)){
    tmp = res[ res$Var2 == i,]
    pr = c(pr,t.test(value~group,data = tmp)$p.value)
}

pr[is.na(pr)] = 1

names(pr) = unique(res$Var2)
res$Var2 = factor(res$Var2,levels = names(pr))

col = scales::brewer_pal(palette = "Set1")(9)

library(ggplot2)
p = ggplot(data = res)+geom_boxplot(aes(x = Var2,y = value,color = group))+
geom_point(aes(x = Var2,y = value,color = group),position = position_jitterdodge(0.1))+
    theme_bw()+
    scale_color_manual(values = c(col[1],col[2]))+
    labs(x ="Cell type",y="The proportion of each cell type in each patient",color = " ")+
    theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
          legend.position = "top",
          axis.text.x = element_text(face = "bold",size = 16),
          axis.text.y = element_text(face = "bold",size = 15),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18),
          axis.title = element_text(size = 18))

 for(i in  1:length(pr)){
     tmp = res[res$Var2 == names(pr)[i],]
     gene_max = max(tmp$value)
     p = p + annotate("text",i-0.5,gene_max,label = paste("p.value:",round(pr[i],4)),size = 5,hjust = 0) 
 }

options(repr.plot.height = 10, repr.plot.width = 30)
p

res2 = vector()
for(i in names(pr)){
    a = mean(res[ res$Var2 == i & res$group =="with_thrombus","value"])
    b = mean(res[ res$Var2 == i & res$group =="without_thrombus","value"])
    res2 = c(res2,(a-b)/(a+b) )
}
pr[is.na(pr)] = 1
res3 = as.data.frame(cbind(pr,res2))
colnames(res3) = c("p.value","FC")
res3$FC[is.na(res3$FC)] = 0
res3$p.value = -log10(res3$p.value)
res3$label = rownames(res3)
res3$color = res3$label

res3$color[ res3$p.value < -log10(0.05)] = "No.significance"

options(repr.plot.height = 6, repr.plot.width = 12)
ggplot(data = res3)+
# geom_point(aes(x = FC,y = p.value))+
ylab("-log10(p.value)")+
theme_classic()+
# ylim(0,2.5)+
xlim(-1.5,1.5)+
geom_hline(aes(yintercept=-log10(0.05)),linetype = "dashed")+
geom_vline(aes(xintercept=0),linetype = "dashed")+
geom_label(aes(x = FC,y = p.value,label = label,fill = color,size = p.value),color = "white",show.legend = F)+
scale_fill_manual(values = c(col[6],col[4],col[5],"grey"))+
guides(size=guide_legend(label = 0))+
theme(axis.title = element_text(size = 18,face = "bold"),
      legend.position = "top",
      axis.text = element_text(size = 12,face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 18,face = "bold"))

res3$VS = "with_TT VS without TT in CT"
WvsWOinCT = res3

# AT (with TT) vs AT (without TT)
sub_data = data[,data$Sample_Tag == "ccRcc_AT"]
tmp_res = vector()
subtmp = sub_data[,!sub_data$Patient %in% names(table(sub_data$Patient)[table(sub_data$Patient) < 50])]
group = list("CAF" = unique(data$All_fifth_cluster))

for(i in names(group)){
   

    res = reshape2::melt(table(subtmp$Patient,subtmp$All_fifth_cluster)[,group[[i]]]/rowSums(table(subtmp$Patient,subtmp$All_fifth_cluster)[,group[[i]]]))
    res$Var1 = as.character(res$Var1)
    res$Var2 = as.character(res$Var2)

    tmp = reshape2::melt(table(subtmp$Patient,subtmp$newmeta))
    tmp = tmp[tmp$value != 0,]
    tmp$value = NULL

    tmp$Var1 = as.character(tmp$Var1)
    tmp$Var2 = as.character(tmp$Var2)

    rownames(tmp) = tmp$Var1

    res$group = tmp[ res$Var1,"Var2"]
    
    tmp_res = rbind(tmp_res,res)
}

res = tmp_res
res[is.na(res$value),"value"] = 0

pr = vector()
for(i in unique(res$Var2)){
    tmp = res[ res$Var2 == i,]
    pr = c(pr,t.test(value~group,data = tmp)$p.value)
}

pr[is.na(pr)] = 1

names(pr) = unique(res$Var2)
res$Var2 = factor(res$Var2,levels = names(pr))

col = scales::brewer_pal(palette = "Set1")(9)

library(ggplot2)
p = ggplot(data = res)+geom_boxplot(aes(x = Var2,y = value,color = group))+
geom_point(aes(x = Var2,y = value,color = group),position = position_jitterdodge(0.1))+
    theme_bw()+
    scale_color_manual(values = c(col[1],col[2]))+
    labs(x ="Cell type",y="The proportion of each cell type in each patient",color = " ")+
    theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
          legend.position = "top",
          axis.text.x = element_text(face = "bold",size = 16),
          axis.text.y = element_text(face = "bold",size = 15),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18),
          axis.title = element_text(size = 18))

 for(i in  1:length(pr)){
     tmp = res[res$Var2 == names(pr)[i],]
     gene_max = max(tmp$value)
     p = p + annotate("text",i-0.5,gene_max,label = paste("p.value:",round(pr[i],4)),size = 5,hjust = 0) 
 }

options(repr.plot.height = 10, repr.plot.width = 30)
p

res2 = vector()
for(i in names(pr)){
    a = mean(res[ res$Var2 == i & res$group =="with_thrombus","value"])
    b = mean(res[ res$Var2 == i & res$group =="without_thrombus","value"])
    res2 = c(res2,(a-b)/(a+b) )
}
pr[is.na(pr)] = 1
res3 = as.data.frame(cbind(pr,res2))
colnames(res3) = c("p.value","FC")
res3$FC[is.na(res3$FC)] = 0
res3$p.value = -log10(res3$p.value)
res3$label = rownames(res3)
res3$color = res3$label

res3$color[ res3$p.value < -log10(0.05)] = "No.significance"

options(repr.plot.height = 6, repr.plot.width = 12)
ggplot(data = res3)+
# geom_point(aes(x = FC,y = p.value))+
ylab("-log10(p.value)")+
theme_classic()+
# ylim(0,2.5)+
xlim(-1.5,1.5)+
geom_hline(aes(yintercept=-log10(0.05)),linetype = "dashed")+
geom_vline(aes(xintercept=0),linetype = "dashed")+
geom_label(aes(x = FC,y = p.value,label = label,fill = color,size = p.value),color = "white",show.legend = F)+
scale_fill_manual(values = c(col[4],col[5],"grey"))+
guides(size=guide_legend(label = 0))+
theme(axis.title = element_text(size = 18,face = "bold"),
      legend.position = "top",
      axis.text = element_text(size = 12,face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 18,face = "bold"))

res3$VS = "with_TT VS without TT in AT"
WvsWOinAT = res3

# PT vs TT from patients with TT
sub_data = data[,data$newmeta == "with_thrombus" & data$Sample_Tag != "ccRcc_AT" & !data$Patient  %in% c("without_thrombus_Patient06","without_thrombus_Patient07")]
subtmp_CT = data[,data$newmeta == "with_thrombus" & data$Sample_Tag == "ccRcc_CT" & !data$Patient  %in% c("without_thrombus_Patient06","without_thrombus_Patient07")]
subtmp_TT = data[,data$newmeta == "with_thrombus" & data$Sample_Tag == "ccRcc_TT" & !data$Patient  %in% c("without_thrombus_Patient06","without_thrombus_Patient07")]

table_meta = table(sub_data$Patient,sub_data$Sample_Tag)
paired_patients = rownames(table_meta)[rowSums(table_meta > 50) == 2]

tmp_res = vector()
group = list("CAF" = unique(data$All_fifth_cluster))

for(i in names(group)){
   
    for(j in c("CT","TT")){
            if(j == "CT"){
                res = reshape2::melt(table(subtmp_CT$Patient,subtmp_CT$All_fifth_cluster)[,group[[i]]]/rowSums(table(subtmp_CT$Patient,subtmp_CT$All_fifth_cluster)[,group[[i]]]))
                res$group = "CT"
            }else{
                res = reshape2::melt(table(subtmp_TT$Patient,subtmp_TT$All_fifth_cluster)[,group[[i]]]/rowSums(table(subtmp_TT$Patient,subtmp_TT$All_fifth_cluster)[,group[[i]]]))
                res$group = "TT"
            }
            res$Var1 = as.character(res$Var1)
            res$Var2 = as.character(res$Var2)

            tmp_res = rbind(tmp_res,res)
    }

}

res = tmp_res
res[is.na(res$value),"value"] = 0
res = res[ res$Var1 %in% paired_patients,]

pr = vector()
for(i in unique(res$Var2)){
    tmp = res[ res$Var2 == i,]
    tmp_TT = tmp[tmp$group == "TT","value"]
    tmp_CT = tmp[tmp$group == "CT","value"]
    pr = c(pr,t.test(tmp_TT,tmp_CT,paired = TRUE,)$p.value)
}

pr[is.na(pr)] = 1

names(pr) = unique(res$Var2)
res$Var2 = factor(res$Var2,levels = names(pr))

col = scales::brewer_pal(palette = "Set1")(9)

library(ggplot2)
p = ggplot(data = res)+geom_boxplot(aes(x = Var2,y = value,color = group))+
geom_point(aes(x = Var2,y = value,color = group),position = position_jitterdodge(0.1))+
    theme_bw()+
    scale_color_manual(values = c(col[1],col[2]))+
    labs(x ="Cell type",y="The proportion of each cell type in each patient",color = " ")+
    theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
          legend.position = "top",
          axis.text.x = element_text(face = "bold",size = 16),
          axis.text.y = element_text(face = "bold",size = 15),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18),
          axis.title = element_text(size = 18))

 for(i in  1:length(pr)){
     tmp = res[res$Var2 == names(pr)[i],]
     gene_max = max(tmp$value)
     p = p + annotate("text",i-0.5,gene_max,label = paste("p.value:",round(pr[i],4)),size = 5,hjust = 0) 
 }

options(repr.plot.height = 10, repr.plot.width = 30)
p

res2 = vector()
for(i in names(pr)){
    a = mean(res[ res$Var2 == i & res$group =="TT","value"])
    b = mean(res[ res$Var2 == i & res$group =="CT","value"])
    res2 = c(res2 ,(a-b)/(a+b))
}
pr[is.na(pr)] = 1
res3 = as.data.frame(cbind(pr,res2))
colnames(res3) = c("p.value","FC")
res3$FC[is.na(res3$FC)] = 0
res3$p.value = -log10(res3$p.value)
res3$label = rownames(res3)
res3$color = res3$label

res3$color[ res3$p.value < -log10(0.05)] = "No.significance"

options(repr.plot.height = 6, repr.plot.width = 12)
ggplot(data = res3)+
# geom_point(aes(x = FC,y = p.value))+
ylab("-log10(p.value)")+
theme_classic()+
ylim(0,2.5)+
xlim(-1.1,1.1)+
geom_hline(aes(yintercept=-log10(0.05)),linetype = "dashed")+
geom_vline(aes(xintercept=0),linetype = "dashed")+
geom_label(aes(x = FC,y = p.value,label = label,fill = color,size = p.value),color = "white",show.legend = F)+
scale_fill_manual(values = c("grey"))+
guides(size=guide_legend(label = 0))+
# ggtitle(label = "Primary Tumor VS Tumor Thrombus")+
theme(axis.title = element_text(size = 18,face = "bold"),plot.title = element_text(size = 28,hjust = 0.5,face = "bold"),
      legend.position = "top",
      axis.text = element_text(size = 12,face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 18,face = "bold"))

res3$VS = "CT VS  TT in Patient with TT"
CTvsTTinWTT = res3

# AT vs TT from patients with TT
sub_data = data[,data$newmeta == "with_thrombus" & data$Sample_Tag != "ccRcc_CT" & !data$Patient  %in% c("without_thrombus_Patient06","without_thrombus_Patient07")]
subtmp_AT = data[,data$newmeta == "with_thrombus" & data$Sample_Tag == "ccRcc_AT" & !data$Patient  %in% c("without_thrombus_Patient06","without_thrombus_Patient07")]
subtmp_TT = data[,data$newmeta == "with_thrombus" & data$Sample_Tag == "ccRcc_TT" & !data$Patient  %in% c("without_thrombus_Patient06","without_thrombus_Patient07")]

table_meta = table(sub_data$Patient,sub_data$Sample_Tag)
paired_patients = rownames(table_meta)[rowSums(table_meta > 50) == 2]

tmp_res = vector()
group = list("CAF" = unique(data$All_fifth_cluster))

for(i in names(group)){
   
    for(j in c("AT","TT")){
            if(j == "AT"){
                res = reshape2::melt(table(subtmp_AT$Patient,subtmp_AT$All_fifth_cluster)[,group[[i]]]/rowSums(table(subtmp_AT$Patient,subtmp_AT$All_fifth_cluster)[,group[[i]]]))
                res$group = "AT"
            }else{
                res = reshape2::melt(table(subtmp_TT$Patient,subtmp_TT$All_fifth_cluster)[,group[[i]]]/rowSums(table(subtmp_TT$Patient,subtmp_TT$All_fifth_cluster)[,group[[i]]]))
                res$group = "TT"
            }
            res$Var1 = as.character(res$Var1)
            res$Var2 = as.character(res$Var2)

            tmp_res = rbind(tmp_res,res)
    }

}

res = tmp_res
res[is.na(res$value),"value"] = 0
res = res[ res$Var1 %in% paired_patients,]

pr = vector()
for(i in unique(res$Var2)){
    tmp = res[ res$Var2 == i,]
    tmp_TT = tmp[tmp$group == "TT","value"]
    tmp_AT = tmp[tmp$group == "AT","value"]
    pr = c(pr,t.test(tmp_TT,tmp_AT,paired = TRUE,)$p.value)
}

pr[is.na(pr)] = 1

names(pr) = unique(res$Var2)
res$Var2 = factor(res$Var2,levels = names(pr))

col = scales::brewer_pal(palette = "Set1")(9)

library(ggplot2)
p = ggplot(data = res)+geom_boxplot(aes(x = Var2,y = value,color = group))+
geom_point(aes(x = Var2,y = value,color = group),position = position_jitterdodge(0.1))+
    theme_bw()+
    scale_color_manual(values = c(col[1],col[2]))+
    labs(x ="Cell type",y="The proportion of each cell type in each patient",color = " ")+
    theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
          legend.position = "top",
          axis.text.x = element_text(face = "bold",size = 16),
          axis.text.y = element_text(face = "bold",size = 15),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18),
          axis.title = element_text(size = 18))

 for(i in  1:length(pr)){
     tmp = res[res$Var2 == names(pr)[i],]
     gene_max = max(tmp$value)
     p = p + annotate("text",i-0.5,gene_max,label = paste("p.value:",round(pr[i],4)),size = 5,hjust = 0) 
 }

options(repr.plot.height = 10, repr.plot.width = 30)
p

res2 = vector()
for(i in names(pr)){
    a = mean(res[ res$Var2 == i & res$group =="TT","value"])
    b = mean(res[ res$Var2 == i & res$group =="AT","value"])
    res2 = c(res2 ,(a-b)/(a+b))
}
pr[is.na(pr)] = 1
res3 = as.data.frame(cbind(pr,res2))
colnames(res3) = c("p.value","FC")
res3$FC[is.na(res3$FC)] = 0
res3$p.value = -log10(res3$p.value)
res3$label = rownames(res3)
res3$color = res3$label

res3$color[ res3$p.value < -log10(0.05)] = "No.significance"

n = nrow(res3)
col = colorRampPalette(brewer.pal(9,"Set1"))(n)
col = c(col[seq(1,n,by = 2)],col[seq(2,n,by = 2)])

options(repr.plot.height = 6, repr.plot.width = 12)
ggplot(data = res3)+
# geom_point(aes(x = FC,y = p.value))+
ylab("-log10(p.value)")+
theme_classic()+
# ylim(0,2.5)+
xlim(-1.1,1.1)+
geom_hline(aes(yintercept=-log10(0.05)),linetype = "dashed")+
geom_vline(aes(xintercept=0),linetype = "dashed")+
geom_label(aes(x = FC,y = p.value,label = label,fill = color,size = p.value),color = "white",show.legend = F)+
scale_fill_manual(values = c(col[which(sort(res3$label) %in% res3$color)][1:5],"grey",col[which(sort(res3$label) %in% res3$color)][6]))+
# ggtitle(label = "Adjacent Tissue VS Tumor Thrombus")+
theme(axis.title = element_text(size = 18,face = "bold"),plot.title = element_text(size = 28,hjust = 0.5,face = "bold"),
      legend.position = "top",
      axis.text = element_text(size = 12,face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 18,face = "bold"))

res3$VS = "AT VS  TT in Patient with TT"
ATvsTTinWTT = res3


# AT vs PT from patients with TT
sub_data = data[,data$newmeta == "with_thrombus" & data$Sample_Tag != "ccRcc_TT" & !data$Patient  %in% c("without_thrombus_Patient06","without_thrombus_Patient07")]
subtmp_AT = data[,data$newmeta == "with_thrombus" & data$Sample_Tag == "ccRcc_AT" & !data$Patient  %in% c("without_thrombus_Patient06","without_thrombus_Patient07")]
subtmp_CT = data[,data$newmeta == "with_thrombus" & data$Sample_Tag == "ccRcc_CT" & !data$Patient  %in% c("without_thrombus_Patient06","without_thrombus_Patient07")]

table_meta = table(sub_data$Patient,sub_data$Sample_Tag)
paired_patients = rownames(table_meta)[rowSums(table_meta > 50) == 2]
tmp_res = vector()
group = list("CAF" = unique(data$All_fifth_cluster))

for(i in names(group)){
   
    for(j in c("AT","CT")){
            if(j == "AT"){
                res = reshape2::melt(table(subtmp_AT$Patient,subtmp_AT$All_fifth_cluster)[,group[[i]]]/rowSums(table(subtmp_AT$Patient,subtmp_AT$All_fifth_cluster)[,group[[i]]]))
                res$group = "AT"
            }else{
                res = reshape2::melt(table(subtmp_CT$Patient,subtmp_CT$All_fifth_cluster)[,group[[i]]]/rowSums(table(subtmp_CT$Patient,subtmp_CT$All_fifth_cluster)[,group[[i]]]))
                res$group = "CT"
            }
            res$Var1 = as.character(res$Var1)
            res$Var2 = as.character(res$Var2)

            tmp_res = rbind(tmp_res,res)
    }

}

res = tmp_res
res[is.na(res$value),"value"] = 0
res = res[ res$Var1 %in% paired_patients,]

pr = vector()
for(i in unique(res$Var2)){
    tmp = res[ res$Var2 == i,]
    tmp_CT = tmp[tmp$group == "CT","value"]
    tmp_AT = tmp[tmp$group == "AT","value"]
    pr = c(pr,t.test(tmp_CT,tmp_AT,paired = TRUE,)$p.value)
}

pr[is.na(pr)] = 1

names(pr) = unique(res$Var2)
res$Var2 = factor(res$Var2,levels = names(pr))

col = scales::brewer_pal(palette = "Set1")(9)

library(ggplot2)
p = ggplot(data = res)+geom_boxplot(aes(x = Var2,y = value,color = group))+
geom_point(aes(x = Var2,y = value,color = group),position = position_jitterdodge(0.1))+
    theme_bw()+
    scale_color_manual(values = c(col[1],col[2]))+
    labs(x ="Cell type",y="The proportion of each cell type in each patient",color = " ")+
    theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
          legend.position = "top",
          axis.text.x = element_text(face = "bold",size = 16),
          axis.text.y = element_text(face = "bold",size = 15),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18),
          axis.title = element_text(size = 18))

 for(i in  1:length(pr)){
     tmp = res[res$Var2 == names(pr)[i],]
     gene_max = max(tmp$value)
     p = p + annotate("text",i-0.5,gene_max,label = paste("p.value:",round(pr[i],4)),size = 5,hjust = 0) 
 }

options(repr.plot.height = 10, repr.plot.width = 30)
p

res2 = vector()
for(i in names(pr)){
    a = mean(res[ res$Var2 == i & res$group =="CT","value"])
    b = mean(res[ res$Var2 == i & res$group =="AT","value"])
    res2 = c(res2 ,(a-b)/(a+b))
}
pr[is.na(pr)] = 1
res3 = as.data.frame(cbind(pr,res2))
colnames(res3) = c("p.value","FC")
res3$FC[is.na(res3$FC)] = 0
res3$p.value = -log10(res3$p.value)
res3$label = rownames(res3)
res3$color = res3$label

res3$color[ res3$p.value < -log10(0.05)] = "No.significance"

n = nrow(res3)
col = colorRampPalette(brewer.pal(9,"Set1"))(n)
col = c(col[seq(1,n,by = 2)],col[seq(2,n,by = 2)])

options(repr.plot.height = 6, repr.plot.width = 12)
ggplot(data = res3)+
# geom_point(aes(x = FC,y = p.value))+
ylab("-log10(p.value)")+
theme_classic()+
# ylim(0,2.5)+
xlim(-1.1,1.1)+
geom_hline(aes(yintercept=-log10(0.05)),linetype = "dashed")+
geom_vline(aes(xintercept=0),linetype = "dashed")+
geom_label(aes(x = FC,y = p.value,label = label,fill = color,size = p.value),color = "white",show.legend = F)+
scale_fill_manual(values = c(col[which(sort(res3$label) %in% res3$color)][1:5],"grey",col[which(sort(res3$label) %in% res3$color)][6]))+
# ggtitle(label = "Adjacent Tissue VS Tumor Thrombus")+
theme(axis.title = element_text(size = 18,face = "bold"),plot.title = element_text(size = 28,hjust = 0.5,face = "bold"),
      legend.position = "top",
      axis.text = element_text(size = 12,face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 18,face = "bold"))

res3$VS = "AT VS  CT in Patient with TT"
ATvsCTinWTT = res3



VS = rbind(ATvsCTinWTT,ATvsTTinWTT,CTvsTTinWTT,WvsWOinAT,WvsWOinCT)
VS$p.value[ VS$VS %in% c("AT VS  CT in Patient with TT","AT VS  TT in Patient with TT","CT VS  TT in Patient with TT")] = -1 * VS$p.value[ VS$VS %in% c("AT VS  CT in Patient with TT","AT VS  TT in Patient with TT","CT VS  TT in Patient with TT")]
VS1 = VS[VS$VS %in% c("with_TT VS without TT in CT","with_TT VS without TT in AT"),]

VS1$FC_label = ifelse(VS1$FC > 0,"Enrich in with_TT","Enrich in without_TT")
VS1$FC2 = abs(VS1$FC)
VS1$label = factor(VS1$label,levels = rev(levels(data$All_fifth_cluster)))


VS2 = VS[!VS$VS %in% c("with_TT VS without TT in CT","with_TT VS without TT in AT"),]

VS2$FC_label = ifelse(VS2$FC > 0 & VS2$VS == "AT VS  CT in Patient with TT","Enrich in CT",
                      ifelse(VS2$FC < 0 & VS2$VS == "AT VS  CT in Patient with TT","Enrich in AT",
                            ifelse(VS2$FC > 0 & VS2$VS == "AT VS  TT in Patient with TT","Enrich in TT",
                                   ifelse(VS2$FC < 0 & VS2$VS == "AT VS  TT in Patient with TT","Enrich in AT",
                                          ifelse(VS2$FC > 0 & VS2$VS == "CT VS  TT in Patient with TT","Enrich in TT","Enrich in CT")
                                                ))))
VS2$FC2 = abs(VS2$FC)
VS2$label = factor(VS2$label,levels = rev(levels(data$All_fifth_cluster)))

VS1$p.value = 10^(-VS1$p.value)
VS2$p.value = 10^(VS2$p.value)
VS = rbind(VS1,VS2)
VS$fdr_P = p.adjust(VS$p.value,method = "BH")

VS1 = VS[1:24,]
VS2 = VS[25:60,]

options(repr.plot.height = 9, repr.plot.width = 10)
p1 = ggplot(VS1)+geom_point(aes(x = label,y = -log10(fdr_P),color = VS,shape = FC_label,size = FC2),position = position_jitter(width = 0.01))+
ylab("-log10(FDR)")+
ggtitle(label = "with TT VS without TT")+
geom_hline(aes(yintercept=-log10(0.2)),linetype = "dashed")+
scale_color_manual(values = scales::brewer_pal(palette = "Set1")(9)[1:2])+
coord_flip()+
theme_classic()+
theme(axis.text.y = element_text(hjust = 0.5,size = 18,face = "bold"),
      legend.text = element_text(size = 14),
      legend.title = element_text(face = "bold",size = 16),
      plot.title = element_text(hjust = 0.5,size = 22,face = "bold"),
      axis.title.x = element_text(face = "bold",size = 18),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 14,face = "bold"))+
guides(color = guide_legend(order = 1,title = "Comparison",ncol = 1,override.aes = list(size=5)),
         shape = guide_legend(order = 2,title = "Enrich",ncol = 1,override.aes = list(size=5)),
         size = guide_legend(order = 3,title = "Fold Change",ncol = 1))


options(repr.plot.height = 9, repr.plot.width = 10)
p2 = ggplot(VS2)+geom_point(aes(x = label,y = log10(fdr_P),color = VS,shape = FC_label,size = FC2),position = position_jitter(width = 0.01))+
ylab("-log10(FDR)")+
ggtitle(label = "AT VS PT VS TT")+
geom_hline(aes(yintercept=log10(0.2)),linetype = "dashed")+
scale_color_manual(values = scales::brewer_pal(palette = "Set1")(9)[3:5])+
coord_flip()+
theme_classic()+
scale_x_discrete(position = "top")+
scale_y_continuous(labels = c(5,4,3,2,1,0),breaks = c(-5,-4,-3,-2,-1,0))+
theme(axis.text.y = element_blank(),
      axis.text.x = element_text(size = 14,face = "bold"),
      legend.box.just = "right",
      legend.title = element_text(hjust = 1,face = "bold",size = 16),
      legend.text = element_text(size = 14),
      legend.position = "left",
      plot.title = element_text(hjust = 0.5,size = 22,face = "bold"),
      axis.title.x = element_text(face = "bold",size = 18),
      axis.title.y = element_blank())+
      guides(color = guide_legend(order = 1,title = "Comparison",ncol = 1,override.aes = list(size=5)),
         shape = guide_legend(order = 2,title = "Enrich",ncol = 1,override.aes = list(size=5)),
         size = guide_legend(order = 3,title = "Fold Change",ncol = 1))

library(patchwork)
options(repr.plot.height = 5, repr.plot.width = 20)
p2 + p1

pdf("Paired_comparison_Figure_3C_FDR.pdf",height = 5,width = 20)
p2 + p1
dev.off()

write.table(VS,file = "Paired_comparison_adjusted_P_Figure_3C.txt",sep = "\t",col.names = T,row.names = F)