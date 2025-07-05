#!/usr/bin/env python
# coding: utf-8

######### R kernal ###########
library(Seurat)
library(harmony)
library(RColorBrewer)
library(ggplot2)

data = readRDS("./data/Mesenchymal_clean.rds")
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

name = colnames(data)

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
npcs = 60
data2 = RunPCA(object = data2,features = VariableFeatures(data2),npcs = npcs)
data2 = RunHarmony(object = data2,group.by.vars = "Patient", plot_convergence = TRUE,max.iter.harmony = 5,dims.use = 1:npcs,theta = 2,lambda =1,reduction = "pca")
data2 = RunHarmony(object = data2,group.by.vars = "Sample_Tag", plot_convergence = TRUE,reduction = "harmony",dims.use = 1:npcs,theta = 2,lambda = 1,max.iter.harmony = 2,verbose = T)
data2@active.ident = data2$All_third_cluster

data2 = RunUMAP(object = data2,reduction = "harmony",dims = c(1:npcs))

n = length(unique(data2$All_fifth_cluster))
col = colorRampPalette(brewer.pal(9,"Set1"))(n)
col = c(col[seq(1,n,by = 2)],col[seq(2,n,by = 2)])

names(col) = levels(data2$All_fifth_cluster)

options(repr.plot.height = 12, repr.plot.width = 14)
DimPlot(object = data2,group.by = "All_fifth_cluster",pt.size = 1,raster = FALSE,label = F,cols = col,label.box = F)+
NoAxes()+
theme(plot.title = element_blank(),legend.position = "top",legend.text = element_text(size = 18,face = "bold"))+
guides(color=guide_legend(override.aes = list(size = 10,text = 1)))

data2$All_fifth_cluster = as.character(data2$All_fifth_cluster)

options(repr.plot.height = 20, repr.plot.width = 20)
DimPlot(data2,split.by = "Patient",ncol = 5)
data2@assays$RNA@scale.data = as.matrix(data2@assays$RNA@counts)

library(SeuratDisk)
SaveH5Seurat(data2,"./RNA_velocity/Discovery_Mesenchymal_cell.h5Seurat",)
Convert("./RNA_velocity/Discovery_Mesenchymal_cell.h5Seurat",dest = "./RNA_velocity/Discovery_Mesenchymal_cell.h5ad")


Detial_info =  data.frame(
                    p7 = c("with_thrombus_Patient07","RNA_velocity/Patient/P7","P7_filtered_cell_barcode.txt","P7_CAF_velocity_meta.txt","/databak/MJC/ccRCC/20201103/result","P7.BAM"),
                    p8 = c("with_thrombus_Patient08","RNA_velocity/Patient/P8","P8_filtered_cell_barcode.txt","P8_CAF_velocity_meta.txt","/databak/MJC/ccRCC/20201118/patient8/result","P8.BAM"),
                    p9 = c("with_thrombus_Patient09","RNA_velocity/Patient/P9","P9_filtered_cell_barcode.txt","P9_CAF_velocity_meta.txt","/databak/MJC/ccRCC/20201118/patient9/result","P9.BAM"),
                    p10 = c("with_thrombus_Patient10","RNA_velocity/Patient/P10","P10_filtered_cell_barcode.txt","P10_CAF_velocity_meta.txt","/databak/MJC/ccRCC/20201124/patient10/result","P10.BAM"),
                    p11 = c("with_thrombus_Patient11","RNA_velocity/Patient/P11","P11_filtered_cell_barcode.txt","P11_CAF_velocity_meta.txt","/databak/MJC/ccRCC/20210103/patient11/result","P11.BAM"),    
                    p12 = c("with_thrombus_Patient12","RNA_velocity/Patient/P12","P12_filtered_cell_barcode.txt","P12_CAF_velocity_meta.txt","/databak/MJC/ccRCC/20210113/patient12/result_patient12","P12.BAM"),
                    p13 = c("with_thrombus_Patient13","RNA_velocity/Patient/P13","P13_filtered_cell_barcode.txt","P13_CAF_velocity_meta.txt","/databak/MJC/ccRCC/20210225/patient13/result_patient13","P13.BAM"),
                    p14 = c("with_thrombus_Patient14","RNA_velocity/Patient/P14","P14_filtered_cell_barcode.txt","P14_CAF_velocity_meta.txt","/databak/MJC/ccRCC/20210225/patient14/result_patient14","P14.BAM"),
                    p15 = c("with_thrombus_Patient15","RNA_velocity/Patient/P15","P15_filtered_cell_barcode.txt","P15_CAF_velocity_meta.txt","/databak/MJC/ccRCC/20210314/patient15/result_patient15","P15.BAM"),
                    p16 = c("with_thrombus_Patient16","RNA_velocity/Patient/P16","P16_filtered_cell_barcode.txt","P16_CAF_velocity_meta.txt","/databak/MJC/ccRCC/20210422/patient16/result_patient16_exchange","P16.BAM"),
                    p17 = c("with_thrombus_Patient17","RNA_velocity/Patient/P17","P17_filtered_cell_barcode.txt","P17_CAF_velocity_meta.txt","/databak/MJC/ccRCC/20210518/patient17/result_patient17","P17.BAM"),
                    p18 = c("with_thrombus_Patient18","RNA_velocity/Patient/P18","P18_filtered_cell_barcode.txt","P18_CAF_velocity_meta.txt","/databak/MJC/ccRCC/20210707/patient18/patient18","P18.BAM"),

                    p01_02 = c("without_thrombus_Patient01|without_thrombus_Patient02","RNA_velocity/Patient/P01_02","P01_02_filtered_cell_barcode.txt","P01_02_CAF_velocity_meta.txt","/databak/MJC/ccRCC/20201218/p01-02/result","P01_02.BAM"),
                    p03_04 = c("without_thrombus_Patient03|without_thrombus_Patient04","RNA_velocity/Patient/P03_04","P03_04_filtered_cell_barcode.txt","P03_04_CAF_velocity_meta.txt","/databak/MJC/ccRCC/20210112/p03-04/result_new","P03_04.BAM"),
                    p05_06 = c("without_thrombus_Patient05|without_thrombus_Patient06","RNA_velocity/Patient/P05_06","P05_06_filtered_cell_barcode.txt","P05_06_CAF_velocity_meta.txt","/databak/MJC/ccRCC/20210208/p05-06new/result_resequence","P05_06.BAM"),
                    p07_08 = c("without_thrombus_Patient07|without_thrombus_Patient08","RNA_velocity/Patient/P07_08","P07_08_filtered_cell_barcode.txt","P07_08_CAF_velocity_meta.txt","/databak/MJC/ccRCC/20210422/p07-08/result_p0708","P07_08.BAM"),
                    p09_010 = c("without_thrombus_Patient09|without_thrombus_Patient10","RNA_velocity/Patient/P09_010","P09_010_filtered_cell_barcode.txt","P09_010_CAF_velocity_meta.txt","/databak/MJC/ccRCC/20210422/p09-10/result_p0910","P09_010.BAM")

                )

meta = cbind(data2@meta.data,data2@reductions$umap@cell.embeddings)
meta$tmp = NULL

write.table(meta,"./RNA_velocity/Discovery_meta.txt",sep = "\t",col.names = T,row.names = T,quote = F)

Run_RNA_velocity = function(data,Detial_info_patient){
    
    #本地创建文件夹
    if( !dir.exists(Detial_info_patient[2]) ){
    dir.create(Detial_info_patient[2])
    }
    
    start_time <- proc.time()
    if(grepl("\\|",Detial_info_patient[1])){
        
        #barcode
        patient_names = strsplit(Detial_info_patient[1],"\\|")[[1]]
        Patient_barcode = do.call(rbind,strsplit(colnames(data)[ data$Patient %in% patient_names],"_"))[,5]
        message(print(patient_names))
        message(print(length(Patient_barcode)))
        write.table(Patient_barcode,paste(Detial_info_patient[2],"/",Detial_info_patient[3],sep=""),col.names = F,row.names = F,quote = F)
        message(paste("SAVE",paste(Detial_info_patient[2],"/",Detial_info_patient[3],sep="")))
        
        #meta data
        tmp_meta = data[,data$Patient %in% patient_names]@meta.data
        write.table(cbind(tmp_meta,data[,data$Patient %in% patient_names]@reductions$umap@cell.embeddings),paste(Detial_info_patient[2],"/",Detial_info_patient[4],sep=""),col.names = T,row.names = F)
        message(paste("SAVE",paste(Detial_info_patient[2],"/",Detial_info_patient[4],sep="")))
    
    }else{
        patient_names = Detial_info_patient[1]
        Patient_barcode = do.call(rbind,strsplit(colnames(data)[ data$Patient %in% patient_names ],"_"))[,5]
        message(print(patient_names))
        message(print(length(Patient_barcode)))
        write.table(Patient_barcode,paste(Detial_info_patient[2],"/",Detial_info_patient[3],sep=""),col.names = F,row.names = F,quote = F)
        message(paste("SAVE",paste(Detial_info_patient[2],"/",Detial_info_patient[3],sep="")))
        
        tmp_meta = data[,data$Patient %in% patient_names]@meta.data
        write.table(cbind(tmp_meta,data[,data$Patient %in% patient_names]@reductions$umap@cell.embeddings),paste(Detial_info_patient[2],"/",Detial_info_patient[4],sep=""),col.names = T,row.names = F)
        message(paste("SAVE",paste(Detial_info_patient[2],"/",Detial_info_patient[4],sep="")))
    }
    end_time <- proc.time()
    Execution <- end_time - start_time
    print(paste("Execution time: ",Execution["elapsed"],"seconds"))
    
    
    message(paste("Now sending Pre_BAM.py to",Detial_info_patient[5]))
    code1 = paste("sshpass -p ****** scp ./RNA_velocity/Pre_BAM.py liyang@192.168.99.141:",Detial_info_patient[5],sep = "")
    system(code1)
    message(paste("Now Running Pre_BAM.py in",Detial_info_patient[5],"and save tmp*bam in PYSAM"))
    code2 = paste("sshpass -p ****** ssh liyang@192.168.99.141 'cd ",Detial_info_patient[5],";ls *BAM | xargs /home/liyang/miniconda3/bin/python Pre_BAM.py'",sep="")
    system(code2)
    end_time <- proc.time()
    Execution <- end_time - start_time
    print(paste("Execution time: ",Execution["elapsed"]/60,"minunes"))
    
    
    message(paste("Now merging tmp*bam to",Detial_info_patient[6]))
    code3 = paste("sshpass -p ****** ssh liyang@192.168.99.141 'cd ",Detial_info_patient[5],"/PYSAM;samtools merge -@ 32 ",Detial_info_patient[6]," tmp_*BAM'",sep="")
    system(code3)
    message(paste("Now removing tmp*bam"))
    code4 = paste("sshpass -p ****** ssh liyang@192.168.99.141 'cd ",Detial_info_patient[5],"/PYSAM;rm tmp*BAM'",sep="")
    system(code4)
    end_time <- proc.time()
    Execution <- end_time - start_time
    print(paste("Execution time: ",Execution["elapsed"]/60,"minunes"))
    
    
    message(paste("Now sending ",Detial_info_patient[2],"/",Detial_info_patient[3]," to ",Detial_info_patient[5],"/PYSAM",sep = ""))
    code5 = paste("sshpass -p ****** scp ",Detial_info_patient[2],"/",Detial_info_patient[3]," liyang@192.168.99.141:",Detial_info_patient[5],"/PYSAM;",sep="")
    system(code5)
    message(paste("Now Running velocyto"))
    code6 = paste("sshpass -p ****** ssh liyang@192.168.99.141 'cd ",Detial_info_patient[5],"/PYSAM;/home/liyang/miniconda3/bin/velocyto  run -b ",Detial_info_patient[3]," -o./ -m /data2/liyang/Reference/GENOME/GRCh38/hg38_repeat_rmsk.gtf ",Detial_info_patient[6]," /data2/liyang/Reference/GENOME/GRCh38/Homo_sapiens.GRCh38.101.gtf'",sep="")
    system(code6)
    message(paste("Now sending loom file to",Detial_info_patient[2]))
    code7 = paste("sshpass -p ****** scp liyang@192.168.99.141:",Detial_info_patient[5],"/PYSAM/*.loom ",Detial_info_patient[2],sep="")
    system(code7)
    
    end_time <- proc.time()
    Execution <- end_time - start_time
    print(paste("Execution time: ",Execution["elapsed"]/60,"minunes"))
}


Run_RNA_velocity(data,Detial_info$p7)
Run_RNA_velocity(data,Detial_info$p8)
Run_RNA_velocity(data,Detial_info$p9)
Run_RNA_velocity(data,Detial_info$p10)
Run_RNA_velocity(data,Detial_info$p11)
Run_RNA_velocity(data,Detial_info$p13)
Run_RNA_velocity(data,Detial_info$p14)
Run_RNA_velocity(data,Detial_info$p15)
Run_RNA_velocity(data,Detial_info$p16)
Run_RNA_velocity(data,Detial_info$p17)
Run_RNA_velocity(data,Detial_info$p18)

Run_RNA_velocity(data,Detial_info$p01_02)
Run_RNA_velocity(data,Detial_info$p03_04)
Run_RNA_velocity(data,Detial_info$p05_06)
Run_RNA_velocity(data,Detial_info$p07_08)
Run_RNA_velocity(data,Detial_info$p09_010)


# # velocyto 结果合并

######### python kernal ###########

import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import loompy
import velocyto as vcy
import logging
from sklearn.svm import SVR
from sklearn.linear_model import LinearRegression
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.interpolate import interp1d

logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)
get_ipython().run_line_magic('matplotlib', 'inline')
plt.rcParams['pdf.fonttype'] = 42


# plotting utility functions
def despline():
    ax1 = plt.gca()
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')
    
def minimal_xticks(start, end):
    end_ = np.around(end, -int(np.log10(end))+1)
    xlims = np.linspace(start, end_, 5)
    xlims_tx = [""]*len(xlims)
    xlims_tx[0], xlims_tx[-1] = f"{xlims[0]:.0f}", f"{xlims[-1]:.02f}"
    plt.xticks(xlims, xlims_tx)

    
def minimal_yticks(start, end):
    end_ = np.around(end, -int(np.log10(end))+1)
    ylims = np.linspace(start, end_, 5)
    ylims_tx = [""]*len(ylims)
    ylims_tx[0], ylims_tx[-1] = f"{ylims[0]:.0f}", f"{ylims[-1]:.02f}"
    plt.yticks(ylims, ylims_tx)


vlm = loompy.combine(["./RNA_velocity/Patient/P01_02/P01_02_OATW5.loom",
                      "./RNA_velocity/Patient/P03_04/P03_04_QDAD9.loom",
                      "RNA_velocity/Patient/P05_06/P05_06_KFE5J.loom",
                      "RNA_velocity/Patient/P07_08/P07_08_HH6E5.loom",
                      "RNA_velocity/Patient/P09_010/P09_010_L0H5Q.loom",
                      "./RNA_velocity/Patient/P7/P7_B6DMP.loom",
                      "./RNA_velocity/Patient/P8/P8_TRPWU.loom",
                      "./RNA_velocity/Patient/P9/P9_M3ZLF.loom",
                      "./RNA_velocity/Patient/P10/P10_N5O7Z.loom",
                      "./RNA_velocity/Patient/P11/P11_QZC3R.loom",
                      "./RNA_velocity/Patient/P12/P12_2ARH6.loom",
                      "./RNA_velocity/Patient/P13/P13_0E42N.loom",
                      "./RNA_velocity/Patient/P14/P14_8VPHS.loom",
                      "./RNA_velocity/Patient/P15/P15_JO5W1.loom",
                      "./RNA_velocity/Patient/P16/P16_UYZ2S.loom",
                      "./RNA_velocity/Patient/P17/P17_NVSWM.loom",
                      "./RNA_velocity/Patient/P18/P18_A6XTL.loom"
                     ],
                     "./RNA_velocity/Patient/Discovery_Mesenchymal.loom"
                    )

import scvelo as scv
import pandas as pd



IDMAP =     {"with_thrombus_Patient07":"P7",
             "with_thrombus_Patient08":"P8",
             "with_thrombus_Patient09":"P9",
             "with_thrombus_Patient10":"P10",
             "with_thrombus_Patient11":"P11",
             "with_thrombus_Patient12":"P12",
             "with_thrombus_Patient13":"P13",
             "with_thrombus_Patient14":"P14",
             "with_thrombus_Patient15":"P15",
             "with_thrombus_Patient16":"P16",
             "with_thrombus_Patient17":"P17",
             "with_thrombus_Patient18":"P18",

             "without_thrombus_Patient01":"P01",
             "without_thrombus_Patient02":"P01",
             "without_thrombus_Patient03":"P03",
             "without_thrombus_Patient04":"P03",
             "without_thrombus_Patient05":"P05",
             "without_thrombus_Patient06":"P05",
             "without_thrombus_Patient07":"P07",
             "without_thrombus_Patient08":"P07",
             "without_thrombus_Patient09":"P09",
             "without_thrombus_Patient10":"P09"
            }

scv.set_figure_params(dpi = 150)

data = scv.read("./RNA_velocity/Discovery_Mesenchymal_cell.h5ad", cache=True,)

data.obs["All_fifth_cluster"] = pd.Categorical(data.obs["All_fifth_cluster"],categories=["Mesangial cell","vSMC 1","vSMC 2","Cycling Pericyte","IFI+ Pericyte","Pericyte 1","Pericyte 2","CYSLTR2+ Fibroblast","FAP+ Fibroblast","FBLN5+ Fibroblast","MFAP5+ Fibroblast","ADH1B+ Fibroblast"])

colors = ['#E41A1C','#419486','#91569A','#FFAD12','#B6742A','#DD87B4','#66628D','#5A9D5A','#D96D3B','#F6EF32','#D26D7A','#999999']
scv.pl.umap(data,color="All_fifth_cluster",legend_loc = "right margin",palette=colors)
# scv.pl.tsne(data,color="subset",legend_loc = "right margin")
scv.pl.pca(data,color="All_fifth_cluster",legend_loc = "right margin")
data.obs_names = data.obs.apply(lambda x: IDMAP[x["Patient"]]+"_"+x["Cell_Index"].split("_ ")[1],axis = 1)


data.varm['PCs'] = data.varm['HARMONY'].copy()
data.obsm['X_pca'] = data.obsm['X_harmony'].copy()

ldata = scv.read("./RNA_velocity/Patient/Discovery_Mesenchymal.loom", cache=True)

ldata.obs_names = list(map(lambda x:x.split("_")[0]+"_"+x.split(":")[1][:-1],ldata.obs_names))

adata = scv.utils.merge(data, ldata)

scv.pl.proportions(adata,groupby="All_third_cluster")
scv.pl.proportions(adata,groupby="Patient")

scv.pp.filter_genes(adata, min_shared_counts=20,min_cells = 5)

scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
scv.pp.log1p(adata)
scv.pp.neighbors(adata,n_pcs=60,n_neighbors=30)

scv.tl.recover_dynamics(adata,n_jobs=24)
scv.tl.velocity(adata, mode='dynamical')

scv.tl.velocity_graph(adata,approx = None,n_jobs = 24)

import matplotlib.pyplot as plt

scv.pl.velocity_embedding_stream(adata,dpi=150,
                                 title="RNA Velocity",
                                 palette=colors,
                                 basis='umap',
                                 color='All_fifth_cluster',
                                 legend_loc='best',
                                 density=3,
                                 size=14,
                                 alpha=1,
                                 arrow_color="black",
                                 arrow_size=1
                                )

