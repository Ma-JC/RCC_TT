#!/usr/bin/env python
# coding: utf-8

########### Cell niche analysis ######################
import pandas as pd
import seaborn
import pickle
import scanpy as sc
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import os


os.environ['PYTHONHASHSEED'] = '100'
np.random.seed(1)
random.seed(1)


Sfile = open("/home/ncpsb/new_mnt/ccRCC_data/Spatial_transcriptome/All_ST_stlearn_stSMEclust.pkl","rb")
All_Data_stSME_clus = pickle.load(Sfile)
Sfile.close()


import pickle
Sfile = open("/home/ncpsb/new_mnt/ccRCC_data/Spatial_transcriptome/data/Spatial_Data_List_cell2location.1.pkl","rb")
Data_interface = pickle.load(Sfile)
Sfile.close()

Sfile = open("/home/ncpsb/new_mnt/ccRCC_data/Spatial_transcriptome/data/Spatial_Core_cell2location.1.pkl","rb")
Data_core = pickle.load(Sfile)
Sfile.close()

Sfile = open("/home/ncpsb/new_mnt/ccRCC_data/Spatial_transcriptome/data/Spatial_immunity_cell2location.1.pkl","rb")
Data_immunity = pickle.load(Sfile)
Sfile.close()

Sfile = open("/home/ncpsb/new_mnt/ccRCC_data/Spatial_transcriptome/data/Spatial_four_TT_cell2location.1.pkl","rb")
Data_TT = pickle.load(Sfile)
Sfile.close()

Sfile = open("/home/ncpsb/new_mnt/ccRCC_data/Spatial_transcriptome/data/Spatial_TTpatient_four_PT_cell2location.1.pkl","rb")
Data_PT = pickle.load(Sfile)
Sfile.close()



All_Data = Data_interface+Data_core+Data_immunity+Data_TT+Data_PT


All_Data2 = []
for i,data in enumerate(All_Data_stSME_clus):
    data1 = data.copy()
    data2 = All_Data[i].copy()
    if data1.obs["sample"].cat.categories == data2.obs["sample"].cat.categories:
        print(data1.obs["sample"].cat.categories,data2.obs["sample"].cat.categories,i)
        data2 = data2[data1.obs_names,:]
        data2.obs.loc[:,"louvain"] = data1.obs["louvain"]
        All_Data2.append(data2)

for i,tmp_data in enumerate(All_Data2):
    print(i)
    tmp_data = tmp_data.copy()
    cell_names = list(map(lambda x:x.split("meanscell_abundance_w_sf_")[1],tmp_data.obsm["means_cell_abundance_w_sf"].columns))
    for j,category in enumerate(tmp_data.obs["louvain"].cat.categories):

        tmp_mat = tmp_data.obs.loc[tmp_data.obs["louvain"] == category,cell_names].sum(0)
        tmp_mat = tmp_mat/tmp_mat.sum()
        symbol = tmp_mat.index
        spatial_domain_name = tmp_data.obs["sample"].cat.categories + "_domain_" + str(category)
        if i == 0 and j == 0:
            ref = pd.DataFrame(tmp_mat.transpose().tolist(),index = symbol,columns=spatial_domain_name)
        else:
            ref = pd.merge(ref,pd.DataFrame(tmp_mat.transpose().tolist(),index = symbol,columns=spatial_domain_name),
                           left_index=True,right_index=True,how = "outer")

pseudo_bulk_domain = sc.AnnData(ref.transpose())

pseudo_bulk_domain.obs["sample"] = [i.split("_domain_")[0] for i in pseudo_bulk_domain.obs_names]
patient_meta = pd.read_csv("/home/ncpsb/new_mnt/ccRCC_data/Spatial_transcriptome/data/ST_meta.csv",index_col=0)
patient_meta = patient_meta.loc[[i.split("_domain_")[0] for i in pseudo_bulk_domain.obs_names],["T_stage_2","T_stage_detial","batch","TT1","TT2","sample_type"]]

for i in patient_meta.columns:
    pseudo_bulk_domain.obs.loc[:,i] = patient_meta.loc[:,i].tolist()
    
pseudo_bulk_domain = pseudo_bulk_domain[ pseudo_bulk_domain.obs["batch"].isin(["immunity","interface","core","PT","TT"]),:]

pseudo_bulk_domain.raw = pseudo_bulk_domain.copy()

meta_data = pseudo_bulk_domain.obs
data_mat = pseudo_bulk_domain.X

import harmonypy as hm
ho = hm.run_harmony(data_mat, meta_data, "batch",max_iter_harmony=5,theta=2)
ho2 = hm.run_harmony(ho.Z_corr.T, meta_data, "sample_type",max_iter_harmony=5,theta=2)
ho3 = hm.run_harmony(ho2.Z_corr.T, meta_data, "sample",max_iter_harmony=5,theta=2)
pseudo_bulk_domain.X = ho3.Z_corr.T

sc.pp.neighbors(pseudo_bulk_domain, n_pcs =0,n_neighbors=10)
sc.tl.umap(pseudo_bulk_domain)



sc.set_figure_params(fontsize = 12)
# Run clustering with leiden
sc.tl.leiden(pseudo_bulk_domain, resolution=0.1)
pseudo_bulk_domain.obs["cluster1"] = pseudo_bulk_domain.obs["leiden"]
sc.tl.leiden(pseudo_bulk_domain, resolution=0.5)
pseudo_bulk_domain.obs["cluster2"] = pseudo_bulk_domain.obs["leiden"]
sc.tl.leiden(pseudo_bulk_domain, resolution=3)
pseudo_bulk_domain.obs["cluster3"] = pseudo_bulk_domain.obs["leiden"]

# Plotting UMAP
sc.pl.umap(pseudo_bulk_domain, color=["cluster1"],size=10,legend_loc="on data",legend_fontsize=8)
sc.pl.umap(pseudo_bulk_domain, color=["cluster2"],size=10,legend_loc="on data",legend_fontsize=8)
sc.pl.umap(pseudo_bulk_domain, color=["cluster3"],size=10,legend_loc="on data",legend_fontsize=8)
sc.pl.umap(pseudo_bulk_domain, color=["sample_type"],size=10)
sc.pl.umap(pseudo_bulk_domain, color=["batch"],size=10)
sc.pl.umap(pseudo_bulk_domain, color=["sample"],size=10)
sc.pl.umap(pseudo_bulk_domain, color=["T_stage_2"],size=10)
sc.pl.umap(pseudo_bulk_domain, color=["TT1"],size=10)

map1 = dict()
for i in pseudo_bulk_domain.obs["cluster1"].index:
    if i.split("_domain_")[0] not in map1.keys():
        map1[i.split("_domain_")[0]] = dict()
        map1[i.split("_domain_")[0]][i.split("_domain_")[1]] = pseudo_bulk_domain.obs["cluster1"][i]
    else:
        map1[i.split("_domain_")[0]][i.split("_domain_")[1]] = pseudo_bulk_domain.obs["cluster1"][i]
        
map2 = dict()
for i in pseudo_bulk_domain.obs["cluster2"].index:
    if i.split("_domain_")[0] not in map2.keys():
        map2[i.split("_domain_")[0]] = dict()
        map2[i.split("_domain_")[0]][i.split("_domain_")[1]] = pseudo_bulk_domain.obs["cluster2"][i]
    else:
        map2[i.split("_domain_")[0]][i.split("_domain_")[1]] = pseudo_bulk_domain.obs["cluster2"][i]
        
map3 = dict()
for i in pseudo_bulk_domain.obs["cluster3"].index:
    if i.split("_domain_")[0] not in map3.keys():
        map3[i.split("_domain_")[0]] = dict()
        map3[i.split("_domain_")[0]][i.split("_domain_")[1]] = pseudo_bulk_domain.obs["cluster3"][i]
    else:
        map3[i.split("_domain_")[0]][i.split("_domain_")[1]] = pseudo_bulk_domain.obs["cluster3"][i]


sc.set_figure_params(fontsize = 12)

n = 0
for i,data in enumerate(All_Data2):
    data = data.copy()
    sample_name = data.obs["sample"].cat.categories
    if list(sample_name)[0] in list(map2.keys()):
        map_tmp1 = map1[sample_name.to_list()[0]]
        map_tmp2 = map2[sample_name.to_list()[0]]
        map_tmp3 = map3[sample_name.to_list()[0]]
#         map_tmp4 = map4[sample_name.to_list()[0]]
        data.obs["cluster1"] = [map_tmp1[i] for i in data.obs["louvain"]]
        data.obs["cluster2"] = [map_tmp2[i] for i in data.obs["louvain"]]
        data.obs["cluster3"] = [map_tmp3[i] for i in data.obs["louvain"]]
#         data.obs["cluster4"] = [map_tmp4[i] for i in data.obs["louvain"]]
        print(sample_name)
#         sc.pl.spatial(data,color = ["cluster1","cluster2","cluster3"])
        mat_c = data.obs[list(map(lambda x:x.split("meanscell_abundance_w_sf_")[1],data.obsm["means_cell_abundance_w_sf"].columns))+["cluster3"]].groupby("cluster3").sum()
        mat_c = mat_c.transpose()
        mat_c = mat_c/mat_c.sum(0)
        mat_c.columns = [ str(i) + "|" + sample_name.tolist()[0] for i in mat_c.columns.tolist()]
        
        if n == 0:
            mat_c_all = mat_c
            n += 1
        else:
            mat_c_all = pd.merge(mat_c_all,mat_c,left_index=True,right_index=True)


n = 0
for i in pseudo_bulk_domain.obs["cluster3"].cat.categories:
    m = [x.split("|")[0] == i for x in mat_c_all.columns]
    mat_inte = mat_c_all.loc[:,m].sum(1)/mat_c_all.loc[:,m].sum(1).sum()
    symbol = mat_inte.index
    spatial_domain_name = ["Domain_"+i]
    if n == 0:
        ref2 = pd.DataFrame(mat_inte.tolist(),index = symbol,columns=spatial_domain_name)
        n += 1
    else:
        ref2 = pd.merge(ref2,pd.DataFrame(mat_inte.tolist(),index = symbol,columns=spatial_domain_name),
                           left_index=True,right_index=True,how = "outer")

color_dict = {'Fibroblast': 'lime',
             'PT-like Cancer cell': 'black',
             'EMT-like Cancer cell': 'blue',
             'Stress-like Cancer cell': 'springgreen',
             'Immunoregulatory Cancer cell': 'yellow',
             'Epithelial cell': 'violet',
             'B cell': 'cyan',
             'CD4+ T cell': 'orange',
             'Macrophage': 'blueviolet',
             'Pericyte': 'red',
             'Endothelial cell 5(Cancer)': 'darkgreen',
             'CD8+ T cell': 'magenta',
             'CD4+ Treg cell': 'tan',
             'Cycling Myeloid cell': 'dimgrey',
             'Cycling Cancer cell': 'darkgrey',
             'Cycling T cell': 'lightgrey',
             'Dendritic cell': 'greenyellow',
             'Endothelial cell 1(Glomerular Capillaries)': 'sienna',
             'Endothelial cell 2(Arterioles)': 'lightcoral',
             'Endothelial cell 4(Vasa Recta)': 'rosybrown',
             'Endothelial cell 6(Cancer)': 'bisque',
             'Mast cell': 'lightsteelblue',
             'Mesangial cell': 'gold',
             'Monocyte': 'pink',
             'Neutrophil': 'darkkhaki',
             'Plasma cell': 'olivedrab',
             'vSMC': 'purple'}


pseudo_bulk_domain.X = data_mat


min_x = np.min(pseudo_bulk_domain.obsm["X_umap"][:,0])
max_x = np.max(pseudo_bulk_domain.obsm["X_umap"][:,0])
min_y = np.min(pseudo_bulk_domain.obsm["X_umap"][:,1])
max_y = np.max(pseudo_bulk_domain.obsm["X_umap"][:,1])


range_x = max_x - min_x
range_y = max_y - min_y

mapped_x = 0 + (pseudo_bulk_domain.obsm["X_umap"][:,0] - min_x) * (12 - 0) / range_x
mapped_y = 0 + (pseudo_bulk_domain.obsm["X_umap"][:,1] - min_y) * (12 - 0) / range_y

import matplotlib.cm as cm
from matplotlib.lines import Line2D
top_num = 3
label_set = []

pseudo_bulk_domain.X = data_mat

for i in range(pseudo_bulk_domain.shape[0]):
#     m = np.argsort(pseudo_bulk_domain.X[i,:])[-top_num:]
    m = pseudo_bulk_domain.X[i,:] >= 0.08
    frac = pseudo_bulk_domain.X[i,:][m]
    labels = pseudo_bulk_domain.var_names[m]
    label_set += labels.tolist()
    
label_set = set(label_set)    

fig, ax = plt.subplots(dpi=50)
for i in range(pseudo_bulk_domain.shape[0]):
#     m = np.argsort(pseudo_bulk_domain.X[i,:])[-top_num:]
    m = pseudo_bulk_domain.X[i,:] >= 0.08
    frac = pseudo_bulk_domain.X[i,:][m]
    labels = pseudo_bulk_domain.var_names[m]
    colors = [color_dict[i] for i in labels]
    x = mapped_x[i]
    y = mapped_y[i]

    ax.pie(frac, startangle=90, counterclock=False,colors= colors,
           radius=0.1, center=(x,y))

      

ax.spines['left'].set_position('zero')
ax.spines['left'].set_linewidth(0.5)
ax.spines['left'].set_color('gray')

ax.spines['bottom'].set_position('zero')
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['bottom'].set_color('gray')

ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')

ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color_dict[label], markersize=28, label=label) for label in sorted(list(label_set))]
plt.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(2.2, 2.2),ncol =1,borderpad=0.3,fontsize=28)
plt.show()


pseudo_bulk_domain.X = ho3.Z_corr.T



import pickle
Sfile = open("/home/ncpsb/new_mnt/ccRCC_data/Spatial_transcriptome/data/pseudo_bulk_domain.CancerCell.pkl","rb")
pseudo_bulk_domain,data_mat,ho3 = pickle.load(Sfile)
Sfile.close()


tmp_data = pseudo_bulk_domain[pseudo_bulk_domain.obs["cluster3"].isin(["3","5","2","11","9","14","28","4","25"]),:]
sc.pp.neighbors(tmp_data, n_pcs =0,n_neighbors=10)
sc.tl.umap(tmp_data)


sc.tl.leiden(tmp_data, resolution=3)
tmp_data.obs["tmp_cluster"] = tmp_data.obs["leiden"]

sc.set_figure_params(fontsize = 12)
sc.pl.umap(tmp_data, color=["cluster3"],size=28,legend_loc="on data",legend_fontsize=8)


fig=sc.pl.umap(tmp_data, color=["tmp_cluster"],size=100,legend_loc="on data",legend_fontsize=8,show=False,return_fig=True)
ax=fig.axes[0]
ax.set_title("")
ax.set_xlabel('UMAP1', fontsize=12,fontproperties=FontProperties(weight='bold'))
ax.set_ylabel('UMAP2', fontsize=12,fontproperties=FontProperties(weight='bold'))
plt.show()


tmp_data.X = data_mat[pseudo_bulk_domain.obs["cluster3"].isin(["3","5","2","11","9","14","28","4","25"]),:]

tmp_min_x = np.min(tmp_data.obsm["X_umap"][:,0])
tmp_max_x = np.max(tmp_data.obsm["X_umap"][:,0])
tmp_min_y = np.min(tmp_data.obsm["X_umap"][:,1])
tmp_max_y = np.max(tmp_data.obsm["X_umap"][:,1])


tmp_range_x = tmp_max_x - tmp_min_x
tmp_range_y = tmp_max_y - tmp_min_y

tmp_mapped_x = 0 + (tmp_data.obsm["X_umap"][:,0] - tmp_min_x) * (10 - 0) / tmp_range_x
tmp_mapped_y = 0 + (tmp_data.obsm["X_umap"][:,1] - tmp_min_y) * (10 - 0) / tmp_range_y

import matplotlib.cm as cm
from matplotlib.lines import Line2D
top_num = 3
label_set = []

for i in range(tmp_data.shape[0]):
#     m = np.argsort(tmp_data.X[i,:])[-top_num:]
    m = tmp_data.X[i,:] >= 0.08
    frac = tmp_data.X[i,:][m]
    labels = tmp_data.var_names[m]
    label_set += labels.tolist()
    
label_set = set(label_set)    
   color_dict[j] = colors(i)
    

fig, ax = plt.subplots(dpi=50)
for i in range(tmp_data.shape[0]):
#     m = np.argsort(tmp_data.X[i,:])[-top_num:]
    m = tmp_data.X[i,:] >= 0.08
    frac = tmp_data.X[i,:][m]
    labels = tmp_data.var_names[m]
    colors = [color_dict[i] for i in labels]
    x = tmp_mapped_x[i]
    y = tmp_mapped_y[i]

    ax.pie(frac, startangle=90, counterclock=False,colors= colors,
           radius=0.1, center=(x,y))

      

ax.spines['left'].set_position('zero')
ax.spines['left'].set_linewidth(0.5)
ax.spines['left'].set_color('gray')

ax.spines['bottom'].set_position('zero')
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['bottom'].set_color('gray')

ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')

ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color_dict[label], markersize=22, label=label) for label in sorted(list(label_set))]
plt.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(3, 1.7),fontsize=22)
plt.show()


tmp_data.obs["New_Domain"] = None
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["4"])] = "Cell_Niche_1"
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["13","8"])] = "Cell_Niche_2"
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["15"])] = "Cell_Niche_3"
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["3","18","1","17"])] = "Cell_Niche_4"

# tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["18","15","17"])] = "Domain_4"
# tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["1"])] = "Domain_5"
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["5"])] = "Cell_Niche_5"
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["7","19"])] = "Cell_Niche_6"
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["2","16","9"])] = "Cell_Niche_7"
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["6","12"])] = "Cell_Niche_8"
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["14","11"])] = "Cell_Niche_9"
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["0"])] = "Cell_Niche_10"
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["10"])] = "Cell_Niche_11"

fig=sc.pl.umap(tmp_data, color=["New_Domain"],size=100,legend_loc="on data",legend_fontsize=12,show=False,return_fig=True)
ax=fig.axes[0]
ax.set_title("")
ax.set_xlabel('UMAP1', fontsize=12,fontproperties=FontProperties(weight='bold'))
ax.set_ylabel('UMAP2', fontsize=12,fontproperties=FontProperties(weight='bold'))
plt.show()


C1 = tmp_data.obs.copy()

tmp_data = pseudo_bulk_domain[~pseudo_bulk_domain.obs["cluster3"].isin(["3","5","2","11","9","14","28","4","25"]),:]
sc.pp.neighbors(tmp_data, n_pcs =0,n_neighbors=10)
sc.tl.umap(tmp_data)

sc.tl.leiden(tmp_data, resolution=3)
tmp_data.obs["tmp_cluster"] = tmp_data.obs["leiden"]

sc.set_figure_params(fontsize = 12)
sc.pl.umap(tmp_data, color=["cluster3"],size=28,legend_loc="on data",legend_fontsize=8)

fig=sc.pl.umap(tmp_data, color=["tmp_cluster"],size=100,legend_loc="on data",legend_fontsize=12,show=False,return_fig=True)
ax=fig.axes[0]
ax.set_title("")
ax.set_xlabel('UMAP1', fontsize=12,fontproperties=FontProperties(weight='bold'))
ax.set_ylabel('UMAP2', fontsize=12,fontproperties=FontProperties(weight='bold'))
plt.show()

tmp_data.X = data_mat[~pseudo_bulk_domain.obs["cluster3"].isin(["3","5","2","11","9","14","28","4","25"]),:]

tmp_min_x = np.min(tmp_data.obsm["X_umap"][:,0])
tmp_max_x = np.max(tmp_data.obsm["X_umap"][:,0])
tmp_min_y = np.min(tmp_data.obsm["X_umap"][:,1])
tmp_max_y = np.max(tmp_data.obsm["X_umap"][:,1])


tmp_range_x = tmp_max_x - tmp_min_x
tmp_range_y = tmp_max_y - tmp_min_y

tmp_mapped_x = 0 + (tmp_data.obsm["X_umap"][:,0] - tmp_min_x) * (10 - 0) / tmp_range_x
tmp_mapped_y = 0 + (tmp_data.obsm["X_umap"][:,1] - tmp_min_y) * (10 - 0) / tmp_range_y

import matplotlib.cm as cm
from matplotlib.lines import Line2D
top_num = 3
label_set = []

for i in range(tmp_data.shape[0]):
#     m = np.argsort(tmp_data.X[i,:])[-top_num:]
    m = tmp_data.X[i,:] >= 0.08
    frac = tmp_data.X[i,:][m]
    labels = tmp_data.var_names[m]
    label_set += labels.tolist()
    
label_set = set(label_set)    
      

fig, ax = plt.subplots(dpi=50)
for i in range(tmp_data.shape[0]):
#     m = np.argsort(tmp_data.X[i,:])[-top_num:]
    m = tmp_data.X[i,:] >= 0.08
    frac = tmp_data.X[i,:][m]
    labels = tmp_data.var_names[m]
    colors = [color_dict[i] for i in labels]
    x = tmp_mapped_x[i]
    y = tmp_mapped_y[i]

    ax.pie(frac, startangle=90, counterclock=False,colors= colors,
           radius=0.1, center=(x,y))

      

ax.spines['left'].set_position('zero')
ax.spines['left'].set_linewidth(0.5)
ax.spines['left'].set_color('gray')

ax.spines['bottom'].set_position('zero')
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['bottom'].set_color('gray')

ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')

ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color_dict[label], markersize=18, label=label) for label in sorted(list(label_set))]
plt.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(3.5, 3),fontsize=18)
plt.show()


tmp_data.obs["New_Domain"] = None
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["11","10","9"])] = "Cell_Niche_12"
# tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["9"])] = "Domain_14"
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["14","24","7"])] = "Cell_Niche_13"
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["3"])] = "Cell_Niche_14"
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["13","17"])] = "Cell_Niche_15"
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["5"])] = "Cell_Niche_16"
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["4","6","23","8"])] = "Cell_Niche_17"
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["2"])] = "Cell_Niche_18"
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["1","12","18"])] = "Cell_Niche_19"
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["0","16","21"])] = "Cell_Niche_20"
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["19"])] = "Cell_Niche_21"
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["15","25"])] = "Cell_Niche_22"
tmp_data.obs["New_Domain"][tmp_data.obs["tmp_cluster"].isin(["22","20"])] = "Cell_Niche_23"

fig=sc.pl.umap(tmp_data, color=["New_Domain"],size=100,legend_loc="on data",legend_fontsize=12,show=False,return_fig=True)
ax=fig.axes[0]
ax.set_title("")
ax.set_xlabel('UMAP1', fontsize=12,fontproperties=FontProperties(weight='bold'))
ax.set_ylabel('UMAP2', fontsize=12,fontproperties=FontProperties(weight='bold'))
plt.show()

C2 = tmp_data.obs.copy()
pseudo_bulk_domain.obs["New_Domain"] = pd.concat([C1,C2])["New_Domain"][pseudo_bulk_domain.obs_names]

color_map_domain = {
             'Cell_Niche_1': 'red',
             'Cell_Niche_2': 'greenyellow',
             'Cell_Niche_3': 'dimgrey',
             'Cell_Niche_4': 'yellow',
             'Cell_Niche_5': 'magenta',
             'Cell_Niche_6': 'springgreen',
             'Cell_Niche_7': 'sienna',
             'Cell_Niche_8': 'blueviolet',
             'Cell_Niche_9': 'rosybrown',
             'Cell_Niche_10': 'tan',
             'Cell_Niche_11': 'pink',
             'Cell_Niche_12': 'orange',
             'Cell_Niche_13': 'violet',
             'Cell_Niche_14': 'lightsteelblue',
             'Cell_Niche_15': 'cyan',
             'Cell_Niche_16': 'purple',
             'Cell_Niche_17': 'black',
             'Cell_Niche_18': 'bisque',
             'Cell_Niche_19': 'darkgreen',
             'Cell_Niche_20': 'blue',
             'Cell_Niche_21': 'lightcoral',
             'Cell_Niche_22': 'darkkhaki',
             'Cell_Niche_23': 'olivedrab',
             'Cell_Niche_24': 'gold'}

pseudo_bulk_domain.uns["New_Domain_colors"]=list(color_map_domain.values())


fig=sc.pl.umap(pseudo_bulk_domain, color=["New_Domain"],size=100,legend_fontsize=8,return_fig=True)
ax=fig.axes[0]
# ax.legend_.set_ncols(1)
# legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color_dict[label], markersize=8, label=label) for label in pseudo_bulk_domain.obs["Domain"].cat.categories]
ax.legend( loc='upper center', bbox_to_anchor=(0.5, -0.1),ncol = 4,borderpad=1,fontsize=8)
ax.set_title("")
ax.set_xlabel('UMAP1', fontsize=12,fontproperties=FontProperties(weight='bold'))
ax.set_ylabel('UMAP2', fontsize=12,fontproperties=FontProperties(weight='bold'))
plt.show()


map3 = dict()
for i in pseudo_bulk_domain.obs["New_Domain"].index:
    if i.split("_domain_")[0] not in map3.keys():
        map3[i.split("_domain_")[0]] = dict()
        map3[i.split("_domain_")[0]][i.split("_domain_")[1]] = pseudo_bulk_domain.obs["New_Domain"][i]
    else:
        map3[i.split("_domain_")[0]][i.split("_domain_")[1]] = pseudo_bulk_domain.obs["New_Domain"][i]


sc.set_figure_params(fontsize = 12)

n = 0
for i,data in enumerate(All_Data2):
    data = data.copy()
    sample_name = data.obs["sample"].cat.categories
    if list(sample_name)[0] in list(map3.keys()):
        map_tmp3 = map3[sample_name.to_list()[0]]
        data.obs["New_Domain"] = [map_tmp3[i] for i in data.obs["louvain"]]
#         data.obs["New_Domain"] = data.obs["New_Domain"].astype("category")
#         data.uns["New_Domain_colors"] = [ color_map_domain[i] for i in data.obs["New_Domain"].cat.categories]
#         data.obs["cluster4"] = [map_tmp4[i] for i in data.obs["louvain"]]
        print(sample_name)
        sc.pl.spatial(data,color = ["New_Domain"],title="color_random")
        data.uns["New_Domain_colors"] = [ color_map_domain[i] for i in data.obs["New_Domain"].cat.categories]
        sc.pl.spatial(data,color = ["New_Domain"])
        mat_c = data.obs[list(map(lambda x:x.split("meanscell_abundance_w_sf_")[1],data.obsm["means_cell_abundance_w_sf"].columns))+["New_Domain"]].groupby("New_Domain").sum()
        mat_c = mat_c.transpose()
        mat_c = mat_c/mat_c.sum(0)
        mat_c.columns = [ str(i) + "|" + sample_name.tolist()[0] for i in mat_c.columns.tolist()]
        
        if n == 0:
            mat_c_all = mat_c
            n += 1
        else:
            mat_c_all = pd.merge(mat_c_all,mat_c,left_index=True,right_index=True)

n = 0
for i in pseudo_bulk_domain.obs["New_Domain"].cat.categories:
    m = [x.split("|")[0] == i for x in mat_c_all.columns]
    mat_inte = mat_c_all.loc[:,m].sum(1)/mat_c_all.loc[:,m].sum(1).sum()
    symbol = mat_inte.index
    spatial_domain_name = [i]
    if n == 0:
        ref2 = pd.DataFrame(mat_inte.tolist(),index = symbol,columns=spatial_domain_name)
        n += 1
    else:
        ref2 = pd.merge(ref2,pd.DataFrame(mat_inte.tolist(),index = symbol,columns=spatial_domain_name),
                           left_index=True,right_index=True,how = "outer")


# CN_stack_plot
ref3 = ref2.copy()
ref3[ ref3 < 0.08 ] = 0
ref3 = ref3.loc[ref3.sum(1)!=0,:]
ref3 = ref3.T.div(ref3.T.sum(axis=1), axis=0)

fig,ax = plt.subplots(figsize=(5,3),dpi=200)
ax = ref3.plot(kind='bar', stacked=True,ax=ax,color=[color_dict[col] for col in ref3.columns])

ax.grid(True, zorder=0)
ax.set_axisbelow(True)    
plt.xlabel('Index')
plt.ylabel('Value')
plt.title('Stacked Bar Plot')
plt.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig("CN_stack_plot.pdf", format='pdf')
plt.show()


from collections import Counter
import matplotlib.pyplot as plt

Tac = Counter(pseudo_bulk_domain.obs.loc[~pseudo_bulk_domain.obs.duplicated(subset=["sample","New_Domain"]),:].copy()["New_Domain"])
Tac = sorted(Tac.items(),  key=lambda d: d[1], reverse=True)
Tac = dict(Tac)

color_map_domain = dict(zip(pseudo_bulk_domain.obs["New_Domain"].cat.categories,pseudo_bulk_domain.uns["New_Domain_colors"]))

sc.set_figure_params(fontsize=12)

fig,ax = plt.subplots(figsize=(8,3),dpi=200)
ax.bar(Tac.keys(),Tac.values(),width=0.9,color=[color_map_domain[i] for i in Tac.keys()],align="center")
ax.set_axisbelow(True)    
plt.xticks(rotation=90)
plt.ylim(0,48)

ax_pie = plt.gca().twiny().twinx()

for i,domain_name in enumerate(Tac.keys()):
    m = ref2[domain_name][ref2[domain_name]>=0.08].sort_values(ascending=False)
    labels = m.index
    ax_pie.pie(m,colors=[color_dict[i] for i in labels],wedgeprops=dict(width=1),radius=2,center=(i*5,Tac[domain_name]))
    ax_pie.axis('equal')

plt.show()

result = []
ps = []
for i in ref2.index:
    for j in ref2.index:
        if j not in ps:
            if ref2.loc[[i,j],:].min().max() >= 0.08:
                result.append([i,j,ref2.loc[[i,j],:].min().max(),ref2.loc[[i,j],:].min().idxmax()])
    ps.append(i)
    
result2 = pd.DataFrame(result,columns=["Cell1","Cell2","weight","Cell_Niche"])


import igraph as ig


connections = result
graph = ig.Graph.TupleList(connections, directed=False, edge_attrs=["weight", "type"])
edge_weights = graph.es["weight"]
edge_types = graph.es["type"]

edge_widths = [weight*20 for weight in edge_weights]
edge_colors = color_map_domain  

edge_colors_list = [edge_colors[edge_type] for edge_type in edge_types]

vertex_labels = graph.vs["name"]

visual_style = {
    "vertex_size": 20,
    "vertex_label": vertex_labels,
    "edge_width": edge_widths,
    "edge_color": edge_colors_list,
    "edge_label_dist": 1.5,
    "layout": graph.layout("kk"),
    "bbox": (600*1.2, 500*1.2),
    "margin": 100
    
}

plot = ig.plot(graph, **visual_style,)
plot


batch_colors = {"interface":"red","core":"blue","immunity":"yellow","PT":"green","TT":"black"}

map_sample = dict()
tmp_data = pseudo_bulk_domain.obs.loc[~pseudo_bulk_domain.obs.duplicated(subset=["sample","batch"]),:].copy()
for i in tmp_data.index:
    map_sample[tmp_data.loc[i,"sample"]] = batch_colors[tmp_data.loc[i,"batch"]]


from sklearn.metrics.pairwise import cosine_similarity
sc.set_figure_params(dpi=150,figsize=[15,15],fontsize=5)

pseudo_bulk_domain.obs["values"] = 1
tmp_data = pseudo_bulk_domain.obs.loc[~pseudo_bulk_domain.obs.duplicated(subset=["sample","New_Domain"]),:].copy()
sample_to_domain = tmp_data.pivot(index="sample",columns="New_Domain",values="values").fillna(0)
sample_corr = pd.DataFrame(cosine_similarity(sample_to_domain),index=sample_to_domain.index,columns=sample_to_domain.index)

b = seaborn.clustermap(sample_corr,method="ward",col_colors=[map_sample[i] for i in sample_corr.columns])
seaborn.clustermap(sample_to_domain.loc[sample_to_domain.index[b.dendrogram_row.reordered_ind],:],standard_scale="none",method="ward",row_cluster=False)

from collections import Counter
for i,data in enumerate(All_Data2):
    data = data.copy()
    sample_name = data.obs["sample"].cat.categories
    if list(sample_name)[0] in list(map3.keys()):
        map_tmp3= map3[sample_name.to_list()[0]]
        data.obs["New_Domain"] = [map_tmp3[i] for i in data.obs["louvain"]]
        cm = Counter(data.obs["New_Domain"])
        cmt = sum(cm.values())
        for i in cm.keys():
            cm[i] = cm[i]/cmt
        tmp_data.loc[ tmp_data["sample"] == sample_name.tolist()[0],"values"] = tmp_data.loc[ tmp_data["sample"] == sample_name.tolist()[0],:].apply(lambda x:cm[x["New_Domain"]],axis=1)

sample_to_domain = tmp_data.pivot(index="sample",columns="New_Domain",values="values").fillna(0)
sample_corr = pd.DataFrame(cosine_similarity(sample_to_domain),index=sample_to_domain.index,columns=sample_to_domain.index)
b = seaborn.clustermap(sample_corr,method="ward",col_colors=[map_sample[i] for i in sample_corr.columns])
seaborn.clustermap(sample_to_domain.loc[sample_to_domain.index[b.dendrogram_row.reordered_ind],:],method="ward",row_cluster=False)

from collections import Counter
for i,data in enumerate(All_Data2):
    data = data.copy()
    sample_name = data.obs["sample"].cat.categories
    if list(sample_name)[0] in list(map3.keys()):
        map_tmp3= map3[sample_name.to_list()[0]]
        data.obs["New_Domain"] = [map_tmp3[i] for i in data.obs["louvain"]]
        cm = Counter(data.obs["New_Domain"])
        tmp_data.loc[ tmp_data["sample"] == sample_name.tolist()[0],"values"] = tmp_data.loc[ tmp_data["sample"] == sample_name.tolist()[0],:].apply(lambda x:cm[x["New_Domain"]],axis=1)

sample_to_domain = tmp_data.pivot(index="sample",columns="New_Domain",values="values").fillna(0)



import scipy.stats as stats
result = []
for i in sample_to_domain.columns:
    x = (sample_to_domain.transpose()).loc[i,pseudo_bulk_domain.obs["sample"][(pseudo_bulk_domain.obs["TT1"] == "with_TT") & (pseudo_bulk_domain.obs["batch"].isin(["PT","immunity","core"]))].unique()]
    y = (sample_to_domain.transpose()).loc[i,pseudo_bulk_domain.obs["sample"][(pseudo_bulk_domain.obs["TT1"] == "without_TT")& (pseudo_bulk_domain.obs["batch"].isin(["PT","immunity","core"]))].unique()]
    result.append([i,np.log(x.mean()+1)-np.log(y.mean()+1),stats.wilcoxon(x.tolist(),y.tolist())[1]])

comparison_result = pd.DataFrame(result)
comparison_result.to_csv("CN_comparison_result.csv")

i = "Cell_Niche_15"
x = (sample_to_domain.transpose()).loc[i,pseudo_bulk_domain.obs["sample"][(pseudo_bulk_domain.obs["TT1"] == "with_TT") & (pseudo_bulk_domain.obs["batch"].isin(["PT","immunity","core"]))].unique()]
y = (sample_to_domain.transpose()).loc[i,pseudo_bulk_domain.obs["sample"][(pseudo_bulk_domain.obs["TT1"] == "without_TT")& (pseudo_bulk_domain.obs["batch"].isin(["PT","immunity","core"]))].unique()]
lm_data = pd.concat([pd.concat([x,y]),
                       pd.Series(x.index.to_list()+y.index.to_list(),
                                 index=x.index.to_list()+y.index.to_list()
                                ),
                       pd.Series(["with_TT"]*8 + ["without_TT"]*8 ,index=x.index.to_list()+y.index.to_list()),
                       pd.Series(["core"]+["immunity"]*3+["PT"]*4+["immunity"]*8,index=x.index.to_list()+y.index.to_list())
                      ],axis=1)

lm_data.columns = ["CN_area","Section_ID","TT_status","data_source"]
metadata = pd.read_csv("/home/ncpsb/new_mnt/ccRCC_data/Spatial_transcriptome/data/ST_meta.csv",index_col=0)
lm_data = pd.concat([lm_data,metadata.loc[lm_data.index,:]],axis=1)
lm_data['TT_status'] = pd.Categorical(lm_data['TT_status'], categories=['without_TT', 'with_TT'])
lm_data['data_source'] = pd.Categorical(lm_data['data_source'], categories=['PT', 'core','immunity'])
lm_data['sample_type'] = pd.Categorical(lm_data['sample_type'], categories=['Frozen', 'FFPE'])
lm_data

import pandas as pd
import statsmodels.formula.api as smf

model = smf.ols('CN_area ~ TT_status + data_source + sample_type', data=lm_data).fit()
print(model.pvalues)
print(model.summary())


coef = model.params
conf = model.conf_int()
conf.columns = ['lower', 'upper']
pvals = model.pvalues

forest_data = pd.DataFrame({
    'coef': coef,
    'lower': conf['lower'],
    'upper': conf['upper'],
    'pval': pvals
})
forest_data['variable'] = forest_data.index

forest_data = forest_data[forest_data['variable'] != 'Intercept']

plt.figure(figsize=(4, 3))
plt.errorbar(forest_data['coef'], forest_data['variable'],
             xerr=[forest_data['coef'] - forest_data['lower'], 
                   forest_data['upper'] - forest_data['coef']],
             fmt='o', capsize=4)

for i, row in forest_data.iterrows():
    plt.text(row['upper'] + 20, row['variable'], f"P={row['pval']:.3f}", va='center')

plt.axvline(x=0, color='gray', linestyle='--')  
plt.xlabel('Coefficient')
plt.title('Forest Plot of Multiple Linear Regression Coefficients')
plt.tight_layout()
x_min = forest_data['lower'].min()
x_max = forest_data['upper'].max()

plt.xlim(x_min-200, x_max + 200)

plt.savefig("forest_plot_Multiple_Linear_Regression.pdf", format='pdf')
plt.show()


#Figure 2F
colors = ['#0000FF', '#FF0000']

plt.figure(figsize=(3, 3), dpi=200)

# Violin plot
parts = plt.violinplot([y, x], showmeans=False, showmedians=False, showextrema=False)

for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor(colors[i])
    pc.set_edgecolor('black')
    pc.set_alpha(0.7)


bp = plt.boxplot([y, x], showfliers=False, widths=0.1, patch_artist=True)
for i, box in enumerate(bp['boxes']):
    box.set_facecolor(colors[i])
    box.set_edgecolor('black')
    box.set_alpha(0.7)


plt.text(1.5, max(max(x), max(y)) + 1, f'P value={p_value:.3f}',
         fontsize=10, ha='center', va='center', color='black')


plt.xticks([1, 2], ["without_TT", "with_TT"], fontsize=12)
plt.ylabel("The area of Cell Niche 15\n(spot numbers)", fontsize=12)


plt.tight_layout()
plt.savefig("Figure 2F.pdf", format='pdf')
plt.show()
