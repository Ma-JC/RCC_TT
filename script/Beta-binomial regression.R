library(Seurat)
library(harmony)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(ggpointdensity)
library(viridis)

data = readRDS("~/mnt/Analysis_for_each_Patient/02.Cell annotation and remove doublets in all patients/discovery_datasets/data/AllData/merge_postprocess.2.rds")
data$All_second_cluster = as.character(data$All_second_cluster)
data$All_third_cluster = data$All_second_cluster

Mesen = readRDS("~/mnt/Analysis_for_each_Patient/02.Cell annotation and remove doublets in all patients/discovery_datasets/data/CellData/Mesenchymal.rds")
data$All_third_cluster[ paste("Mesenchymal_cell",colnames(Mesen),sep = "_") ] = as.character(Mesen$All_third_cluster)

sub_data = data[,data$Sample_Tag %in% c("ccRcc_CT") & data$All_third_cluster != "Noise" ]

tmp_res = vector()

group = list("CAF" = setdiff(unique(Mesen$All_third_cluster),"Noise"),
    "EC" = c('Endothelial cell 5(Cancer)','Endothelial cell 6(Cancer)','Endothelial cell 2(Arterioles)','Endothelial cell 4(Vasa Recta)','Endothelial cell 3(Vein)','Endothelial cell 1(Glomerular Capillaries)'),
    "Epi_cancer"=c("Epithelial cell","Cancer cell"),
    "immune" = vector())

group[["immune"]] = setdiff(unique(data$All_third_cluster),c("Noise",group[["CAF"]],group[["EC"]],group[["Epi_cancer"]]))

for(i in names(group)){
   

    res = reshape2::melt(table(sub_data$Patient,sub_data$All_third_cluster)[,group[[i]]]/rowSums(table(sub_data$Patient,sub_data$All_third_cluster)[,group[[i]]]))
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

tmp2 = dcast(res,formula = Var1~Var2,fill="value")
tmp2$Var1 = NULL
tmp2 = as.data.frame(apply(tmp2,2,as.numeric))
tmp3 = scale(tmp2)
tmp3 = as.data.frame(apply(tmp2,2,function(x){(x-min(x))/(max(x)-min(x))}))
fibro_subtype = c("FAP+ Fibroblast","CYSLTR2+ Fibroblast","MFAP5+ Fibroblast","FBLN5+ Fibroblast","ADH1B+ Fibroblast")
BN_MATRIX = cbind(Total_count = rowSums(table(sub_data$Patient,sub_data$All_third_cluster)[,group[["immune"]]]),
                  'NK cell' = table(sub_data$Patient,sub_data$All_third_cluster)[,c("NK cell")],
                  tmp3[,fibro_subtype])

BN_MATRIX = as.data.frame(BN_MATRIX)


BN_MATRIX$NonNK <- BN_MATRIX$Total_count - BN_MATRIX$`NK cell`
BN_MATRIX$NK_count = BN_MATRIX$`NK cell`

library(glmmTMB)

model = glmmTMB(cbind(NK_count, Total_count - NK_count) ~ `FAP+ Fibroblast`+`CYSLTR2+ Fibroblast`+`ADH1B+ Fibroblast`+`FBLN5+ Fibroblast`+`MFAP5+ Fibroblast`,
        family = betabinomial(),
        data = BN_MATRIX)

summary(model)

library(ggplot2)
library(dplyr)

# 提取模型结果
coefs <- summary(model)$coefficients$cond
df <- data.frame(
  term = rownames(coefs),
  estimate = coefs[, "Estimate"],
  std_error = coefs[, "Std. Error"],
  p_value = coefs[, "Pr(>|z|)"]
)

# 计算置信区间
df <- df %>%
  mutate(
    conf_low = estimate - 1.96 * std_error,
    conf_high = estimate + 1.96 * std_error
  ) %>%
  filter(term != "(Intercept)")

# 手动定义分组（根据你的变量名调整）
df <- df %>%
  mutate(
    group = case_when(
      term %in% c("FAP+ Fibroblast", "CYSLTR2+ Fibroblast") ~ "Group 1: Fibroblast Types",
      term %in% c("ADH1B+ Fibroblast", "FBLN5+ Fibroblast", "MFAP5+ Fibroblast") ~ "Group 2: Other Fibroblasts",
      TRUE ~ "Other"
    ),
    is_reference = FALSE
  )

# 给每组第一个变量设为reference（举例）
df <- df %>%
  group_by(group) %>%
  mutate(is_reference = ifelse(row_number() == 1, TRUE, FALSE)) %>%
  ungroup()

# 调整因子顺序，按estimate降序排列
df <- df %>%
  arrange(desc(estimate)) %>%
  mutate(term = factor(term, levels = term))

# 绘图
options(repr.plot.height = 5, repr.plot.width = 7)
p = ggplot(df, aes(x = estimate, y = term)) +
  geom_point(size = 5,color="#1F77B4") +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high),size=1, height = 0.3,color="#1F77B4") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  # “reference”标签
#   geom_text(data = df %>% filter(is_reference), aes(label = "reference", x = max(conf_high) + 0.3), hjust = 0, size = 3) +
  # p值标签，只给非reference且有p值的
  geom_text(data = df, aes(label = paste0("p = ", signif(p_value, 3)), x = conf_high + 0.3), hjust = 0, size = 6) +
#   facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.4))) + # 右边留白
  labs(
    x = "Coefficient",
    y = NULL,
    title = "Beta-binomial regression for NK cell abundance"
  ) +
  theme_bw() +
  theme(
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor = element_blank(),
    strip.text.y = element_text(angle = 0, face = "bold"),
    plot.title = element_text(hjust = 0.5,face = "bold",size = 15),
    axis.text.y = element_text(size = 14,face = "bold"),
    axis.text.x = element_text(size = 12,face = "bold"),
      axis.title.x = element_text(size = 12,face = "bold")
  )
p

pdf(file = "Beta_binomial_NK.pdf",width = 7,height = 5)
p
dev.off()
