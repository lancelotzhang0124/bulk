#已经做好了表达矩阵matrix1和matrix2
################预处理##################
# 去除表达量过低的基因
# 只选择纯数字列
numeric_data1 <- matrix1[, -1]
numeric_data2 <- matrix2[, -1]
# 计算行均值，并过滤基于行均值
filtered_rows1 <- rowMeans(numeric_data1) > 1
filtered_rows2 <- rowMeans(numeric_data2) > 1
# 使用这些行选择原始数据
matrix1 <- matrix1[filtered_rows1, ]
matrix2 <- matrix2[filtered_rows2, ]

saveRDS(matrix1, file="matrix1.rds")
saveRDS(matrix2, file="matrix2.rds")

#############构建分组向量##############
library("stringr")

group_list=c(rep('control',2),rep('USP11KO',2))
## 强制限定顺序
group_list <- factor(group_list,levels = c("control","USP11KO"),ordered = T)

#############差异分析################

# 为每个样本转换其每个基因的FPKM到TPM

FPKMtoTPM <- function(matrix1)
{
  exp(log(matrix1) - log(sum(matrix1)) + log(1e6))
}

tpm_M1 <- apply(matrix1,2,FPKMtoTPM)
tpm_M1[1:3,]
colSums(tpm_M1)

tpm_M2 <- apply(matrix2,2,FPKMtoTPM)
tpm_M2[1:3,]
colSums(tpm_M2)

library(limma)
library(edgeR)
library(ggplot2)

exprSet <- tpm_M1
#检测数据分布
pdf("test.pdf")
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
dev.off()
#exprSet=normalizeBetweenArrays(exprSet) #标准化
#exprSet <- log2(exprSet+1) #正态分布

#差异分析
DEG_Analysis <- function(exprSet, group_list) {
  
  #数据预处理
  exprSet=normalizeBetweenArrays(exprSet) #标准化
  exprSet <- log2(exprSet+1) #正态分布
  # 设定小数位数
  options(digits = 4)
  # 加载必要的库
  library(limma)
  
  # 设计矩阵
  design = model.matrix(~factor(group_list))
  
  # 应用线性模型
  fit = lmFit(exprSet, design)
  fit = eBayes(fit)
  
  # 获取差异表达的基因
  deg = topTable(fit, coef=2, adjust='BH', number=Inf)
  
  # 返回结果和函数
  list(deg=deg)
}
zero_variance_genes <- apply(exprSet, 1, var) == 0
names(zero_variance_genes)[zero_variance_genes]


results_M1 = DEG_Analysis(tpm_M1, group_list = group_list)
results_M2 = DEG_Analysis(tpm_M2, group_list = group_list)

bp = function(g) {
  library(ggpubr)
  df = data.frame(gene=g, stage=group_list)
  
  pdf("bp_boxplot.pdf")
  # 绘制盒状图
  p = ggboxplot(df, x = "stage", y = "gene", color = "stage", palette = "jco", add = "jitter")
  
  # 添加p值
  p + stat_compare_means()
  dev.off()
}

deg=topTable(fit,coef=2,adjust='BH',number = Inf)
head(deg) 

# 增加一列到results_M1，表示-log10(p-value)
results_M1$deg$neg_log10_pvalue <- -log10(results_M1$deg$P.Value)
results_M2$deg$neg_log10_pvalue <- -log10(results_M2$deg$P.Value)


# 定义阈值
alpha <- 0.05
lfc_threshold <- 1

log2FC_cutoff <- 1
padj_cutoff <- 0.01


significant_genes_M1 <- results_M1$deg[abs(results_M1$deg$logFC) > log2FC_cutoff & results_M1$deg$P.Value < padj_cutoff, ]
significant_genes_M2 <- results_M2$deg[abs(results_M2$deg$logFC) > log2FC_cutoff & results_M2$deg$P.Value < padj_cutoff, ]


# 绘制火山图
volcano <- ggplot(results_M1$deg, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = (P.Value < padj_cutoff & abs(logFC) > log2FC_cutoff)), size = 1) + 
  geom_point(data = significant_genes_M1, aes(x = logFC, y = -log10(P.Value)), color = "blue", size = 2) +
  geom_text_repel(data = significant_genes_M1, aes(label = rownames(significant_genes_M1)), nudge_y = 0.5, size = 3) +
  scale_color_manual(values = c("grey", "blue")) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), linetype = "dashed", color = "black") +
  labs(title = "DEG in E13.5", x = "Log2 Fold Change", y = "-Log10 p-value", color = "Significant") +
  theme_minimal()

# 只有在火山图创建成功后，才执行保存操作
ggsave(filename = "volcano_plot_M1.pdf", plot = volcano, path = 'Figures/', dpi = 300)

volcano <- ggplot(results_M2$deg, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = (P.Value < padj_cutoff & abs(logFC) > log2FC_cutoff)), size = 1) + 
  geom_point(data = significant_genes, aes(x = logFC, y = -log10(P.Value)), color = "blue", size = 2) +
  geom_text_repel(data = significant_genes_M2, aes(label = rownames(significant_genes)), nudge_y = 0.5, size = 3) +
  scale_color_manual(values = c("grey", "blue")) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), linetype = "dashed", color = "black") +
  labs(title = "DEG in E18.5", x = "Log2 Fold Change", y = "-Log10 p-value", color = "Significant") +
  theme_minimal()

# 只有在火山图创建成功后，才执行保存操作
ggsave(filename = "volcano_plot_M2.pdf", plot = volcano, path = 'Figures/', dpi = 300)



