# 描述

`Version0.1.0`

本包主要目的对`TCGA`上面的`miRNA`基因进行差异基因分析。主要使用手段有：

1. 使用`edgeR`包推荐方法进行分析
2. 使用`edgeR`包推荐classic analysis method进行分析
3. 使用`DESeq2`包推荐方法进行分析
4. 使用`limma`包推荐方法进行分析
5. 使用ggplot绘制火山图
6. 使用ggplot绘制热图
7. 使用ggplot绘制PCA图（主成分分析图）
8. 使用edgeR绘制BCV图（生物学差异图）

`version0.2.0`

本包主要目的对`TCGA`上面的`miRNA`基因进行差异基因分析(集成)。主要使用手段有：

1. 使用`edgeR`包推荐方法进行分析
2. 使用`edgeR`包推荐classic analysis method进行分析
3. 使用`DESeq2`包推荐方法进行分析
4. 使用`limma`包推荐方法进行分析
5. 使用ggplot绘制火山图
6. 使用ggplot绘制热图
7. 使用ggplot绘制PCA图（主成分分析图）
8. 使用edgeR绘制BCV图（生物学差异图）

# 使用方法

+ main() : 主要使用`edgeR`包进行经典基因差异分析。

+ draw_pca(dat,group_list,path) : 绘制PCA分析图

+ draw_heatmap(dat,group_list,path):绘制热图

+ envirment() : 安装及检测本包所需要的环境

+ autoAnalyssisSingleCancer(cancer_type) : 单癌症差异基因分析(包括DESeq2,edgeR,limma)

+ autoAnalyssisByDESeq2AndedgeR(cancer_type) : 单癌症差异基因分析(包括DESeq2,edgeR)

+ runByDEseq2AndEdgeRmain(cancer_type) : 单癌症差异基因分析(包括DESeq2,edgeR)

## 参数说明

+ path 为解压后的数据的总文件

```R
path = "C:\\Users\\tanzicai\\Downloads\\gdc_download_20210915_122538.396858"
```

+ dat 列名时样本名，行名时探针名

```
dat = log(原始样本数据+1)
```

+ group_list

```
经因子化后的分组信息表
```

+ path 储存路径
+ cancer_type 

```
癌症类型:TCGA-KICH
```

## 安装方法

```R
install.packages("devtools")
library(devtools)
install_github("tanzicai/miRNA_TCGA_Analysis/miRNADiff")
library(miRNADiff)
```



