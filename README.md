# 描述

本包主要目的对`TCGA`上面的`miRNA`基因进行差异基因分析。主要使用手段有：

1. 使用`edgeR`包推荐方法进行分析
2. 使用`edgeR`包推荐classic analysis method进行分析
3. 使用`DESeq2`包推荐方法进行分析
4. 使用`limma`包推荐方法进行分析
5. 使用ggplot绘制火山图
6. 使用ggplot绘制热图
7. 使用ggplot绘制PCA图（主成分分析图）
8. 使用edgeR绘制BCV图（生物学差异图）

# 使用方法

主要函数说明

## main函数

主要使用`edgeR`包进行经典基因差异分析。

### 参数说明

```R
main(path)
```

+ path 为解压后的数据的总文件

```R
path = "C:\\Users\\tanzicai\\Downloads\\gdc_download_20210915_122538.396858"
```

+ 注意需要把'metadata.cart.2021-09-15.json'文件放入总文件夹的同级文件夹`

```R
C:\\Users\\tanzicai\\Downloads
```

+ 注意保留解压后的'MANIFEST.txt'

## main2函数


主要使用`edgeR`、`DESeq2`、`limma`进行包进行经典基因差异分析。
```R
main()
```

## 筛选条件

+ 默认筛选参数 PValue = 0.05

```R
main(path)
```

+ 可选参数  PValue

```R
main(path,pValue)
```

```R
main(path,p = 0.5)
```

## 安装方法

```R
install.packages("devtools")
library(devtools)
install_github("tanzicai/miRNA_TCGA_Analysis/miRNADiff")
library(miRNADiff)
```



