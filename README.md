## 描述

本包主要目的是对'TCGA'上面的'miRNA'基因进行差异基因分析的，主要使用了'edgeR'包进行分析，并对差异基因生成火山图。

## 使用方法

### 参数说明

+ path 为解压后的数据的总文件

```R
path = "C:\\Users\\tanzicai\\Downloads\\gdc_download_20210915_122538.396858"
```

+ 注意需要把'metadata.cart.2021-09-15.json'文件放入总文件夹的同级文件夹`

```R
C:\\Users\\tanzicai\\Downloads
```

+ 注意保留解压后的'MANIFEST.txt'


### 开始分析

+ 调用主函数'main'开始基因分析

```R
main(path)
```

## 主要使用的R包

+ rjson
+ BiocManager
+ edgeR
+ statmod
+ limma

