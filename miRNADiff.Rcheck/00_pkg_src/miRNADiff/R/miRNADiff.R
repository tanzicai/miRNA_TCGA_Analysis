#' @title miRNADiff

#' @description use edgeR to analyze the different genes

#' @param path

#' @return NULL

#' @examples main("C:\\Users\\tanzicai\\OneDrive\\文档\\cc_test\\test")

#' @export main


##' tips warming Keep the MANIFEST.txt file

##' run the "main" function to start analysis

#' @functionName munzip

#' @author tanzicai

#' @purpose move all files to a folder

#' @createtime 2021.9.17

#' @time 2021.9.20

#' @input work path

#' @return none

#' @param create folder named "temp" and "source"

#  move all files from source forlder

#  delete the file called "annotations.txt"

munzip <- function( path ){

  dir.create("./temp")

  dir.create("./temp/source")

  filePath = dir(path = path ,full.names = T)

  filePath = filePath[1:length(filePath)]

  print(length(filePath))

  for(wd in filePath){

    files = dir(path = wd)

    fromFilePath = paste0(wd,"/",files)

    toFilePath = paste0("./temp/source/",files)

    file.copy(fromFilePath,toFilePath)
  }

  TRUE
}



#' @functionName combin

#' @author tanzicai

#' @purpose read all files and merge a table

#' @createtime 2021.9.17

#' @time 2021.9.20

#' @input temp files' source path

#' @return merged the table DATA

#' @param create read all the files in the source forlder
#  merge the files and create a table
#  return the table
combin <- function(path){

  countFilePath = dir(path = "./temp/source",pattern = "*.txt")

  counts_merge = NULL

  for (i in countFilePath) {

    x=read.delim(paste0("./temp/source/",i),col.names = c("miRNA_ID",substr(i,1,9),substr(i,1,9),substr(i,1,9)))

    if(is.null(counts_merge)){

      counts_merge = x

    }else{

      counts_merge =merge(counts_merge,x,by = "miRNA_ID")

    }

  }

  write.csv(counts_merge,file = "./temp/源文件合并.csv")

  rownames(counts_merge)=counts_merge$miRNA_ID

  DATA = counts_merge

  save(DATA,file = "./DATA.Rdata")

  TRUE
}















#' @functionName convertName

#' @author tanzicai

#' @purpose Comment file extraction

#' @createtime 2021.9.17

#' @time 2021.9.20

#' @input work space path

#' @return void

#' @param read json file
#  read json's file name / TCGA code and create a table
#return the table



convertName <- function(path){

  if(require("rjson")){

    print("成功载入包:rjson")

  } else {

    print("不存在这个包，正在尝试安装")

    install.packages("rjson")

    if(require("rjson")){

      print("成功安装并载入包:rjson")

    } else {

      stop("安装失败")

    }
  }

  library(rjson)

  metaJSON = fromJSON(file = paste0("../",dir(path = "../",pattern = "*.json",full.names = F)))

  jsonDataInfo = data.frame(fileName = c(),TCGA_Barcode = c())

  for (i in 1:length(metaJSON)) {

    TCGA_barcode = metaJSON[[i]][["associated_entities"]][[1]][["entity_submitter_id"]]

    fileName = metaJSON[[i]][["file_name"]]

    jsonDataInfo = rbind(jsonDataInfo,data.frame(fileName = substr(fileName,1,9),TCGA_barcode = substr(TCGA_barcode,1,15)))

  }

  rownames(jsonDataInfo) = jsonDataInfo[,1]

  write.csv(jsonDataInfo,file = "../temp/文件名与TCGA编号对照表.csv")


  DATA = read.csv("../temp/源文件合并.csv")

  DATA = DATA[-1]

  fileNameToTCGA_Barcode = jsonDataInfo[-1]

  jsonDataInfo = jsonDataInfo[order(jsonDataInfo$fileName),]

  data = DATA

  DATA = DATA[-1]

  colnames(DATA) = jsonDataInfo$TCGA_barcode

  DATA = cbind(data[1],DATA)

  write.csv(DATA,file = "../temp/重命名_列为TCGA编号.csv")

  TRUE
}














#' @functionName edgRAnalyze

#' @author tanzicai

#' @purpose use edgR carries on difference analysis

#' @createtime 2021.9.17

#' @time 2021.9.20

#' @input work space path

#' @return Boolean

#' @param read json file
#  read json's file name / TCGA code and create a table

edgRAnalyze <- function(DATA,pValue){


  if(require("BiocManager")){

    print("成功载入包:BiocManager")

  } else {

    print("不存在这个包，正在尝试安装")

    install.packages("BiocManager")

    if(require("BiocManager")){

      print("成功安装并载入包")

    } else {

      stop("安装失败")

    }
  }
  if(require("edgeR")){

    print("成功载入包:edgeR")

  } else {

    print("不存在这个包，正在尝试安装")

    BiocManager::install("edgeR")

    BiocManager::install("limma")
  }

  if(require("edgeR")){

    print("成功安装并载入包")

  } else {

    stop("安装失败")

  }



  library(edgeR)

  library(limma)

  #将第一列换成行名

  DATA = read.csv("../temp/重命名_列为TCGA编号.csv")

  DATA = DATA[-1]

  row.names(DATA) <- DATA[, 1]

  DATA <- DATA[, -1]

  #分组信息

  group_list <- ifelse(as.numeric(substr(colnames(DATA),14,15))<10,"tumor","normal")

  group_list <- factor(group_list,levels = c("normal","tumor"))

  print(table(group_list))

  dge <- DGEList(counts=DATA,group=group_list)



  dge$samples$lib.size <- colSums(dge$counts)

  dge <- calcNormFactors(dge)

  design <- model.matrix(~0+group_list)

  rownames(design)<-colnames(dge)

  colnames(design)<-levels(group_list)

  dge <- estimateGLMCommonDisp(dge,design)

  dge <- estimateGLMTrendedDisp(dge, design)

  dge <- estimateGLMTagwiseDisp(dge, design)

  fit <- glmFit(dge, design)

  fit2 <- glmLRT(fit, contrast=c(-1,1))

  DEG2=topTags(fit2, n=nrow(DATA))

  DEG2=as.data.frame(DEG2)

  logFC_cutoff2 <- with(DEG2,mean(abs(logFC)) + 2*sd(abs(logFC)))

  DEG2$change = as.factor(ifelse(DEG2$PValue < pValue & abs(DEG2$logFC) > logFC_cutoff2,ifelse(DEG2$logFC > logFC_cutoff2 ,'UP','DOWN'),'NOT'))

  print("本次差异基因分析结果为")

  print(table(DEG2$change))

  edgeR_DEG <- DEG2

  dir.create("../out_put")

  write.csv(edgeR_DEG,"../out_put/差异分析结果.csv")

  TRUE

}








#' @functionName edgRAnalyze2

#' @author tanzicai

#' @purpose use edgR carries on difference analysis(second way)

#' @createtime 2021.9.17

#' @time 2021.9.20

#' @input work space path

#' @return Boolean

#' @param use edgR carries on difference analysis(second way)

edgRAnalyze2 <- function(DATA,pValue){


  if(require("BiocManager")){

    print("成功载入包:BiocManager")

  } else {

    print("不存在这个包，正在尝试安装")

    install.packages("BiocManager")

    if(require("BiocManager")){

      print("成功安装并载入包")

    } else {

      stop("安装失败")

    }
  }
  if(require("edgeR")){

    print("成功载入包:edgeR")

  } else {

    print("不存在这个包，正在尝试安装")

    BiocManager::install("edgeR")

    BiocManager::install("limma")
  }

  if(require("edgeR")){

    print("成功安装并载入包")

  } else {

    stop("安装失败")

  }

  if(require("statmod")){

    print("成功载入包:statmod")

  } else {

    print("不存在这个包，正在尝试安装")

    install.packages("statmod")

  }

  if(require("edgeR")){

    print("成功安装并载入包:")

  } else {

    stop("安装失败")

  }


  library(edgeR)

  library(limma)

  library(statmod)

  y<-DGEList(counts = DATA[,2:length(DATA)],genes = DATA[,1])

  o <- order(rowSums(y$counts), decreasing=TRUE)

  y <- y[o,]

  d <- duplicated(y$genes$genes)

  y <- y[!d,]

  y$samples$lib.size <- colSums(y$counts)

  #重命名基因表达矩阵的行名

  rownames(y$counts)=y$genes$genes

  y<-calcNormFactors(y)

  plotMDS(y)

  group_list <- ifelse(as.numeric(substr(colnames(DATA),14,15))<10,"tumor","normal")

  Samples <- factor(group_list,levels = c("normal","tumor"))

  design=model.matrix(~Samples)

  rownames(design)=colnames(y)

  y1=estimateDisp(y,design,robust = T)

  plotBCV(y1)

  fit<-glmFit(y1,design)

  et1=glmLRT(fit)

  print("本次差异基因分析结果为")

  o <- order(et1$table$PValue)

  print(summary(decideTests(et1,p.value = pValue)))

  plotMD(et1)

  abline(h=c(-1, 1), col="blue")

  write.csv(et1,file = "../out_put/edgeR2差异分析结果.csv")

  TRUE
}












#' @functionName edgRAnalyze3

#' @author tanzicai

#' @purpose use edgR carries on difference analysis(second way)

#' @createtime 2021.9.17

#' @time 2021.9.20

#' @input work space path

#' @return Boolean

#' @param use edgR carries on difference analysis(third way)

edgRAnalyze3 <- function(DATA,pValue){


  if(require("BiocManager")){

    print("成功载入包:BiocManager")

  } else {

    print("不存在这个包，正在尝试安装")

    install.packages("BiocManager")

    if(require("BiocManager")){

      print("成功安装并载入包")

    } else {

      stop("安装失败")

    }
  }
  if(require("edgeR")){

    print("成功载入包:edgeR")

  } else {

    print("不存在这个包，正在尝试安装")

    BiocManager::install("edgeR")

    BiocManager::install("limma")
  }

  if(require("edgeR")){

    print("成功安装并载入包")

  } else {

    stop("安装失败")

  }

  if(require("statmod")){

    print("成功载入包:statmod")

  } else {

    print("不存在这个包，正在尝试安装")

    install.packages("statmod")

  }

  if(require("edgeR")){

    print("成功安装并载入包:")

  } else {

    stop("安装失败")

  }


  library(edgeR)

  library(limma)

  library(statmod)

  y<-DGEList(counts = DATA[,2:length(DATA)],genes = DATA[,1])

  o <- order(rowSums(y$counts), decreasing=TRUE)

  y <- y[o,]

  d <- duplicated(y$genes$genes)

  y <- y[!d,]

  y$samples$lib.size <- colSums(y$counts)

  #重命名基因表达矩阵的行名

  rownames(y$counts)=y$genes$genes

  y<-calcNormFactors(y)

  plotMDS(y)

  group_list <- ifelse(as.numeric(substr(colnames(DATA),14,15))<10,"tumor","normal")

  Samples <- factor(group_list,levels = c("normal","tumor"))

  design=model.matrix(~Samples)

  rownames(design)=colnames(y)

  y1=estimateDisp(y,design,robust = T)

  plotBCV(y1)

  Samples = Samples[-1]

  y1$samples$group=Samples

  et=exactTest(y1)

  write.csv(et,file = "../out_put/edgeR3差异分析结果.csv")

  print("本次差异基因分析结果为")

  o <- order(et$table$PValue)

  print(summary(decideTests(et,p.value = pValue)))

  plotMD(et)

  abline(h=c(-1, 1), col="blue")

  write.csv(et,file = "../out_put/edgeR2差异分析结果.csv")

  TRUE

}















#' @functionName main

#' @author tanzicai

#' @purpose main function to start the analyze programe

#' @createtime 2021.9.17

#' @time 2021.9.20

#' @input work space path

#' @return void

#' @param start munzip
#' @param start combin
#' @param start convertName
#' @param start edgRAnalyze
#' @param start edgRAnalyze2

main <- function(path , p = NULL){

  if(!is.null(p))
  {

    pValue = p

  }else{

    pValue = 0.05

  }

  setwd(path)

  print("start to analyze the defferent genes in miRNA")

  print("step 1 : move all files to a folder......")

  if(munzip(path) == TRUE )print("successful")

  print("step 2 ：read all files and merge a table.......")

  if(combin(path) == TRUE )print("successful")

  print("step 3 : Comment file extraction......")

  if(convertName(path) == TRUE )print("successful")

  DATA = read.csv("../temp/重命名_列为TCGA编号.csv")

  DATA = DATA[,-1]

  print("step 4 : use first way of edgR carries on difference analysis.....")

  if(edgRAnalyze(DATA,pValue) == TRUE )print("successful")

  print("step 5 : use second way of edgR carries on difference analysis.....")

  if(edgRAnalyze2(DATA,pValue) == TRUE )print("successful")

  print("step 6 : use third way of edgR carries on difference analysis.....")

  if(edgRAnalyze3(DATA,pValue) == TRUE )print("successful")

}





#' @functionName draw_pca

#' @author tanzicai

#' @purpose draw pca and save it

#' @createtime 2021.9.17

#' @time 2021.9.20

#' @input dat,grouplist,save_path

#' @return void

#' @param draw pca


draw_pca<-function(dat,group_list,path){
  ## 下面是画PCA的必须操作，需要看说明书。

  dat=t(dat)#画PCA图时要求是行名时样本名，列名时探针名，因此此时需要转换

  dat=as.data.frame(dat)#将matrix转换为data.frame

  dat=cbind(dat,group_list) #cbind横向追加，即将分组信息追加到最后一列

  library("FactoMineR")#画主成分分析图需要加载这两个包

  library("factoextra")

  dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)#现在dat最后一列是group_list，需要重新赋值给一个dat.pca,这个矩阵是不含有分组信息的

  fviz_pca_ind(dat.pca,

               geom.ind = "point", # show points only (nbut not "text")

               col.ind = dat$group_list, # color by groups

               #palette = c("#00AFBB", "#E7B800"),

               addEllipses = TRUE, # Concentration ellipses

               legend.title = "Groups"
  )

  ggsave(path)

}








#' @functionName draw_heatmap

#' @author tanzicai

#' @purpose draw heat map and save it

#' @createtime 2021.9.17

#' @time 2021.9.20

#' @input dat,grouplist,save_path

#' @return void

#' @param draw heat map

#' @param save heat map



draw_heatmap<-function(dat,group_list,path){

  cg=names(tail(sort(apply(dat,1,sd)),1000))#apply按行（'1'是按行取，'2'是按列取）取每一行的方差，从小到大排序，取最大的1000个

  library(pheatmap)

  pheatmap(dat[cg,],show_colnames =F,show_rownames = F) #对那些提取出来的1000个基因所在的每一行取出，组合起来为一个新的表达矩阵

  n=t(scale(t(dat[cg,]))) # 'scale'可以对log-ratio数值进行归一化

  n[n>2]=2

  n[n< -2]= -2

  n[1:4,1:4]

  pheatmap(n,show_colnames =F,show_rownames = F)

  ac=data.frame(g=group_list)

  rownames(ac)=colnames(n) #把ac的行名给到n的列名，即对每一个探针标记上分组信息

  pheatmap(n,show_colnames =F,show_rownames = F,
           annotation_col=ac,filename = path)


}







#' @functionName envirment

#' @author tanzicai

#' @purpose download the envirment of miRNA anakysis

#' @createtime 2021.9.17

#' @time 2021.9.20

#' @input dat,grouplist,save_path

#' @return void

#' @param download 正在加载：stringr

#' @param download 正在加载：patchwork

#' @param download 正在加载：DESeq2

#' @param download 正在加载：edgeR

#' @param download 正在加载：stringr

#' @param download TCGAbiolinks


envirment <- function(){

  print("正在安装所需要的包")

  print("正在加载：stringr")
  if(!require(stringr))
    BiocManager::install('stringr')

  print("正在加载：patchwork")
  if(!require(patchwork))
    BiocManager::install('patchwork')
  print("正在加载：cowplot")

  if(!require(cowplot))
    BiocManager::install('cowplot')
  print("正在加载：DESeq2")

  if(!require(DESeq2))
    BiocManager::install('DESeq2')
  print("正在加载：edgeR")

  if(!require(edgeR))
    BiocManager::install('edgeR')
  print("正在加载：limma")

  if(!require(limma))
    BiocManager::install('limma')

  if(!require(TCGAbiolinks))
    BiocManager::install('TCGAbiolinks')
  library(TCGAbiolinks)
}


#' @functionName autoAnalyssis

#' @author tanzicai

#' @purpose auto analyze miRNA difference

#' @createtime 2021.9.17

#' @time 2021.9.20

#' @input NULL / cancer type

#' @return void

#' @param download the cancer's miRNA data

#' @param conbine the cancer's miRNA data

#' @param use DESeq2 to analyze cancer's miRNA data

#' @param use edgeR to analyze cancer's miRNA data

autoAnalyssisSingleCancer <-function(input_type = NULL){

  print("您正在使用自动化分析函数,请确认你的癌症有癌症组和正常组!!!")

  print("环境监测中:")

  envirment()

  if(!is.null(input_type)){

    print("支持的癌症种类的缩写如下，请选择下载的癌症类型")

    print(TCGAbiolinks:::getGDCprojects()$project_id)

    cancer_type = readline()
  }else{

    cancer_type = input_type

  }

  print("正在请求数据:")

  clinical <- GDCquery_clinic(project = cancer_type, type = "clinical")

  query <- GDCquery(project = cancer_type,

                    data.category = "Transcriptome Profiling",

                    data.type = "miRNA Expression Quantification",

                    workflow.type = "BCGSC miRNA Profiling",

                    data.format = "txt")

  print("正在下载数据")

  GDCdownload(query, method = "api", files.per.chunk = 50)

  print("正在合并数据")

  expdat <- GDCprepare(query = query)


  library(tibble)

  print("重命名行名")

  rownames(expdat) <- NULL

  expdat <- column_to_rownames(expdat,var = "miRNA_ID")

  expdat[1:3,1:3]

  print("行列置换:耗时,请耐心等待")

  exp = t(expdat[,seq(1,ncol(expdat),3)])

  print("部分结果:")

  print(exp[1:4,1:4])

  print("行名分割替换")

  rowName <- str_split(rownames(exp),'_',simplify = T)[,3]

  expr<- apply(expr,2,as.numeric)

  expr<- na.omit(expr)

  print("处理结果:")

  print(dim(expr))

  print("处理低/零表达基因")

  expr <- expr[,apply(expr, 2,function(x){sum(x>1)>10})]

  rownames(expr) <- rowName

  print(dim(expr))

  print("处理结束,部分结果:")

  print(expr[1:4,1:4])

  print("保存数据")

  save(expr,clinical,file = paste0("tcga-",cancer_type,"-download.Rdata"))

  print("处理附带信息:")

  meta <- clinical

  colnames(meta)

  meta <- meta[,c("submitter_id","vital_status",

                  "days_to_death","days_to_last_follow_up",

                  "race",

                  "age_at_diagnosis",

                  "gender" ,

                  "ajcc_pathologic_stage")]

  print("转置样本:")

  expr=t(expr)

  expr[1:4,1:4]

  #分出tumor组和normal组

  print("分组:")


  group_list <- ifelse(as.numeric(str_sub(colnames(expr),14,15))<10,"tumor","normal")

  group_list <- factor(group_list,levels = c("normal","tumor"))


  print("分组信息如下：")

  table(group_list)

  save(expr,group_list,file = paste0("tcga-",cancer_type,"-raw.Rdata"))

  print("使用DESeq2分析中:")

  library(DESeq2)

  colData <- data.frame(row.names =colnames(expr), condition=group_list)


  write.csv(colData,"./分组信息.csv")

  dds <- DESeqDataSetFromMatrix(
    countData = expr,

    colData = colData,

    design = ~ condition)

  dds <- DESeq(dds)

  # 两两比较/

  res <- results(dds, contrast = c("condition",rev(levels(group_list))))

  resOrdered <- res[order(res$pvalue),] # 按照P值排序

  DEG <- as.data.frame(resOrdered)

  head(DEG)

  DEG <- na.omit(DEG)

  #logFC_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)) )
  logFC_cutoff <- 1

  DEG$change = as.factor(ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,ifelse(DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT'))

  head(DEG)

  DESeq2_DEG <- DEG

  print("使用DESeq2分析结果如下")

  print(table(DESeq2_DEG$change))

  write.csv(DESeq2_DEG,"./DESeq2_DEG分析结果.csv")

  print("使用edgeR方法中:")

  library(edgeR)

  dge <- DGEList(counts=expr,group=group_list)

  dge$samples$lib.size <- colSums(dge$counts)

  dge <- calcNormFactors(dge)

  design <- model.matrix(~0+group_list)

  rownames(design)<-colnames(dge)

  colnames(design)<-levels(group_list)

  dge <- estimateGLMCommonDisp(dge,design)

  dge <- estimateGLMTrendedDisp(dge, design)

  dge <- estimateGLMTagwiseDisp(dge, design)

  fit <- glmFit(dge, design)

  fit2 <- glmLRT(fit, contrast=c(-1,1))

  DEG=topTags(fit2, n=nrow(expr))

  DEG=as.data.frame(DEG)

  logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )

  logFC_cutoff <- 1

  DEG$change = as.factor(

    ifelse(DEG$PValue < 0.05 & abs(DEG$logFC) > logFC_cutoff,

           ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')

  )

  head(DEG)

  print("使用edgeR分析结果：")

  print(table(DEG$change))

  edgeR_DEG <- DEG

  write.csv(edgeR_DEG,"./edgeR_DEG分析结果.csv")

  print("使用limma分析中:")

  library(limma)

  design <- model.matrix(~0+group_list)

  colnames(design)=levels(group_list)

  rownames(design)=colnames(expr)


  dge <- DGEList(counts=expr)

  dge <- calcNormFactors(dge)

  logCPM <- cpm(dge, log=TRUE, prior.count=3)


  v <- voom(dge,design, normalize="quantile")

  fit <- lmFit(v, design)


  constrasts = paste(rev(levels(group_list)),collapse = "-")

  cont.matrix <- makeContrasts(contrasts=constrasts,levels = design)

  fit2=contrasts.fit(fit,cont.matrix)

  fit2=eBayes(fit2)

  DEG = topTable(fit2, coef=constrasts, n=Inf)

  DEG = na.omit(DEG)
  #logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )

  logFC_cutoff <- 1

  DEG$change = as.factor(

    ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,

           ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
  )


  head(DEG)

  limma_voom_DEG <- DEG

  write.csv(limma_voom_DEG,"./limma_voom_DEG分析结果.csv")

  print("使用limma分析：")

  print(table(limma_voom_DEG$change))

  save(DESeq2_DEG,edgeR_DEG,limma_voom_DEG,group_list,file = "DEG.Rdata")

  print("计算画图数据中:")

  cg1 = rownames(DESeq2_DEG)[DESeq2_DEG$change !="NOT"]

  cg2 = rownames(edgeR_DEG)[edgeR_DEG$change !="NOT"]

  cg3 = rownames(limma_voom_DEG)[limma_voom_DEG$change !="NOT"]

  logFC_cutoff <- 1

  print("画图中:")

  dat = log(expr+1)

  print("PCA图(主成分分析图)绘制中:")

  p1 = draw_pca(dat,group_list,"./pca.png")

  print("heatmap(热图)绘制中:")

  h1 = draw_heatmap(expr[cg1,],group_list,"./heatmap_DESeq2.png")

  h2 = draw_heatmap(expr[cg2,],group_list,"./heatmap_edgeR.png")

  h3 = draw_heatmap(expr[cg3,],group_list,"./heatmap_limma.png")

  print("volcano(火山图)数据导出中:")

  data = DESeq2_DEG[,c(2,5,7)]

  write.csv(data,"./DESeq2_DEG_volcano.csv")

  data = edgeR_DEG[,c(1,4,6)]

  write.csv(data,"./edgeR_DEG_volcano.csv")

  data = limma_voom_DEG[,c(1,4,7)]

  write.csv(data,"./limma_DEG_volcano.csv")

}



#' @functionName autoAnalyssisByDESeq2AndedgeR

#' @author tanzicai

#' @purpose auto analyze miRNA difference by DESeq2 And edgeR

#' @createtime 2021.9.17

#' @time 2021.9.20

#' @input  cancer type

#' @return void

#' @param download the cancer's miRNA data

#' @param conbine the cancer's miRNA data

#' @param use DESeq2 to analyze cancer's miRNA data

#' @param use edgeR to analyze cancer's miRNA data

autoAnalyssisByDESeq2AndedgeR<-function(input_type){

  print("开始分析:")

  cancer_type = input_type

  clinical <- GDCquery_clinic(project = cancer_type, type = "clinical")

  print("正在请求数据")

  query <- GDCquery(project = cancer_type,

                    data.category = "Transcriptome Profiling",

                    data.type = "miRNA Expression Quantification",

                    workflow.type = "BCGSC miRNA Profiling",

                    data.format = "txt")

  print("请求数据成功")

  print("正在下载数据")

  GDCdownload(query, method = "api", files.per.chunk = 50)

  print("下载数据成功")

  print("正在合并数据")

  expdat <- GDCprepare(query = query)


  library(tibble)

  print("重命名行名")

  rownames(expdat) <- NULL

  expdat <- column_to_rownames(expdat,var = "miRNA_ID")

  expdat[1:3,1:3]

  print("行列置换")

  exp = t(expdat[,seq(1,ncol(expdat),3)])

  exp[1:4,1:4]

  expr=exp

  print("行名分割替换")

  rowName <- str_split(rownames(exp),'_',simplify = T)[,3]

  expr<- apply(expr,2,as.numeric)

  expr<- na.omit(expr)

  dim(expr)

  expr <- expr[,apply(expr, 2,function(x){sum(x>1)>10})]

  rownames(expr) <- rowName

  dim(expr)

  expr[1:4,1:4]

  print("合并数据成功")

  print("正在保存数据")

  #保存数据

  save(expr,clinical,file = paste0("tcga-",cancer_type,"-download.Rdata"))

  meta <- clinical

  colnames(meta)

  meta <- meta[,c("submitter_id","vital_status",

                  "days_to_death","days_to_last_follow_up",

                  "race",

                  "age_at_diagnosis",

                  "gender" ,

                  "ajcc_pathologic_stage")]


  print("转置样本")

  expr=t(expr)

  expr[1:4,1:4]

  print("正在分组:")

  #分出tumor组和normal组

  group_list <- ifelse(as.numeric(str_sub(colnames(expr),14,15))<10,"tumor","normal")

  group_list <- factor(group_list,levels = c("normal","tumor"))


  print("分组信息如下：")

  table(group_list)

  save(expr,group_list,file = paste0("tcga-",cancer_type,"-raw.Rdata"))

  print("使用DESeq2分析中")

  library(DESeq2)

  colData <- data.frame(row.names =colnames(expr), condition=group_list)


  write.csv(colData,"./分组信息.csv")

  dds <- DESeqDataSetFromMatrix(
    countData = expr,

    colData = colData,

    design = ~ condition)

  dds <- DESeq(dds)

  # 两两比较/

  res <- results(dds, contrast = c("condition",rev(levels(group_list))))

  resOrdered <- res[order(res$pvalue),] # 按照P值排序

  DEG <- as.data.frame(resOrdered)

  head(DEG)

  DEG <- na.omit(DEG)

  #logFC_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)) )
  logFC_cutoff <- 1

  DEG$change = as.factor(ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,ifelse(DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT'))

  head(DEG)

  DESeq2_DEG <- DEG

  print("使用DESeq2分析结果如下")

  print(table(DESeq2_DEG$change))

  write.csv(DESeq2_DEG,"./DESeq2_DEG分析结果.csv")

  print("使用edgeR方法分析中:")

  library(edgeR)

  dge <- DGEList(counts=expr,group=group_list)

  dge$samples$lib.size <- colSums(dge$counts)

  dge <- calcNormFactors(dge)

  design <- model.matrix(~0+group_list)

  rownames(design)<-colnames(dge)

  colnames(design)<-levels(group_list)

  dge <- estimateGLMCommonDisp(dge,design)

  dge <- estimateGLMTrendedDisp(dge, design)

  dge <- estimateGLMTagwiseDisp(dge, design)

  fit <- glmFit(dge, design)

  fit2 <- glmLRT(fit, contrast=c(-1,1))

  DEG=topTags(fit2, n=nrow(expr))

  DEG=as.data.frame(DEG)

  logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )

  logFC_cutoff <- 1

  DEG$change = as.factor(

    ifelse(DEG$PValue < 0.05 & abs(DEG$logFC) > logFC_cutoff,

           ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')

  )

  head(DEG)

  print("使用edgeR分析结果：")
  print(table(DEG$change))

  edgeR_DEG <- DEG

  write.csv(edgeR_DEG,"./edgeR_DEG分析结果.csv")


  imfor = c(cancer_type,dim(expr)[2],dim(expr)[1], table(group_list)[1],table(group_list)[2],table(DESeq2_DEG$change)[1],table(DESeq2_DEG$change)[2],table(DESeq2_DEG$change)[3],table(DEG$change)[1],table(DEG$change)[2],table(DEG$change)[3])

  imfor

}







#' @functionName convertName2

#' @author tanzicai

#' @purpose auto analyze miRNA difference by DESeq2 And edgeR

#' @createtime 2021.9.17

#' @time 2021.9.20

#' @input  cancer type

#' @return void

#' @param download the cancer's miRNA data

#' @param conbine the cancer's miRNA data

#' @param use DESeq2 to analyze cancer's miRNA data

#' @param use edgeR to analyze cancer's miRNA data

convertName2 <- function(cancer_type){

  load("./DATA.Rdata")

  metaMatrix.MIR <- gdcParseMetadata(project.id = cancer_type,
                                     data.type  = 'miRNAs',
                                     write.meta = FALSE)

  metaMatrix.MIR = metaMatrix.MIR[c(1,5)]

  for(i in 1:length(metaMatrix.MIR[,1])){

    metaMatrix.MIR[i,1] = substr( metaMatrix.MIR[i,1],1,9)

  }

  row.names(metaMatrix.MIR) = metaMatrix.MIR$file_name

  write.csv(metaMatrix.MIR,file = "./temp/文件名与TCGA编号对照表.csv")

  DATA = DATA[-1]

  jsonDataInfo = metaMatrix.MIR

  jsonDataInfo = jsonDataInfo[order(jsonDataInfo$file_name),]

  colnames(DATA) = jsonDataInfo$submitter_id

  save(DATA,file = "./DATA.Rdata")

  write.csv(DATA,file = "../temp/重命名_列为TCGA编号.csv")

  imfor = c(cancer_type,dim(expr)[2],dim(expr)[1], table(group_list)[1],table(group_list)[2],table(DESeq2_DEG$change)[1],table(DESeq2_DEG$change)[2],table(DESeq2_DEG$change)[3],table(DEG$change)[1],table(DEG$change)[2],table(DEG$change)[3])

  imfor
}














#' @functionName autoAnalyssisByDESeq2AndedgeRmain

#' @author tanzicai

#' @purpose auto analyze miRNA difference by DESeq2 And edgeR

#' @createtime 2021.9.17

#' @time 2021.9.20

#' @input  cancer type

#' @return void

#' @param download the cancer's miRNA data

#' @param conbine the cancer's miRNA data

#' @param use DESeq2 to analyze cancer's miRNA data

#' @param use edgeR to analyze cancer's miRNA data

runByDEseq2AndEdgeRmain <- function(input_type){

  print(paste0("开始分析:",input_type))

  cancer_type = input_type

  print("正在请求数据")

  clinical <- GDCquery_clinic(project = cancer_type, type = "clinical")

  query <- GDCquery(project = cancer_type,

                    data.category = "Transcriptome Profiling",

                    data.type = "miRNA Expression Quantification",

                    workflow.type = "BCGSC miRNA Profiling",

                    data.format = "txt")

  print("请求数据成功")

  print("正在下载数据")

  GDCdownload(query, method = "api", files.per.chunk = 50)

  print("下载数据成功")

  print("正在合并数据")

  expdat <- GDCprepare(query = query)

  expdat = combin_TCGA(cancer_type)

  library(tibble)

  #重命名行名

  rownames(expdat) <- NULL

  expdat <- column_to_rownames(expdat,var = "miRNA_ID")

  #行列置换

  exp = DATA

  exp = t(expdat[,seq(1,ncol(expdat),3)])

  expr=exp

  #行名分割替换

  rowName <- str_split(rownames(exp),'_',simplify = T)[,3]

  expr<- apply(expr,2,as.numeric)

  expr<- na.omit(expr)

  dim(expr)

  print("正在去除0表达基因")

  expr <- expr[,apply(expr, 2,function(x){sum(x>1)>10})]

  rownames(expr) <- rowName

  dim(expr)

  print("去除0表达基因成功")

  print("合并数据成功")

  print("正在保存数据")

  保存数据

  save(expr,clinical,file = paste0("tcga-",cancer_type,"-download.Rdata"))

  meta <- clinical

  colnames(meta)

  meta <- meta[,c("submitter_id","vital_status",

                  "days_to_death","days_to_last_follow_up",

                  "race",

                  "age_at_diagnosis",

                  "gender" ,

                  "ajcc_pathologic_stage")]

  #转置样本
  expr=t(expr)

  #分出tumor组和normal组

  group_list <- ifelse(as.numeric(str_sub(colnames(expr),14,15))<10,"tumor","normal")

  group_list <- factor(group_list,levels = c("normal","tumor"))


  print("分组信息如下：")

  table(group_list)

  save(expr,group_list,file = paste0("tcga-",cancer_type,"-raw.Rdata"))

  # 使用DESeq2

  library(DESeq2)

  colData <- data.frame(row.names =colnames(expr), condition=group_list)


  write.csv(colData,"./分组信息.csv")

  dds <- DESeqDataSetFromMatrix(
    countData = expr,

    colData = colData,

    design = ~ condition)

  dds <- DESeq(dds)

  # 两两比较/

  res <- results(dds, contrast = c("condition",rev(levels(group_list))))

  resOrdered <- res[order(res$pvalue),] # 按照P值排序

  DEG <- as.data.frame(resOrdered)

  head(DEG)

  DEG <- na.omit(DEG)

  #logFC_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)) )
  logFC_cutoff <- 1

  DEG$change = as.factor(ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,ifelse(DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT'))

  head(DEG)

  DESeq2_DEG <- DEG

  print("使用DESeq2分析结果如下")

  print(table(DESeq2_DEG$change))

  write.csv(DESeq2_DEG,"./DESeq2_DEG分析结果.csv")

  # 使用edgeR方法


  dge <- DGEList(counts=expr,group=group_list)

  dge$samples$lib.size <- colSums(dge$counts)

  dge <- calcNormFactors(dge)

  design <- model.matrix(~0+group_list)

  rownames(design)<-colnames(dge)

  colnames(design)<-levels(group_list)

  dge <- estimateGLMCommonDisp(dge,design)

  dge <- estimateGLMTrendedDisp(dge, design)

  dge <- estimateGLMTagwiseDisp(dge, design)

  fit <- glmFit(dge, design)

  fit2 <- glmLRT(fit, contrast=c(-1,1))

  DEG=topTags(fit2, n=nrow(expr))

  DEG=as.data.frame(DEG)

  logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )

  logFC_cutoff <- 1

  DEG$change = as.factor(

    ifelse(DEG$PValue < 0.05 & abs(DEG$logFC) > logFC_cutoff,

           ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')

  )

  head(DEG)

  print("使用edgeR分析结果：")

  print(table(DEG$change))

  edgeR_DEG <- DEG

  write.csv(edgeR_DEG,"./edgeR_DEG分析结果.csv")

  imfor = c(cancer_type,dim(expr)[2],dim(expr)[1], table(group_list)[1],table(group_list)[2],table(DESeq2_DEG$change)[1],table(DESeq2_DEG$change)[2],table(DESeq2_DEG$change)[3],table(DEG$change)[1],table(DEG$change)[2],table(DEG$change)[3])

  imfor
}



#分析信息返回

#项目名称 cancer_type
# 样本总数 dim(expr)[2]
# 样本基因总数 dim(expr)[1]
#normal table(group_list)[1]
#tumor table(group_list)[2]
#分析结果DESeq2 DOWN table(DESeq2_DEG$change)[1]
#分析结果DESeq2 NOT table(DESeq2_DEG$change)[2]
#分析结果DESeq2 UP table(DESeq2_DEG$change)[3]
#分析结果edgeR DOWN table(DEG$change)[1]
#分析结果edgeR NOT table(DEG$change)[2]
#分析结果edgeR UP table(DEG$change)[3]









#' @functionName combin_TCGA

#' @author tanzicai

#' @purpose auto analyze miRNA difference by DESeq2 And edgeR

#' @createtime 2021.9.17

#' @time 2021.9.20

#' @input  cancer type

#' @return void

#' @param download the cancer's miRNA data

#' @param conbine the cancer's miRNA data

#' @param use DESeq2 to analyze cancer's miRNA data

#' @param use edgeR to analyze cancer's miRNA data



combin_TCGA <-function(cancer_type){

  path = paste0("./GDCdata/",cancer_type,"/harmonized/Transcriptome_Profiling/miRNA_Expression_Quantification")

  print("start to analyze the defferent genes in miRNA")

  print("step 1 : move all files to a folder......")

  if(munzip(path) == TRUE )print("successful")

  print("step 2 ：read all files and merge a table.......")

  if(combin(path) == TRUE )print("successful")

  load("./DATA.Rdata")

  DATA = convertName2(cancer_type)
}







#' @functionName autoAnalyssisByDESeq2AndedgeR

#' @author tanzicai

#' @purpose auto analyze miRNA difference by DESeq2 And edgeR

#' @createtime 2021.9.17

#' @time 2021.9.20

#' @input  cancer type

#' @return void

#' @param download the cancer's miRNA data

#' @param conbine the cancer's miRNA data

#' @param use DESeq2 to analyze cancer's miRNA data

#' @param use edgeR to analyze cancer's miRNA data



aotoAllCancerAnalysis <- function(){

  cancer_list = TCGAbiolinks:::getGDCprojects()$project_id

  table = gregexpr("TCGA",cancer_list)

  needAnalyze = NULL

  for(i in 1:length(cancer_list)){

    if(table[[i]][1] == 1){

      needAnalyze = c(needAnalyze,cancer_list[i])

    }

  }

  save(needAnalyze,file = "./cancer_infomation.Rdata")

  cancer_infomation = data.frame(row.names = c("项目名称","样本总数","样本非0基因总数","normal","tumor","DESeq2_DOWN","DESeq2_NOT","DESeq2_UP","edgeR_DOWN","edgeR_NOT","edgeR_UP"))

  for(i in 1:length(needAnalyze)){

    load("./cancer_infomation.Rdata")

    input_type =  needAnalyze[i]

    dir.create(input_type)

    setwd(paste0("./",input_type))

    imfo = try{
      runByDEseq2AndEdgeRmain(input_type)
    }

    setwd("../")

    cancer_infomation = cbind(cancer_infomation,imfo)

    save(cancer_infomation,needAnalyze,"./cancer_infomation.Rdata")

    rm(list = ls())

  }

}






#' @functionName autoAnalyssisByDESeq2AndedgeRTest

#' @author tanzicai

#' @purpose auto analyze miRNA difference by DESeq2 And edgeR

#' @createtime 2021.9.17

#' @time 2021.9.20

#' @input  cancer type

#' @return void

#' @param download the cancer's miRNA data

#' @param conbine the cancer's miRNA data

#' @param use DESeq2 to analyze cancer's miRNA data

#' @param use edgeR to analyze cancer's miRNA data

aotoAllCancerAnalysisTest <- function(){
  # setwd("/home/TCGA_Analysis")

  cancer_infomation = data.frame(row.names = c("项目名称","样本总数","样本非0基因总数","normal","tumor","DESeq2_DOWN","DESeq2_NOT","DESeq2_UP","edgeR_DOWN","edgeR_NOT","edgeR_UP"))

  cancer = c("TCGA-KICH","TCGA-KICH")

  for(i in 1:length(cancer)){

    cancer = c("TCGA-KICH","TCGA-KICH")

    input_type =  cancer[i]

    dir.create(input_type)

    setwd(paste0("./",input_type))

    imfo = runByDEseq2AndEdgeR(input_type)

    setwd("../")

    if(i == 2)load("./cancer_infomation.Rdata")

    cancer_infomation = cbind(cancer_infomation,imfo)

    save(cancer_infomation,file = "./cancer_infomation.Rdata")

    rm(list = ls())

  }

}
## TCGA-HNSC
## TCGA-KICH
