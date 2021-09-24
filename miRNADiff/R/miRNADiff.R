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

  dir.create("../temp")

  dir.create("../temp/source")

  filePath = dir(path = path ,full.names = T)

  filePath = filePath[1:length(filePath)-1]

  for(wd in filePath){

    files = dir(path = wd,pattern = "txt$")

    fromFilePath = paste0(wd,"/",files)

    toFilePath = paste0("../temp/source/",files)

    file.copy(fromFilePath,toFilePath)
  }

  unlink("../temp/source/annotations.txt")

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

  countFilePath = dir(path = "../temp/source",pattern = "*.txt")

  counts_merge = NULL

  for (i in countFilePath) {

    x=read.delim(paste0("../temp/source/",i),col.names = c("miRNA_ID",substr(i,1,9),1,2))

    x=x[,1:length(x)-1]

    x=x[,1:length(x)-1]

    if(is.null(counts_merge)){

      counts_merge = x

    }else{

      counts_merge =merge(counts_merge,x,by = "miRNA_ID")

    }

  }

  write.csv(counts_merge,file = "../temp/源文件合并.csv")

  rownames(counts_merge)=counts_merge$miRNA_ID

  DATA = counts_merge[,-1]

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

  table(group_list)

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


main2 <-function(){
  rm(list = ls())
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

  devtools::install_github('kevinblighe/EnhancedVolcano')


  if(!require(TCGAbiolinks))
    BiocManager::install('TCGAbiolinks')
  library(TCGAbiolinks)

  devtools::install_github('kevinblighe/EnhancedVolcano')
  library(EnhancedVolcano)


  print("支持的癌症种类的缩写如下，请选择下载的癌症类型")
  print(TCGAbiolinks:::getGDCprojects()$project_id)
  TCGA-KICHcancer_type = readline()

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

  #重命名行名

  rownames(expdat) <- NULL

  expdat <- column_to_rownames(expdat,var = "miRNA_ID")

  expdat[1:3,1:3]

  #行列置换

  exp = t(expdat[,seq(1,ncol(expdat),3)])

  exp[1:4,1:4]

  expr=exp

  #行名分割替换

  rowName <- str_split(rownames(exp),'_',simplify = T)[,3]

  expr<- apply(expr,2,as.numeric)

  expr<- na.omit(expr)

  dim(expr)

  expr <- expr[,apply(expr, 2,function(x){sum(x>1)>10})]

  rownames(expr) <- rowName

  dim(expr)

  expr[1:4,1:4]

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

  #转置样本
  expr=t(expr)

  expr[1:4,1:4]

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

  table(DESeq2_DEG$change)

  write.csv(DESeq2_DEG,"./DESeq2_DEG分析结果.csv")

  # 使用edgeR方法

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
  table(DEG$change)

  edgeR_DEG <- DEG

  write.csv(edgeR_DEG,"./edgeR_DEG分析结果.csv")

  #limma

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

  table(limma_voom_DEG$change)

  save(DESeq2_DEG,edgeR_DEG,limma_voom_DEG,group_list,file = "DEG.Rdata")

  cg1 = rownames(DESeq2_DEG)[DESeq2_DEG$change !="NOT"]

  cg2 = rownames(edgeR_DEG)[edgeR_DEG$change !="NOT"]

  cg3 = rownames(limma_voom_DEG)[limma_voom_DEG$change !="NOT"]

  logFC_cutoff <- 1


  dat = log(expr+1)

  p1 = draw_pca(dat,group_list,"pc.png")

  h1 = draw_heatmap(expr[cg1,],group_list,"heatmap_DESeq2.png")

  h2 = draw_heatmap(expr[cg2,],group_list,"heatmap_edgeR.png")

  h3 = draw_heatmap(expr[cg3,],group_list,"heatmap_limma.png")

 # v1 = EnhancedVolcano(DESeq2_DEG[,c(2,5,7)],
 #                      lab = rownames(DESeq2_DEG[,c(2,5,7)]),
 #                      x = 'log2FoldChange',
  #                     y = 'padj',
  #                     xlim = c(-6,7),
  #                     xlab = bquote(~Log[2]~ 'fold change'),
  #                     pCutoff = 10e-14,
   #                    FCcutoff = 2
  #)

 # v2 =  EnhancedVolcano(edgeR_DEG[,c(1,4)],
  #                      lab = rownames(edgeR_DEG[,c(1,4)]),
   #                     x = 'logFC',
  #                      y = 'PValue  ',
  #                      xlim = c(-6,7),
  #                      xlab = bquote(~Log[2]~ 'fold change'),
  #                      pCutoff = 10e-14,
  #                      FCcutoff = 2
 # )

 # v3 = EnhancedVolcano(limma_voom_DEG[,c(1,4,7)],
  #                     lab = rownames(limma_voom_DEG[,c(1,4,7)]),
  #                     x = 'logFC',
  #                     y = 'P.Value  ',
   #                    xlim = c(-6,7),
  #                     xlab = bquote(~Log[2]~ 'fold change'),
  #                     pCutoff = 10e-14,
  #                     FCcutoff = 2
  #)


 # BiocManager::install("ggplot")
 # library(ggplot)

#  van = ggplot(DESeq2_DEG[,c(2,5,7)],aes(log2FoldChange,padj))
 # van + geom_point(aes())+
 #   ylim(0,20) + xlim(-5,5) +
 #   scale_color_manual(values=c("blue","grey", "red"))+
 #   geom_vline(xintercept = c(-1, 1), lty = 2,colour="#000000")+ #增加虚线
 #   geom_hline(yintercept = c(1), lty = 2,colour="#000000")+
  #  theme(
  #    axis.text=element_text(size=20),
  #    axis.title=element_text(size=20)
  #  )
 #
}


