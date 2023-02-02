#writer:slm
.libPaths(c("/home/shenlm/R/x86_64-pc-linux-gnu-library/4.1","/opt/R/4.1.2/lib/R/library"))
#matrix
fpkm2tpm <- function(data1){
  #data1格式:1列为名,余为fpkm
  data1_name <- as.data.frame(data1[,1])
  colnames(data1_name) <- colnames(data1)[1]
  fpkm <- data1[,-1]
  fpkm <- exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
  data1 <- cbind(data1_name,fpkm)
  return(data1)
}
de_point <- function(data1,col){
  #1位:列表 2位:要de_point的列
  library(stringi)
  for (i in 1:nrow(data1)) {
    data1[i,col] <- str_sub(data1[i,col],1,str_locate(data1[i,col],"[.]")[1]-1)
  }
  return(data1)
}
de_1 <- function(data1,col){
  library(stringi)
  library(stringr)
  for (i in 1:nrow(data1)) {
    if (str_detect(data1[i,col],"[|]")) {
      data1[i,col] <- str_sub(data1[i,col],1,str_locate(data1[i,col],"[|]")[1]-1)
    }else if (str_detect(data1[i,col],"[|]")==F) {
      data1[i,col] <- data1[i,col]
    }
  }
  return(data1)
}
de_1_af <- function(data1,col){
  library(stringi)
  library(stringr)
  for (i in 1:nrow(data1)) {
    if (str_detect(data1[i,col],"[|]")) {
      data1[i,col] <- str_sub(data1[i,col],str_locate(data1[i,col],"[|]")[1]+1,str_length(data1[i,col]))
    }else if (str_detect(data1[i,col],"[|]")==F) {
      data1[i,col] <- data1[i,col]
    }
  }
  return(data1)
}
de_rowall_zero <- function(data1){
  #data1格式:1列为名,余为表达值
  for (i in 2:ncol(data1)) {
    data1[,i] <- as.numeric(data1[,i])
  }
  data2 <- data1[which(rowSums(data1[2:ncol(data1)])>0),]
  return(data2)
}
tss2loc <- function(species,typ,data1){
  #function:根据dom的染色体位置,将其转换为gene
  #1位:物种;2位:类型(ENSEMBL/SYMBOL);3位:转换为df的dom
  options(bedtools.path="/data1/shenluemou/biosoft/anaconda3/bin/")
  library(bedtoolsr)
  if (species=="mouse") {
    if (typ=="ENSEMBL") {
      ref <- read.csv("/data1/shenluemou/database/gene_local/R/mouse_grcm38_ens_geneloc.gtf",header = F)
    }else if (typ=="SYMBOL") {
      ref <- read.csv("/data1/shenluemou/database/gene_local/R/mouse_grcm38_geneloc.gtf",header = F)
    }
  }else if (species=="macacaf") {
    if (typ=="ENSEMBL") {
      ref <- read.csv("/data1/shenluemou/database/gene_local/R/macaca_fac_ens_geneloc.gtf",header = F)
    }else if (typ=="SYMBOL") {
      ref <- read.csv("/data1/shenluemou/database/gene_local/R/macaca_fac_sys_geneloc.gtf",header = F)
    }
  }else if (species=="human") {
    if (typ=="ENSEMBL") {
      ref <- read.csv("/data1/shenluemou/database/gene_local/R/human_hg38_ens_geneloc.gtf",header = F) 
    }else if (typ=="SYMBOL") {
      ref <- read.csv("/data1/shenluemou/database/gene_local/R/human_hg38_sys_geneloc.gtf",header = F)
    }
  }
  tmp1 <-select_(data1,"seqnames","start","end")
  tmp2 <- bt.intersect(a = tmp1,b = ref,wa = T,wb = T)
  tmp3 <- select_(tmp2,"V1","V2","V7")
  colnames(tmp3) <- c("seqnames","dominant_ctss","gene")
  return(tmp3)
}
genesymol_all <- function(species,type,data1){
  library(clusterProfiler)
  library(AnnotationDbi)
  library(dplyr)
  colnames(data1) <- "gene"
  if (species=="mouse") {
    library(org.Mm.eg.db)
    org <- org.Mm.eg.db
  }else if (species=="human") {
    library(org.Hs.eg.db)
    org <- org.Hs.eg.db
  }else if (species=="macacaf") {
    org <- loadDb("/data1/luofc/hanx/database/Enrichment/GO.db/org.Macaca_fascicularis.eg.sqlite")
  }
  if (type=="ENSEMBL") {
    tmp1 <- bitr(data1$gene,fromType = "ENSEMBL",
                 toType = c("ENTREZID","SYMBOL"),org,drop = TRUE)
  }else if (type=="ENTREZID") {
    tmp1 <- bitr(data1$gene,fromType = "ENTREZID",
                 toType = c("ENSEMBL","SYMBOL"),org,drop = TRUE)
  }else if (type=="SYMBOL") {
    if (species=="macacaf") {
      tmp1 <- bitr(data1$gene,fromType = "SYMBOL",
                   toType = c("ENTREZID"),org,drop = TRUE)
    }else if (species!="macacaf") {
      tmp1 <- bitr(data1$gene,fromType = "SYMBOL",
                   toType = c("ENSEMBL","ENTREZID"),org,drop = TRUE)     
    }
  }
  if (species=="macacaf") {
    data2 <- select_(tmp1,"SYMBOL","ENTREZID")
  }else if (species!="macacaf") {
    data2 <- select_(tmp1,"ENSEMBL","SYMBOL","ENTREZID")
  }
  return(data2)
}
enrich_go <- function(species,list1){
  library(AnnotationDbi)
  library(clusterProfiler)
  if (species=="mouse") {
    library(org.Mm.eg.db)
    org <- org.Mm.eg.db
  }else if (species=="human") {
    library(org.Hs.eg.db)
    org <- org.Hs.eg.db
  }else if (species=="macacaf") {
    org <- loadDb("/data1/luofc/hanx/database/Enrichment/GO.db/org.Macaca_fascicularis.eg.sqlite")
  }
  go_re <- enrichGO(gene=list1,
                    OrgDb = org,
                    keyType = "ENTREZID",
                    ont="ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = T)
  return(go_re)
}
enrich_kegg <- function(species,list1){
  library(clusterProfiler)
  if (species=="mouse") {
    org <- "mmu"
  }else if (species=="human") {
    org <- "hsa"
  }
  re <- enrichKEGG(list1, 
             organism=org, 
             pvalueCutoff=0.01, 
             pAdjustMethod="BH",
             keyType="kegg")
  return(re)
}
tolog10 <- function(data1){
  data1_name <- as.data.frame(data1[,1])
  colnames(data1_name) <- colnames(data1)[1]
  data1 <- data1[,-1]
  data1 <- log10(data1+1)
  data1 <- cbind(data1_name,data1)
  return(data1)
}
de_rowatleast <- function(data1,num){
  #function:选出至少有某列大于某值的行
  #data1格式:第一行为标记，其余为matrix
  library(dplyr)
  data1_name <- as.data.frame(data1[,1])
  colnames(data1_name) <- colnames(data1)[1]
  data1 <- as.data.frame(data1[,-1])
  data1_name <- data1_name[apply(data1>num, 1, any),]
  data1 <- data1[apply(data1>num, 1, any),]
  data1 <- cbind(data1_name,data1)
  return(data1)
}
de_rowall <- function(data1,num){
  #function:选出所有列大于某值的行
  #data1格式:第一行为标记，其余为matrix
  library(dplyr)
  data1_name <- as.data.frame(data1[,1])
  colnames(data1_name) <- colnames(data1)[1]
  data1 <- as.data.frame(data1[,-1])
  data1_name <- data1_name[apply(data1>num, 1, all),]
  data1 <- data1[apply(data1>num, 1, all),]
  data1 <- cbind(data1_name,data1)
  return(data1)
}
select_col <- function(data1,colnum){
  #function:选择列表中的某一或几列单独成表格
  for (i in 1:length(colnum)) {
    tmp_num <- colnum[i]
    tmp <- as.data.frame(data1[,tmp_num])
    colnames(tmp) <- colnames(data1)[tmp_num]
    if (i==1) {
      data2 <- tmp
    }else if (i>1) {
      data2 <- cbind(data2,tmp)
    }
  }
  return(data2)
}
cluster_re <- function(data,num){
  library(stats)
  #----data include : "gene_id" and sample
  colnames(data)[1] <- "gene_id"
  rownames(data) <- make.unique(data$gene_id)
  data$gene_id <- NULL
  data <- data[apply(data>1,1,any),]
  data_cl <- data %>% t() %>% scale(center = FALSE,scale = TRUE) %>%
    t() %>% hcut(method="pearson")
  x <- c(1.7, 2.7, 3.7, -3.7, -5.7,3.7, 5.7, 4.8, 10.3, -12.9, 13.8, 12.3)
  y <- c(1.2, 2.3, 3.2, -3.5, -3.2, 2.1,4.7, .8, 1.2, 11.5, 1.3, 3.2)
  plot(x, y, cex = 1, pch = 3,xlab ="x", ylab ="y",col ="black")
  data_cluster <- rect.hclust(data_cl,num)
  return(data_cluster)
}
edger_single_compare <- function(data1,pv,fc){
  #用于edger的单样本差异分析
  #data1的格式为1列为基因,23列分别为控制组和实验组的count
  library(edgeR)
  group <- 1:2
  data1[,2] <- as.double(data1[,2])
  data1[,3] <- as.double(data1[,3])
  y <- DGEList(counts = data1[,2:3],genes = data1[,1],group = group)
  keep <- rowSums(cpm(y)>1) >= 1
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  y_bcv <- y
  bcv <- 0.1
  et <- exactTest(y_bcv, dispersion = bcv ^ 2)
  gene1 <- decideTestsDGE(et, p.value = pv, lfc = fc)  
  colnames(gene1) <- "Signifi"
  results <- cbind(y$genes,y$counts,et$table,gene1)
  colnames(results)[4] <- "log2FoldChange"
  colnames(results)[6] <- "padj"
  colnames(results)[7] <- "threshold"
  results <- results[order(results$padj),]
  rownames(results) <- make.unique(results$genes)
  results <- results[,-1]
  results$threshold <- as.character(results$threshold)
  results$threshold <- ifelse(results$threshold=="1","up",ifelse(results$threshold=="-1",
                                                                 "down", "no_dif"))
  re_up <- results[which(results$threshold=="up"),]
  re_down <- results[which(results$threshold=="down"),]
  re_change <- results[which(results$threshold!="no_dif"),]
  all_res <- list(results,re_up,re_down,re_change)
  names(all_res) <- c("all","up","down","change")
  return(all_res)
}
deseq2_multi_compare <- function(df,group,pv,fc){
  #使用DESeq2进行多平行的差异分析,group为2-1
  #df的格式为1列为基因名标记,其余为表达矩阵
  library(DESeq2)
  rownames(df) <- make.unique(df[,1])
  df <- df[,-1]
  df <- df[rowMeans(df)>1,]
  colData <- data.frame(row.names = colnames(df),group)
  dds <- DESeqDataSetFromMatrix(countData = df, 
                                colData = colData, 
                                design = ~ group)
  dds1 <- DESeq(dds, fitType = 'mean', 
                minReplicatesForReplace = 7, 
                parallel = FALSE) 
  res <- results(dds1)
  res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  res1$threshold <- ifelse(res1$padj>pv,"no",ifelse(abs(res1$log2FoldChange)<2,"no",
                                                          ifelse(res1$log2FoldChange>2,"up","down")))
  res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
  res1_up<- res1[which(res1$log2FoldChange >= fc & res1$padj < pv),]
  res1_down<- res1[which(res1$log2FoldChange <= -fc & res1$padj < pv),]
  res1_change <- rbind(res1_up,res1_down)
  all_res <- list(res1,res1_up,res1_down,res1_change)
  names(all_res) <- c("all","up","down","change")
  return(all_res)
}
macacaf_sys_tran <- function(df1,type){
  db <- read.csv("/data1/shenluemou/database/gene_local/R/macacaf_ens_sys.csv",
                 header = F)  
  colnames(db) <- c("ens","sys")
  if (type=="ENSEMBL") {
      colnames(df1) <- "ens"
      data1 <- merge(df1,db,by="ens")
      colnames(data1)[2] <- "sys"
  }else if (type=="SYMBOL") {
      colnames(df1) <- "sys"
      data1 <- merge(df1,db,by="sys")
      colnames(data1)[2] <- "ens"
  }
  data2 <- select_(data1,"ens","sys")
  return(data2)
}

#plot
plot_clu2heatmap <- function(list1,raw_matrix,rown){
  colnames(raw_matrix)[1] <- "gene"
  for (i in 1:length(list1)) {
    tmp <- list1[[i]]
    colnames(tmp) <- "gene"
    tmp2 <- merge(tmp,raw_matrix,by="gene")
    if (i==1) {
      all <- tmp2
    }else if (i>1) {
      all <- rbind(all,tmp2)
    }
  }
  p <- heat_map(all,F,F,rown)
  return(p)
}
venn_plot <- function(list1,numsize){
  library(venn)
  cross <- venn(list1,
                zcolor = "style",
                opacity = 0.3,
                box = F,
                ilcs = numsize,
                sncs = 1)
  return(cross)
}
violin_plot <- function(data1,xname,yname){
  #data1格式:第一列为gene,其它列为表达matrix格式
  library(reshape2)
  colnames(data1)[1] <- "gene"
  data2 <- melt(data1)
  a <- ggplot(data2,aes(x=variable,y=value))+
    geom_violin(aes(color=variable),trim=T)+
    geom_boxplot(width=0.01)+
    labs(x=xname,y=yname)
  return(a)
}
tpm2pca_plot_group <- function(data1,group){
  library(factoextra)
  data1 <- de_rowall_zero(data1)
  rownames(data1) <- make.unique(data1[,1])
  data1 <- data1[,-1]
  data1 <- t(data1)
  data.pca <- prcomp(data1, scale = T)
  p <- fviz_pca_ind(data.pca, col.ind=group, 
                    mean.point=F,  # 去除分组的中心点
                    label = "none", # 隐藏样本标签
                    addEllipses = T, # 添加边界线
                    legend.title="Groups",
                    pointsize = 5,
                    #ellipse.type="confidence", # 绘制置信椭圆 
                    ellipse.level=0.9)
  return(p)
}
heat_map <- function(data1,clu_col,clu_row,rown){
  #heatmap，第一列为标记
  #1:列表;2.是否列聚类;3.是否行聚类;4.是否显示行名
  library(pheatmap)
  rownames(data1) <- make.unique(data1[,1])
  data1 <- data1[,-1]
  if (ncol(data1)<=2) {
    data1 <- t(data1) %>% scale(center = FALSE) %>% t()
  }else if (ncol(data1)>2) {
    data1 <- t(data1) %>% scale() %>% t()
  }
  p <- pheatmap::pheatmap(data1,cluster_rows = clu_row,
                          cluster_cols = clu_col,border=FALSE,
                          show_rownames = rown)
}
plot_vac <- function(Dat,pv,fc){
  library(ggplot2)
  library(ggrepel)
  Gene <- rownames(Dat[1:20,])
  ggplot(Dat,aes(x=log2FoldChange,y=-log10(padj),color=threshold))+
    geom_point()+
    scale_color_manual(values=c("#00008B","#808080","#DC143C"))+#确定点的颜色
    geom_text_repel(
      data = Dat[1:20,],
      aes(label = Gene))+#添加关注的点的基因名
    theme_bw()+#修改图片背景
    theme(
      legend.title = element_blank(),#不显示图例标题
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank()
      #不显示图例标题 
    )+
    ylab('-log10 (p-adj)')+#修改y轴名称
    xlab('log2 (FoldChange)')+#修改x轴名称
    geom_vline(xintercept=c(-fc,fc),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange
    geom_hline(yintercept = -log10(pv),lty=3,col="black",lwd=0.5)
}
plot_MA <- function(df,fdr,fc){
  library(ggpubr)
  p <- ggmaplot(df, fdr=fdr, fc = fc, size = 0.4,
           palette = c("#B31B21", "#1465AC", "darkgray"),
           legend = "top", top = 20,
           font.label = c("bold", 11),
           font.legend = "bold",
           font.main = "bold",
           ggtheme = ggplot2::theme_minimal())+
    theme_bw()+#修改图片背景
    theme(
      legend.title = element_blank(),#不显示图例标题
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank()
      #不显示图例标题 
    )
  return(p)
}
plot_go <- function(S4,outdir,name){
  df <- as.data.frame(S4)
  n1 <- nrow(df[which(df[,1]=="BP"),])
  n2 <- nrow(df[which(df[,1]=="CC"),])
  n3 <- nrow(df[which(df[,1]=="MF"),])
  n <- max(n1,n2,n3)
  p_dot <- dotplot(S4,split="ONTOLOGY")+facet_grid(ONTOLOGY~.,scale = "free")
  p_bar <- barplot(S4,split="ONTOLOGY")+facet_grid(ONTOLOGY~.,scale = "free")
  if (n>10) {
    ggsave(paste0(outdir,name,"_dot.png"),p_dot,width = 10,height = 15)
    ggsave(paste0(outdir,name,"_bar.png"),p_bar,width = 10,height = 15)
  }else if (n<10) {
    ggsave(paste0(outdir,name,"_dot.png"),p_dot,width = 10,height = 1.5*n)
    ggsave(paste0(outdir,name,"_bar.png"),p_bar,width = 10,height = 1.5*n)
  }
}
plot_kegg <- function(S4,outdir,name){
  df <- as.data.frame(S4)
  n <- nrow(df)
  p_dot <- dotplot(S4)
  p_bar <- barplot(S4)
  if (n >= 10) {
    ggsave(paste0(outdir,name,"_dot.png"),p_dot,width = 10,height = 10)
    ggsave(paste0(outdir,name,"_bar.png"),p_bar,width = 10,height = 10)
  }else if (n < 10) {
    ggsave(paste0(outdir,name,"_dot.png"),p_dot,width = 10,height = 1*n)
    ggsave(paste0(outdir,name,"_bar.png"),p_bar,width = 10,height = 1*n)
  }
}







