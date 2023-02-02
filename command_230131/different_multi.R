#writer:slm
.libPaths(c("/home/shenlm/R/x86_64-pc-linux-gnu-library/4.1","/opt/R/4.1.2/lib/R/library"))
#function
source("universal_function.R")
library(getopt)
library(ggplot2)
library(ggbiplot)
spec = matrix(c('matrix','m',1,"character",
                'fpkm','f',1,"character",
                'count','c',1,"character",
                'outdir','o',1,"character",
                'species','s',1,"character"),
              byrow = TRUE,ncol = 4)
opt=getopt(spec)

meta <- read.csv(opt$matrix,header = T)
fpkm <- read.csv(opt$fpkm,header = T)
count <- read.csv(opt$count,header = T)
out_dir <- opt$outdir
species <- opt$species

#--------matrix--------
count_ens <- de_1(count,1)
count_sys <- de_1_af(count,1)
meta$pare <- "1"
for (i in 1:nrow(meta)){
  now_sam <- meta[i,1]
  fir_sam <- meta[1,1]
  if (now_sam != fir_sam) {
    meta[i,3] <- "2"
  }
}

coln <- as.data.frame(colnames(count_ens[,-1]))
colnames(coln) <- "Group"
all <- merge(coln,meta,by="Group")
group <- all$pare
group2 <- all$Sample
group3 <- unique(group2)
p_name <- paste0(group3[1],"-",group3[2])

#---------different_analysis---------
result_ens <- deseq2_multi_compare(count_ens,group,0.01,2)
result_ens$all <- na.omit(result_ens$all)
val_ens <- plot_vac(result_ens$all,0.01,2)
ggsave(paste0(out_dir,p_name,"_ens_volcano.png"),val_ens,width = 10,height = 10)
write.csv(result_ens$all,paste0(out_dir,p_name,"_ens_diff.csv"))

result_sym <- deseq2_multi_compare(count_sys,group,0.01,2)
result_sym$all <- na.omit(result_sym$all)
val_sym <- plot_vac(result_sym$all,0.01,2)
ggsave(paste0(out_dir,p_name,"_sym_volcano.png"),val_sym,width = 10,height = 10)
write.csv(result_sym$all,paste0(out_dir,p_name,"_sym_diff.csv"))

#----------MA_plot---------------
library(ggpubr)
ma_ens <- plot_MA(result_ens$all,0.05,2)
ggsave(paste0(out_dir,p_name,"_ens_MA.png"),ma_ens,width = 10,height = 10,bg="white")
ma_sym <- plot_MA(result_sym$all,0.05,2)
ggsave(paste0(out_dir,p_name,"_sym_MA.png"),ma_sym,width = 10,height = 10,bg="white")

#----------PCA_plot-------------
tpm <- fpkm2tpm(fpkm)
pca <- tpm2pca_plot_group(tpm,group2)
ggsave(paste0(out_dir,p_name,"_pca.png"),pca,width = 10,height = 10,bg="white")

#---------GO--------------------
list_change <- as.data.frame(rownames(result_ens$change))
colnames(list_change) <- c("ens")
list_change <- de_point(list_change,1)
list_change2 <- genesymol_all(species,"ENSEMBL",list_change)
go_change <- enrich_go(species,list_change2$ENTREZID)
df_go_change <- as.data.frame(go_change)
write.csv(df_go_change,paste0(out_dir,p_name,"_change_go.csv"))
plot_go(go_change,out_dir,paste0(p_name,"_change_go"))

list_up <- as.data.frame(rownames(result_ens$up))
colnames(list_up) <- c("ens")
list_up <- de_point(list_up,1)
list_up2 <- genesymol_all(species,"ENSEMBL",list_up)
go_up <- enrich_go(species,list_up2$ENTREZID)
df_go_up <- as.data.frame(go_up)
write.csv(df_go_up,paste0(out_dir,p_name,"_up_go.csv"))
plot_go(go_up,out_dir,paste0(p_name,"_up_go"))

list_down <- as.data.frame(rownames(result_ens$down))
colnames(list_down) <- c("ens")
list_down <- de_point(list_down,1)
list_down2 <- genesymol_all(species,"ENSEMBL",list_down)
go_down <- enrich_go(species,list_down2$ENTREZID)
df_go_down <- as.data.frame(go_down)
write.csv(df_go_down,paste0(out_dir,p_name,"_down_go.csv"))
plot_go(go_down,out_dir,paste0(p_name,"_down_go"))


#-----KEGG_analysis-------------
kegg_change <- enrich_kegg(species,list_change2$ENTREZID)
kegg_up <- enrich_kegg(species,list_up2$ENTREZID)
kegg_down <- enrich_kegg(species,list_down2$ENTREZID)
kn_change <- nrow(as.data.frame(kegg_change))
kn_up<- nrow(as.data.frame(kegg_up))
kn_down<- nrow(as.data.frame(kegg_down))
if (kn_change>0) {
  plot_kegg(kegg_change,out_dir,paste0(p_name,"_kegg_change"))
  write.csv(as.data.frame(kegg_change),paste0(out_dir,p_name,"_kegg_change.csv"))
}
if (kn_up>0) {
  plot_kegg(kegg_up,out_dir,paste0(p_name,"_kegg_up"))
  write.csv(as.data.frame(kegg_up),paste0(out_dir,p_name,"_kegg_up.csv"))
}
if (kn_down>0) {
  plot_kegg(kegg_down,out_dir,paste0(p_name,"_kegg_down"))
  write.csv(as.data.frame(kegg_down),paste0(out_dir,p_name,"_kegg_down.csv"))
}


