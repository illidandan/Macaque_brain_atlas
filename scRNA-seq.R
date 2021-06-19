res0.8



## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj,
                          feature,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = rel(1), angle = 0),
          axis.text.y = element_text(size = rel(1)),
          plot.margin = plot.margin )
  return(p)
}
## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}
## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  # change the y-axis tick to only max value
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x +
                            scale_y_continuous(breaks = c(y)) +
                            expand_limits(y = y))
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}
StackedVlnPlot(obj = mes_sub, group.by = "type2", features = rev(my_features)) + theme(axis.text.x = element_text(angle=90, hjust=1))


my_features <- c("PKHD1L1","ITLN1", "MSLN", 
                 "IL7R", "SKAP1", "PTPRC", "CD3D", "TRAC", "NKG7")
my_features<-c("PTPRC","IKZF1","SKAP1","TRAC","THEMIS","IL7R","EFNA5","ITLN1","MSLN")

StackedVlnPlot(obj = mes_sub, group.by = "type2", features = rev(my_features)) + theme(axis.text.x = element_text(angle=90, hjust=1))


mes_sub<-subset(mes,subset=organ=="Visceral_adipose")
mes_sub<-subset(mes_sub,idents = c("0","1","2","5","10","11"))
mes_sub$type2<-"Meso"
mes_sub$type2[which(mes_sub$seurat_clusters=="5")]<-"Immune"
mes_sub$type2[which(mes_sub$seurat_clusters=="10")]<-"Immune"
mes_sub$type2[which(mes_sub$seurat_clusters=="11")]<-"Meso"



SK$mutype<-"Type I"
SK@meta.data$mutype[which(SK@meta.data$seurat_clusters=="0")]<-"Type IIa"
SK@meta.data$mutype[which(SK@meta.data$seurat_clusters=="11")]<-"MTJ"
SK@meta.data$mutype[which(SK@meta.data$seurat_clusters=="13")]<-"NMJ"
SK@meta.data$mutype[which(SK@meta.data$seurat_clusters=="10")]<-"FAP"
SK@meta.data$mutype[which(SK@meta.data$seurat_clusters=="12")]<-"Type IIa"

SK@meta.data$mutype[which(SK@meta.data$seurat_clusters=="3")]<-"Type IIa"
SK@meta.data$mutype[which(SK@meta.data$seurat_clusters=="4")]<-"Type IIa"
SK@meta.data$mutype[which(SK@meta.data$seurat_clusters=="1")]<-"Type IIb"
SK@meta.data$mutype[which(SK@meta.data$seurat_clusters=="2")]<-"Type IIb"
SK@meta.data$mutype[which(SK@meta.data$seurat_clusters=="5")]<-"Type IIa"





gene_FAPs<-c("ARHGAP24","EGFR","EGF2","LVRN")
gene_Myoten_junc<-c("CD22A1","SLC24A2","LAMA2","NAV3","ERBB4")
gene_neuron_junc<-c("NCAM1")


gene<-c("MYH7","MYH2","GPD2","LVRN","NAV3")
setwd("D:\\项目\\猴脑\\body atlas\\Figure\\Fig2")
SK<-readRDS("New.SK.mutype.RDS")
SK<-readRDS("Seurat_commoncellSK_new_0202.integrated.Find.RDS")
libgene<-c("MYH7","MYH2","GPD2","LVRN","NAV3")
FeaturePlot(SK,features = gene)
FeaturePlot(SK,features = gene,cols=inferno(n=10),ncol=3)

Idents(SK)<-SK$mutype
DEG<-FindAllMarkers(SK,logfc.threshold = 0.5,only.pos = T,min.pct = 0.25)
library(dplyr)
library(viridis)
library(ggplot2)
DotPlot(SK,features = unique(gene_top15))+scale_color_viridis(direction=-1)+theme(axis.text.x = element_text(size=8,angle=90, hjust=1, vjust=1))

ggmeta<-SK@meta.data
p1<-ggplot(ggmeta,aes(x=organ,fill=mutype))+geom_bar(stat="count",position="fill")+ theme_bw()

SK$label<-paste0(SK$mutype,"_",SK$organ)
Idents(SK)<-SK$label
avg.all.cells<- AverageExpression(SK,slot="data",assays = "RNA")
tmp<-as.data.frame(avg.all.cells$RNA)
#co<-cor(tmp,method="spearman")
#cop<-cor(tmp,method="pearson")
DEG_label<-FindAllMarkers(SK,min.pct = 0.5,logfc.threshold=0.69315,only.pos = T)
pos<-which(rownames(tmp) %in% DEG_label$gene)
tmp_cut<-tmp[pos,]
cop_cut<-cor(tmp_cut,method="pearson")
co_cut<-cor(tmp_cut,method="spearman")
pheatmap(co_cut,show_colnames = T,show_rownames = T,fontsize = 20,
         cluster_rows = T,cluster_cols = T,scale="none")
pheatmap(cop_cut,show_colnames = T,show_rownames = T,fontsize = 20,
         cluster_rows = T,cluster_cols = T,scale="none")






co_cut_cut<-co_cut[c(-5,-11,-7,-12,-16),c(-5,-11,-7,-12,-16)]






library(org.Hs.eg.db)
library(clusterProfiler)
  filter_gene<-DEG_T1
  filter_gene<-filter_gene[filter_gene$p_val_adj<= 0.01,]
  filter_gene<-filter_gene[filter_gene$pct.1 >= 0.25,]
  filter_gene$FC<-exp(filter_gene$avg_logFC)
  filter_gene<-filter_gene[filter_gene$FC >= 2.5 ,]
  c=unique(filter_gene$cluster)
  n=length(unique(filter_gene$cluster))
  gene=list()
  j =1
  for(x in c){
    sub=subset(filter_gene,filter_gene$cluster==x)
    test = bitr(sub$Gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
    gene[[j]]=as.character(test$ENTREZID)
    j=j+1
  }
  names(gene)=c
  compGO_SV_2.5 <- compareCluster(geneCluster   = gene,
                           fun           = "enrichGO",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH",
                           OrgDb = org.Hs.eg.db,
                           ont = "BP",
                           readable = TRUE)
  pdf(paste0(gsub("_DEG.list.xls","",id_list[i]),"Go.pdf"),width=12,height=12)
  p_SV_2.5<-dotplot(compGO_SV_2.5,showCategory = 10)
  p_SV_2.5
  dev.off()
  out=compGO_SV_2.5@compareClusterResult
  write.table(out,"SKmutype_Filtered_GOList.xls",row.names = FALSE,quote = FALSE,sep = "\t")
  #write.table(filter_gene,paste0(organ,"_",id_list[i],"_Filtered_GeneList.xls"),row.names = FALSE,quote = FALSE,sep = "\t")
}



DEG<-DEG_T2a
DEG$FC<-exp(DEG$avg_logFC)
DEG<-DEG[DEG$p_val_adj<0.01,]
DEG<-DEG[DEG$avg_logFC>0.5 | DEG$avg_logFC<(-0.5),]
DEG$cluster="small"
DEG$cluster[which(DEG$avg_logFC < 0)]<-"large"
DEG$gene<-rownames(DEG)
DEG$FC<-exp(DEG$avg_logFC)
DEG<-DEG[order(DEG$FC,decreasing = T),]



library(org.Hs.eg.db)
library(clusterProfiler)
filter_gene<-DEG
#filter_gene<-filter_gene[filter_gene$p_val_adj<= 0.01,]
#filter_gene<-filter_gene[filter_gene$pct.1 >= 0.25,]
#filter_gene$FC<-exp(filter_gene$avg_logFC)
#filter_gene<-filter_gene[filter_gene$FC >= 2,]
c=unique(filter_gene$cluster)
n=length(unique(filter_gene$cluster))
gene=list()
j =1

for(x in c){
  sub=subset(filter_gene,filter_gene$cluster==x)
  test = bitr(sub$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  gene[[j]]=as.character(test$ENTREZID)
  j=j+1
}

names(gene)=c
compGO <- compareCluster(geneCluster   = gene,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",
                         readable = TRUE)
pdf(paste0(gsub("_DEG.list.xls","",id_list[i]),"Go.pdf"),width=12,height=12)
p1<-dotplot(compGO,showCategory = 10)
print(p1)
dev.off()
out=compGO@compareClusterResult
write.table(out,file="Tmes_Filtered_GOList.xls",row.names = FALSE,quote = FALSE,sep = "\t")
#write.table(filter_gene,paste0(organ,"_",id_list[i],"_Filtered_GeneList.xls"),row.names = FALSE,quote = FALSE,sep = "\t")
}



gene<-c("PKHD1L1","PTPRC","HOXD4","CLU","PAX8")
gene<-c("MYH7","MYH2","GPD2","LVRN","NAV3","ETV5")
gene<-c("LGR5","IL7R","CD3G","CD3D","THEMIS","NKG7","TRAC","CD4")
gene<-c("LVRN","NAV3","ETV5")


gene<-read.table("top10SK.marker.xls",sep="\t",header=T,row.names=1)
gene<-as.vector(unique(gene$gene))
gene<-c(gene,"LVRN","ETV5","MUSK")




gene<-c("MYH7","MYH2","GPD2","LVRN","NAV3","ETV5","DCN","COL22A1","MUSK","AMPD1","PRKAG3","BICD1","TNNT1","LPL","PHLDB2","APOD","SCN7A","MYLK3","PAX7")

for (i in c(18:19,1)){

#gene<-c(")
#pdf(paste0("order_F/",gene[i],"_SK.pdf"))
png(paste0("order_F_new/",gene[i],"_SK.png"),width=600,height=600)
p1<-FeaturePlot(SK,features = gene[i],pt.size=0.5,order=F)+
  theme_bw() + 
  theme(panel.grid =element_blank()) +  
  theme(axis.text = element_blank()) +  
  theme(axis.ticks = element_blank()) +   
  theme(panel.border = element_blank()) +
  labs(x = "", y = "", title = "")+
  NoLegend()+scale_color_viridis(direction=1)

print(p1)
dev.off()
}

top30 <- DEG %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
mes$ov_sub<-"none"
mes$ov_sub[which(mes@meta.data$seurat_clusters=="6")]<-"HOXD4+ Ovary Mesothelial Cell"
mes$ov_sub[which(mes@meta.data$seurat_clusters=="3")]<-"PAX8+ Ovary Mesothelial Cell"
mes$ov_sub[which(mes@meta.data$seurat_clusters=="4")]<-"PAX8+ Ovary Mesothelial Cell"
levels(mes)


Idents(mes)<-mes$ov_sub

DotPlot(subset(sub,subset = type !="none"),features=as.vector(top15$gene))+theme(axis.text.x = element_text(size=10,angle=90, hjust=1, vjust=1))+scale_color_viridis(direction=-1)+coord_flip()

gene_for_vln<-c("CD3D","THEMIS","TRAC","NKG7","IL7R","CD4","MSLN","ITN1","EFNA5")
VlnPlot(sub,features = top15$gene,ncol=3)



gene_for_vln1<-c("MSLN","ITN1","EFNA5","CD3D","THEMIS","TRAC","NKG7","IL7R","CD4")
gene_for_vln2<-c("CD44","CD24","LGR5","MYC","MYB","PAX8","HOXD4","LIMA1","PLCB1","EPCAM")

mes@meta.data$mutype<-"none"
mes@meta.data$mutype[which(mes@meta.data$seurat_clusters=="0")]<-"ba_Mesothelial"
mes@meta.data$mutype[which(mes@meta.data$seurat_clusters=="1")]<-"ba_Mesothelial"
mes@meta.data$mutype[which(mes@meta.data$seurat_clusters=="2")]<-"ba_Mesothelial"
mes@meta.data$mutype[which(mes@meta.data$seurat_clusters=="11")]<-"ba_Mesothelial"
mes@meta.data$mutype[which(mes@meta.data$seurat_clusters=="5")]<-"a_Immune like mesothelial"
mes@meta.data$mutype[which(mes@meta.data$seurat_clusters=="10")]<-"a_Immune like mesothelial"

mes@meta.data$organ_2<-mes@meta.data$organ
mes@meta.data$organ_2[which(mes@meta.data$seurat_clusters=="0")]<-"Va_Meso"
mes@meta.data$organ_2[which(mes@meta.data$seurat_clusters=="1")]<-"Va_Meso"
mes@meta.data$organ_2[which(mes@meta.data$seurat_clusters=="2")]<-"Va_Meso"
mes@meta.data$organ_2[which(mes@meta.data$seurat_clusters=="11")]<-"Va_Meso"
mes@meta.data$organ_2[which(mes@meta.data$seurat_clusters=="5")]<-"Va_IM"
mes@meta.data$organ_2[which(mes@meta.data$seurat_clusters=="10")]<-"Va_IM"
mes@meta.data$organ_2[which(mes@meta.data$seurat_clusters=="3")]<-"Ovary Surface Epi"
mes@meta.data$organ_2[which(mes@meta.data$seurat_clusters=="4")]<-"Ovary Pro_like Epi"
mes@meta.data$organ_2[which(mes@meta.data$seurat_clusters=="6")]<-"Ovary Meso"

DimPlot(mes,group.by = "organ_2",label=T,label.size = 8)
DimPlot(mes,group.by = "organ_2")+NoLegend()

  






DimPlot(mes,group.by = "type",cols=c("blue","grey","red"))


sub<-subset(mes,subset=ov_sub!="none")

StackedVlnPlot(obj = mes, group.by = "ov_sub", features = gene_for_vln2) + theme(axis.text.x = element_text(angle=90, hjust=1))


ave<-AverageExpression(sub)

tmp<-ave$RNA
tmp_cut<-tmp[gene_for_vln2,]

library(pheatmap)
pheatmap(tmp_cut,show_colnames = T,show_rownames = T,fontsize = 20,
         cluster_rows = F,cluster_cols = F,scale = "row")

pheatmap(tmp_cut,scale="row",show_colnames = T,show_rownames = T,fontsize = 20,
         cluster_rows = F,cluster_cols = F)


my_features <- c("PDE1A", "ABI3BP", "MARCH1", "PAX8", "CD44",
                 "EPCAM", "MYC", "LGR5", "HOXD4", "MSLN", "ITLN1")
sub<-subset(mes,idents = c("3","4","6"))
StackedVlnPlot(obj = sub, group.by = "seurat_clusters", features = my_features) + theme(axis.text.x = element_text(angle=90, hjust=1))



SK<-readRDS("New.SK.mutype.RDS")
SK$label<-paste0(SK$organ,"_",SK$mutype)
a<-table(SK$label)
a<-a[a>200]
Idents(SK)<-SK$label
ave<-AverageExpression(SK,assays = "RNA",slot = "data")
tmp<-ave$RNA[,names(a)]
co<-cor(tmp,method = "spearman")
pheatmap(co,show_colnames = T,show_rownames = T,fontsize = 20,
         cluster_rows = T,cluster_cols = T,scale="none")
Idents(SK)<-SK$label

T1<-subset(SK,idents = c("Abdominal_wall_Type I","Diaphragm_Type I","Tongue_Type I"))
T2a<-subset(SK,idents = c("Esophagus_Type IIa","Diaphragm_Type IIa","Tongue_Type IIa"))
T2b<-subset(SK,idents = c("Diaphragm_Type IIb","Abdominal_wall_Type IIb"))

DEG_T1<-FindAllMarkers(T1,logfc.threshold = 0.5,min.pct = 0.25,return.thresh = 0.01,only.pos = T)
DEG_T2a<-FindAllMarkers(T2a,logfc.threshold = 0.5,min.pct = 0.25,return.thresh = 0.01,only.pos = T)
DEG_T2b<-FindAllMarkers(T2b,logfc.threshold = 0.5,min.pct = 0.25,return.thresh = 0.01,only.pos = T)
DEG_T1$FC<-exp(DEG_T1$avg_logFC)
DEG_T2a$FC<-exp(DEG_T2a$avg_logFC)
DEG_T2b$FC<-exp(DEG_T2b$avg_logFC)

write.table(DEG_T1,file="DEG_T1_new.xls",quote=F,sep="\t")
write.table(DEG_T2a,file="DEG_T2a_new.xls",quote=F,sep="\t")
write.table(DEG_T2b,file="DEG_T2b_new.xls",quote=F,sep="\t")


library(dplyr)
library(viridis)

top10_T1 <- DEG_T1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10_T2a <- DEG_T2a %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10_T2b <- DEG_T2b %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

p1<-DotPlot(T1,features = unique(as.vector(top10_T1$gene)))+scale_color_viridis(direction=-1)+theme(axis.text.x = element_text(size=8,angle=90, hjust=1, vjust=1))

p2<-DotPlot(T2a,features = unique(as.vector(top10_T2a$gene)))+scale_color_viridis(direction=-1)+theme(axis.text.x = element_text(size=8,angle=90, hjust=1, vjust=1))

p3<-DotPlot(T2b,features = unique(as.vector(top10_T2b$gene)))+scale_color_viridis(direction=-1)+theme(axis.text.x = element_text(size=8,angle=90, hjust=1, vjust=1))




DEG_mutype<-FindAllMarkers(SK,logfc.threshold = 0.5,min.pct = 0.25,only.pos = T,return.thresh = 0.01)
DEG_mutype$FC<-exp(DEG_mutype$avg_logFC)

write.table(DEG_mutype,file="DEG_mutype.new.xls",quote=F,sep="\t")

SK@meta.data$mutype[which(SK@meta.data$mutype=="Type2a Fast muscle")]<-"Type IIa Myonuclei"
SK@meta.data$mutype[which(SK@meta.data$integrated_snn_res.0.8=="4")]<-"Type IIa Myonuclei"
DimPlot(SK,label=T,label.size=8,group.by = "mutype")

ave_SK_mutype<-AverageExpression(SK,assays = "RNA")
  
top20_mutype <- DEG_mutype %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
gene_mutype<-unique(as.vector(top20_mutype$gene))
tmp<-ave_SK_mutype$RNA[gene_mutype,]
library(pheatmap)
p1<-pheatmap(tmp,scale="row",show_colnames = T,show_rownames = T,fontsize = 10,
             cluster_rows = F,cluster_cols = F,border_color = "NA",color = colorRampPalette(c("#0000ff", "#ffffff", "#ff0000"))(256))
p1<-p1+coord_flip()


scanpy_color<-c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B","#FEE500","#8A9FD1","#C06CAB","#D8A767","#90D5E4","#89C75F","#F37B7D")

pdf(paste0("MUSK","_SK.pdf"))
p1
dev.off()

pdf(paste0("COL22A1","_SK.pdf"))
p2
dev.off()

pdf(paste0("DCN","_SK.pdf"))
p3
dev.off()

p1<-FeaturePlot(SK,features = "MUSK",pt.size=0.8,order=F)+
  theme_bw() + 
  theme(panel.grid =element_blank()) +  
  theme(axis.text = element_blank()) +  
  theme(axis.ticks = element_blank()) +   
  theme(panel.border = element_blank()) +
  labs(x = "", y = "", title = "")+
  NoLegend()+scale_color_viridis(direction=1)


gene_WNT<-c("FZD1","FZD2","FZD3","FZD4","FZD5","FZD6","FZD7","FZD8","FZD9","FZD10","LRP5","LRP6","MUSK","PTK7","ROR1","ROR2","RYK","WNT","WNT1","WNT2","WNT2B","WNT3","WNT3A","WNT4","WNT5A","WNT5B","WNT6","WNT7A","WNT7B","WNT8A","WNT8B","WNT9A","WNT9B","WNT10A","WNT16","LGR4","LGR5","LGR6","RSPO1","RSPO2","RSPO3","RSPO4")

pos<-which(gene_WNT %in% rownames(mes))
gene_WNT<-gene_WNT[pos]


for (i in 1:length(gene)){
  
  #gene<-c(")
  png(paste0("18.196_libs/",gene[i],"_pbmc.png"),width=600,height=600)
  
  p1<-FeaturePlot(pbmc,features = gene[i],pt.size=0.4,cols=c("lightgrey","red"))+
    theme_bw() + 
    theme(panel.grid =element_blank()) +  
    theme(axis.text = element_blank()) +  
    theme(axis.ticks = element_blank()) +   
    theme(panel.border = element_blank()) +
    labs(x = "", y = "", title = "")+
    NoLegend()
  
  print(p1)
  dev.off()
}

DotPlot(mes,features = gene_WNT)



ave<-AverageExpression(mes)
tmp<-ave$RNA
tmp_cut<-tmp[gene_WNT,]
library(pheatmap)
pheatmap(tmp_cut,show_colnames = T,show_rownames = T,fontsize = 20,
         cluster_rows = F,cluster_cols = F)



kidney@meta.data$Celltype = "Ascending LOH"
kidney@meta.data$Celltype[which(kidney@meta.data$seurat_clusters=="1")]<-"Proximal tubule S1"
kidney@meta.data$Celltype[which(kidney@meta.data$seurat_clusters=="2")]<-"Proximal tubule S3"
kidney@meta.data$Celltype[which(kidney@meta.data$seurat_clusters=="3")]<-"Principal"
kidney@meta.data$Celltype[which(kidney@meta.data$seurat_clusters=="4")]<-"Connecting tubule"
kidney@meta.data$Celltype[which(kidney@meta.data$seurat_clusters=="5")]<-"Distal convoluted tubule"
kidney@meta.data$Celltype[which(kidney@meta.data$seurat_clusters=="7")]<-"Proximal tubule S3"
kidney@meta.data$Celltype[which(kidney@meta.data$seurat_clusters=="9")]<-"Descending LOH"
kidney@meta.data$Celltype[which(kidney@meta.data$seurat_clusters=="10")]<-"ENDO"
kidney@meta.data$Celltype[which(kidney@meta.data$seurat_clusters=="11")]<-"Proximal tubule S3"
kidney@meta.data$Celltype[which(kidney@meta.data$seurat_clusters=="12")]<-"Principal"
kidney@meta.data$Celltype[which(kidney@meta.data$seurat_clusters=="13")]<-"Podocyte"
kidney@meta.data$Celltype[which(kidney@meta.data$seurat_clusters=="15")]<-"Intercalated"
kidney@meta.data$Celltype[which(kidney@meta.data$seurat_clusters=="16")]<-"Intercalated"
kidney@meta.data$Celltype[which(kidney@meta.data$seurat_clusters=="17")]<-"Principal"
kidney@meta.data$Celltype[which(kidney@meta.data$seurat_clusters=="18")]<-"Intercalated"
kidney@meta.data$Celltype[which(kidney@meta.data$seurat_clusters=="19")]<-"Myofibro"
kidney@meta.data$Celltype[which(kidney@meta.data$seurat_clusters=="20")]<-"Principal"
kidney@meta.data$Celltype[which(kidney@meta.data$seurat_clusters=="21")]<-"Descending LOH"
kidney@meta.data$Celltype[which(kidney@meta.data$seurat_clusters=="22")]<-"mDC"


gene<-c("LPIN1", "COL5A3", "ACSL1")
gene_pro<-c("CXCL14", "APOD", "CD34","CFD")
gene_immune<-c("CD163","MRC1","PTPRC","CSF1R")
gene_pro2<-c("PDGFRA","ITGB1","CD34","ATXN1")
gene_beg<-c("UCP1","ADIPOQ","CD36","PPARG")
gene_beg2<-c("TNFRSF9","TMEM26","TBX1")
gene_brown<-c("ABCD2","BLCAP","PFKFB3","TFF2","IFI44","CALR3","XAF1","GSDMB","EMP2")
gene_white<-c("DUSP26","STEAP4","ADAM28","F11","C4BPA","AOX1","CD52","CLDN4","FBXO2")


all<-c(gene,gene_pro,gene_immune)
FeaturePlot(AD ,features = gene)
FeaturePlot(AD ,features = gene_pro2)
FeaturePlot(AD ,features = gene_beg)
DotPlot(AD ,features = all)

library(dplyr)

top3_ST <- ST_DEG %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
top3_MA <- MA_DEG %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)


AD@meta.data$Celltype = "Adipocyte progenitor"
AD@meta.data$Celltype[which(AD@meta.data$seurat_clusters=="0")] = "Adipocyte"
AD@meta.data$Celltype[which(AD@meta.data$seurat_clusters=="5")] = "Adipocyte"
AD@meta.data$Celltype[which(AD@meta.data$seurat_clusters=="8")] = "Adipocyte"
AD@meta.data$Celltype[which(AD@meta.data$seurat_clusters=="9")] = "Adipocyte"
AD@meta.data$Celltype[which(AD@meta.data$seurat_clusters=="10")] = "Adipocyte"
AD@meta.data$Celltype[which(AD@meta.data$seurat_clusters=="12")] = "Adipocyte"
AD@meta.data$Celltype[which(AD@meta.data$seurat_clusters=="14")] = "Adipocyte"
AD@meta.data$Celltype[which(AD@meta.data$seurat_clusters=="6")] = "Immune"



gene<-c("NR1H3","MYC","KLF4","ESR1","RELA","BCL6")
FeaturePlot(AD,features = gene)






######Fig2-S2##############
setwd("D:/项目/猴脑/body atlas/Figure/Fig2/Commoncell")
rds<-readRDS("Seurat_commoncellSmooth_200.integrated.Find.RDS")
deg<-read.table("commoncell_Smooth_200_cut.DEG.xls",sep="\t",header=T)
library(dplyr)
library(viridis)
library(ggplot2)
top3<-deg %>% group_by(cluster) %>% top_n(n=3,wt =avg_logFC)
Idents(rds)<-rds$organ
p1<-DotPlot(rds,features = unique(as.vector(top3$gene)))+scale_color_viridis(direction=-1)+theme(axis.text.x = element_text(size=8,angle=90, hjust=1, vjust=1))


pdf("SM.dot.pdf",height=5,width=15)
p1
dev.off()


pdf("SM.umap.pdf",height=8,width=8)
DimPlot(rds,group.by = "organ",label=F)+NoLegend()
DimPlot(rds,group.by = "organ",label=T)+NoLegend()
dev.off()


SK$mustype2<-SK$mustype
SK$mustype2[which(SK$mutype == "TypeIIa")]<-"Type II"
SK$mustype2[which(SK$mutype == "TypeIIb")]<-"Type II"


library(org.Hs.eg.db)
library(clusterProfiler)

####GO term####
setwd("D:\\项目\\猴脑\\body atlas\\Figure\\Fig3")
DEG<-readRDS("pg.Folli.DEG.RDS")
DEG$FC<-exp(DEG$avg_logFC)
a_h<-DEG[DEG$FC>=2 & DEG$p_val_adj<=0.01,]
deg<-a_h
#if(nrow(a_h) > 0){
  eg = bitr(rownames(a_h), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  id<-as.character(eg[,2])
  
  ego_bp_h <- enrichGO(gene = id,
                       OrgDb=org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)
  #write.table(as.data.frame( ego_bp_h@result), file=paste("D:\\项目\\杂项\\冠状病毒\\北京患者\\ww\\6.BJT_BZJ_ZDD_D2_vs_3healthy\\",n,"D2.up.GO.xls",sep = ""),sep="\t",quote=F)
  #pdf(paste("D:\\项目\\杂项\\冠状病毒\\北京患者\\ww\\6.BJT_BZJ_ZDD_D2_vs_3healthy\\",n,"D2.up.GO.pdf"),width=15,height=15)
  a<-barplot(ego_bp_h, showCategory=20)
  print(a)
  dev.off()
  
  
  
  ego_mf_h <- enrichGO(gene = id,
                       OrgDb=org.Hs.eg.db,
                       ont = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)
  
  library(ggplot2)
  tmp<-head(x =ego_bp_h,n=20)
  library(forcats)
  library(dplyr)
  tmp$q_va<- (-1)*(log10(tmp$qvalue))
  tmp$Description = with(tmp, reorder(Description, q_va, median))

  p_bp<-ggplot(tmp,aes(x=Description,y=q_va,fill=-(qvalue)))+
        geom_bar(stat="identity")+
        coord_flip() +
        scale_fill_gradient(low="blue",high="red")+
         theme_bw()
   p_bp     
   
   
   
   
   
   
   
   
   
   
   
   
   library(pheatmap)
   
   gene=c("LGR5","LGR6","RSPO1","RSPO2","RSPO3","RSPO4","WNT1","WNT2","WNT2B","WNT3","WNT3A","WNT4","WNT5A","WNT5B","WNT6","WNT7A",
          "WNT7B","WNT8A","WNT8B","WNT9A","WNT9B","WNT10A","WNT10B","WNT11","WNT16",
          "FZD1","FZD2","FZD3","FZD4","FZD5","FZD6","FZD7",
          "FZD8","FZD9","FZD10")
   
   ###Cor celltype
   ave_cell=read.table("/Users/xiaoyuwei/Desktop/1.Research/1.Monkey_Atlas_RNA/000000.Figure-V1/SNP/Celltype_raw_mean_mat.tsv",header = T,sep = "\t")
   rownames(ave_cell)=ave_cell[,1]
   ave_cell=ave_cell[,-1]
   ave_cell_t=t(ave_cell)
   ave_cell_t_sub=ave_cell_t[na.omit(match(gene,rownames(ave_cell_t))),]
   pheatmap(cor(ave_cell_t_sub),fontsize_row = 6,fontsize_col = 6)
   ###Cor organ
   ave_organ=read.table("/Users/xiaoyuwei/Desktop/1.Research/1.Monkey_Atlas_RNA/000000.Figure-V1/SNP/Celltype_raw_mean_mat.tsv",header = T,sep = "\t")
   rownames(ave_organ)=ave_organ[,1]
   ave_organ=ave_organ[,-1]
   ave_organ_t=t(ave_organ)
   ave_organ_t_sub=ave_organ_t[na.omit(match(gene,rownames(ave_organ_t))),]
   pheatmap(cor(ave_organ_t_sub),fontsize_row = 6,fontsize_col = 6)
   
   
   
   
   
   SK$mutype<-"Type I"
   SK@meta.data$mutype[which(SK@meta.data$seurat_clusters=="0")]<-"Type IIb"
   SK@meta.data$mutype[which(SK@meta.data$seurat_clusters=="12")]<-"MTJ"
   SK@meta.data$mutype[which(SK@meta.data$seurat_clusters=="14")]<-"NMJ"
   SK@meta.data$mutype[which(SK@meta.data$seurat_clusters=="11")]<-"FAP"
   SK@meta.data$mutype[which(SK@meta.data$seurat_clusters=="1")]<-"Type IIb"
   SK@meta.data$mutype[which(SK@meta.data$seurat_clusters=="3")]<-"Type IIa"
   SK@meta.data$mutype[which(SK@meta.data$seurat_clusters=="4")]<-"Type IIa"
   SK@meta.data$mutype[which(SK@meta.data$seurat_clusters=="15")]<-"Type IIb"
   SK@meta.data$mutype[which(SK@meta.data$seurat_clusters=="13")]<-"Satellite"
   SK@meta.data$mutype[which(SK@meta.data$seurat_clusters=="5")]<-"Type IIa"
   SK@meta.data$mutype[which(SK@meta.data$seurat_clusters=="2")]<-"Type IIa"
   SK@meta.data$mutype[which(SK@meta.data$seurat_clusters=="10")]<-"Type IIa"
   
   DimPlot(SK,group.by = "mutype",label=T)
   c<-c(brewer.pal(8, "Dark2"),brewer.pal(8, "Dark2"))
c<-c[-6]


ggmeta<-SK@meta.data
p1<-ggplot(ggmeta,aes(x=organ,fill=mutype))+geom_bar(stat="count",position="fill")+ theme_bw()+scale_fill_manual(values = c)



DEG<-read.table("DEG_mutype.new.xls",header=T,sep="\t",row.names=1)
top20_mutype <- DEG %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
tmp<-tmp[as.character(top20_mutype$gene),]
p1<-pheatmap(tmp,scale="row",show_colnames = T,show_rownames = T,fontsize = 10,
             cluster_rows = F,cluster_cols = F,border_color = "NA",color = colorRampPalette(c("#0000ff", "#ffffff", "#ff0000"))(256))





p1<-FeaturePlot(Human,features = "LGR6",pt.size=0.5,order=T)+
  theme_bw() + 
  theme(panel.grid =element_blank()) +  
  theme(axis.text = element_blank()) +  
  theme(axis.ticks = element_blank()) +   
  theme(panel.border = element_blank()) +
  labs(x = "", y = "", title = "")+
  NoLegend()+scale_color_viridis(direction=1)

p2<-FeaturePlot(Mouse,features = "LGR6",pt.size=0.5,order=T)+
  theme_bw() + 
  theme(panel.grid =element_blank()) +  
  theme(axis.text = element_blank()) +  
  theme(axis.ticks = element_blank()) +   
  theme(panel.border = element_blank()) +
  labs(x = "", y = "", title = "")+
  NoLegend()+scale_color_viridis(direction=1)

p3<-FeaturePlot(Macaca,features = "LGR6",pt.size=0.5,order=T)+
  theme_bw() + 
  theme(panel.grid =element_blank()) +  
  theme(axis.text = element_blank()) +  
  theme(axis.ticks = element_blank()) +   
  theme(panel.border = element_blank()) +
  labs(x = "", y = "", title = "")+
  NoLegend()+scale_color_viridis(direction=1)


heart$CT<-"Cardiomyocyte"
heart$CT[which(heart$seurat_clusters == "2")]<-"Peri"
heart$CT[which(heart$seurat_clusters == "5")]<-"Stromal"
heart$CT[which(heart$seurat_clusters == "6")]<-"Endo"
heart$CT[which(heart$seurat_clusters == "8")]<-"Stromal"
heart$CT[which(heart$seurat_clusters == "9")]<-"Macrophage"
heart$CT[which(heart$seurat_clusters == "12")]<-"Endo"
heart$CT[which(heart$seurat_clusters == "13")]<-"Smooth muscle"
heart$CT[which(heart$seurat_clusters == "14")]<-"Endo"
heart$CT[which(heart$seurat_clusters == "15")]<-"Unknown"




Msub<-FindVariableFeatures(Msub,nfeatures = 3000)
Msub<-ScaleData(Msub)
Msub<-RunPCA(Msub)
Msub<-RunUMAP(Msub, reduction = "pca",dims = 1:20)
Msub<-FindNeighbors(Msub,dims = 1:20)
Msub<- FindClusters(Msub,resolution = 0.3)
Msub$Celltype<-"RP+ Mono 1"
Msub$Celltype[which(Msub$seurat_clusters == "2")] <- "RP+ Mono 2"
Msub$Celltype[which(Msub$seurat_clusters == "0")] <- "AHR+ Mono"
Msub$Celltype[which(Msub$seurat_clusters == "1")] <- "HLA+ Mono"
Msub$Celltype[which(Msub$seurat_clusters == "4")] <- "CD69+ Mono"
saveRDS(Msub,file="anno_Msub.RDS")


DEG_T1<-read.csv("allexlayer_v1.csv",header=T,row.names=1)
DEG_T2<-read.csv("allexlayer_pms.csv",header=T,row.names=1)





library(ggplot2)
library(RColorBrewer)
mycol <- brewer.pal(11,"RdYlBu")

setwd("D:\\项目\\疫苗\\V2方法学文章")
a<-read.table("merge.all.clean.xls",sep="\t",header=T,row.names=1)
a$pval<-"a:0"
a$pval[which(a$donor == ">=3")] <- "b:>=3"
a$label<-paste0(a$value,"_",a$X)
plot <- ggplot(data = a, mapping = aes_string(x = 'label', y = 'variable'))+
geom_point(mapping = aes_string(size = 'pval', color = "pval_group"))+
facet_grid(cols = vars(time),scales = "free")+coord_flip()+
scale_color_gradientn(colors=rev(mycol))+
theme_bw()+theme(axis.text.y=element_text(size=15),axis.text.x=element_text(angle=90, vjust=0.5, hjust=1,size=15))
