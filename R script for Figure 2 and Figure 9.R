rm(list = ls())
#setting work directory
mainDir <-"C:/Bioinformatics/000 EPC rt/radiosensitivity"
subDir <-"Step07.diff"
outputDir <- file.path(mainDir, subDir)
setwd(outputDir)
getwd()

suppressMessages(library(ggpubr));suppressMessages(library(ggthemes))

load("Step07.02.GSE45670.DEG.limma.rda")
#write.csv(DEG.limma,file = "dataset for Figure 2.csv")

deg.data <- DEG.limma
head(deg.data)
deg.data$logP <- -log10(deg.data$P.Value)

symbol<-rownames(deg.data)
deg.data<-cbind(symbol,deg.data)
deg.data$Label=""

deg.data <- deg.data[order(deg.data$P.Value),]

up.genes <- head(deg.data$symbol[which(deg.data$change == "UP")],10)
down.genes <- head(deg.data$symbol[which(deg.data$change =="DOWN")],10)

deg.top10.genes <- c(as.character(up.genes),as.character(down.genes))
deg.data$Label[match(deg.top10.genes,deg.data$symbol)] <- deg.top10.genes

#plotting
suppressMessages(library(ggrepel))
suppressMessages(library(ggplot2))
df = data.frame('logFC'=deg.data$logFC,
                'FDR'=deg.data$logP,
                'Gene'=rownames(deg.data),
                'change'=deg.data$change,
                'label'=deg.data$Label,
                stringsAsFactors = F)
colnames(df)

wn_volcano = function(df,logFC_cutoff=1,P.Value_cutoff=0.05){
  require(RColorBrewer,quietly = T,warn.conflicts =F)
  require(ggplot2,quietly = T,warn.conflicts =F)
  volcano1=ggplot(data = df, aes(x = logFC, y = FDR))
  volcano1=volcano1+geom_point(alpha=0.8, size=2, aes(color=change))  +  
    scale_color_manual(values=c('#00AFBB', '#999999', '#FC4E07')) +
    geom_vline(xintercept=c(-logFC_cutoff,logFC_cutoff),lty=3,col="azure4",lwd=1)+  
    geom_hline(yintercept = -log10(P.Value_cutoff),lty=3,col="azure4",lwd=1)+  
    ylab('-log10(FDR)')+
    xlab('log2(FoldChange)')+
    theme_bw()
  
  #p1
  volcano<-volcano1 + geom_text_repel(data=df, aes(label= label), 
                                      color="black", family="Times New Roman",
                                      size=2.5, fontface="italic", 
                                      arrow = arrow(ends="first", length = unit(0.01, "npc")),
                                      box.padding = 0.2,
                                      point.padding = 0.3, 
                                      segment.color = 'black', 
                                      segment.size = 0.3, 
                                      max.overlaps = getOption("ggrepel.max.overlaps", default = 10),
                                      segment.alpha=0.8,
                                      force = 1, 
                                      max.iter = 3e3)
  return(volcano)
}

##call the function
volcano = wn_volcano(df=df,logFC_cutoff=logFC_cutoff,P.Value_cutoff=P.Value_cutoff)
volcano
class(volcano)
##################heatmap#################
load("C:/Bioinformatics/000 EPC rt/radiosensitivity/Step07.diff/Step07.04.diffLncRNAExp.rda")
write.csv(diffLncRNAExp,file = "heatmapdata for Figure 2.csv")


#install.packages("pheatmap")
#install.packages("ggplot2") 
wn_pheatmap = function(heat.df,Group){
  suppressMessages(library(pheatmap)) 
  suppressMessages(library(ggplot2)) 
  annotation_col = data.frame(Group)
  rownames(annotation_col)
  colnames(heat.df)
  rownames(annotation_col) <- colnames(heat.df)
 
  pheatmap <- pheatmap(heat.df,scale="row",
                       annotation_col = annotation_col,
                       border="white", 
                       cluster_cols = F,treeheight_col = 25,
                       cluster_rows = T,treeheight_row = 25,
                       show_rownames = T, 
                       show_colnames = T,
                       legend = T, 
                       fontsize_row = 5,fontsize_col = 5, 
                       clustering_distance_rows = "correlation", 
                       clustering_method="single", 
                       angle_col = 45, 
                       cellwidth = 8,cellheight = 8, 
                       cutree_cols = 6, cutree_rows =5
                       ) 
  
  return(pheatmap)
}
pheatmap = wn_pheatmap(heat.df=diffLncRNAExp,Group=group_list)

print(pheatmap)


if (!require("ggplotify", character.only=T, quietly=T)) {
  install.packages("ggplotify")
  library("ggplotify", character.only=T)
}
if (!require("patchwork", character.only=T, quietly=T)) {
  install.packages("patchwork")
  library("patchwork", character.only=T)
}
if (!require("cowplot", character.only=T, quietly=T)) {
  install.packages("cowplot")
  library("cowplot", character.only=T)
}

pheatmap.grid <- grid2grob(print(pheatmap))
class(pheatmap)

save(volcano,pheatmap,pheatmap.grid,file="Step07.05.picdata.rda")

p3<-plot_grid(pheatmap.grid,volcano,
              labels = "AUTO",
              rel_widths = c(1.5, 1),
              rel_heights=c(1.5,1),
              ncol=2)
print(p3)

#save plotting
phname<-"FIGURE 2  Differential expression results of immune-related lncRNAs in GSE45670. (A) Heatmap, (B) Volcano map"
phtype<-"tiff"
ggsave(filename = paste0(phname,".",phtype),p3,width =10, 
       height = 5.5, dpi = 300, units = "in", device=phtype)
