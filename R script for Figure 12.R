############# GO KEGG ###################

rm(list = ls())
#Create a folder, set the working directory
mainDir <-"C:/Bioinformatics/000 EPC rt/GSE53625"
subDir <-"Step21 GSE45670 irgenes GO.KEGG"
outputDir <- file.path(mainDir, subDir)

if (!dir.exists(outputDir)){
  dir.create(outputDir)
} else {
  print("Dir already exists!")
}

setwd(mainDir)

file=paste0("Step20.11.DEG.limma.diff.rda")
filefrom<-paste0("Step20 GSE45670 irgenes DEA/",file)
fileto<-paste0("Step21 GSE45670 irgenes GO.KEGG/",file)

file.copy(filefrom, 
          fileto, 
          overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE) # 文件复制

setwd(outputDir)
getwd()

#BiocManager::install("ReactomePA")
library(ReactomePA)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("reactome.db")
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
#install.packages("ggridges")
library("ggridges")
library("DO.db")

load(file)

genelist_input <- DEG.limma.diff
genelist_input$Gene<-rownames(genelist_input)
class(genelist_input$Gene)
genename <- rownames(genelist_input) 
gene_map <- select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID"))
colnames(gene_map)[1]<-"Gene"

rt<-inner_join(gene_map,genelist_input,by = "Gene")
rt<-rt[,-1]
rt<-na.omit(rt)
rt$logFC<-sort(rt$logFC,decreasing = T)

geneList = rt[,2]
names(geneList) = as.character(rt[,1])
geneList

save(geneList,rt,file="Step21.01.GO.KEGG.geneList.Rda")

geneFC=rt$logFC
gene=rt$ENTREZID
names(geneFC)=gene
class(gene)


# BP  enrichGO 
BPenrich <- enrichGO(gene = gene, 
                     OrgDb = org.Hs.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff = 1, 
                     qvalueCutoff = 1,
                     readable = TRUE)
head(BPenrich,2)
summary.BP<-summary(BPenrich)
write.csv(summary(BPenrich),"Step21.02.GO.BPenrich.csv",row.names =FALSE)
save(summary.BP,file = "Step21.02.GO.BPenrich.Rda")

tiff(file="Step21.02.GO.BPenrich.dotplot.tiff",width =20,height =30,units ="cm",compression="lzw",bg="white",res=300)
dotplot(BPenrich, showCategory = 47)
dev.off()

tiff(file="Step21.02.GO.BPenrich.barplot.tiff",width =20,height =30,units ="cm",compression="lzw",bg="white",res=300)
barplot(BPenrich, showCategory = 20)
dev.off()

tiff(file="Step21.02.GO.BPenrich.plotGOgraph.tiff",width =20,height =30,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(BPenrich)
dev.off()

# CC enrichGO
CCenrich <- enrichGO(gene = gene, 
                     OrgDb = org.Hs.eg.db, 
                     ont = "CC", 
                     pAdjustMethod = "BH",
                     pvalueCutoff = 1, 
                     qvalueCutoff = 1,
                     readable = TRUE)
head(CCenrich,2)
summary.CC<-summary(CCenrich)
write.csv(summary(CCenrich),"Step21.03.GO.CCenrich.csv",row.names =FALSE)
save(summary.CC,file = "Step21.03.GO.CCenrich.Rda")

tiff(file="Step21.03.GO.CCenrich.dotplot.tiff",width =20,height =30,units ="cm",compression="lzw",bg="white",res=300)
dotplot(CCenrich, showCategory = 47)
dev.off()

tiff(file="Step21.03.GO.CCenrich.barplot.tiff",width =20,height =30,units ="cm",compression="lzw",bg="white",res=300)
barplot(CCenrich, showCategory = 20)
dev.off()

tiff(file="Step21.03.GO.CCenrich.plotGOgraph.tiff",width =20,height =30,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(CCenrich)
dev.off()

# MF enrichGO
MFenrich <- enrichGO(gene = gene, 
                     OrgDb = org.Hs.eg.db, 
                     ont = "MF", 
                     pAdjustMethod = "BH",
                     pvalueCutoff = 1, 
                     qvalueCutoff = 1,
                     readable = TRUE)
head(MFenrich,2)
summary.MF<-summary(MFenrich)
write.csv(summary(MFenrich),"Step21.04.GO.MFenrich.csv",row.names =FALSE)
save(summary.MF,file = "Step21.04.GO.MFenrich.Rda")

tiff(file="Step21.04.GO.MFenrich.dotplot.tiff",width =20,height =30,units ="cm",compression="lzw",bg="white",res=300)
dotplot(MFenrich, showCategory = 20)
dev.off()

tiff(file="Step21.04.GO.MFenrich.barplot.tiff",width =20,height =30,units ="cm",compression="lzw",bg="white",res=300)
barplot(MFenrich, showCategory = 20)
dev.off()

tiff(file="Step21.04.GO.MFenrich.plotGOgraph.tiff",width =20,height =30,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(MFenrich)
dev.off()
gene
genename


KEGGenrich <- enrichKEGG(gene = gene, 
                         organism = "hsa", 
                         pvalueCutoff =1, 
                         qvalueCutoff =1,
                         use_internal_data = T)
save(KEGGenrich, file ="Step21.05.KEGGenrich.Rda")
write.csv(summary(KEGGenrich),"Step21.05.KEGG.enrich.csv",row.names =FALSE)
#KEGG bar
tiff(file="Step21.05.KEGGenrich.barplot.tiff",width = 20,height = 20,units ="cm",compression="lzw",bg="white",res=300)
barplot(KEGGenrich, drop = TRUE, showCategory = 12)
dev.off()

#KEGG  dot
tiff(file="Step21.05.KEGGenrich.dotplot.tiff",width = 20,height = 20,units ="cm",compression="lzw",bg="white",res=300)
dotplot(KEGGenrich)
dev.off()

#Merge plotting

BP.plotdata<- BPenrich@result
BP.plotdata <- BP.plotdata[BP.plotdata$p.adjust<0.05,]
library(ggplot2)
BP.bar<-ggplot(BP.plotdata, aes(x=Description, y=Count,fill=p.adjust)) + 
  geom_bar(stat="identity") + 
  scale_fill_gradient(low = 'blue',high='red')+
  scale_x_discrete(name ="pathway names") +
  scale_y_continuous(name ="Count") +
  coord_flip() + 
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("GO biological process enrichment ")
print(BP.bar)

CC.plotdata<- CCenrich@result
CC.plotdata <- CC.plotdata[CC.plotdata$p.adjust<0.3,]
library(ggplot2)
CC.bar<-ggplot(CC.plotdata, aes(x=Description, y=Count,fill=p.adjust)) + 
  geom_bar(stat="identity") + 
  scale_fill_gradient(low = 'blue',high='red')+
  scale_x_discrete(name ="pathway names") +
  scale_y_continuous(name ="Count") +
  coord_flip() + 
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("GO cellular component enrichment ")
print(CC.bar)


MF.plotdata<- MFenrich@result
MF.plotdata <- MF.plotdata[MF.plotdata$p.adjust<0.05,]
library(ggplot2)
MF.bar<-ggplot(MF.plotdata, aes(x=Description, y=Count,fill=p.adjust)) + 
  geom_bar(stat="identity") + 
  scale_fill_gradient(low = 'blue',high='red')+
  scale_x_discrete(name ="pathway names") +
  scale_y_continuous(name ="Count") +
  coord_flip() + 
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("GO molecular function enrichment ")
print(MF.bar)

KEGG.plotdata<- KEGGenrich@result
KEGG.plotdata <- KEGG.plotdata[KEGG.plotdata$pvalue<0.05,]
library(ggplot2)
KEGG.bar<-ggplot(KEGG.plotdata, aes(x=Description, y=Count,fill=pvalue)) + 
  geom_bar(stat="identity") + 
  scale_fill_gradient(low = 'blue',high='red')+
  scale_x_discrete(name ="pathway names") +
  scale_y_continuous(name ="Count") +
  coord_flip() + 
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("KEGG Pathway Enrichment")
print(KEGG.bar)

library(cowplot)
p3<-plot_grid(BP.bar, CC.bar, MF.bar, KEGG.bar,labels = "AUTO")
class(p3)
ggsave(filename = "Step21.06.Pathway Enrichment.png",p3,width =30, 
       height = 10, dpi = 300, units = "in", device='png')
ggsave(filename = "Step21.06.Pathway Enrichment.tiff",p3,width =30, 
       height = 10, dpi = 300, units = "in", device='tiff')


save(BPenrich,CCenrich,MFenrich,KEGGenrich,file = "Step21.06.Pathway Enrichment.Rda")


#pathview
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager") 
BiocManager::install("pathview")
library("pathview")
keggxls=KEGG_gseresult
for(i in keggxls$ID){
  pv.out <- pathview(gene.data = geneFC, pathway.id = i, species = "hsa", out.suffix = "pathview")
}



