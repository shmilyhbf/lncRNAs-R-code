setwd("C:/Bioinformatics/000 EPC rt/GSE53625/Step19.immuneCor/CIBERSORT")#设置三个文件所在的文件夹
rm(list = ls()) 

#data preparation
load("Figure 8 data.Rdata")
rowname<-rownames(exprSet.annotated.radiosensity)
data<-data.frame(rowname,exprSet.annotated.radiosensity)
library(tidyverse)
colnames(data)[1] <- 'Gene symbol'
write.table(data,file ="data.txt" ,sep ='\t',row.names=F,col.names=T)

library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)

source('Cibersort.R')#Download from Cibersort official website

LM22.file <- "LM22.txt"#Download from Cibersort official website
TCGA_exp.file <- "data.txt"

TCGA_TME.results <- CIBERSORT(LM22.file ,TCGA_exp.file, perm =1000, QN = F)  

save(TCGA_TME.results,group_list, file = "CIBERSORT.Rda")

load("CIBERSORT.Rda")
group_list
Group<-factor(as.character(group_list),levels=c('2','1'),labels=c("responder","non_responder"))
table(Group) # Normal 43   Tumor 43 

## 3. ploting
# 3.1 Coarse processing of data
TME_data <- as.data.frame(TCGA_TME.results[,1:22])
TME_data$group <- Group
TME_data$sample <- row.names(TME_data)

# 3.2 data melt
TME_New = melt(TME_data)

colnames(TME_New)=c("Group","Sample","Celltype","Composition")  #设置行名
head(TME_New)

# 3.3 Plot by median proportion of immune cells (optional)
plot_order = TME_New[TME_New$Group=="non_responder",] %>%
  group_by(Celltype) %>% 
  summarise(m = median(Composition)) %>% 
  arrange(desc(m)) %>% 
  pull(Celltype)

TME_New$Celltype = factor(TME_New$Celltype,levels = plot_order)


# 3.3 plotting
if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) }

box_TME <- ggplot(TME_New, aes(x = Celltype, y = Composition))+ 
  labs(y="Cell composition",x= NULL,title = "TME Cell composition")+  
  geom_boxplot(aes(fill = Group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group =  Group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)

box_TME;ggsave("EPC_TME.pdf",box_TME,height=15,width=25,unit="cm")


##4.1. Extract the top 20 immune cells used in the literature
# 4.1  Extract the top 20 immune cells
TCGA_TME_four = as.data.frame(TCGA_TME.results[,1:20])
head(TCGA_TME_four,3)

# 4.2 Classification of immune cells based on literature
immCell_four_type <- read.table("Cibersort_four_types.txt", header = T, row.names = NULL, sep = "\t")
colname<- colnames(TCGA_TME_four)
colname1<-gsub('\\.', ' ', colname)
colnames(TCGA_TME_four)<-colname1
####
colnames(TCGA_TME_four) == immCell_four_type$Immune.cells #T

##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

head(immCell_four_type)
# 4.3 
TCGA_TME_four$group = Group
TCGA_TME_four$sample <- row.names(TCGA_TME_four)
TME_four_new = melt(TCGA_TME_four)

## Using group, sample as id variables

colnames(TME_four_new) = c("Group","Sample","Immune.cells","Composition")

TCGA_TME_four_new2 = left_join(TME_four_new, immCell_four_type, by = "Immune.cells") %>% 
  group_by(Sample,Group,Types) %>%
  drop_na(Types)%>% 
  summarize(Sum = sum(Composition))

## `summarise()` regrouping output by 'Sample', 'Group' (override with `.groups` argument)

# plotting
box_four_immtypes <- ggplot(TCGA_TME_four_new2, aes(x = Group, y = Sum))+ 
  labs(y="Cell composition",x= NULL)+  #,title = "TCGA"
  geom_boxplot(aes(fill = Group),position=position_dodge(0.5),width=0.5,size=0.4,
               outlier.alpha = 1, outlier.size = 0.5)+ 
  theme_bw() + mytheme + 
  scale_fill_manual(values = c("#1CB4B8","#EB7369"))+ 
  scale_y_continuous(labels = scales::percent)+
  facet_wrap(~ Types,scales = "free",ncol = 4) + 
  stat_compare_means(aes(group =  Group),
                     label = "p.format",
                     method = "wilcox.test",
                     size = 3.5,
                     hide.ns = T)
box_four_immtypes;ggsave("Cibersort_four_immune_cell_types.pdf",
                         box_four_immtypes ,height= 12,width=25,unit="cm")

#Merge the two graphs
class(box_four_immtypes)
library(ggplotify)
library(patchwork)
library(cowplot)
p1 <- grid2grob(print(forst1))
p2 <- grid2grob(print(forst2))

p3<-plot_grid(box_TME,box_four_immtypes,labels = "AUTO",ncol = 1)
class(p3)
ggsave(filename = " infiltrating immune cell types.tiff",p3,width =12, 
       height = 15, dpi = 300, units = "in", device='tiff')












