
rm(list = ls())
#Create a folder, set the working directory
mainDir <-"C:/Bioinformatics/000 EPC rt/GSE53625"
subDir <-"Step20 GSE45670 irgenes DEA"
outputDir <- file.path(mainDir, subDir)

if (!dir.exists(outputDir)){
  dir.create(outputDir)
} else {
  print("Dir already exists!")
}

setwd(mainDir)

setwd(outputDir)
getwd()

#
load("Step03.02.GSE45670.exprSet.annotated.radiosensity.Rdata")

lncRNAlist<-c("LINC01121","FAM167A-AS1","ADAMTS9-AS2",
              "MGC12916","MIR124-2HG","FAM167A-AS1")


lncRNAexp.diff<-exprSet.annotated.radiosensity[which(rownames(exprSet.annotated.radiosensity) %in% lncRNAlist),]
#Screening tumor samples
GSE45670.lncRNAexp.diff<-lncRNAexp.diff

#Load immune-related genes
load("0irgene.rda")


save(exprSet.annotated.radiosensity,lncRNAlist,
     group_list,irgene,
     file = "Figure 10 data.rda")



GSE45670.irgene.mRNAexp<-exprSet.annotated.radiosensity[which(rownames(exprSet.annotated.radiosensity) %in% irgene$irgene),]

#GSE45670.irgene.mRNAexp<-GSE45670.irgene.mRNAexp[,GSE45670.Group=="non_responder"]

library(limma) 
#载入lncRNA表达文件,并对数据进行处理
lncRNAExpMarix<-GSE45670.lncRNAexp.diff

class(lncRNAExpMarix)
exp <-lncRNAExpMarix;dim(exp)
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames);dim(data)
data=avereps(data)
dim(data)
data=data[rowMeans(data)>0.5,];dim(data);class(data)
save(data, file = "Step20.01.GSE45670.data.rda")
#

#
#Load immune gene expression files and process data
immGeneExp<-GSE45670.irgene.mRNAexp
class(immGeneExp)

dimnames=list(rownames(immGeneExp),colnames(immGeneExp))
immuneGene=matrix(as.numeric(as.matrix(immGeneExp)),
                  nrow=nrow(immGeneExp),dimnames=dimnames)
immuneGene=avereps(immuneGene)
immuneGene=immuneGene[rowMeans(immuneGene)>0.5,]#提取行均值大于0.5的条目
dim(immuneGene)

save(data,immuneGene, file = "Step20.02.GSE45670.data.immuneGene.rda")

#Correlation test

memory.limit(102400)
dim(immuneGene)
dim(data)
immuneGene<-immuneGene
lncRNA<-data

library(dplyr)
samplenum<-ceiling(dim(lncRNA)[1]/100)
samplenum
grp<-as.factor(rep(1:samplenum, each = 100,len =nrow(lncRNA)))
data <- cbind(grp,as.data.frame(lncRNA))
data_split <- split(data, data$grp)
class(data_split)

save(data_split,file = "Step20.03.data_split.rda")


load("Step20.03.data_split.rda")

outTablist<-list()

for (k in 1:length(data_split)){
  tmp <- as.matrix(data_split[[k]][,-1]) 
  corX = pvalueX =NULL
  lncRNAi = immuneGenej  =NULL
  outTab = NULL
  for(i in row.names(tmp)){
    if(sd(tmp[i,])>0.5){
      for(j in row.names(immuneGene)){
        x=as.numeric(tmp[i,])
        y=as.numeric(immuneGene[j,])
        corT=cor.test(x,y)
        cor=corT$estimate
        pvalue=corT$p.value
        lncRNAi = c(lncRNAi,i)
        immuneGenej = c(immuneGenej,j)
        corX = c(corX,cor)
        pvalueX= c(pvalueX,pvalue)
      }
    }
  }
  outTab = data.frame(immuneGenej,lncRNAi,corX,pvalueX)
  colnames(outTab)<-c("immuneGene","lncRNA","cor","pvalue")
  outTablist[[k]]<-outTab
  nameoutput=paste0("Step20.04.outTab.",k,".rda", sep = "")
  save(outTab, file = nameoutput)
  save(outTablist, file = "Step20.05.outTablist.rda")
}

load("Step20.05.outTablist.rda")
length(outTablist)

library("plyr") 
corResults.all<-do.call(rbind.fill,outTablist)

corResults.all$Regulation[corResults.all$cor>0]<-"postive"
corResults.all$Regulation[corResults.all$cor<0]<-"negative"
table(corResults.all$Regulation)
save(corResults.all,file = "Step20.06.GSE45670.corResults.all.rda")

load("Step20.06.GSE45670.corResults.all.rda")

library(limma)
corFilter.negtive=-0.4
pvalueFilter.negtive=0.01
corFilter.positive=0.4
pvalueFilter.positive=0.01

negtivedf<-subset(corResults.all, cor <= corFilter.negtive & pvalue<=pvalueFilter.negtive)
dim(negtivedf)

positivedf<-subset(corResults.all, cor >= (corFilter.positive) & pvalue<=pvalueFilter.positive)
dim(positivedf)

corResults<-rbind(negtivedf,positivedf)

save(negtivedf,positivedf,corResults,file = "Step20.07.GSE45670.corResults.rda")

GSE45670.imugene=unique(as.vector(corResults[,"immuneGene"]))

GSE45670.irGene.exp=exprSet.annotated.radiosensity[GSE45670.imugene,]

save(GSE45670.irGene.exp,GSE45670.imugene,group_list,file = "Step20.08.GSE45670.irGene.exp.rda")

expSet<-as.data.frame(GSE45670.irGene.exp)
expSet <- expSet[rowMeans(expSet)>0,] 
dim(expSet)

suppressMessages(library("stringr"))
class(group_list)
group_list
Group<-as.factor(ifelse(group_list==2,0,1))#注意1才是对照组
Group <- factor(Group,labels = c("responder","non_responder"))
table(Group)

###################=====limma DEA=====######################
suppressMessages(library(reshape2))
names.expSet<-rownames(expSet)
expSet.tmp<-cbind(names.expSet,expSet)
exp_L = melt(expSet.tmp)
head(exp_L)
colnames(exp_L)=c('symbol','sample','value')
head(exp_L)


# get group info
exp_L$group = rep(group_list,each=nrow(expSet))
head(exp_L)

# ggplot2 plotting 
library(ggplot2)
p = ggplot(exp_L,
           aes(x=sample,y=value,fill=group))+geom_boxplot()
print(p)

p=ggplot(exp_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
p=p+stat_summary(fun="mean",geom="point",shape=23,size=3,fill="red")
p=p+theme_set(theme_set(theme_bw(base_size=20)))
p=p+theme(text=element_text(face='bold'),axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank())
print(p)

#batch infect
library(limma) 
expSet = normalizeBetweenArrays(expSet)

#############Check sample grouping information#################

#hclust plotting
head(expSet)
group_list0<-group_list[group_list=="non_responder"]
group_list1<-group_list[group_list=="responder"]
clname<-c(paste(group_list0,1:17,sep=''),paste(group_list1,1:11,sep=''))

colnames(expSet) = clname
head(expSet)
#nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
# clust
hc=hclust(dist(t(expSet)))
par(mar=c(5,5,5,10)) 
# plotting
plot(as.dendrogram(hc), nodePar = nodePar,  horiz = TRUE)

###############PCA##################
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
if(! require("ggfortify")) install.packages("ggfortify")
#install.packages("ggfortify")

df=as.data.frame(t(expSet))
dim(df)
dim(expSet)

exp[1:6,1:6]
df[1:6,1:6]

df$group=group_list 
autoplot(prcomp(df[,1:(ncol(df)-1)] ), data=df,colour = 'group')
save(expSet,group_list,file = "step20.09.changecolname.Rdata")


##################数据准备阶段###############################
suppressMessages(library(limma))
suppressMessages(library("stringr"))

####log2 #####
log2expSet=log2(expSet)
exp<-normalizeBetweenArrays(log2expSet)
par(mfrow=c(1,2))
n.sample<-ncol(log2expSet)
cols <- rainbow(n.sample*1.2)
boxplot(data.frame(exp),col=cols,main="expression value",las=2) ## 画箱式图，比较数据分布情况

####grouping matrix
#group_list <-as.data.frame(group_list)
class(group_list)
group_list
group_list<-as.factor(ifelse(group_list==2,0,1))#注意1才是对照组
group_list <- factor(group_list,labels = c("responder","non_responder"))
table(group_list)
group_list

design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exp)
design
####3.0Difference Comparison Matrix
contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
contrast.matrix 
#####4.0 fit
##step1
fit <- lmFit(exp,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2)  

tempOutput = topTable(fit2, coef=1, n=Inf)
DEG.limma = na.omit(tempOutput) 
head(DEG.limma)

##Select differentially expressed genes
#或者logFC_cutoff <- with(DEG.limma,mean(abs(logFC)) + 2*sd(abs(logFC)) )
logFC_cutoff<-1
logFC_cutoff
P.Value_cutoff=0.05
k1 = (DEG.limma$P.Value < 0.05)&(DEG.limma$logFC < -logFC_cutoff)
k2 = (DEG.limma$P.Value < 0.05)&(DEG.limma$logFC > logFC_cutoff)
DEG.limma$change = as.factor(ifelse(k1,"DOWN",ifelse(k2,"UP","NOT")))
table(DEG.limma$change)
head(DEG.limma)
save(DEG.limma,logFC_cutoff,P.Value_cutoff,file = paste0("Step20.10.irgene.DEG.limma.rda"))

library("tidyverse")
DEG.limma.diff<-DEG.limma%>% filter(change %in% c("UP","DOWN"))
save(DEG.limma.diff, file = "Step20.11.DEG.limma.diff.rda")

####Output differentially expressed genes
diffirgenelist=as.vector(rownames(DEG.limma.diff))
log2.exp.normlize<-as.data.frame(exp)
diffirgene=log2.exp.normlize[diffirgenelist,]
save(diffirgene,group_list,diffirgenelist, file = "Step20.11.diffirgene.rda")

#####################Plot histogram of differentially expressed genes########################
lncRNAlist<-diffirgenelist
group_list<-group_list
sample_list<-colnames(diffirgene)

save(lncRNAlist,group_list,sample_list, file = "Step20.12.grp.Rda")


#选择差异的lncRNA
lncRNAexp.diff<-diffirgene

#转置表达矩阵
lncRNAexp.trans<-t(lncRNAexp.diff)
lncRNAexp.sig<-data.frame(group_list,lncRNAexp.trans)
colnames(lncRNAexp.sig)<-gsub("\\.","-",colnames(lncRNAexp.sig))

#整理数据，使之适合面板箱型图
col<-colnames(lncRNAexp.sig)
nrow.lncRNAexp.sig<-nrow(lncRNAexp.sig)
rownames(lncRNAexp.sig)<-NULL
length(col)
exp.lncRNA<-data.frame()
#函数形式：rep(x, time = , length = , each = ,)
for (i in (2:length(col))){#
  exp.temp<-NULL
  gene<-rep(col[i],nrow.lncRNAexp.sig)
  exp.temp<-subset(lncRNAexp.sig,select = c(1,i))
  colnames(exp.temp)<-c("Group","exp")
  exp.temp<-data.frame(exp.temp,gene)
  exp.lncRNA<-rbind(exp.lncRNA,exp.temp)
}

save(exp.lncRNA,file = "Step20.13.exp.lncRNA.Rda")


library("dplyr")
library("ggplot2")
library("ggpubr")
dev.off() 
p <- ggplot(data=exp.lncRNA,aes(x=Group,y=exp))
p+geom_boxplot(aes(fill=Group))+facet_wrap(~gene)+ 
  labs(x="Group", y = "Gene expression")+
  stat_compare_means(method = "t.test",aes(label = paste0("p =",as.character(sprintf("%0.3f", as.numeric(..p.format..))))))
ggsave(paste0("Differential expression analysis of immune-related genes",".tiff"),width = 15,height = 10)
