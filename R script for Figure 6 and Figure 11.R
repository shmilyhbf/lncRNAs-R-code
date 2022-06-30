
#clear data
rm(list = ls())

#Create a folder, set the working directory
mainDir <-"C:/Bioinformatics/000 EPC rt/GSE53625"
subDir <-"Step14.difflncRNA.plot"
outputDir <- file.path(mainDir, subDir)

if (!dir.exists(outputDir)){
  dir.create(outputDir)
} else {
  print("Dir already exists!")
}

setwd(mainDir)

file=paste0("immuneLncRNAexp.rda")
filefrom<-paste0("Step05.irlncRNAselect/",file)
fileto<-paste0("Step14.difflncRNA.plot/",file)

file.copy(filefrom, 
          fileto, 
          overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE) # 文件复制

setwd(outputDir)
getwd()

load(file)

load("Step14.01.grp.Rda")
load("Step12.03.risk.df.rda")

lncRNAlist<-c("LINC01121","FAM167A-AS1","ADAMTS9-AS2",
              "MGC12916","MIR124-2HG","FAM167A-AS1")
group_risk<-risk.df$risk
sample_list<-rownames(risk.df)

save(lncRNAlist,group_risk,sample_list, file = "Step14.01.grp.Rda")

#lncRNAs selection
lncRNAexp.diff<-lncRNAexp.select[which(rownames(lncRNAexp.select) %in% lncRNAlist),which(colnames(lncRNAexp.select) %in% sample_list)]

#Transpose Expression Matrix
lncRNAexp.trans<-t(lncRNAexp.diff)
lncRNAexp.sig<-data.frame(group_risk,lncRNAexp.trans)
colnames(lncRNAexp.sig)<-gsub("\\.","-",colnames(lncRNAexp.sig))

#Organize data to fit in a panel boxplot
#Grouped by risk

col<-colnames(lncRNAexp.sig)
nrow.lncRNAexp.sig<-nrow(lncRNAexp.sig)
rownames(lncRNAexp.sig)<-NULL
length(col)
exp.lncRNA<-data.frame()

for (i in (2:length(col))){
  exp.temp<-NULL
  gene<-rep(col[i],nrow.lncRNAexp.sig)
  exp.temp<-subset(lncRNAexp.sig,select = c(1,i))
  colnames(exp.temp)<-c("Group","exp")
  exp.temp<-data.frame(exp.temp,gene)
  exp.lncRNA<-rbind(exp.lncRNA,exp.temp)
}
save(exp.lncRNA,file = "Step14.02.01.exp.lncRNA.Rda")


library("dplyr")
library("ggplot2")
library("ggpubr")
p <- ggplot(data=exp.lncRNA,aes(x=Group,y=exp))
print(p)
p+geom_boxplot(aes(fill=Group))+facet_wrap(~gene)+ 
  labs(x="Different group of riskscore", y = "Gene expression")+
  stat_compare_means(method = "t.test",aes(label = paste0("p =",as.character(sprintf("%0.3f", as.numeric(..p.format..))))))
ggsave(paste0("Fig6 different expression of ir-lncRNAs",".tiff"),width = 10,height = 10)




