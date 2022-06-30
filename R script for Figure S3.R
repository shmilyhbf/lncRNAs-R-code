#Clinical characteristics and riskscore correlation analysis

rm(list = ls())
memory.limit(102400)
#Create a folder, set the working directory
mainDir <-"C:/Bioinformatics/000 EPC rt/GSE53625"
subDir <-"Step17.cliCor"
outputDir <- file.path(mainDir, subDir)

if (!dir.exists(outputDir)){
  dir.create(outputDir)
} else {
  print("Dir already exists!")
}

setwd(mainDir)

#copy files
file.copy("Step12.ROC/Step12.03.risk.df.rda", 
          "Step17.cliCor/Step12.03.risk.df.rda", 
          overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE) 

file.copy("Step14.indep/Step14.01.cli.factors.rda", 
          "Step17.cliCor/Step14.01.cli.factors.rda", 
          overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE) 

setwd(outputDir)
getwd()

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")


#引用包
library(limma)
library(ggpubr)

#load risk files
load("Step12.03.risk.df.rda")
risk=risk.df[,c("futime","fustat","riskScore")]

#load clinical data files
load("Step14.01.cli.factors.rda")

cli<-cli.factors
#Consolidate data
samSample=intersect(row.names(risk), row.names(cli))
risk=risk[samSample,"riskScore",drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk, cli)


library(tidyverse)

boxplot<-list()
m<-0
for(clinical in colnames(rt[,2:ncol(rt)])){
  m=m+1
  data=rt[c("riskScore", clinical)]
  colnames(data)=c("riskScore", "clinical")
  data=data[(data[,"clinical"]!="unknow"),]
  data<-data %>% drop_na()
  #set comparison group
  group=levels(factor(data$clinical))
  data$clinical=factor(data$clinical, levels=group)
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  #draw a boxplot
  boxplot[[m]]=ggboxplot(data, x="clinical", y="riskScore", color="clinical",
                    xlab=clinical,
                    ylab="Risk score",
                    legend.title=clinical,
                    add = "jitter")+ 
    stat_compare_means(comparisons = my_comparisons)
}
p3<-plot_grid(boxplot[[1]], boxplot[[2]], 
              boxplot[[3]], boxplot[[4]], 
              boxplot[[5]], boxplot[[6]], 
              boxplot[[7]], boxplot[[8]], 
              boxplot[[9]], 
              labels = "AUTO")

class(p3)
ggsave(filename = "clinical.riskScore.png",p3,width =20, 
       height = 15, dpi = 300, units = "in", device='png')


