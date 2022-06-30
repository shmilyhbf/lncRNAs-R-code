#survival analysis

rm(list = ls())
memory.limit(102400)
mainDir <-"C:/Bioinformatics/000 EPC rt/GSE53625"
subDir <-"Step13.survival"
outputDir <- file.path(mainDir, subDir)

if (!dir.exists(outputDir)){
  dir.create(outputDir)
} else {
  print("Dir already exists!")
}

setwd(mainDir)

file.copy("Step12.ROC/Step12.03.risk.df.rda", 
          "Step13.survival/Step12.03.risk.df.rda", 
          overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE) # 文件复制

setwd(outputDir)
getwd()



#install.packages("survival")
#install.packages("survminer")


library(survival)
library(survminer)
load("Step12.03.risk.df.rda")
inputdata<-risk.df

bioSurvival=function(inputFile=null,outFile=null){
  rt<-inputFile
  diff=survdiff(Surv(futime, fustat) ~risk, data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%0.3f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=TRUE,
                     pval=pValue,
                     pval.size=6,
                     palette=c("red", "blue"),
                     legend.title="Risk",
                     legend.labs=c("High risk", "Low risk"),
                     xlab="Time(years)",
                     break.time.by = 12,
                     risk.table=TRUE,
                     risk.table.title="",	  
                     risk.table.height=.25)
  pdf(file=outFile,onefile = FALSE,width = 6.5,height =5.5)
  print(surPlot)
  dev.off()
}
bioSurvival(inputFile=inputdata,outFile="Step13.01.survival.pdf")
