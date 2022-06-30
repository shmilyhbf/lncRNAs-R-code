###############ROC#####################

rm(list = ls())
memory.limit(102400)
#Create a folder, set the working directory
mainDir <-"C:/Bioinformatics/000 EPC rt/GSE53625"
subDir <-"Step12.ROC"
outputDir <- file.path(mainDir, subDir)

if (!dir.exists(outputDir)){
  dir.create(outputDir)
} else {
  print("Dir already exists!")
}

setwd(mainDir)

file.copy("Step11.model/Step11.07.riskOut.rda", "Step12.ROC/Step11.07.riskOut.rda", 
          overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE) 

setwd(outputDir)
getwd()


library(survivalROC)  
load("Step11.07.riskOut.rda")
rt<-riskOut
class(rt$futime)
class(rt$fustat)
class(rt$riskScore)
#ROC曲线绘制

predictTime=1*12 
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker=rt$riskScore, predict.time =predictTime, method="KM")
pdf(file="Step12.01.ROC.pdf", width=5.5, height=5.5)
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="black", 
     xlab="False positive rate", ylab="True positive rate",
     lwd = 2, cex.main=1.2, cex.lab=1.2, cex.axis=1.2, font=1.2)
polygon(x=c(0,roc$FP,1,0),y=c(0,roc$TP,1,0),col="#24B35D",border=NA)
text(0.85, 0.1, paste0("AUC=",sprintf("%.3f",roc$AUC)), cex=1.2)
segments(0,0,1,1,lty=2)
dev.off()

#Get the best cutoff
predictTime=1*12 
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker=rt$riskScore, predict.time =predictTime, method="KM")
sum=roc$TP-roc$FP
cutOp=roc$cut.values[which.max(sum)]
cutTP=roc$TP[which.max(sum)]
cutFP=roc$FP[which.max(sum)]
pdf(file="Step12.02.ROC.cutoff.pdf",width=5.5,height=5.5)
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="black", 
     xlab="False positive rate", ylab="True positive rate",
     lwd = 2, cex.main=1.2, cex.lab=1.2, cex.axis=1.2, font=1.2)
polygon(x=c(0,roc$FP,1,0),y=c(0,roc$TP,1,0),col="#24B35D",border=NA)
segments(0,0,1,1,lty=2)
points(cutFP,cutTP, pch=20, col="red",cex=1.5)
text(cutFP+0.15,cutTP-0.05,paste0("Cutoff:",sprintf("%0.3f",cutOp)))
text(0.85, 0.1, paste0("AUC=",sprintf("%.3f",roc$AUC)), cex=1.2)
dev.off()


#Grouped by best  cutoff
risk=as.vector(ifelse(rt$riskScore>cutOp,"high","low"))
risk.df=cbind(rt, risk)
class(risk.df)

save(risk.df,file = "Step12.03.risk.df.rda")

######ROC curves for multiple times######
rocCol=c("red", "green", "blue")
aucText=c()
pdf(file="Step12.04.ROC.multiTime.pdf",width=6,height=6)
#1 year ROC
predictTime=1*12
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, 
                status=rt$fustat, 
                marker=rt$riskScore, 
                predict.time=predictTime, 
                method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText,paste0("one year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
abline(0,1)
#2 year ROC
predictTime=2*12
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker=rt$riskScore, predict.time =predictTime, method="KM")
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)
aucText=c(aucText,paste0("two year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
#3 year ROC
predictTime=3*12
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker=rt$riskScore, predict.time =predictTime, method="KM")
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],lwd = 2)
aucText=c(aucText,paste0("three year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
#draw legend
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
dev.off()


