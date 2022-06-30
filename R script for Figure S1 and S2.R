
rm(list = ls())
memory.limit(102400)
#Create a folder, set the working directory
mainDir <-"C:/Bioinformatics/000 EPC rt/GSE53625"
subDir <-"Step11.model"
outputDir <- file.path(mainDir, subDir)

if (!dir.exists(outputDir)){
  dir.create(outputDir)
} else {
  print("Dir already exists!")
}

setwd(mainDir)

#copy files
file.copy("Step10.uniCox/Step10.04.uniCoxdf.select.rda", "Step11.model/Step10.04.uniCoxdf.select.rda", 
          overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE) 
file.copy("Step10.uniCox/Step10.05.surSigExp.rda", "Step11.model/Step10.05.surSigExp.rda", 
          overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)

file.copy("Step09.mergeTime/Step09.01.survdata.expselect.cliselect.rda", 
          "Step11.model/Step09.01.survdata.expselect.cliselect.rda", 
          overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)



setwd(outputDir)
getwd()


library(survival)
library(survminer)
library(glmnet)

#load files
load("Step10.05.surSigExp.rda")

class(surSigExp)
class(surSigExp$fustat)
surSigExp$fustat<-as.numeric(surSigExp$fustat)-1


library(glmnet) 
library(rms)
library(VIM) 
library(survival)
#load dataset
dt <-surSigExp

str(dt) 
aggr(dt,prop=T,numbers=T)
dt <- na.omit(dt)

###################step 01 Data curation######################

for(i in names(dt)[c(3:ncol(dt))]) {dt[,i] <- as.factor(dt[,i])}
str(dt)
x <- data.matrix(dt[,c(3:ncol(dt))])
y <- data.matrix(Surv(dt$futime,dt$fustat))


#################step 02 Lasso Regression Analysis##########################
fit.cox <-glmnet(x,y,family = "cox",alpha = 1, nlambda = 50)
dev.new()
tiff("Step11.01.lasso.lambda.tiff")
plot(fit.cox, xvar = "lambda", label = TRUE)
dev.off()
print(fit.cox)

fit.cv <- cv.glmnet(x,y,family="cox", alpha=1,nfolds=10)

#save picture of LASSO lamda
tiff("Step11.02.lasso.cvfit.tiff")
plot(fit.cv)
abline(v=log(c(fit.cv$lambda.min,fit.cv$lambda.1se)),lty="dashed")
dev.off()
print(fit.cv)

coef <-coef(fit.cox, s=c(fit.cox$lambda[6],0.1))
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoGene <- na.omit(lassoGene)
lassoGene=c("futime","fustat",lassoGene)

lassoSigExp=surSigExp[,lassoGene]


#save data
save(lassoSigExp,file = "Step11.03.lassoSigExp.rda")

##############Step 3: Built a Cox model#################
load("Step11.03.lassoSigExp.rda")
dim(lassoSigExp)
datainput<-lassoSigExp[,c(1:ncol(lassoSigExp))]

multiCox=coxph(Surv(futime, fustat) ~ ., data = datainput,x=T,y=T)
multiCox=step(multiCox, direction="both")
multiCoxSum=summary(multiCox)

save(multiCox,multiCoxSum,file = "Step11.04.multiCox.rda")


#output pm for the model
resultdf.cox=data.frame()
resultdf.cox=cbind(b=multiCoxSum$coefficients[,"coef"],
  se=multiCoxSum$coefficients[,"se(coef)"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])

row.names(resultdf.cox)=gsub("`","",row.names(resultdf.cox))

#savedata
#save(resultdf.cox,file = "Step11.04.resultdf.cox.rda")
saveRDS(resultdf.cox,file = "Step11.05.resultdf.cox.rds")
write.csv(resultdf.cox,file = "Step11.05.resultdf.cox.csv")


#rm(list = ls())
load("Step11.03.lassoSigExp.rda")
datainput<-lassoSigExp

load(file = "Step11.04.multiCox.rda")
riskScore=predict(multiCox, type="risk", newdata=datainput)
multiCoxSum=summary(multiCox)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`", "", coxGene)
outCol=c("futime", "fustat", coxGene)
riskOut=cbind(datainput[,outCol], riskScore)

#save data
save(riskOut,file = "Step11.07.riskOut.rda")

load("Step10.04.uniCoxdf.select.rda")
uniCoxdata<-uniCoxdf.select
rownames(uniCoxdata) <-uniCoxdata[,1]

uniRT<-uniCoxdata[,-1]
class(uniRT[,1])
df.unicox<-as.data.frame(lapply(uniRT,as.numeric))
rownames(df.unicox)<-rownames(uniCoxdata)

resultdf.unicox=df.unicox[coxGene,]
class(resultdf.unicox)
saveRDS(resultdf.unicox,file = "Step11.08.resultdf.unicox.rds")
