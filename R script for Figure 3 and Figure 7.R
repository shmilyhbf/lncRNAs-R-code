rm(list = ls())
memory.limit(102400)

mainDir <-"C:/Bioinformatics/000 EPC rt/GSE53625"
subDir <-"Step10.uniCox"
outputDir <- file.path(mainDir, subDir)

if (!dir.exists(outputDir)){
  dir.create(outputDir)
} else {
  print("Dir already exists!")
}

setwd(mainDir)
setwd(outputDir)
getwd()


library(plyr)
library(dplyr)
suppressMessages(library(forestplot))



data.prepare= function(dataset=data.read){
  names(dataset)[1]<-"Variable"
  input.df<-dataset  %>% plyr::rename(c(HR="estimate",HR.95L="lowerCI",HR.95H="upperCI"))%>% 
    dplyr::select(Variable, estimate, lowerCI,upperCI,pvalue)
  effectindex="HR"
  if(effectindex=="HR") {
    effectname<-"HR with 95% CI"
  } else if(effectindex=="OR") {
    effectname<-"OR with 95% CI"
  } else if(effectindex=="RR" ) {
    effectname<-"RR with 95% CI"
  }
  Variable <- input.df[,"Variable"]
  
  hr <- sprintf("%.3f",input.df$"estimate")
  hrLow  <- sprintf("%.3f",input.df$"lowerCI")
  hrHigh <- sprintf("%.3f",input.df$"upperCI")
  effect.indicator <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  NA.check<-grepl("NA",effect.indicator,ignore.case=T)
  effect.indicator[NA.check==TRUE]<-NA
  Pvalue <- ifelse(input.df$pvalue<0.001, "<0.001", sprintf("%.3f", input.df$pvalue))
  labelmatrix<-cbind(c("Variable",Variable),
                     c(effectname,effect.indicator),
                     c("P value",Pvalue))
  return(list(input.df,labelmatrix))
}  

#  
my_forest= function(input.df=input.df,labelmatrix=labelmatrix,effectindex="HR"){  
  forst<-forestplot(labeltext=labelmatrix,
                    mean=c(NA,signif(input.df$"estimate",4)),
                    lower=c(NA,signif(input.df$"lowerCI",4)),
                    upper=c(NA,signif(input.df$"upperCI",4)),
                    #title = "Forestplot",
                    clip = c(-4, 4),
                    hrzl_lines = list("1" = gpar(lty=1, lwd=2,columns=1:4, col="black"),
                                      "2" = gpar(lty=2, lwd=2,columns=1:4, col="black"),
                                      "5" = gpar(lty=1, lwd=2,columns=1:4, col="black")),
                    graphwidth = unit(0.4,"npc"),
                    align=c("l","c","c","c"),
                    col=fpColors(box='#458B00',
                                 summary='#8B008B',
                                 lines = 'black',
                                 zero = '#7AC5CD'),
                    zero=1,
                    boxsize=0.5,
                    mar=unit(rep(1.25, times = 4), "cm"), 
                    colgap=unit(5,'mm'),
                    lineheight=unit(9,'mm'),line.margin = 0.08,
                    graph.pos=2,
                    txt_gp=fpTxtGp(label=gpar(cex=0.8), ticks=gpar(cex=0.8), xlab=gpar(cex = 0.8), title=gpar(cex = 0.8)),
                    xlab=effectindex,
                    xticks = c(-1.0,0, 1,2, 4.0),
                    lty.ci = "solid",
                    ci.vertices=T,
                    lwd.zero=1.5, lwd.ci=2,lwd.xaxis =1)
  return(forst)
}

#load data
data.unicox.clin <- read.csv("uniCox.clin.csv", header = T)

#data processing
data.prepare.clin<-data.prepare(dataset=data.unicox.clin)
#get last line
linelast <-nrow(data.prepare.clin[[2]])+1
print(linelast)
#
forst.uniCOX.clin<-my_forest(input.df=data.prepare.clin[[1]],
                             labelmatrix=data.prepare.clin[[2]],
                             effectindex="HR")
print(forst.uniCOX.clin)

save(data.unicox.clin,data.prepare.clin,forst.uniCOX.clin,
     file = "unicox.clin.Rda")
#save plot
library(ggplotify)
library(cowplot)
library(ggplot2)
forst.uniCOX.clin.grid <- grid2grob(print(forst.uniCOX.clin))

phname<-"forst.uniCOX.clin"
phtype<-"tiff"
savefile<-forst.uniCOX.clin.grid 
ggsave(filename = paste0(phname,".",phtype),savefile,width =8, 
       height = 6, dpi = 300, units = "in", device=phtype)




data.multicox.clin <- read.csv("multiCox.clin.csv", header = T)
data.prepare.clin<-data.prepare(dataset=data.multicox.clin)
linelast <-nrow(data.prepare.clin[[2]])+1
print(linelast)
forst.multicox.clin<-my_forest(input.df=data.prepare.clin[[1]],
                               labelmatrix=data.prepare.clin[[2]],
                               effectindex="HR")
print(forst.multicox.clin)
save(data.multicox.clin,data.prepare.clin,forst.multicox.clin,
     file = "multicox.clin.Rda")
forst.multicox.clin.grid  <- grid2grob(print(forst.multicox.clin))

phname<-"forst.multicox.clin"
phtype<-"tiff"
savefile<-forst.multicox.clin.grid
ggsave(filename = paste0(phname,".",phtype),savefile,width =8, 
       height = 5, dpi = 300, units = "in", device=phtype)


#merge plotting
pmerge.clin<-plot_grid(forst.uniCOX.clin.grid, 
                       forst.multicox.clin.grid,
                       rel_heights=c(1.8,1),
                       labels = "AUTO",ncol = 1)

phname<-"forst.merge.clin"
phtype<-"pdf"
savefile<-pmerge.clin
ggsave(filename = paste0(phname,".",phtype),savefile,width =6, 
       height = 9, dpi = 300, units = "in", device=phtype)

#
data.unicox.lncRNAs <- read.csv("unicox.lncRNAs.csv", header = T)
data.prepare.lncRNAs<-data.prepare(dataset=data.unicox.lncRNAs)
linelast <-nrow(data.prepare.lncRNAs[[2]])+1
print(linelast)
forst.uniCOX.lncRNAs<-my_forest(input.df=data.prepare.lncRNAs[[1]],
                                labelmatrix=data.prepare.lncRNAs[[2]],
                                effectindex="HR")
print(forst.uniCOX.lncRNAs)

save(data.unicox.lncRNAs,data.prepare.lncRNAs,forst.uniCOX.lncRNAs,
     file = "unicox.lncRNAs.Rda")
forst.uniCOX.lncRNAs.grid  <- grid2grob(print(forst.uniCOX.lncRNAs))

phname<-"forst.uniCOX.lncRNAs"
phtype<-"tiff"
savefile<-forst.uniCOX.lncRNAs.grid
ggsave(filename = paste0(phname,".",phtype),savefile,width =8, 
       height = 5, dpi = 300, units = "in", device=phtype)


data.multicox.lncRNAs <- read.csv("multicox.lncrna.csv", header = T)
data.prepare.lncRNAs<-data.prepare(dataset=data.multicox.lncRNAs)
linelast <-nrow(data.prepare.lncRNAs[[2]])+1
print(linelast)
forst.multicox.lncRNAs<-my_forest(input.df=data.prepare.lncRNAs[[1]],
                                  labelmatrix=data.prepare.lncRNAs[[2]],
                                  effectindex="HR")
print(forst.multicox.lncRNAs)

save(data.multicox.lncRNAs,data.prepare.lncRNAs,forst.multicox.lncRNAs,
     file = "multicox.lncRNAs.Rda")
forst.multicox.lncRNAs.grid  <- grid2grob(print(forst.multicox.lncRNAs))

phname<-"forst.multicox.lncRNAs"
phtype<-"tiff"
savefile<-forst.multicox.lncRNAs.grid
ggsave(filename = paste0(phname,".",phtype),savefile,width =8, 
       height = 5, dpi = 300, units = "in", device=phtype)

#merge plotting
pmerge.lncRNAs<-plot_grid(forst.uniCOX.lncRNAs.grid, 
                          forst.multicox.lncRNAs.grid,
                          rel_heights=c(1.8,1),
                          labels = "AUTO",ncol = 1)

phname<-"forst.merge.lncRNAs"
phtype<-"pdf"
savefile<-pmerge.lncRNAs
ggsave(filename = paste0(phname,".",phtype),savefile,width =8, 
       height = 6, dpi = 300, units = "in", device=phtype)

save(pmerge.lncRNAs,pmerge.clin,
     file = "pmerge.Rda")

