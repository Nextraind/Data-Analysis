

library(limma)
corFilter=0.4           #相关系数过滤标准
pvalueFilter=0.001       #p值过滤标准
#读取lncRNA表达文件,并对数据进行处理
rt = read.table("TCGA.normalize.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]
lncRNA=data
#读取免疫基因表达文件,并对数据进行处理
rt = read.table("TCGA.Cysteine.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
immuneGene=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
immuneGene=avereps(immuneGene)
immuneGene=immuneGene[rowMeans(immuneGene)>0.5,]
#相关性检验
outTab=data.frame()
for(i in row.names(lncRNA)){
  if(sd(lncRNA[i,])>0.5){
    for(j in row.names(immuneGene)){
      x=as.numeric(lncRNA[i,])
      y=as.numeric(immuneGene[j,])
      corT=cor.test(x,y)
      cor=corT$estimate
      pvalue=corT$p.value
      if((cor>corFilter) & (pvalue<pvalueFilter)){
        outTab=rbind(outTab,cbind(immuneGene=j,lncRNA=i,cor,pvalue,Regulation="postive"))
      }
      if((cor< -corFilter) & (pvalue<pvalueFilter)){
        outTab=rbind(outTab,cbind(immuneGene=j,lncRNA=i,cor,pvalue,Regulation="negative"))
      }
    }
  }
}

#输出相关性结果
write.table(file="corResult.txt",outTab,sep="\t",quote=F,row.names=F)

#输出免疫lncRNA表达量
immuneLncRNA=unique(as.vector(outTab[,"lncRNA"]))
immuneLncRNAexp=data[immuneLncRNA,]
immuneLncRNAexp=rbind(ID=colnames(immuneLncRNAexp), immuneLncRNAexp)
write.table(immuneLncRNAexp,file="ReCysteineExp.txt",sep="\t",quote=F,col.names=F)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)
corFilter=0.3           #相关系数过滤标准
pvalueFilter=0.001       #p值过滤标准
#读取lncRNA表达文件,并对数据进行处理
rt = read.table("TCGA.normalize.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]
lncRNA=data
#读取免疫基因表达文件,并对数据进行处理
rt = read.table("TCGA.Cysteine.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
immuneGene=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
immuneGene=avereps(immuneGene)
immuneGene=immuneGene[rowMeans(immuneGene)>0.5,]
#相关性检验
outTab=data.frame()
for(i in row.names(lncRNA)){
  if(sd(lncRNA[i,])>0.5){
    for(j in row.names(immuneGene)){
      x=as.numeric(lncRNA[i,])
      y=as.numeric(immuneGene[j,])
      corT=cor.test(x,y)
      cor=corT$estimate
      pvalue=corT$p.value
      if((cor>corFilter) & (pvalue<pvalueFilter)){
        outTab=rbind(outTab,cbind(immuneGene=j,lncRNA=i,cor,pvalue,Regulation="postive"))
      }
      if((cor< -corFilter) & (pvalue<pvalueFilter)){
        outTab=rbind(outTab,cbind(immuneGene=j,lncRNA=i,cor,pvalue,Regulation="negative"))
      }
    }
  }
}

#输出相关性结果
write.table(file="corResult.txt",outTab,sep="\t",quote=F,row.names=F)

#输出免疫lncRNA表达量
immuneLncRNA=unique(as.vector(outTab[,"lncRNA"]))
immuneLncRNAexp=data[immuneLncRNA,]
immuneLncRNAexp=rbind(ID=colnames(immuneLncRNAexp), immuneLncRNAexp)
write.table(immuneLncRNAexp,file="ReCysteineExp.txt",sep="\t",quote=F,col.names=F)


#引用包
library("survival")
library("caret")
library(glmnet)
library(survminer)
library(timeROC)

coxPfilter=0.05        #cox方法显著性过滤标准
setwd("C:\\biowolf\\NMF\\22.model")      #设置工作目录
rt=read.table("TCGA.expTime.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件
rt$futime[rt$futime<=0]=1
rt$futime=rt$futime/365

#对数据进行分组，构建模型
n=20   #分组的次数
for(i in 1:n){
  #############对数据进行分组#############
  inTrain=createDataPartition(y=rt[,3], p=0.7, list=F)
  train=rt[inTrain,]
  test=rt[-inTrain,]
  trainOut=cbind(id=row.names(train),train)
  testOut=cbind(id=row.names(test),test)
  
  #单因素cox分析
  outUniTab=data.frame()
  sigGenes=c("futime","fustat")
  for(i in colnames(train[,3:ncol(train)])){
    #cox分析
    cox <- coxph(Surv(futime, fustat) ~ train[,i], data = train)
    coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
    
    #保留显著性基因
    if(coxP<coxPfilter){
      sigGenes=c(sigGenes,i)
      outUniTab=rbind(outUniTab,
                      cbind(id=i,
                            HR=coxSummary$conf.int[,"exp(coef)"],
                            HR.95L=coxSummary$conf.int[,"lower .95"],
                            HR.95H=coxSummary$conf.int[,"upper .95"],
                            pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
      )
    }
  }
  uniSigExp=train[,sigGenes]
  uniSigExpOut=cbind(id=row.names(uniSigExp),uniSigExp)
  if(length(sigGenes)<5){next}
  
  #lasso回归
  x=as.matrix(uniSigExp[,c(3:ncol(uniSigExp))])
  y=data.matrix(Surv(uniSigExp$futime,uniSigExp$fustat))
  fit <- glmnet(x, y, family = "cox", maxit = 1000)
  cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
  coef <- coef(fit, s = cvfit$lambda.min)
  index <- which(coef != 0)
  actCoef <- coef[index]
  lassoGene=row.names(coef)[index]
  lassoSigExp=uniSigExp[,c("futime", "fustat", lassoGene)]
  lassoSigExpOut=cbind(id=row.names(lassoSigExp), lassoSigExp)
  geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
  if(nrow(geneCoef)<2){next}
  
  #############构建COX模型#############
  multiCox <- coxph(Surv(futime, fustat) ~ ., data = lassoSigExp)
  multiCox=step(multiCox,direction = "both")
  multiCoxSum=summary(multiCox)
  
  #输出模型公式
  outMultiTab=data.frame()
  outMultiTab=cbind(
    coef=multiCoxSum$coefficients[,"coef"],
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
  
  #输出train组风险文件
  riskScore=predict(multiCox,type="risk",newdata=train)          #利用train得到模型预测train样品风险
  coxGene=rownames(multiCoxSum$coefficients)
  coxGene=gsub("`","",coxGene)
  outCol=c("futime","fustat",coxGene)
  medianTrainRisk=median(riskScore)
  risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
  trainRiskOut=cbind(id=rownames(cbind(train[,outCol],riskScore,risk)),cbind(train[,outCol],riskScore,Risk=risk))
  
  #输出test组风险文件
  riskScoreTest=predict(multiCox,type="risk",newdata=test)        #利用train得到模型预测test样品风险
  riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
  testRiskOut=cbind(id=rownames(cbind(test[,outCol],riskScoreTest,riskTest)),cbind(test[,outCol],riskScore=riskScoreTest,Risk=riskTest))
  
  #输出GEO的风险值
  GEO=read.table("GEO.expTime.txt", header=T, sep="\t", check.names=F, row.names=1)
  geoScore=predict(multiCox, type="risk", newdata=GEO)
  geoRisk=as.vector(ifelse(geoScore>medianTrainRisk, "high", "low"))
  GEO=cbind(GEO[,outCol], riskScore=as.vector(geoScore), Risk=geoRisk)
  geoRiskOut=cbind(id=rownames(GEO), GEO)
  
  #生存差异pvalue	
  diff=survdiff(Surv(futime, fustat) ~Risk,data = trainRiskOut)
  pValue=1-pchisq(diff$chisq, df=1)
  diffTest=survdiff(Surv(futime, fustat) ~Risk,data = testRiskOut)
  pValueTest=1-pchisq(diffTest$chisq, df=1)
  diffGEO=survdiff(Surv(futime, fustat) ~Risk, data=GEO)
  pValueGEO=1-pchisq(diffGEO$chisq, df=1)
  
  #ROC曲线下面积
  predictTime=1    #预测时间
  roc=timeROC(T=train$futime, delta=train$fustat,
              marker=riskScore, cause=1,
              times=c(predictTime), ROC=TRUE)
  rocTest=timeROC(T=test$futime, delta=test$fustat,
                  marker=riskScoreTest, cause=1,
                  times=c(predictTime), ROC=TRUE)	
  
  if((pValue<0.01) & (roc$AUC[2]>0.58) & (pValueTest<0.1) & (rocTest$AUC[2]>0.63) & (pValueGEO<0.049)){
    #输出分组结果
    write.table(trainOut,file="data.train.txt",sep="\t",quote=F,row.names=F)
    write.table(testOut,file="data.test.txt",sep="\t",quote=F,row.names=F)
    #输出单因素结果
    write.table(outUniTab,file="uni.trainCox.txt",sep="\t",row.names=F,quote=F)
    write.table(uniSigExpOut,file="uni.SigExp.txt",sep="\t",row.names=F,quote=F)
    #lasso结果
    write.table(lassoSigExpOut,file="lasso.SigExp.txt",sep="\t",row.names=F,quote=F)
    pdf("lasso.lambda.pdf")
    plot(fit, xvar = "lambda", label = TRUE)
    dev.off()
    pdf("lasso.cvfit.pdf")
    plot(cvfit)
    abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
    dev.off()
    #输出多因素结果
    outMultiTab=outMultiTab[,1:2]
    write.table(outMultiTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)
    write.table(trainRiskOut,file="risk.TCGAtrain.txt",sep="\t",quote=F,row.names=F)
    write.table(testRiskOut,file="risk.TCGAtest.txt",sep="\t",quote=F,row.names=F)
    write.table(geoRiskOut,file="risk.GEO.txt",sep="\t",quote=F,row.names=F)
    #所有样品的风险值
    allRiskOut=rbind(trainRiskOut, testRiskOut)
    write.table(allRiskOut,file="risk.TCGAall.txt",sep="\t",quote=F,row.names=F)
    break
  }
}

library(pheatmap)     
setwd("D:\\BaiduNetdiskDownload\\159macrophage\\macrophage\\9.RISK_heatmap")      #???ù???Ŀ¼

bioRiskPlot=function(inputFile=null, riskScoreFile=null, survStatFile=null, heatmapFile=null){
  rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    #??ȡ?????ļ?
  rt=rt[order(rt$riskScore),]      #???շ??մ??ֶ???Ʒ????
  
  #???Ʒ???????
  riskClass=rt[,"Risk"]
  lowLength=length(riskClass[riskClass=="low"])
  highLength=length(riskClass[riskClass=="high"])
  lowMax=max(rt$riskScore[riskClass=="low"])
  line=rt[,"riskScore"]
  line[line>10]=10
  pdf(file=riskScoreFile, width=7, height=4)
  plot(line, type="p", pch=20,
       xlab="Patients (increasing risk socre)", ylab="Risk score",
       col=c(rep("#182E4F",lowLength),rep("#FF681F",highLength)) )
  abline(h=lowMax,v=lowLength,lty=2)
  legend("topleft", c("High risk", "Low Risk"),bty="n",pch=19,col=c("#FF681F","#182E4F"),cex=1.2)
  dev.off()
  
  #????????״̬ͼ
  color=as.vector(rt$fustat)
  color[color==1]="#FF681F"
  color[color==0]="#182E4F"
  pdf(file=survStatFile, width=7, height=4)
  plot(rt$futime, pch=19,
       xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
       col=color)
  legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("#FF681F","#182E4F"),cex=1.2)
  abline(v=lowLength,lty=2)
  dev.off()
  
  #???Ʒ?????ͼ
  rt1=rt[c(3:(ncol(rt)-2))]
  rt1=t(rt1)
  annotation=data.frame(type=rt[,ncol(rt)])
  rownames(annotation)=rownames(rt)
  pdf(file=heatmapFile, width=7, height=4)
  pheatmap(rt1, 
           annotation=annotation, 
           cluster_cols = F,
           cluster_rows = F,
           show_colnames = F,
           scale="row",
           color = colorRampPalette(c(rep("#182E4F",5), "white", rep("#FF681F",5)))(50),
           fontsize_col=3,
           fontsize=7,
           fontsize_row=8)
  dev.off()
}

#???ú????????Ʒ???????
bioRiskPlot(inputFile="risk.TCGAall.txt",
            riskScoreFile="allriskScore.pdf",
            survStatFile="allsurvStat.pdf",
            heatmapFile="allheatmap.pdf")

bioRiskPlot(inputFile="risk.TCGAtest.txt",
            riskScoreFile="testriskScore.pdf",
            survStatFile="testsurvStat.pdf",
            heatmapFile="testheatmap.pdf")

bioRiskPlot(inputFile="risk.TCGAtrain.txt",
            riskScoreFile="trainriskScore.pdf",
            survStatFile="trainsurvStat.pdf",
            heatmapFile="trainheatmap.pdf")
bioRiskPlot(inputFile="risk.GEO.txt",
            riskScoreFile="GEOriskScore.pdf",
            survStatFile="GEOsurvStat.pdf",
            heatmapFile="GEOheatmap.pdf")



#install.packages("survivalROC")


library(survivalROC)         #引用包
riskFile="risk.txt"          #风险文件
cliFile="clinical.txt"       #临床数据文件
setwd("D:\\biowolf\\irlncRNA\\18.cliROC")     #修改工作目录

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime","fustat","riskScore")]

#读取临床数据文件
cli=read.table(cliFile,sep="\t",header=T,check.names=F,row.names=1)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

#绘制risk打分的ROC曲线
rocCol=rainbow(ncol(rt)-2)
aucText=c()
pdf(file="cliROC.pdf", width=6, height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=risk$futime, status=risk$fustat, marker=risk$riskScore, predict.time=1, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
	xlab="False positive rate", ylab="True positive rate",
  	lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText,paste0("risk score"," (AUC=",sprintf("%.3f",roc$AUC),")"))
abline(0,1)

#绘制其他临床性状的ROC曲线
j=1
for(i in colnames(rt[,4:ncol(rt)])){
	roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,i], predict.time =1, method="KM")
	j=j+1
	lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 2)
	aucText=c(aucText,paste0(i," (AUC=",sprintf("%.3f",roc$AUC),")"))
}
legend("bottomright", aucText, lwd=2, bty="n", col=rocCol)
dev.off()



library(survival)
library(survminer)
setwd("C:\\biowolf\\NMF\\24.survival")       #设置工作目录

#绘制生存曲线函数
bioSurvival=function(inputFile=null, outFile=null){
  #读取输入文件
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  #比较高低风险组生存差异，得到显著性p值
  diff=survdiff(Surv(futime, fustat) ~Risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt)
  
  #绘制生存曲线
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=T,
                     pval=pValue,
                     pval.size=6,
                     legend.title="Risk",
                     legend.labs=c("High risk", "Low risk"),
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette=c("#182E4F", "#FF681F"),
                     risk.table=TRUE,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  pdf(file=outFile,onefile = FALSE,width = 9.5,height =5.5)
  print(surPlot)
  dev.off()
}

#调用函数，绘制生存曲线
bioSurvival(inputFile="risk.TCGAtrain.txt", outFile="surv.TCGAtrain.pdf")
bioSurvival(inputFile="risk.TCGAtest.txt", outFile="surv.TCGAtest.pdf")
bioSurvival(inputFile="risk.TCGAall.txt", outFile="surv.TCGAall.pdf")
bioSurvival(inputFile="risk.GEO.txt", outFile="surv.GEO.pdf")



#install.packages("survival")
#install.packages("survminer")
#install.packages("timeROC")


#引用包
library(survival)
library(survminer)
library(timeROC)
setwd("C:\\biowolf\\NMF\\25.ROC")      #设置工作目录

#定义绘制ROC曲线函数
bioROC=function(inputFile=null, rocFile=null){
	#读取输入文件
	rt=read.table(inputFile, header=T, sep="\t", check.names=F)
	#ROC曲线
	ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
	               marker=rt$riskScore,cause=1,
	               weighting='aalen',
	               times=c(1,3,5),ROC=TRUE)
	pdf(file=rocFile,width=5,height=5)
	plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
	plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
	plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
	legend('bottomright',
	        c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
	          paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
	          paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
	        col=c("green",'blue','red'),lwd=2,bty = 'n')
	dev.off()
}

#调用函数,绘制ROC曲线
bioROC(inputFile="risk.TCGAtrain.txt", rocFile="ROC.TCGAtrain.pdf")
bioROC(inputFile="risk.TCGAtest.txt", rocFile="ROC.TCGAtest.pdf")
bioROC(inputFile="risk.TCGAall.txt", rocFile="ROC.TCGAall.pdf")
bioROC(inputFile="risk.GEO.txt", rocFile="ROC.GEO.pdf")

#install.packages("survival")
#install.packages("survminer")
#install.packages("regplot")
#install.packages("rms")


#引用包
library(survival)
library(regplot)
library(rms)
library(survminer)

riskFile="risk.all.txt"     #风险文件
cliFile="clinical.txt"      #临床数据文件
setwd("C:\\biowolf\\Anoikis\\33.Nomo")     #设置工作目录

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
cli$Age=as.numeric(cli$Age)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1[,c("futime", "fustat", "risk")], cli)

#绘制列线图
res.cox=coxph(Surv(futime, fustat) ~ . , data = rt)
nom1=regplot(res.cox,
             plots = c("density", "boxes"),
             clickable=F,
             title="",
             points=TRUE,
             droplines=TRUE,
             observation=rt[2,],
             rank="sd",
             failtime = c(1,3,5),
             prfail = F)
dev.copy2pdf(file="Nomo.pdf", width=8, height=6, out.type="pdf")

#输出列线图的风险得分
nomoRisk=predict(res.cox, data=rt, type="risk")
rt=cbind(risk1, Nomogram=nomoRisk)
outTab=rbind(ID=colnames(rt), rt)
write.table(outTab, file="nomoRisk.txt", sep="\t", col.names=F, quote=F)

#校准曲线
pdf(file="calibration.pdf", width=5, height=5)
#1年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=1)
cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="green", sub=F)
#3年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=3)
cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=1.5, col="blue", sub=F, add=T)
#5年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=5)
cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=1.5, col="red", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("green","blue","red"), lwd=1.5, bty = 'n')
dev.off()

#累计风险曲线
nomoRisk=ifelse(rt$Nomogram>median(rt$Nomogram), "High", "Low")
fit=survfit(Surv(futime, fustat) ~ nomoRisk, data=rt)
gg=ggsurvplot(fit,
              conf.int = T,
              risk.table.col="strata",
              ggtheme = theme_bw(),
              #palette = "lancet",
              fun = 'cumhaz')
pdf(file="cumulative.pdf", width=5, height=4.8, onefile=F)
print(gg)
dev.off()



#install.packages('survival')


library(survival)       #???冒?
setwd("D:\\NO5\\1111\\24duliyuhou")     #???霉???目录

############????森??图????############
bioForest=function(coxFile=null, forestFile=null, forestCol=null){
	#??取?????募?
	rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
	gene <- rownames(rt)
	hr <- sprintf("%.3f",rt$"HR")
	hrLow  <- sprintf("%.3f",rt$"HR.95L")
	hrHigh <- sprintf("%.3f",rt$"HR.95H")
	Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
	pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
		
	#????图??
	pdf(file=forestFile, width=6.5, height=4.5)
	n <- nrow(rt)
	nRow <- n+1
	ylim <- c(1,nRow)
	layout(matrix(c(1,2),nc=2),width=c(3,2.5))
		
	#????森??图???叩??俅???息
	xlim = c(0,3)
	par(mar=c(4,2.5,2,1))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
	text.cex=0.8
	text(0,n:1,gene,adj=0,cex=text.cex)
	text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
	text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
		
	#?????冶叩?森??图
	par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
	xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
	arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=3)
	abline(v=1, col="black", lty=2, lwd=2)
	boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
	points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=2)
	axis(1)
	dev.off()
}
############????森??图????############

#??????立预??????????
indep=function(riskFile=null,cliFile=null,uniOutFile=null,multiOutFile=null,uniForest=null,multiForest=null){
	risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #??取?????募?
	cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)      #??取?俅??募?
	
	#???莺喜?
	sameSample=intersect(row.names(cli),row.names(risk))
	risk=risk[sameSample,]
	cli=cli[sameSample,]
	rt=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])
	
	#?????囟?立预??????
	uniTab=data.frame()
	for(i in colnames(rt[,3:ncol(rt)])){
		 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
		 coxSummary = summary(cox)
		 uniTab=rbind(uniTab,
		              cbind(id=i,
		              HR=coxSummary$conf.int[,"exp(coef)"],
		              HR.95L=coxSummary$conf.int[,"lower .95"],
		              HR.95H=coxSummary$conf.int[,"upper .95"],
		              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
		              )
	}
	write.table(uniTab,file=uniOutFile,sep="\t",row.names=F,quote=F)
	bioForest(coxFile=uniOutFile, forestFile=uniForest, forestCol="green")

uniTab$pvalue <- as.character(uniTab$pvalue)
uniTab$id<- as.character(uniTab$id)
uniTab$HR<- as.character(uniTab$HR)
uniTab$HR.95L<- as.character(uniTab$HR.95L)
uniTab$HR.95H<- as.character(uniTab$HR.95H)

	#?????囟?立预??????
	uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<1,]
	rt1=rt[,c("futime", "fustat", as.vector(uniTab[,"id"]))]
	multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
	multiCoxSum=summary(multiCox)
	multiTab=data.frame()
	multiTab=cbind(
	             HR=multiCoxSum$conf.int[,"exp(coef)"],
	             HR.95L=multiCoxSum$conf.int[,"lower .95"],
	             HR.95H=multiCoxSum$conf.int[,"upper .95"],
	             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
	multiTab=cbind(id=row.names(multiTab),multiTab)
	write.table(multiTab,file=multiOutFile,sep="\t",row.names=F,quote=F)
	bioForest(coxFile=multiOutFile, forestFile=multiForest, forestCol="red")
}


indep(riskFile="risk.all.txt",
      cliFile="clinical.txt",
      uniOutFile="uniCox.txt",
      multiOutFile="multiCox.txt",
      uniForest="uniForest.pdf",
      multiForest="multiForest.pdf")


#引用包
library(survival)
library(survminer)

clusterFile="cluster.txt"      #分型结果文件
cliFile="time.txt"             #生存数据文件
setwd("C:\\biowolf\\NMF\\17.clusterSur")      #设置工作目录

#读取输入文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365

#数据合并
sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], cluster[sameSample,,drop=F])

#生存差异统计
length=length(levels(factor(rt$Cluster)))
diff=survdiff(Surv(futime, fustat) ~ Cluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ Cluster, data = rt)
#print(surv_median(fit))

#绘制生存曲线
bioCol=c("#182E4F", "#FF681F")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="Cluster",
                   legend.labs=levels(factor(rt[,"Cluster"])),
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   ylab="Overall survival",
                   break.time.by = 1,
                   palette = bioCol,
                   #surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.3)
pdf(file="survival.pdf",onefile = FALSE,width=10,height=6)
print(surPlot)
dev.off()


#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

library(stringr)

str_replace_all(s$id,'\\.','-')

?order
#???ð?
library(limma)
library(estimate)
inputFile="TPM_symbol.txt"       #?????????ļ?
setwd("D:\\BaiduNetdiskDownload\\BBB\\26.MY_DAFEN_dafen")       #???ù???Ŀ¼

#??ȡ?ļ?,?????????ļ?????????
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#???????????ľ????ļ?
out=rbind(ID=colnames(data),data)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)

#????estimate??
filterCommonGenes(input.f="uniq.symbol.txt", 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct")

#????ÿ????Ʒ?Ĵ???
scores=read.table("estimateScore.gct", skip=2, header=T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.", "\\-", rownames(scores))
out=rbind(ID=colnames(scores), scores)
write.table(out, file="TMEscores.txt", sep="\t", quote=F, col.names=F)






#引用包
options(stringsAsFactors=F)
library(limma)
library(ggpubr)
library(reshape2)

riskFile="cluster.txt"          #风险文件
scoreFile="ssgseaOut.txt"         #ssgsea结果文件
setwd("D:\\BaiduNetdiskDownload\\BBB\\49.ssgsea")         
data=read.table(scoreFile,sep="\t",header=T,check.names=F,row.names=1)       #读取ssGSEA结果文件

#去除正常样品
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(data))
data=avereps(t(data))

risk=read.table(riskFile,header=T,sep="\t",row.names=1,check.names=F)
risk$v1=risk$Cluster
colnames(risk)=c("V1","Cluster")
#合并数据
sameSample=intersect(row.names(data),row.names(risk))
data=data[sameSample,]
risk=risk[sameSample,]
colnames(risk)
rt=cbind(data,risk[,c("Cluster","Cluster")])
rt=rt[,-(ncol(rt)-1)]

#对免疫细胞绘制箱线图
immCell=c("aDCs","B_cells","CD8+_T_cells","DCs","iDCs","Macrophages",
          "Mast_cells","Neutrophils","NK_cells","pDCs","T_helper_cells",
          "Tfh","Th1_cells","Th2_cells","TIL","Treg")
rt1=rt[,c(immCell,"Cluster.1")]
data=melt(rt1,id.vars=c("Cluster.1"))
colnames(data)=c("Cluster","Type","Score")
data$Cluster=factor(data$Cluster, levels=c("C1","C2"))
p=ggboxplot(data, x="Type", y="Score", color = "Cluster",
            ylab="Score",add = "none",xlab="",palette = c("#ffbb78","#9467bd") )
p=p+rotate_x_text(50)
pdf(file="immCell.boxplot.pdf",width=7,height=6)            
p+stat_compare_means(aes(group=Cluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()

#对免疫相关功能绘制箱线图
immFunction=c("APC_co_inhibition","APC_co_stimulation","CCR",
              "Check-point","Cytolytic_activity","HLA","Inflammation-promoting",
              "MHC_class_I","Parainflammation","T_cell_co-inhibition",
              "T_cell_co-stimulation","Type_I_IFN_Reponse","Type_II_IFN_Reponse")
rt1=rt[,c(immFunction,"Cluster.1")]
data=melt(rt1,id.vars=c("Cluster.1"))
colnames(data)=c("Cluster","Type","Score")
data$Cluster=factor(data$Cluster, levels=c("C1","C2"))
p=ggboxplot(data, x="Type", y="Score", color = "Cluster",
            ylab="Score",add = "none",xlab="",palette = c("#ffbb78","#9467bd") )
p=p+rotate_x_text(50)
pdf(file="immFunction.boxplot.pdf",width=7,height=6)            
p+stat_compare_means(aes(group=Cluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()



#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")


#???ð?
library(limma)
library(reshape2)
library(ggpubr)

riskFile="risk.all.txt"            #?????ļ?
immFile="CIBERSORT-Results.txt"     #????ϸ???????????ļ?
pFilter=0.05                        #????ϸ???????????Ĺ???????
setwd("D:\\BaiduNetdiskDownload\\BBB\\32.MY_JINRUN_plot")      #???ù???Ŀ¼

#??ȡ????ϸ???????ļ??????????ݽ???????
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<pFilter,]
data=as.matrix(immune[,1:(ncol(immune)-3)])

#ɾ????????Ʒ
group=sapply(strsplit(row.names(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
data=avereps(data)

#??ȡ?????ļ?
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(risk))
rt=cbind(data[sameSample,,drop=F], risk[sameSample,"risk",drop=F])
data=rt[order(rt$risk, decreasing=T),]

#??????ת????ggplot2?????ļ?
data=melt(data, id.vars=c("risk"))
colnames(data)=c("Risk", "Immune", "Expression")
#????????ͼ
group=levels(factor(data$Risk))
data$Risk=factor(data$Risk, levels=c("low","high"))
bioCol=c("#4DBBD5B2","#E64B35B2","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="Risk",
				  xlab="",
				  ylab="Fraction",
				  legend.title="Risk",
				  width=0.8,
				  palette=bioCol)+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", " ns")), label="p.signif")
#????ͼƬ
pdf(file="immune.diff.pdf", width=8, height=6)
print(boxplot)
dev.off()



#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")

rm(list = ls())
#???ð?
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
expFile="fpkm_symbol.txt"      #?????????ļ?
riskFile="cluster.txt"       #?????????ļ?
geneFile="gene.txt"       #???߼??????Ļ????ļ?
setwd("D:\\BaiduNetdiskDownload\\BBB\\50.MY_JCD")     #???ù???Ŀ¼

#??ȡ?????????ļ?,???????ݽ??д???
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#??ȡ?????ļ?
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data),as.vector(gene[,1]))
data=t(data[sameGene,])
data=log2(data+1)

#ɾ????????Ʒ
group=sapply(strsplit(row.names(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",row.names(data))
data=avereps(data)

#?ϲ?????
risk=read.table(riskFile, sep="\t", header=T, check.names=F, row.names=1)
risk$v1=risk$Cluster
sameSample=intersect(row.names(data),row.names(risk))
rt1=cbind(data[sameSample,],risk[sameSample,])
colnames(rt1)
rt1=rt1[,c(sameGene,"Cluster")]

#??ȡ?????????Ļ???
sigGene=c()
for(i in colnames(rt1)[1:(ncol(rt1)-1)]){
  if(sd(rt1[,i])<0.001){next}
  wilcoxTest=wilcox.test(rt1[,i] ~ rt1[,"Cluster"])
  pvalue=wilcoxTest$p.value
  if(wilcoxTest$p.value<0.05){
    sigGene=c(sigGene, i)
  }
}
sigGene=c(sigGene, "Cluster")
rt1=rt1[,sigGene]

#??????ת????ggplot2?????ļ?
rt1=melt(rt1,id.vars=c("Cluster"))
colnames(rt1)=c("Cluster","Gene","Expression")

#???ñȽ???
group=levels(factor(rt1$Cluster))
rt1$Cluster=factor(rt1$Cluster, levels=c("C1","C2"))
comp=combn(group,2)
my_comparisons=list()
for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}

#????????ͼ
boxplot=ggboxplot(rt1, x="Gene", y="Expression", fill="Cluster",
                  xlab="",
                  ylab="Gene expression",
                  legend.title="Cluster",
                  width=0.8,
                  palette = c("#ffbb78","#9467bd") )+
  rotate_x_text(50)+
  stat_compare_means(aes(group=Cluster),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")

#????ͼƬ
pdf(file="checkpoint.diff.pdf", width=10, height=5)
print(boxplot)
dev.off()

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("scatterplot3d")


#引用包
library(limma)
library(scatterplot3d)
setwd("C:\\Users\\lexb\\Desktop\\DRG\\25.PCA")      #设置工作目录

#定义PCA分析的函数
myPCA=function(input=null,output=null){
	#读取表达数据文件
	rt=read.table(input, header=T, sep="\t", check.names=F)
	rt=as.matrix(rt)
	rownames(rt)=rt[,1]
	exp=rt[,2:ncol(rt)]
	dimnames=list(rownames(exp),colnames(exp))
	data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
	data=avereps(data)
	data=data[rowMeans(data)>0.5,]
	
	#删除正常样品
	type=sapply(strsplit(colnames(data),"\\-"),"[",4)
	type=sapply(strsplit(type,""),"[",1)
	type=gsub("2","1",type)
	data=t(data[,type==0])
	rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(data))
		
	#读取risk风险文件
	risk=read.table("risk.all.txt", header=T, sep="\t", row.names=1, check.names=F)
	sameSample=intersect(rownames(data),rownames(risk))
	data=data[sameSample,]
	risk=risk[sameSample,]
	group=as.vector(risk[,"risk"])
		
	#PCA分析
	data.class <- rownames(data)
	data.pca <- prcomp(data, scale. = TRUE)
	pcaPredict=predict(data.pca)

	#绘制PCA图形
	color=ifelse(group=="low",4,2)
	pdf(file=output, width=7, height=7)
	par(oma=c(1,1,2.5,1))
	s3d=scatterplot3d(pcaPredict[,1:3], pch = 16, color=color, angle=35)
	legend("top", legend = c("Low risk","High risk"),pch = 16, inset = -0.2, box.col="white", xpd = TRUE, horiz = TRUE,col=c(4,2))
	dev.off()
}

######绘制所有基因的PCA图，将04节课symbol.txt复制到当前目录
myPCA(input="symbol.txt", output="PCA.allGene.pdf")
######绘制双硫死亡基因的PCA图，将09节课disulfidptosisExp.txt复制到当前目录
myPCA(input="disulfidptosisExp.txt", output="PCA.disulfidptosisGene.pdf")
######绘制双硫死亡lncRNA的PCA图，将09节课disulfidptosisLncExp.txt复制到当前目录
myPCA(input="disulfidptosisLncExp.txt", output="PCA.disulfidptosisLncRNA.pdf")


######读取风险文件,绘制模型lncRNA的PCA图，将14节课risk.all.txt复制到当前目录
risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F, row.names=1)
data=risk[,3:(ncol(risk)-2)]
group=as.vector(risk[,"risk"])
		
#PCA分析
data.class <- rownames(data)
data.pca <- prcomp(data, scale. = TRUE)
pcaPredict=predict(data.pca)

#可视化
color=ifelse(group=="low",4,2)
pdf(file="PCA.riskLnc.pdf", width=6.5, height=6)
par(oma=c(1,1,2.5,1))
s3d=scatterplot3d(pcaPredict[,1:3], pch = 16, color=color, angle=35)
legend("top", legend = c("Low risk","High risk"),pch = 16, inset = -0.2, box.col="white", xpd = TRUE, horiz = TRUE,col=c(4,2))
dev.off()

#if (!require("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install("maftools")


library(maftools)      #引用包
setwd("C:\\Users\\lexb\\Desktop\\macrophage\\35.maftools")     #设置工作目录

#读取风险文件，得到瀑布图的注释信息
risk=read.table("risk.TCGA.txt", header=T, sep="\t", check.names=F)
outTab=risk[,c(1, ncol(risk))]
colnames(outTab)=c("Tumor_Sample_Barcode", "Risk")
write.table(outTab, file="ann.txt", sep="\t", quote=F, row.names=F)

#读取基因突变的文件
geneNum=20     #设置展示基因的数目
geneMut=read.table("geneMut.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneMut)[1:geneNum]

#定义注释的颜色
ann_colors=list()
col=c("#0088FF", "#FF5555")
names(col)=c("low", "high")
ann_colors[["Risk"]]=col

#绘制低风险组瀑布图
pdf(file="low.pdf", width=6, height=6)
maf=read.maf(maf="low.maf", clinicalData="ann.txt")
oncoplot(maf=maf, clinicalFeatures="Risk", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()

#绘制高风险组瀑布图
pdf(file="high.pdf", width=6, height=6)
maf=read.maf(maf="high.maf", clinicalData="ann.txt")
oncoplot(maf=maf, clinicalFeatures="Risk", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()



#install.packages("survival")
#install.packages("survminer")


#???冒?
library(survival)
library(survminer)
tmbFile="TMB.txt"            #????突?涓???募?
riskFile="risk.all.txt"      #?????募?
setwd("D:\\NO66\\33TMBsur")       #?薷墓???目录

#??取?????募?
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #??取?????募?
tmb=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)      #??取TMB?????募?

#?喜?????
sameSample=intersect(row.names(tmb), row.names(risk))
tmb=tmb[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(risk, tmb)

#??取????突?涓??????cutoff
res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("TMB"))
cutoff=as.numeric(res.cut$cutpoint[1])
tmbType=ifelse(data[,"TMB"]<=cutoff, "L-TMB", "H-TMB")
scoreType=ifelse(data$risk=="low", "low risk", "high risk")
mergeType=paste0(tmbType, "+", scoreType)

#????????????????
bioSurvival=function(surData=null, outFile=null){
	diff=survdiff(Surv(futime, fustat) ~ group, data=surData)
	length=length(levels(factor(surData[,"group"])))
	pValue=1-pchisq(diff$chisq, df=length-1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ group, data = surData)
	#print(surv_median(fit))
	
	#????????????
	bioCol=c("#FF0000","#0066FF","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
	bioCol=bioCol[1:length]
	surPlot=ggsurvplot(fit, 
			           data=surData,
			           conf.int=F,
			           pval=pValue,
			           pval.size=6,
			           legend.title="",
			           legend.labs=levels(factor(surData[,"group"])),
			           font.legend=10,
			           legend = c(0.8, 0.8),
			           xlab="Time(years)",
			           break.time.by = 1,
			           palette = bioCol,
			           surv.median.line = "hv",
			           risk.table=F,
			           cumevents=F,
			           risk.table.height=.25)
	#????图??
	pdf(file=outFile, onefile = FALSE, width=5.5, height=4.8)
	print(surPlot)
	dev.off()
}

#????????突?涓?傻?????????
data$group=tmbType
bioSurvival(surData=data, outFile="TMB.survival.pdf")

#????????突?涓??联?细叩头??盏?????????
data$group=mergeType
bioSurvival(surData=data, outFile="TMB-risk.survival.pdf")



#install.packages("ggpubr")


#引用包
library(ggpubr)
library(reshape2)

tmbFile="TMB.txt"             #肿瘤突变负荷文件
riskFile="risk.all.txt"       #风险文件
cluFile="cluster.txt"     #基因分型文件
setwd("D:\\NO66\\32TMBdiff")       #设置工作目录

#读取输入文件
tmb=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)       #读取TMB数据文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)     #读取风险文件
clu=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)       #读取基因分型类文件

#合并数据
tmb=as.matrix(tmb)
tmb[tmb>quantile(tmb,0.975)]=quantile(tmb,0.975)
sameSample=intersect(row.names(tmb), row.names(risk))
tmb=tmb[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
#rownames(clu)=gsub("(.*?)\\_(.*?)", "\\2", rownames(clu))
clu=clu[sameSample,,drop=F]
data=cbind(risk, tmb, clu)
data=data[,c("riskScore", "risk", "Cluster", "TMB")]

#设置比较组
data$risk=factor(data$risk, levels=c("low", "high"))
risk=levels(factor(data$risk))
comp=combn(risk, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#设置颜色
bioCol=c("#4DBBD5B2","#E64B35B2","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(risk)]

#绘制箱线图
boxplot=ggboxplot(data, x="risk", y="TMB", fill="risk",
		          xlab="",
		          ylab="Tumor Burden Mutation",
		          legend.title="Risk",
		          palette = bioCol )+ 
	    stat_compare_means(comparisons = my_comparisons)
pdf(file="boxplot.pdf",width=5,height=4.5)
print(boxplot)
dev.off()

#相关性图形
length=length(levels(factor(data$Cluster)))
bioCol=c("#E18727","#21854F","#CC66CC","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
p1=ggplot(data, aes(riskScore, TMB)) + 
		  xlab("Risk score")+ylab("Tumor Burden Mutation")+
		  geom_point(aes(colour=Cluster))+
		  scale_color_manual(values=bioCol[1:length])+ 
		  geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
		  stat_cor(method = 'spearman', aes(x =riskScore, y =TMB))
#相关性图形
pdf(file="cor.pdf", width=6, height=4.5)
print(p1)
dev.off()


#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")


#引用包
library(limma)
library(ggplot2)
library(ggpubr)

pFilter=0.001                      #pvalue过滤条件
riskFile="risk.TCGA.txt"           #风险文件
drugFile="DrugPredictions.csv"     #药物敏感性文件
setwd("C:\\Users\\lexb\\Desktop\\macrophage\\38.boxplot")     #设置工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取药物敏感性文件
senstivity=read.csv(drugFile, header=T, sep=",", check.names=F, row.names=1)
colnames(senstivity)=gsub("(.*)\\_(\\d+)", "\\1", colnames(senstivity))

#数据合并
sameSample=intersect(row.names(risk), row.names(senstivity))
risk=risk[sameSample, "Risk",drop=F]
senstivity=senstivity[sameSample,,drop=F]
rt=cbind(risk, senstivity)

#设置比较组
rt$Risk=factor(rt$Risk, levels=c("low", "high"))
type=levels(factor(rt[,"Risk"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#对药物进行循环, 绘制箱线图
for(drug in colnames(rt)[2:ncol(rt)]){
	rt1=rt[,c(drug, "Risk")]
	colnames(rt1)=c("Drug", "Risk")
	rt1=na.omit(rt1)
	rt1$Drug=log2(rt1$Drug+1)
	#差异分析
	test=wilcox.test(Drug ~ Risk, data=rt1)
	diffPvalue=test$p.value
	#对满足条件的药物绘制箱线图
	if(diffPvalue<pFilter){
		boxplot=ggboxplot(rt1, x="Risk", y="Drug", fill="Risk",
					      xlab="Risk",
					      ylab=paste0(drug, " senstivity"),
					      legend.title="Risk",
					      palette=c("green", "red")
					     )+ 
			stat_compare_means(comparisons=my_comparisons)
		#输出图形
		pdf(file=paste0("drugSenstivity.", drug, ".pdf"), width=5, height=4.5)
		print(boxplot)
		dev.off()
	}
}




#引用包
library(survival)
library(survminer)

riskFile="risk.all.txt"      #风险文件
cliFile="clinical.txt"       #临床数据文件
setwd("C:\\biowolf\\m6aDrug\\23.cliGroupSur")       #设置工作目录

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取临床文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cliName=colnames(cli)[1]

#数据合并
sameSample=intersect(row.names(cli), row.names(risk))
risk=risk[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(futime=risk[,1], fustat=risk[,2], cli, risk[,"risk",drop=F])
colnames(rt)=c("futime", "fustat", "clinical", "Risk")
tab=table(rt[,"clinical"])
tab=tab[tab!=0]

#对每个临床信息里面的每个分组进行循环
for(j in names(tab)){
	rt1=rt[(rt[,"clinical"]==j),]
	tab1=table(rt1[,"Risk"])
	tab1=tab1[tab1!=0]
	labels=names(tab1)
	if(length(labels)!=2){next}
	if((cliName=="age") | (cliName=="Age") | (cliName=="AGE")){
		titleName=paste0("age",j)
	}
	
	#计算高低风险组差异pvalue
	diff=survdiff(Surv(futime, fustat) ~Risk,data = rt1)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	
	#绘制生存曲线
	fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt1)
	surPlot=ggsurvplot(fit, 
			           data=rt1,
			           conf.int=F,
			           pval=pValue,
			           pval.size=6,
			           title=paste0("Patients with ",j),
			           legend.title="Risk",
			           legend.labs=labels,
			           font.legend=12,
			           xlab="Time(years)",
			           break.time.by = 1,
			           palette=c("red", "blue"),
			           risk.table=F,
			       	   risk.table.title="",
			           risk.table.col = "strata",
			           risk.table.height=.25)
	
	#输出图片
	j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
	pdf(file=paste0("survival.",cliName,"_",j,".pdf"), onefile = FALSE,
			       width = 6,        #图片的宽度
			       height =5)        #图片的高度
	print(surPlot)
	dev.off()
}



####首先认识数据结构####
###10XGenomics 的filtered_feature_bc_matrix
##https://www.ncbi.nlm.n[ih.gov/geo/query/acc.cgi?acc=GSE134520
###BD RSEC_MolsPerCell.csv
###SMART-seq RPKM.txt

rm(list = ls())#清空环境
#install.packages("Seurat")
####首先我们先处理10X Genomics的数据
library(Seurat)##version 4.0.5
library(dplyr)
library(future)
library(future.apply)
library(DoubletFinder)
plan("multicore", workers = 3) ###compute cores
options(future.globals.maxSize = 10000 * 1024^2)

####SRR780####
SRR780<-Read10X(data.dir = "CID4530N/")#膀胱癌的单细胞数据，使用read10X读取，其实是一个巨大的矩阵
#有些是read10X读取读取出来的结果其实就是一个表达矩阵max，如果直接是一个矩阵的话可以直接使用fread函数进行读取
SRR780_object<- CreateSeuratObject(counts = SRR780, project = "SRR780", min.cells = 0, min.features = 200)#创建Seurat对象,min.cells一个基因在多少细胞中表达。min.features,是一个细胞中有多少的基因。
x=SRR780_object@meta.data#创建出来的Seurat对象里的meta.data部分储存了
#行名对应是每一个细胞
#ordient,是那组的细胞
#ncount_Feature测到了多少count数
#nFeature_RNA是对应的基因数
SRR780_object[["percent.mt"]] <- PercentageFeatureSet(SRR780_object, pattern = "^MT-")#找出来线粒体基因。记住这个函数
hist(SRR780_object[["percent.mt"]]$percent.mt)
VlnPlot(SRR780_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(SRR780_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SRR780_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))#这两个图是判断一下COUNT数和线粒体基因数的对应情况，
#和count数和基因数的对应情况4409 sample
SRR780_object
##select the cells with 300 genes at least and 4000 at most, the percent of mitochondrion genes is less than 10%
SRR780_val<- subset(SRR780_object, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)#基因大于200并且小于7000且线粒体数要小于20%
#过滤完成4325 samples
#SRR780_val@assays#储存表达矩阵的地方
#SRR780_val@assays[["RNA"]]#储存原始count的地方
#SRR780_val@assays[["RNA"]]@counts#储存count的地方
#SRR780_val@assays[["RNA"]]@data#在count值标准化后储存的地方
#SRR780_val@meta.data#储存每个细胞信息的地方。比如线粒体，count数 对应基因，或者细胞注释
dim(SRR780_val@meta.data)
rownames(SRR780_val@meta.data)
colname<-paste("SRR780_",colnames(SRR780_val),sep="")#把每个细胞名字前面加上对应的文件或者样本名字
colname
SRR780_val<-RenameCells(object = SRR780_val,colname)#对细胞重新命名，防止多样本的时候细胞bulk重复
#rm(list = c("BC2","bc2_object"))
SRR780_val@meta.data[1:5,]#对SRR780_val重新命名后对应的meta.data也进行了重新命名

SRR780_val$tissue_type<-"tormal"#同样的也可以新增一列进行注释
SRR780_val@meta.data[1:5,]

SRR780_val <- NormalizeData(SRR780_val, normalization.method = "LogNormalize", scale.factor = 10000)#进行标准化
SRR780_val <- FindVariableFeatures(SRR780_val, selection.method = "vst", nfeatures = 2000)#找出2000个不同的高变基因
SRR780_val@assays$RNA@var.features#找到2000个差异基因后存储到SRR780_val@assays$RNA@var.features之中
length(SRR780_val@assays$RNA@var.features)
###scaling the data###
#随后进行降维处理
SRR780_val <- ScaleData(SRR780_val,vars.to.regress = c("percent.mt"))#将线粒体基因过滤掉
###perform linear dimensional reduction###
SRR780_val <- RunPCA(SRR780_val, features = VariableFeatures(object = SRR780_val))#根据前面2000个高变基因进行降维
print(SRR780_val[["pca"]], dims = 1:5, nfeatures = 5)
#这里的主成分分析还是储存到了mata.data中
dev.off()
VizDimLoadings(SRR780_val, dims = 1:2, reduction = "pca")


ElbowPlot(SRR780_val,ndims = 50)

####cluster the cells###
SRR780_val <- FindNeighbors(SRR780_val, dims = 1:30)#选择主成分，当拐点很小时候就是最好的拐点
####Run non-linear dimensional reduction (UMAP/tSNE)
##UMAP和tSNE进行可视化聚集降维到2维空间中
SRR780_val <- RunUMAP(SRR780_val, dims = 1:30)
SRR780_val<-RunTSNE(SRR780_val,dims=1:30)

SRR780_val <- FindClusters(SRR780_val, resolution = 0.6)###找出合适的群更细的可以调整resol
#分群信息可以存在mate.data中
DimPlot(SRR780_val, reduction = "umap",label = T)
#?DimPlot
DimPlot(SRR780_val, reduction = "tsne",label = T)

table(Idents(SRR780_val))

####用singleR软件来注释细胞####  有空闲时跑
#devtools::install_github("dviraran/SingleR")
BiocManager::install("SingleR")
BiocManager::install("celldex",force = TRUE)
library(SingleR)
library(celldex)
sc_data<-as.matrix(SRR780_val@assays$RNA@counts)
load("HumanPrimaryCellAtlas_hpca.se_human.RData")
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se$label.main
clusters <- SRR780_val@meta.data$seurat_clusters
colnames(x)
cellpred <- SingleR(test = sc_data, ref = hpca.se, labels = hpca.se$label.main,method = "cluster", clusters = clusters, 
                   assay.type.test = "logcounts", assay.type.ref = "logcounts")
###制作细胞类型的注释文件
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = FALSE)
###保存一下
write.csv(celltype,"celltype_singleR.csv",row.names = FALSE)
##把singler的注释写到metadata中 有两种方法
###方法一
x=SRR780_val@meta.data
SRR780_val@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  SRR780_val@meta.data[which(SRR780_val@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
###因为我把singler的注释加载到metadata中时候，命名的名字叫celltype，所以画图时候，group.by="celltype"
DimPlot(SRR780_val, group.by="celltype", label=T, label.size=5, reduction='tsne')
DimPlot(SRR780_val, group.by="celltype", label=T, label.size=5, reduction='umap')

##鉴定结果展示
p1 = DimPlot(scRNA1, group.by="celltype", label=T, label.size=5, reduction='tsne')
p1
p2 = DimPlot(scRNA1, , label=T, label.size=5, reduction='umap')
p2
p3 = plotc <- p1+p2+ plot_layout(guides = 'collect')
p3 
ggsave("tSNE_celltype.pdf", p1, width=7 ,height=6)
ggsave("UMAP_celltype.pdf", p2, width=7 ,height=6)
ggsave("celltype.pdf", p3, width=10 ,height=5)
ggsave("celltype.png", p3, width=10 ,height=5)
)
######biomarkers###

FeaturePlot(SRR780_val,features = c("EIF4EBP1","USP30","RPL27A",	"CELSR2",	"KLHDC7B",	"NUDCD1"),reduction = "tsne",label = T)




