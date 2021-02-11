set.seed(12345)
source('functions.R')
if(!dir.exists('out'))dir.create('out')

alleles<-dnar::fillDown(unlist(read.csv('dat/monkeyData.csv',nrows=1,header=FALSE,stringsAsFactors=FALSE,check.names=FALSE))[-1])
dates<-unlist(read.csv('dat/monkeyData.csv',nrows=2,header=FALSE,stringsAsFactors=FALSE)[2,])[-1]
dat<-read.csv('dat/monkeyData.csv',stringsAsFactors=FALSE,row.names=1,skip=2,header=FALSE)
monk<-data.frame('rep'=unlist(dat),'virus'=rep(rownames(dat),ncol(dat)),'date'=rep(dates,each=nrow(dat)),'allele'=rep(alleles,each=nrow(dat)),stringsAsFactors=FALSE)
monk<-monk[monk$virus!='SIVden CD1',]
monk$siv<-ifelse(!grepl('SIV',monk$virus),'SIVcpz',ifelse(grepl('SIVgor',monk$virus),'SIVgor','Monkey SIV'))
monk$siv[grepl('Mock',monk$virus)]<-'Mock'
monk$species<-ifelse(monk$siv=='SIVcpz','cpz',ifelse(monk$siv=='SIVgor','gor',sub('SIV(...).*','\\1',monk$virus)))
monk$repNo0<-ifelse(monk$rep==0&!is.na(monk$rep),min(monk$rep[monk$rep>0],na.rm=TRUE),monk$rep)

alleles<-dnar::fillDown(unlist(read.csv('dat/gorillaData.csv',nrows=1,header=FALSE,stringsAsFactors=FALSE)))[-1]
dates<-unlist(read.csv('dat/gorillaData.csv',nrows=2,header=FALSE,stringsAsFactors=FALSE)[2,])[-1]
dat<-read.csv('dat/gorillaData.csv',row.names=1,skip=2,header=FALSE,stringsAsFactors=FALSE)
gor<-data.frame('rep'=unlist(dat),'virus'=rep(rownames(dat),ncol(dat)),'date'=rep(dates,each=nrow(dat)),'allele'=rep(alleles,each=nrow(dat)),stringsAsFactors=FALSE)
gor$siv<-ifelse(!grepl('SIV',gor$virus),'SIVcpz',ifelse(grepl('SIVgor',gor$virus),'SIVgor','Monkey SIV'))
gor$siv[grepl('Mock',gor$virus)]<-'Mock'
gor$species<-ifelse(gor$siv=='SIVcpz','cpz',ifelse(gor$siv=='SIVgor','gor',sub('SIV(...).*','\\1',gor$virus)))

mocks<-c(gor[grepl('Mock',gor$virus)&!is.na(gor$rep),'rep'],monk[monk$virus=='Mock'&!is.na(monk$rep),'rep'])


fitMonk<-dnar::withAs(xx=monk[monk$virus!='Mock'&!is.na(monk$rep),],fitModel(xx$virus,xx$allele,xx$species,xx$repNo0,xx$date,mocks,modAllele,chains=50,nIter=50000,thin=4,logFunc=logit))
fitGor<-dnar::withAs(xx=gor[!grepl('Mock',gor$virus)&!is.na(gor$rep)&!grepl('Bonobo',gor$allele),],fitModel(xx$virus,ifelse(xx$allele=='AHSM N15T','AHSM',xx$allele),xx$species,xx$rep,xx$date,mocks,modAlleleWithMods,ifelse(xx$allele=='AHSM N15T','N15T','Unmodified'),chains=50,nIter=50000,thin=4,logFunc=logit))
fitBono<-dnar::withAs(xx=gor[!grepl('Mock',gor$virus)&!is.na(gor$rep)&grepl('Bonobo|Human',gor$allele),],fitModel(xx$virus,xx$allele,xx$species,xx$rep,xx$date,mocks,modAllele,chains=50,nIter=50000,thin=4,logFunc=logit))

monkDiffs<-getDiffs(fitMonk,TRUE)
monkDiffs$group<-'Other'
gorDiffs<-getDiffs(fitGor,TRUE)
gorDiffs$group<-'Gorilla'
bonoDiffs<-getDiffs(fitBono,TRUE)
bonoDiffs$group<-'Bonobo'

out<-outPrep<-rbind(monkDiffs,gorDiffs,bonoDiffs)[,c('group','virus','all1','all2','p','mean','lower','upper')]
rownames(out)<-NULL
outPrep$censor<-out$mean>log(1000)|((out$upper-out$lower)>5&out$lower>0)
out[outPrep$censor,'upper']<-out[outPrep$censor,'lower']
out[outPrep$censor,'mean']<-out[outPrep$censor,'lower']
sigFig<-function(xx,digits=2)sub('\\.$','',formatC(signif(xx,digits=digits), digits=digits,format='fg',flag="#"))
out$upper<-sigFig(exp(out$upper))
out$lower<-sigFig(exp(out$lower))
out$mean<-sigFig(exp(out$mean))
out$upper[outPrep$upper>log(1000)]<-'>1000*'
out$upper[out$mean==out$lower&outPrep$censor]<-sprintf('>%s*',(out$lower[out$mean==out$lower&outPrep$censor]))
out$mean[out$mean==out$lower&outPrep$censor]<-sprintf('>%s*',out$lower[out$mean==out$lower&outPrep$censor])
out$p<-sigFig(out$p)
out[outPrep$p<.00001,'p']<-'<0.00001'
colnames(out)<-c('Group','Virus','AlleleA','AlleleB','p(fold change < 1)','Estimated fold change','Lower 95% CrI','Upper 95% CrI')
for(ii in unique(out$Group))write.csv(out[out$Group==ii,-1],sprintf('out/%sDiffs.csv',ii),row.names=FALSE)

mat<-as.matrix(fitGor$fit)
stats<-as.data.frame(t(apply(mat[,grep('[mM]od',colnames(mat))],2,crIMean)))
colnames(stats)<-c('lower','upper','estimate','p')
stats<-stats[rownames(stats)!='modSd[1]',]
stats$virus<-ifelse(rownames(stats)=='modMeans[1]','Overall',names(fitGor$virusId)[suppressWarnings(as.numeric(sub('virusMod\\[([0-9]+),.*','\\1',rownames(stats))))])
out<-stats[,c('virus','lower','upper','estimate','p')]
out[,c('lower','upper','estimate')]<-apply(exp(out[,c('lower','upper','estimate')]),2,sigFig)
out[,c('p')]<-apply(out[,c('p'),drop=FALSE],2,sigFig)
out[stats$p<.00001,'p']<-'<0.00001'
out$AlleleA<-'AHSM'
out$AlleleB<-'AHSM N15T'
out<-out[,c('virus','AlleleA','AlleleB','p','estimate','lower','upper')]
colnames(out)<-c('Virus','AlleleA','AlleleB','p(fold change < 1)','Estimated fold change','Lower 95% CrI','Upper 95% CrI')
write.csv(out,'out/N15T.csv',row.names=FALSE)

