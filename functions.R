library('rstan')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

plotFit<-function(fit,means){
  #virLookup, alleleLookup
  mat<-as.matrix(fit$fit)
  diffs<-do.call(rbind,lapply(fit$virusId,function(vir){
    do.call(rbind,lapply(fit$alleleId,function(xx){
      out<-data.frame(do.call(rbind,lapply(fit$alleleId,function(yy){
        crIMean<-function(xx,crI=c(.025,.975))c(quantile(xx,crI),mean(xx<0))
        if(xx==yy)return(c(0,0,.5))
        return(crIMean(mat[,sprintf('alleleVirus[%d,%d]',vir,xx)]-mat[,sprintf('alleleVirus[%d,%d]',vir,yy)]))
      })))
      out$allele1<-xx
      out$allele2<-fit$alleleId
      out$virus<-vir
      return(out)
    }))
  }))
  diffs<-diffs[diffs$allele1<diffs$allele2,]
  diffs$vir<-names(fit$virusId)[diffs$vir]
  diffs$all1<-names(fit$alleleId)[diffs$allele1]
  diffs$all2<-names(fit$alleleId)[diffs$allele2]
  sig01<-diffs[diffs$V3<.005|diffs$V3>.995,]
  par(mar=c(5,7,.1,.6))
  virOrder<-rev(unlist(lapply(virLookup,names)))
  tmp<-means[!grepl('Mock',rownames(means)),]
  tmp2<-tmp[virOrder,alleleLookup]
  rownames(tmp2)<-do.call(c,unname(virLookup))[rownames(tmp2)]
  sig01$row<-structure(1:nrow(tmp2),.Names=names(rownames(tmp2)))[sig01$vir]
  sig01$col1<-structure(1:ncol(tmp2),.Names=colnames(tmp2))[sig01$all1]
  sig01$rowOffset<-ave(sig01$row,sig01$row,FUN=function(xx)xx+seq(-length(xx)/50,length(xx)/65,length.out=length(xx)))+ave(sig01$col1,sig01$row,FUN=function(xx)cumsum(c(FALSE,xx[-1]!=xx[-length(xx)])*.02))
  sig01$col2<-structure(1:ncol(tmp2),.Names=colnames(tmp2))[sig01$all2]
  plotFunc(tmp2)
  abline(h=cumsum(rev(unlist(sapply(virLookup,length))))[-length(virLookup)]+.5,lwd=2)
  segments(sig01$col1,sig01$rowOffset,sig01$col2,sig01$rowOffset,col='#00000099')
  points(c(sig01$col1,sig01$col2),c(sig01$rowOffset,sig01$rowOffset),col='#00000099',pch=21,bg='#00000099',cex=.3)
}

plotAlleles<-function(fit,xRange=quantile(allGroup,c(.005,.995)),isT=TRUE,expFunc=exp,sdMult=FALSE){
  mat<-as.matrix(fit$fit)
  par(mar=c(0,0,0,0))
  nAllele<-max(fit$alleleId)
  nGroup<-max(fit$groupId)
  layout(cbind(0,rbind(matrix(1:(nAllele*nGroup),nrow=nGroup),0)),height=c(rep(1,nGroup),.5),width=c(.5,rep(1,nAllele)))
  allGroup<-mat[,grepl('alleleGroup',colnames(mat))]
  if(sdMult){
    allSd<-mat[,grepl('alleleSd',colnames(mat))]*exp(mat[,grepl('sdMult',colnames(mat))])
    colnames(allSd)<-sub('sdMult','alleleSd',colnames(allSd))
  } else{
    allSd<-mat[,grepl('alleleSd',colnames(mat))]
  }
  if(isT)plusSd<-lapply(structure(colnames(allGroup),.Names=colnames(allGroup)),function(xx)density(LaplacesDemon::rst(nrow(allGroup),allGroup[,xx],allSd[,sub('alleleGroup\\[([0-9]+),.*','alleleSd[\\1]',xx)],fit$tNu)))
  else plusSd<-lapply(structure(colnames(allGroup),.Names=colnames(allGroup)),function(xx)density(rnorm(nrow(allGroup),allGroup[,xx],allSd[,sub('alleleGroup\\[([0-9]+),.*','alleleSd[\\1]',xx)])))
  denses<-apply(allGroup,2,density)
  comboDenses<-apply(mat[,grepl('alleleMeans',colnames(mat))],2,density)
  yLim<-c(0,max(sapply(denses,function(xx)max(xx$y))))
  print(expFunc)
  for(jj in fit$alleleId){
    for(ii in fit$groupId){
      thisCol<-sprintf('alleleGroup[%d,%d]',jj,ii)
      thisDat<-denses[[thisCol]]
      thisVirusIds<-fit$virusId[names(fit$virusGroup)[fit$virusGroup==ii]]
      thisVirus<-mat[,sprintf('alleleVirus[%d,%d]',thisVirusIds,jj),drop=FALSE]
      virDense<-apply(thisVirus,2,density)
      plot(1,1,type='n',xaxt='n',yaxt='n',xlab='',ylab='',xlim=expFunc(xRange),ylim=yLim,log='x',yaxs='i')
      #plot(1,1,type='n',xaxt='n',yaxt='n',xlab='',ylab='',xlim=expFunc(xRange),ylim=c(-yLim[2],yLim[2]),log='x')
      #polygon(expFunc(thisDat$x),thisDat$y,col='#0000FF33')
      #lapply(virDense,function(xx)polygon(expFunc(xx$x),-xx$y,col='#00000022',border='#00000011'))
      lapply(virDense,function(xx)polygon(expFunc(xx$x),xx$y,col='#00000011',border='#00000022'))
      polygon(expFunc(plusSd[[thisCol]]$x),plusSd[[thisCol]]$y,col='#FF000099')
      polygon(expFunc(comboDenses[[jj]]$x),comboDenses[[jj]]$y,lty=2)
      if(jj==1)mtext(names(fit$groupId)[ii],2,xpd=NA,line=.5)
      if(jj==1&ii==ceiling(nGroup/2))mtext('Probability density',2,xpd=NA,cex=1.2,line=3)
      if(ii==1)text(dnar::convertLineToUser(-2,2),mean(par('usr')[3:4]),names(fit$alleleId)[jj],xpd=NA)
      if(ii==nGroup)dnar::logAxis(1,axisVals=c(-10,-8,-6,-4,-2))
      if(ii==nGroup&jj==4)mtext('Predicted proportion positive cells',1,xpd=NA,line=2.25)
    }
  }
}

logit<-function(xx)log(xx)-log(1-xx)

fitModel<-function(virus,allele,group,rep,experiment,mock,mod,modification=rep(1,length(rep)),baseModify=dnar::mostAbundant(modification),chains=50,nIter=5000,tNu=3,thin=2,baseAllele='Human',logFunc=log,...){
  virusId<-structure(1:length(unique(virus)),.Names=sort(unique(virus)))
  alleleId<-structure(1:length(unique(allele)),.Names=c(baseAllele,sort(unique(allele[allele!=baseAllele]))))
  groupId<-structure(1:length(unique(group)),.Names=sort(unique(group)))
  modId<-structure(1:length(unique(modification)),.Names=c(baseModify,sort(unique(modification[modification!=baseModify]))))
  experimentId<-structure(1:length(unique(experiment)),.Names=sort(unique(experiment)))
  groupLookup<-unique(cbind(virus,group))
  if(any(table(groupLookup[,'virus'])>1))stop('Problem assigning group')
  rownames(groupLookup)<-groupLookup[,'virus']
  virusGroup<-groupId[groupLookup[names(virusId),'group']]
  names(virusGroup)<-names(virusId)
  rep[rep==0]<-min(rep[rep>0])/2
  dat<-list(
    nVirus=max(virusId),
    nAllele=max(alleleId),
    nRep=length(rep),
    nGroup=max(groupId),
    nExperiment=max(experimentId),
    rep=logFunc(rep/100),
    allele=alleleId[allele],
    virus=virusId[virus],
    virusGroup=virusGroup,
    experiment=experimentId[experiment],
    backMean=mean((mock/100)),
    backSd=sd((mock/100)),
    nMock=length(mock),
    tNu=tNu,
    mock=mock/100,
    modification=modId[modification],
    nModification=max(modId)
  )
  fit <- sampling(mod, data = dat, iter=nIter, chains=chains,thin=thin,control=list(adapt_delta=.99,max_treedepth=15),...)
  return(list('fit'=fit,'dat'=dat,'virusId'=virusId,'groupId'=groupId,'virusGroup'=virusGroup,'alleleId'=alleleId,'mod'=mod,'tNu'=tNu,'modId'=modId))
}

crIMean<-function(xx,crI=c(.025,.975))c(quantile(xx,crI),mean(xx),mean(xx<0))

plotFunc<-function(tmp){
  breaks<-seq(0,max(means),length.out=101)
  #cols<-hcl.colors(100, "YlOrRd", rev = TRUE)
  cols<-colorRampPalette(c('white','#f7c252','#eb5500','#7d0025'))(100)
  image(1:ncol(tmp),1:nrow(tmp),t(tmp),xaxt='n',yaxt='n',ylab='',xlab='',breaks=breaks,col=cols)
  dnar::slantAxis(1,1:ncol(tmp),colnames(tmp))
  axis(2,1:nrow(tmp),rownames(tmp),las=1)
  par(lheight=.7)
  dnar::insetScale(breaks,cols,main='Percentage of\ncells infected')
  abline(h=2:nrow(means)-.5,v=2:nrow(means)-.5,col='#FFFFFF33',lwd=2)
  box()
}

stanCodeWithMod<-'
  data {
    int<lower=0> nVirus;
    int<lower=0> nAllele;
    int<lower=0> nRep;
    int<lower=0> nGroup;
    int<lower=0> nExperiment;
    int<lower=0> nModification;
    vector[nRep] rep;
    int<lower=0> allele[nRep];
    int<lower=0> virus[nRep];
    int<lower=0> virusGroup[nVirus];
    int<lower=0> experiment[nRep];
    int<lower=0> modification[nRep];
    real backMean;
    real<lower=0> backSd;
    real<lower=0> tNu;
  }
  parameters {
    matrix[nVirus,nAllele] alleleVirus;
    matrix[nVirus,nExperiment-1] virusExperiment;
    matrix[nAllele,nGroup] alleleGroup;
    real<lower=0> sdWithinSpecies;
    real<lower=0> repSd;
    vector<lower=0>[nRep] background;
    real<lower=0> alleleSd;
    vector[nAllele] sdMult;
    vector[nAllele] alleleMeans;
    vector[nVirus] virusMeans;
    real modMeans[nModification-1];
    real modSd[nModification-1];
    matrix[nVirus,nModification-1] virusMod;
  }
  transformed parameters{
    vector[nRep] beta;
    for(ii in 1:nRep){
      beta[ii]=alleleVirus[virus[ii],allele[ii]]+virusMeans[virus[ii]];
      if(experiment[ii]!=1)beta[ii]=beta[ii]+virusExperiment[virus[ii],experiment[ii]-1];
      if(modification[ii]!=1)beta[ii]=beta[ii]+virusMod[virus[ii],modification[ii]-1];
      beta[ii]=logit(fmin(inv_logit(beta[ii])+background[ii],.9999));
    }
  }
  model {
    rep~normal(beta,repSd);
    for(ii in 1:nVirus)alleleVirus[ii,]~normal(alleleGroup[,virusGroup[ii]],sdWithinSpecies);
    for(ii in 1:nAllele)alleleGroup[ii,]~normal(alleleMeans[ii],alleleSd*exp(sdMult[ii]));
    for(ii in 2:nExperiment)virusExperiment[,ii-1]~normal(0,5);
    repSd~gamma(1,.1);
    alleleSd~gamma(1,.1);
    sdMult~normal(0,.693); //log(2)
    background~normal(backMean,backSd);
    alleleMeans~normal(-2,2.5);
    //alleleMeansMean~normal(-2,5);
    virusMeans~normal(0,2);
    sdWithinSpecies~gamma(1,.1);
    modMeans~normal(0,5);
    modSd~gamma(1,.1);
    for(ii in 1:nVirus)virusMod[ii,]~normal(modMeans,modSd);
  }
'
modAlleleWithMods <- stan_model(model_code = stanCodeWithMod)

stanCode<-'
  data {
    int<lower=0> nVirus;
    int<lower=0> nAllele;
    int<lower=0> nRep;
    int<lower=0> nGroup;
    int<lower=0> nExperiment;
    vector[nRep] rep;
    int<lower=0> allele[nRep];
    int<lower=0> virus[nRep];
    int<lower=0> virusGroup[nVirus];
    int<lower=0> experiment[nRep];
    real backMean;
    real<lower=0> backSd;
    real<lower=0> tNu;
  }
  parameters {
    matrix[nVirus,nAllele] alleleVirus;
    matrix[nVirus,nExperiment-1] virusExperiment;
    matrix[nAllele,nGroup] alleleGroup;
    real<lower=0> sdWithinSpecies;
    real<lower=0> repSd;
    vector<lower=0>[nRep] background;
    real<lower=0> alleleSd;
    vector[nAllele] sdMult;
    vector[nAllele] alleleMeans;
    vector[nVirus] virusMeans;
  }
  transformed parameters{
    vector[nRep] beta;
    for(ii in 1:nRep){
      beta[ii]=alleleVirus[virus[ii],allele[ii]]+virusMeans[virus[ii]];
      if(experiment[ii]!=1)beta[ii]=beta[ii]+virusExperiment[virus[ii],experiment[ii]-1];
      beta[ii]=logit(fmin(inv_logit(beta[ii])+background[ii],.9999));
    }
  }
  model {
    rep~normal(beta,repSd);
    for(ii in 1:nVirus)alleleVirus[ii,]~normal(alleleGroup[,virusGroup[ii]],sdWithinSpecies);
    for(ii in 1:nAllele)alleleGroup[ii,]~normal(alleleMeans[ii],alleleSd*exp(sdMult[ii]));
    for(ii in 2:nExperiment)virusExperiment[,ii-1]~normal(0,5);
    repSd~gamma(1,.1);
    alleleSd~gamma(1,.1);
    sdMult~normal(0,.693); //log(2)
    background~normal(backMean,backSd);
    alleleMeans~normal(-2,2.5);
    //alleleMeansMean~normal(-2,5);
    virusMeans~normal(0,2);
    sdWithinSpecies~gamma(1,.1);
  }
'
modAllele <- stan_model(model_code = stanCode)

getDiffs<-function(fit){
  mat<-as.matrix(fit$fit)
  diffs<-do.call(rbind,parallel::mclapply(c('General'=-99,fit$virusId),function(vir){
    do.call(rbind,lapply(fit$alleleId,function(xx){
      out<-data.frame(do.call(rbind,lapply(fit$alleleId,function(yy){
        if(xx==yy)stats<-c(0,0,0,.5)
        else if(vir==-99)stats<-crIMean(mat[,sprintf('alleleMeans[%d]',xx)]-mat[,sprintf('alleleMeans[%d]',yy)])
        else stats<-crIMean(mat[,sprintf('alleleVirus[%d,%d]',vir,xx)]-mat[,sprintf('alleleVirus[%d,%d]',vir,yy)])
        #if(out[3]<0){
          #tmp<-xx;xx<-yy;yy<-tmp
          #stats<-crIMean(mat[,sprintf('alleleVirus[%d,%d]',vir,xx)]-mat[,sprintf('alleleVirus[%d,%d]',vir,yy)])
        #}
        out<-data.frame('allele1'=xx,'allele2'=yy,'mean'=stats[3],'lower'=stats[1],'upper'=stats[2],'p'=stats[4],stringsAsFactors=FALSE)
        return(out)
      })))
      out$vir<-vir
      return(out)
    }))
  },mc.cores=20))
  diffs<-diffs[diffs$mean>0&diffs$allele1!=diffs$allele2,]
  if(nrow(diffs)!=choose(max(fit$alleleId),2)*(1+max(fit$virusId)))stop('Missing rows')
  diffs$virus[diffs$vir!=-99]<-names(fit$virusId)[diffs$vir[diffs$vir!=-99]]
  diffs$virus[diffs$vir==-99]<-'General'
  diffs$all1<-names(fit$alleleId)[diffs$allele1]
  diffs$all2<-names(fit$alleleId)[diffs$allele2]
  return(diffs)
}
