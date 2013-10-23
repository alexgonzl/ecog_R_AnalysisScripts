
library("leaps")
library("lme4")
library("nlme")
library('ez')
library('HH')
library("fork")

dataPath = '/Users/alexgonzalez/Google Drive/Research/ECoG/ERSP_SummarySTATS/data/'
options(contrasts=c("contr.sum","contr.poly"))



info <-read.csv(paste(path,subjs_str,'Subjs_BigMatInfoChanMean',SplitType,'.csv',sep=''),header=F)
names(info) <- c('M','T','Z','BinNum','chanId','subjId', 'blockId', 'ROI','dP','HR','FAR','hem')

sub_id  = as.factor(info$sub_id);
bac     = info$bac;
region  = as.factor(info$region);
levels(region) = c('IPS','SPL','AG')
condition = as.factor(info$condition)
levels(condition) = c('Hits','CRs')
chan_id = as.factor(info$chan_id)
bin_num = as.factor(info$bin_num);
bin_num_numeric = (info$bin_num);
RTtype = as.factor(info$RTtype);
RTs = info$RTs;
hem = as.factor(info$hem);
levels(hem) = c('left','right')

if (baseline==1){
  str1 ='Relative'
} else if (baseline==2){
  str1 = 'SubtractedBaseline'}

if (dataType==1)      {str2 = '' 
} else if (dataType==2) {str2 = 'Pct'}  

if (analyzeType==1) { str3 = 'Amplitude'
} else if (analyzeType==2){str3 = 'logAmp'
} else if (analyzeType==3){str3 = 'Power'
} else if (analyzeType==4){str3 = 'logPower'}

str5 = '_binsize50_sldwin50'

if (RTlocked==TRUE){str6 = 'RTlocked'
} else {str6 =''}  

if (standarized==TRUE){str7 = 'zscored' 
} else {str7 = ''}

str8 = SplitType;

if (erpReject==TRUE){str9 = 'ERPRejected'
} else {str9 ='';}  

for (band in bands){  
  str4 = band_names[band];
  extstr = paste(str1,str2,str3,str4,str5,str6,str7,str8,str9,date_str,sep='');  
  erp   <-read.csv(paste(path,subjs_str,'Subjs_BigMatERSPChanMean',extstr,'.csv',sep=''),header=F)
  
  nregions=3;
  nobs = length(sub_id)
  id1 = rep(1,nobs)
  id1 = as.factor(id1)
  erp = erp[allregions!=4,];
  conditionCont = (condition=='Hits')-(condition=='CRs') 
  nchans = length(levels(chan_id))
  
  data = data.frame(erp,conditionCont,region,chan_id,sub_id,bin_num,RTs,id1,hem, RTtype, bin_num_numeric)
  
  ######################## 
  #if(subjs_str=='all'){
  if (FALSE){
    
    # 1) Test of Region X Condition X Time interaction
    subset = (hem==hem)
    ModelFull = lm(erp ~ region*conditionCont*bin_num+chan_id,data=data,subset=subset);
    ModelRed = lm(erp ~  region*conditionCont*bin_num-region:conditionCont:bin_num+chan_id,data=data,subset=subset);
    anova(ModelRed,ModelFull)
    
    subset = (hem=='left') 
    ModelFull = lm(erp ~ region*conditionCont*bin_num+chan_id,data=data,subset=subset);
    ModelRed = lm(erp ~  region*conditionCont*bin_num-region:conditionCont:bin_num+chan_id,data=data,subset=subset);
    anova(ModelRed,ModelFull)
    
    subset = (hem=='right')
    ModelFull = lm(erp ~ region*conditionCont*bin_num+chan_id,data=data,subset=subset);
    ModelRed = lm(erp ~  region*conditionCont*bin_num-region:conditionCont:bin_num+chan_id,data=data,subset=subset);
    anova(ModelRed,ModelFull)
    
    # 1a) Test of Condition X Time interactions on lefts (separated by region)
    subset = (hem=='left') & (region=='AG')
    ModelFull = lm(erp ~ conditionCont*bin_num+chan_id,data=data,subset=subset);
    ModelRed  = lm(erp ~ conditionCont*bin_num+chan_id-conditionCont:bin_num,data=data,subset=subset);
    anova(ModelRed,ModelFull)
    subset = (hem=='left') & (region=='IPS')
    ModelFull = lm(erp ~ conditionCont*bin_num+chan_id,data=data,subset=subset);
    ModelRed  = lm(erp ~ conditionCont*bin_num+chan_id-conditionCont:bin_num,data=data,subset=subset);
    anova(ModelRed,ModelFull)
    subset = (hem=='left') & (region=='SPL')
    ModelFull = lm(erp ~ conditionCont*bin_num+chan_id,data=data,subset=subset);
    ModelRed  = lm(erp ~ conditionCont*bin_num+chan_id-conditionCont:bin_num,data=data,subset=subset);
    anova(ModelRed,ModelFull)
    
    # 2) Test of Condition X Time X RT interactions on lefts (separated by region)  
    subset = (hem=='left') & (region=='AG')
    ModelFull = lm(erp ~ conditionCont*bin_num*RTtype+chan_id,data=data,subset=subset);
    ModelRed  = lm(erp ~ conditionCont*bin_num*RTtype+chan_id-conditionCont:bin_num:RTtype,data=data,subset=subset);
    anova(ModelRed,ModelFull)
    subset = (hem=='left') & (region=='IPS')
    ModelFull = lm(erp ~ conditionCont*bin_num*RTtype+chan_id,data=data,subset=subset);
    ModelRed  = lm(erp ~ conditionCont*bin_num*RTtype+chan_id-conditionCont:bin_num:RTtype,data=data,subset=subset);
    anova(ModelRed,ModelFull)
    subset = (hem=='left') & (region=='SPL')
    ModelFull = lm(erp ~ conditionCont*bin_num*RTtype+chan_id,data=data,subset=subset);
    ModelRed  = lm(erp ~ conditionCont*bin_num*RTtype+chan_id-conditionCont:bin_num:RTtype,data=data,subset=subset);
    anova(ModelRed,ModelFull)    
    
    # 3a) Test of Condition X Time X Hemisphere interaction separated by region
    subset = (region=='AG')
    ModelFull = lm(erp ~ conditionCont*bin_num*hem+chan_id,data=data,subset=subset);
    ModelRed  = lm(erp ~ conditionCont*bin_num*hem+chan_id-conditionCont:bin_num:hem,data=data,subset=subset);
    anova(ModelRed,ModelFull)
    subset = (region=='IPS')
    ModelFull = lm(erp ~ conditionCont*bin_num*hem+chan_id,data=data,subset=subset);
    ModelRed  = lm(erp ~ conditionCont*bin_num*hem+chan_id-conditionCont:bin_num:hem,data=data,subset=subset);
    anova(ModelRed,ModelFull)
    subset = (region=='SPL')
    ModelFull = lm(erp ~ conditionCont*bin_num*hem+chan_id,data=data,subset=subset);
    ModelRed  = lm(erp ~ conditionCont*bin_num*hem+chan_id-conditionCont:bin_num:hem,data=data,subset=subset);
    anova(ModelRed,ModelFull)    
    
    # 3b) Test of Condition X Time X Hemisphere interaction not separating by region    
    ModelFull = lm(erp ~ conditionCont*bin_num*hem+chan_id,data=data);
    ModelRed  = lm(erp ~ conditionCont*bin_num*hem+chan_id-conditionCont:bin_num:hem,data=data);
    anova(ModelRed,ModelFull)   
    
    # 4) Test for main effect of condition: separated by region.
    subset = (hem=='left') & (region=='AG')
    ModelFull = lm(erp ~ conditionCont*bin_num+chan_id,data=data,subset=subset);
    anova(ModelFull)
    subset = (hem=='left') & (region=='IPS')
    ModelFull = lm(erp ~ conditionCont*bin_num+chan_id,data=data,subset=subset);
    anova(ModelFull)
    subset = (hem=='left') & (region=='SPL')
    ModelFull = lm(erp ~ conditionCont*bin_num+chan_id,data=data,subset=subset);
    anova(ModelFull)
    subset = (hem=='right') & (region=='AG')
    ModelFull = lm(erp ~ conditionCont*bin_num+chan_id,data=data,subset=subset);
    anova(ModelFull)
    subset = (hem=='right') & (region=='IPS')
    ModelFull = lm(erp ~ conditionCont*bin_num+chan_id,data=data,subset=subset);
    anova(ModelFull)
    subset = (hem=='right') & (region=='SPL')
    ModelFull = lm(erp ~ conditionCont*bin_num+chan_id,data=data,subset=subset);
    anova(ModelFull)
    
  }
  
  ################################################################################
  ################################################################################
  # all trials  
  subset = (hem==hemisphere)
  nchans = length(unique(chan_id[subset]));
  modelFull=lm(erp~conditionCont:bin_num:region-1+chan_id,data, subset = subset)
  # main effects
  idx = -(1:nchans);
  x = summary(modelFull)
  x = coef(x);
  Betas  = matrix(x[idx,1],3,24,byrow=T);
  SE     = matrix(x[idx,2],3,24,byrow=T);
  Pvals  = matrix(x[idx,4],3,24,byrow=T);
  
  # set up for interactionEffects
  x = lm(erp~conditionCont:region:bin_num-1+chan_id,data,subset=subset);
  ncofs = length(coef(x));
  # test at bin_num9
  erp_interactPvals = matrix(nrow=nregions,ncol=nbins)
  for (b in 1:nbins) {
    bin_start_id = nchans+nregions*(b-1);
    bin_end_id =   ncofs-bin_start_id-3;
    K <- matrix(c( 1, -1,  0, -1, 0, 1, 0 ,1 ,-1), nrow=3)
    rownames(K) <- c('IPS-SPL','AG-IPS', 'SPL-AG')
    K = cbind(matrix(0,nrow=3,ncol=bin_start_id),K,matrix(0,nrow=3,ncol=bin_end_id))
    y = summary(glht(x,linfct=K))
    
    erp_interactPvals[,b]=y[9]$test[6]$pvalues[1:3]
  }
  
  path2 = '/Users/alexgonzalez/Documents/ECOG/Results/Spectral_Data/group/';
  write.table(Pvals,paste(path2,hemisphere,'SubjsERSP',str4, '_PValFixedEffectsTimeInteractRegionChanMean_',extstr,sep=''),col.names=F,row.names=F,sep=',')
  
  write.table(Betas,paste(path2,hemisphere,'SubjsERSP',str4,'_BetaEstimateFixedEffectsTimeInteractRegionChanMean_',extstr,sep=''),col.names=F,row.names=F,sep=',')
  
  write.table(SE,paste(path2,hemisphere,'SubjsERSP',str4,'_BetaSEFixedEffectsTimeInteractRegionChanMean_',extstr,sep=''),col.names=F,row.names=F,sep=',')
  
  write.table(erp_interactPvals,paste(path2,hemisphere,'SubjsERSP',str4,'_PVal_Interact_FixedEffectsTimeInteractRegionChanMean_',extstr,sep=''),col.names=F,row.names=F,sep=',')
  
  ################################################################################
  ################################################################################
  if (TRUE){
    #fasttrials
    subset = (RTtype==1) & (hem==hemisphere);
    modelFull=lm(erp~conditionCont:bin_num:region-1+chan_id,data,subset=subset);
    # main effects
    idx = -(1:nchans);
    x = summary(modelFull)
    x = coef(x);
    Betas  = matrix(x[idx,1],3,24,byrow=T);
    SE     = matrix(x[idx,2],3,24,byrow=T);
    Pvals  = matrix(x[idx,4],3,24,byrow=T);
    
    # set up for interactionEffects
    x = lm(erp~conditionCont:region:bin_num-1+chan_id,data=data,subset=subset);
    ncofs = length(coef(x));
    erp_interactPvals = matrix(nrow=nregions,ncol=nbins)
    for (b in 1:nbins) {
      bin_start_id = nchans+nregions*(b-1);
      bin_end_id =   ncofs-bin_start_id-3;
      K <- matrix(c( 1, -1,  0, -1, 0, 1, 0 ,1 ,-1), nrow=3)
      rownames(K) <- c('IPS-SPL','AG-IPS', 'SPL-AG')
      K = cbind(matrix(0,nrow=3,ncol=bin_start_id),K,matrix(0,nrow=3,ncol=bin_end_id))
      y = summary(glht(x,linfct=K))
      
      erp_interactPvals[,b]=y[9]$test[6]$pvalues[1:3]
    }
    write.table(Pvals,paste(path2,hemisphere,'SubjsERSP','fastRTs_',str4,'_PValFixedEffectsTimeInteractRegionChanMean_',extstr,sep=''),col.names=F,row.names=F,sep=',')
    
    write.table(Betas,paste(path2,hemisphere,'SubjsERSP','fastRTs_',str4,'_BetaEstimateFixedEffectsTimeInteractRegionChanMean_',extstr,sep=''),col.names=F,row.names=F,sep=',')
    
    write.table(SE,paste(path2,hemisphere,'SubjsERSP','fastRTs_',str4,'_BetaSEFixedEffectsTimeInteractRegionChanMean_',extstr,sep=''),col.names=F,row.names=F,sep=',')
    
    write.table(erp_interactPvals,paste(path2,hemisphere,'SubjsERSP','fastRTs_',str4,'_InteractionPValFixedEffectsTimeInteractRegionChanMean_',extstr,sep=''),col.names=F,row.names=F,sep=',')
    
    ################################################################################
    ################################################################################
    #slow trials
    subset = (RTtype==2) & (hem==hemisphere);
    modelFull=lm(erp~conditionCont:bin_num:region-1+chan_id,data=data,subset=subset);
    # main effects
    idx = -(1:nchans);
    x = summary(modelFull)
    x = coef(x);
    Betas  = matrix(x[idx,1],3,24,byrow=T);
    SE     = matrix(x[idx,2],3,24,byrow=T);
    Pvals  = matrix(x[idx,4],3,24,byrow=T);
    
    # set up for interactionEffects
    x = lm(erp~conditionCont:region:bin_num-1+chan_id,data=data,subset=subset);
    ncofs = length(coef(x));
    erp_interactPvals = matrix(nrow=nregions,ncol=nbins)
    for (b in 1:nbins) {
      bin_start_id = nchans+nregions*(b-1);
      bin_end_id =   ncofs-bin_start_id-3;
      K <- matrix(c( 1, -1,  0, -1, 0, 1, 0 ,1 ,-1), nrow=3)
      rownames(K) <- c('IPS-SPL','AG-IPS', 'SPL-AG')
      K = cbind(matrix(0,nrow=3,ncol=bin_start_id),K,matrix(0,nrow=3,ncol=bin_end_id))
      y = summary(glht(x,linfct=K))
      
      erp_interactPvals[,b]=y[9]$test[6]$pvalues[1:3]
    }
    
    write.table(Pvals,paste(path2,hemisphere,'SubjsERSP','slowRTs_',str4,'_PValFixedEffectsTimeInteractRegionChanMean_',extstr,sep=''),col.names=F,row.names=F,sep=',')
    
    write.table(Betas,paste(path2,hemisphere,'SubjsERSP','slowRTs_',str4,'_BetaEstimateFixedEffectsTimeInteractRegionChanMean_',extstr,sep=''),col.names=F,row.names=F,sep=',')
    
    write.table(SE,paste(path2,hemisphere,'SubjsERSP','slowRTs_',str4,'_BetaSEFixedEffectsTimeInteractRegionChanMean_',extstr,sep=''),col.names=F,row.names=F,sep=',')
    
    write.table(erp_interactPvals,paste(path2,hemisphere,'SubjsERSP','slowRTs_',str4,'_InteractionPValFixedEffectsTimeInteractRegionChanMean_',extstr,sep=''),col.names=F,row.names=F,sep=',')
  }
}