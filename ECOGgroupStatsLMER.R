
library('lme4')
dataPath = '~/Documents/ECOG/Results/Rdata/';
options(contrasts=c("contr.sum","contr.poly"))

hem 		= 'all'; #{'l','r','all'}
lockType 	= 'stim'; #{'stim','RT'}
ref 		= 'nonLPCleasL1TvalCh10';
type 		= 'power';
band 		= 'hgam'; # ignored if erp
binType     = 'Bin';
byBlockStr 	= 'bySubj';

if (lockType == 'stim'){
	timeStr 	= '200msTo900ms'	
}else if (lockType == 'RT'){
	timeStr 	= 'n600msTo100ms'	
}

if (type=='erp'){  
  band 			= '';
  preFix 		= paste(binType,timeStr,hem,'ERPs',sep='')
  baseLine 		= 'sub';
  analysisType 	= 'Amp';
}else if (type=='power'){
  preFix 		= paste(binType,timeStr,hem,'ERSPs',band,sep='')
  baseLine 		= 'sub';
  analysisType 	= 'logPower';
}

fileName = paste(preFix,'Group',lockType,'Lock',baseLine,analysisType,ref,byBlockStr,sep='');

data 		<- read.csv(paste(dataPath,fileName,'.csv',sep=''),header=F)
names(data) <- c('M','blank1','T','blank2','Z','blank3','BinNum','chanId','subjId', 'blockId', 
				 'ROI','dP','HR','FAR','hem','subROI','ITC_Z','ITP_Z',
				 'Zc','ZcH','ZcCRs','ZH','ZCRs')

data$BinNum 		= as.factor(data$BinNum)
data$chanId 		= as.factor(data$chanId)
data$subjId  		= as.factor(data$subjId)
data$blockId  		= as.factor(data$blockId)
data$ROI  			= as.factor(data$ROI);
levels(data$ROI)	= c('IPS','SPL','AG')
data$chanId 		= as.factor(data$chanId)
data$hem 			= as.factor(data$hem);
levels(data$hem)	= c('left','right')
data$subROI  		= as.factor(data$subROI);
levels(data$subROI) = c('','pIPS','aIPS','pSPL','aSPL')

nBins 	= length(unique(data$BinNum));
nChans 	= length(unique(data$chanId));
nSubjs 	= length(unique(data$subjId));
nBlocks = length(unique(data$blockId));
nROIs 	= length(unique(data$ROI));

########################################## 
## linear models ##
# lefts
subL 			= data$hem=='left'
LModel.lm 	= lm(Z ~ 1+ ROI*BinNum+subjId+chanId,data=data,subset=subL);
anova(LModel.lm)

# rights
subR 			= data$hem=='right'
RModel.lm 	= lm(Z ~ 1+ ROI*BinNum+subjId+chanId,data=data,subset=subR);
anova(RModel.lm)

#' all
FModel.lm 		= lm(Z ~ 1+ ROI*BinNum*hem+subjId+chanId,data=data)
anova(FModel.lm)

########################################
#' linear models with mixed effects 
#' lefts
LModel.lmer		= lmer(Z ~ 1 + subjId + BinNum + ROI + (1 + BinNum|chanId) , data=data , subset=subL , REML = F)
LModelInt.lmer	= lmer(Z ~ 1 + subjId + BinNum * ROI + (1 + BinNum|chanId) , data=data , subset=subL , REML = F)
LRT(LModel.lmer,3)
LRT(LModelInt.lmer,4)

#' rights
RModel.lmer	    = lmer(Z ~ 1 + subjId + BinNum + ROI + (1 + BinNum|chanId) , data=data , subset=subR , REML = F)
RModelInt.lmer	= lmer(Z ~ 1 + subjId + BinNum * ROI + (1 + BinNum|chanId) , data=data , subset=subR , REML = F)
LRT(RModel.lmer,3)
LRT(RModelInt.lmer,4)

#' both hemispheres
FModel.lmer  	= lmer(Z ~ 1 + BinNum * ROI * hem + (1 + BinNum|chanId) , data=data , REML = F)
LRT(FModel.lmer,7)





