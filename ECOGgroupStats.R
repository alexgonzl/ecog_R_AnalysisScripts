dataPath = '~/Documents/ECOG/Results/Rdata/';
options(contrasts=c("contr.sum","contr.poly"))

hem 		= 'all'; #{'l','r','all'}
lockType 	= 'stim'; #{'stim','RT'}
ref 		= 'nonLPCleasL1TvalCh10';
type 		= 'power';
band 		= 'hgam'; # ignored if erp
binType     = 'Bin';
timeStr 	= '200msTo900ms'; #{'n600msTo100ms', '200msTo900ms'}
byBlockStr 	= 'bySubj'; # {'byBlock','bySubj'}

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

data$BinNum 		= as.factor(data$BinNum) # modeling time as a factor
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
# lefts
type
band
subset 			= data$hem=='left'
leftFullModel 	= lm(Z ~ 1+ ROI*BinNum+subjId+chanId,data=data,subset=subset);
anova(leftFullModel)

# rights
subset 			= data$hem=='right'
rightFullModel 	= lm(Z ~ 1+ ROI*BinNum+subjId+chanId,data=data,subset=subset);
anova(rightFullModel)

##########################################
# all
type
band
FullModel 		= lm(Z ~ 1+ ROI*BinNum*hem+subjId+chanId,data=data)
anova(FullModel)

##########################################
# lefts sub regions
subset 			= (data$hem=='left' ) & (data$ROI=='IPS' )
model 			= lm(Z ~ 1+ subROI*BinNum+subjId+chanId,data=data,subset=subset);
anova(model)

subset 			= (data$hem=='left' ) & (data$ROI=='SPL' )
model 			= lm(Z ~ 1+ subROI*BinNum+subjId+chanId,data=data,subset=subset);
anova(model)

##########################################
# right sub regions
subset 			= (data$hem=='right' ) & (data$ROI=='IPS' )
model 			= lm(Z ~ 1+ subROI*BinNum+subjId+chanId,data=data,subset=subset);
anova(model)

subset 			= (data$hem=='right' ) & (data$ROI=='SPL' )
model 			= lm(Z ~ 1+ subROI*BinNum+subjId+chanId,data=data,subset=subset);
anova(model)
