
ECoG Stats Report
========================================================

This report uses lmer and linear models to evaluate the statistical validity 
of the data. 

Data details:

```r
lockType = "RT"
type = "erp"
band = ""
```


Additional setup:
-------------------------------------------------------
Load necessary libraries:

```r
library("lme4")
```


Set additional variables:

```r
dataPath = "~/Documents/ECOG/Results/Rdata/"
hem = "all"
ref = "nonLPCleasL1TvalCh10"
binType = "Bin"
byBlockStr = "bySubj"

if (lockType == "stim") {
    timeStr = "200msTo900ms"
} else if (lockType == "RT") {
    timeStr = "n600msTo100ms"
}

if (type == "erp") {
    band = ""
    preFix = paste(binType, timeStr, hem, "ERPs", sep = "")
    baseLine = "sub"
    analysisType = "Amp"
} else if (type == "power") {
    preFix = paste(binType, timeStr, hem, "ERSPs", band, sep = "")
    baseLine = "sub"
    analysisType = "logPower"
}

fileName = paste(preFix, "Group", lockType, "Lock", baseLine, analysisType, 
    ref, byBlockStr, sep = "")
```


Load data:

```r
data <- read.csv(paste(dataPath, fileName, ".csv", sep = ""), header = F)
names(data) <- c("M", "blank1", "T", "blank2", "Z", "blank3", "BinNum", "chanId", 
    "subjId", "blockId", "ROI", "dP", "HR", "FAR", "hem", "subROI", "ITC_Z", 
    "ITP_Z", "Zc", "ZcH", "ZcCRs", "ZH", "ZCRs")

data$BinNum = as.factor(data$BinNum)
data$chanId = as.factor(data$chanId)
data$subjId = as.factor(data$subjId)
data$blockId = as.factor(data$blockId)
data$ROI = as.factor(data$ROI)
levels(data$ROI) = c("IPS", "SPL", "AG")
data$chanId = as.factor(data$chanId)
data$hem = as.factor(data$hem)
levels(data$hem) = c("left", "right")
data$subROI = as.factor(data$subROI)
levels(data$subROI) = c("", "pIPS", "aIPS", "pSPL", "aSPL")

nBins = length(unique(data$BinNum))
nChans = length(unique(data$chanId))
nSubjs = length(unique(data$subjId))
nBlocks = length(unique(data$blockId))
nROIs = length(unique(data$ROI))
```


Simple linear models, evaluated with anova
------------------------------------------------
**For the lefts subjects:**

```r
subL = data$hem == "left"
LModel.lm = lm(Z ~ 1 + ROI * BinNum + subjId + chanId, data = data, subset = subL)
anova(LModel.lm)
```

```
## Analysis of Variance Table
## 
## Response: Z
##             Df Sum Sq Mean Sq F value  Pr(>F)    
## ROI          2      5    2.61    2.41   0.091 .  
## BinNum       6    135   22.51   20.76 < 2e-16 ***
## subjId       3     47   15.72   14.50 5.5e-09 ***
## chanId      64    105    1.63    1.51   0.010 *  
## ROI:BinNum  12     49    4.07    3.76 2.0e-05 ***
## Residuals  402    436    1.08                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


**For the right subjects:**

```r
subR = data$hem == "right"
RModel.lm = lm(Z ~ 1 + ROI * BinNum + subjId + chanId, data = data, subset = subR)
anova(RModel.lm)
```

```
## Analysis of Variance Table
## 
## Response: Z
##             Df Sum Sq Mean Sq F value  Pr(>F)    
## ROI          2    5.3    2.65    3.69 0.02641 *  
## BinNum       6   18.8    3.14    4.36 0.00033 ***
## subjId       2   10.6    5.28    7.34 0.00080 ***
## chanId      41   89.0    2.17    3.01 5.6e-08 ***
## ROI:BinNum  12    7.2    0.60    0.83 0.62110    
## Residuals  258  185.7    0.72                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


**Evaluate hemishphere interatction**

```r
FModel.lm = lm(Z ~ 1 + ROI * BinNum * hem + subjId + chanId, data = data)
anova(FModel.lm)
```

```
## Analysis of Variance Table
## 
## Response: Z
##                 Df Sum Sq Mean Sq F value  Pr(>F)    
## ROI              2      0    0.01    0.02 0.98502    
## BinNum           6     79   13.18   13.99 5.0e-15 ***
## hem              1      5    5.19    5.51 0.01920 *  
## subjId           5     58   11.54   12.25 2.2e-11 ***
## chanId         107    204    1.91    2.02 8.8e-08 ***
## ROI:BinNum      12     36    3.00    3.18 0.00019 ***
## BinNum:hem       6     73   12.09   12.84 9.4e-14 ***
## ROI:BinNum:hem  12     22    1.86    1.97 0.02439 *  
## Residuals      660    621    0.94                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```



Mixed Models evaluated with lmer
-----------------------------------------------------

First we define the `LRT` function that takes a MER model and an index 
for a factor, and returns the Chi square statistic based on the log-likelihood 
of the full and nested models.

```r
LRT <- function(theModel, idx) {
    modFull = theModel
    modelTerms_h <- terms(theModel)
    modelTerms = attr(modelTerms_h, "term.labels")
    updateStr <- paste(".~.-", as.character(modelTerms[idx]))
    modNested = update(theModel, updateStr)
    
    llFull = logLik(modFull)
    llNested = logLik(modNested)
    
    chsqval <- 2 * (llFull[1] - llNested[1])
    dfDiff <- attr(llFull, "df") - attr(llNested, "df")
    pVal <- 1 - pchisq(chsqval, dfDiff)
    
    out <- c(as.character(modelTerms[idx]), pVal, chsqval, dfDiff)
    names(out) <- c("Factor", " P val ", " ChiSq ", " dF ")
    return(out)
}
```


**Mixed model for lefts**

```r
LModel.lmer = lmer(Z ~ 1 + subjId + BinNum + ROI + (1 + BinNum | chanId), data = data, 
    subset = subL, REML = F)
LModelInt.lmer = lmer(Z ~ 1 + subjId + BinNum * ROI + (1 + BinNum | chanId), 
    data = data, subset = subL, REML = F)
LRT(LModel.lmer, 3)
```

```
##                Factor                P val                 ChiSq  
##                 "ROI" "0.00345700044253328"    "11.3347079800637" 
##                   dF  
##                   "2"
```

```r
LRT(LModelInt.lmer, 4)
```

```
##                 Factor                 P val                  ChiSq  
##           "BinNum:ROI" "1.48812617515315e-05"     "44.0659957651958" 
##                    dF  
##                   "12"
```


**Mixed model for rights**

```r
RModel.lmer = lmer(Z ~ 1 + subjId + BinNum + ROI + (1 + BinNum | chanId), data = data, 
    subset = subR, REML = F)
RModelInt.lmer = lmer(Z ~ 1 + subjId + BinNum * ROI + (1 + BinNum | chanId), 
    data = data, subset = subR, REML = F)
LRT(RModel.lmer, 3)
```

```
##              Factor              P val               ChiSq  
##               "ROI" "0.962732058118835"  "0.07596028507146" 
##                 dF  
##                 "2"
```

```r
LRT(RModelInt.lmer, 4)
```

```
##              Factor              P val               ChiSq  
##        "BinNum:ROI" "0.339831404375435"  "13.4118226000512" 
##                 dF  
##                "12"
```


**Mixed model to test hemisphere interaction.**
Note that the use of subject as a factor in this model results in an 
overdetermined matrix. 

```r
FModel.lmer = lmer(Z ~ 1 + BinNum * ROI * hem + (1 + BinNum | chanId), data = data, 
    REML = F)
LRT(FModel.lmer, 7)
```

```
##              Factor              P val               ChiSq  
##    "BinNum:ROI:hem" "0.124322672288091"  "17.7242714300469" 
##                 dF  
##                "12"
```

```r
FModel.lmer = lmer(Z ~ 1 + subjId + BinNum * ROI * hem + (1 + BinNum | chanId), 
    data = data, REML = F)
```

```
## Error: Downdated X'X is not positive definite, 16.
```



