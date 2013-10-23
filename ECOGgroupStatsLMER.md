
ECoG Stats Report
========================================================

This report uses lmer and linear models to evaluate the statistical validity 
of the data. 

Data details:

```r
lockType = "stim"
type = "power"
band = "hgam"
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
## ROI          2    126    63.1   82.18 < 2e-16 ***
## BinNum       6     58     9.7   12.67 4.2e-13 ***
## subjId       3     27     8.9   11.65 2.4e-07 ***
## chanId      64    388     6.1    7.90 < 2e-16 ***
## ROI:BinNum  12     25     2.0    2.67  0.0019 ** 
## Residuals  402    309     0.8                    
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
## ROI          2   10.2    5.11    7.60 0.00062 ***
## BinNum       6    2.1    0.35    0.51 0.79766    
## subjId       2   11.5    5.74    8.54 0.00026 ***
## chanId      41  191.1    4.66    6.93 < 2e-16 ***
## ROI:BinNum  12   10.9    0.91    1.36 0.18729    
## Residuals  258  173.4    0.67                    
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
## ROI              2     46    22.8   31.23 1.1e-13 ***
## BinNum           6     43     7.2    9.93 1.7e-10 ***
## hem              1     49    48.9   67.01 1.4e-15 ***
## subjId           5     39     7.7   10.60 8.3e-10 ***
## chanId         107    659     6.2    8.43 < 2e-16 ***
## ROI:BinNum      12     12     1.0    1.37 0.17643    
## BinNum:hem       6     15     2.5    3.47 0.00220 ** 
## ROI:BinNum:hem  12     25     2.1    2.88 0.00069 ***
## Residuals      660    482     0.7                    
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
##                 Factor                 P val                  ChiSq  
##                  "ROI" "0.000337709604170677"      "15.986648382057" 
##                    dF  
##                    "2"
```

```r
LRT(LModelInt.lmer, 4)
```

```
##               Factor               P val                ChiSq  
##         "BinNum:ROI" "0.0456485033833335"   "21.3377167309745" 
##                  dF  
##                 "12"
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
##               "ROI" "0.285719272594439"   "2.5054910291343" 
##                 dF  
##                 "2"
```

```r
LRT(RModelInt.lmer, 4)
```

```
##              Factor              P val               ChiSq  
##        "BinNum:ROI" "0.199763348045416"  "15.8169781417179" 
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
##               Factor               P val                ChiSq  
##     "BinNum:ROI:hem" "0.0216629206232617"    "23.798517562062" 
##                  dF  
##                 "12"
```

```r
FModel.lmer = lmer(Z ~ 1 + subjId + BinNum * ROI * hem + (1 + BinNum | chanId), 
    data = data, REML = F)
```

```
## Error: Downdated X'X is not positive definite, 16.
```



