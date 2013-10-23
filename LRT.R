LRT <- function(theModel, idx){
	
	modFull = theModel
	modelTerms_h<-terms(theModel)
	modelTerms = attr(modelTerms_h, "term.labels")
	updateStr<-paste(".~.-",as.character(modelTerms[idx]))
	modNested = update(theModel, updateStr)
	
	llFull = logLik(modFull)
	llNested = logLik(modNested)
	
	chsqval<-2*(llFull[1] - llNested[1])
	dfDiff<- attr(llFull,"df") - attr(llNested,"df")
	
	pVal<-1-pchisq(chsqval,dfDiff)
	
	out <- c(as.character(modelTerms[idx]),pVal, chsqval, dfDiff)
	names(out) <- c('Factor', ' P val ' , ' ChiSq ', ' dF ')
	return(out)	
}