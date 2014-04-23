
# wrapper function for ecogStatistics

library('knitr')

# filepaths
fileNameLM   = '~/Documents/ECOG/RScripts/ECOGgroupStatsLM.Rmd'
fileNameLMER = '~/Documents/ECOG/RScripts/ECOGgroupStatsLMER.Rmd'
outPath  = '~/Google Drive/Research/ECoG Manuscript/stats/'
setwd(outPath)

# settings
lockType = 'stim'
dataType = 'power'
band     = 'hgam'
outFileLM  = paste(lockType,dataType,band,'LM',sep='')
outFileLMER  = paste(lockType,dataType,band,'LMER',sep='')

# edit LM analysis script
text = readLines(fileNameLM)
text[10] = paste0('lockType = "',lockType,'"');
text[11] = paste0('type = "'    ,dataType,'"');
text[12] = paste0('band = "'    ,band    ,'"');
writeLines(text,fileNameLM)
knit2html(fileNameLM, output = paste0(outFileLM,'.html'))
 
# edit LMER analysis scripts
text = readLines(fileNameLMER)
text[10] = paste0('lockType = "',lockType,'"');
text[11] = paste0('type = "'    ,dataType,'"');
text[12] = paste0('band = "'    ,band    ,'"');
writeLines(text,fileNameLMER)
knit2html(fileNameLMER, output = paste0(outFileLMER,'.html'))