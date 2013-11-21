
# wrapper function for ecogStatistics

library('knitr')

# filepaths
fileName = '~/Google Drive/Research/ECoG/ecogScripts/ecog_R_AnalysisScripts/ECOGgroupStatsLM.Rmd'
outPath  = '~/Google Drive/Research/ECoG Manuscript/stats/'

# settings
lockType = 'stim'
dataType = 'power'
band     = 'hgam'
outFile  = paste(lockType,dataType,band,'LM',sep='')

# edit analysis script
text = readLines(fileName)
text[10] = paste0('lockType = "',lockType,'"');
text[11] = paste0('type = "'    ,dataType,'"');
text[12] = paste0('band = "'    ,band    ,'"');
writeLines(text,fileName)

knit2html(fileName, output = paste0(outPath,outFile,'.html'))

