
# wrapper function for ecogStatistics

library('knitr')

# filepaths
fileName = '~/Documents/ECOG/RScripts/ECOGgroupStatsLMER2.Rmd'
outPath  = '~/Google Drive/Research/ECoG Manuscript/stats/'

# settings
lockType = 'RT'
dataType = 'power'
band     = 'lgam'
outFile  = paste(lockType,dataType,band,sep='')

# edit analysis script
text = readLines(fileName)
text[10] = paste0('lockType = "',lockType,'"');
text[11] = paste0('type = "'    ,dataType,'"');
text[12] = paste0('band = "'    ,band    ,'"');
writeLines(text,fileName)

knit2html(fileName, output = paste0(outPath,outFile,'.html'))

