Initialize = function()
{
  options(stringsAsFactors = F)
  
  #### Loading required libraries ####
  options(stringsAsFactors = F)
  
  #BiocManager::install('pacman')
  
  pac <- list('DESeq2', 'ggplot2', 'limma', 'pheatmap',
              'VennDiagram','e1071' ,'reshape2', 'ggrepel', 'RColorBrewer', 
              'plyr', 'gridExtra','gridExtra','caTools')
  
  print(paste(pac , lapply(pac, require, character.only = TRUE), sep = ' : '))
  
  pacman::p_load('DESeq2', 'ggplot2', 'limma', 'pheatmap',
                 'VennDiagram','e1071' ,'reshape2', 'ggrepel', 'RColorBrewer',
                 'plyr','gridExtra','gridExtra','caTools')
  
}




#### Functions ####
FindRepetitionIndex = function(aVector)
{
  i = 1
  repetitionIndex = integer(length(aVector))
  repetitionIndex[!duplicated(aVector)] = i
  aVector[!duplicated(aVector)] = NA
  
  while (!all(is.na(aVector))) {
    i = i+1
    repetitionIndex[!duplicated(aVector) & !is.na(aVector)] = i
    aVector[!duplicated(aVector)] = NA
  }
  
  repetitionIndex  
}



CalcSummaryForCellLines = function(selectedCounts)
{
  cellLinesName = unique(gsub(colnames(selectedCounts), pattern = '_.*', replacement =  ''))
  cellLinesSummary = lapply(cellLinesName, function(aCellLineName)
  {
    theCellLineCounts = selectedCounts[, grepl(colnames(selectedCounts), pattern = aCellLineName)]
    summary(rowMeans(theCellLineCounts))
  })
  names(cellLinesSummary) = cellLinesName
  
  cellLinesSummary
}
